"""Does the literature mention a (gene, disease) pair, and does any abstract state the association?

Two checks a pruner can ask of any pair, including one no source wrote (docs/PROTOCOL.md, second
pruning round):

  CODE ONLY   PubMed is searched for papers with the gene symbol in the title or abstract and the
              disease by PubMed's own mapping. The count of papers that mention both, and the lift
              (that count over the count expected if the two were mentioned independently).

  THE READER  For a pair some paper mentions, a model reads the first five records and answers
              `asserts`, with a PMID and a passage copied from that record, or `mentions_only`.
              Code checks the PMID was one it was given and the passage is in that record. The reader
              can only flag a pair the code-only check kept. It cannot restore one.

Every search, record and model answer is kept under data/eval/literature, so a pair is asked once.
NCBI asks callers to name a contact address: set `NCBI_EMAIL`. None is stored in this repository.

Run:  python -m pharos_gap_kg.litcheck --selftest
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from pharos_gap_kg import paths

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "pharos-gap-kg"
RECORDS, CUT, MIN_PASSAGE = 5, 1800, 20
EVERYTHING = "all[sb]"                       # every PubMed record: the denominator of the lift
OPERATORS = {"and", "or", "not"}
PUNCT = re.compile(r"[^\w\s'\-]", re.UNICODE)

BRIEF = """You check whether published abstracts state an association between a human gene and a disease.
You are given a gene symbol, a disease name, and up to five PubMed records (PMID, title, abstract). Decide from these records only, not from what you know.

Answer "asserts" if at least one record states that this gene, its transcript or its protein is associated with this disease or with a named form of it: a causal variant, a risk association, a biomarker, a change in expression or activity in patients or in a model of the disease, a functional role in the disease, or a drug target for it. Copy one sentence or clause from that record, exactly as it is written there, that states it.

Answer "mentions_only" if no record states such an association: the symbol means something else there, the gene and the disease both appear without being related to each other, or the record reports that no association was found.

Reply with one line of JSON and nothing else:
{"verdict": "asserts", "pmid": "<PMID of the record>", "quote": "<text copied exactly from that record>"}
or
{"verdict": "mentions_only"}"""
BRIEF_SHA = hashlib.sha256(BRIEF.encode()).hexdigest()[:12]


def lit_dir() -> Path:
    return paths.EVAL / "literature"


# ---- queries -----------------------------------------------------------------------------------------------

def clean_name(name: str) -> str:
    """A disease name with nothing PubMed would read as query syntax: no punctuation, no bare operators."""
    words = PUNCT.sub(" ", name).lower().split()
    return " ".join(w for w in words if w not in OPERATORS)


def gene_term(symbol: str) -> str:
    return f'"{symbol}"[tiab]'


def disease_term(name: str) -> str:
    return f"({clean_name(name)})"


def pair_term(symbol: str, name: str) -> str:
    return f"{gene_term(symbol)} AND {disease_term(name)}"


def parse_esearch(response: dict) -> dict:
    """count, the first ids, and every phrase PubMed says it could not find."""
    r = response.get("esearchresult", {})
    if "ERROR" in r or "count" not in r:
        return {"count": None, "ids": [], "not_found": [], "error": str(r.get("ERROR", "no count in the response"))[:200]}
    nf = list((r.get("errorlist") or {}).get("phrasesnotfound") or []) + list((r.get("warninglist") or {}).get("quotedphrasesnotfound") or [])
    return {"count": int(r["count"]), "ids": [str(i) for i in r.get("idlist", [])], "not_found": nf, "error": None}


class Pacer:
    """Hands out start times at least `gap` seconds apart, to any number of threads."""

    def __init__(self, gap: float) -> None:
        self.gap, self.lock, self.next = gap, threading.Lock(), 0.0

    def wait(self) -> None:
        with self.lock:
            now = time.monotonic()
            start = max(self.next, now)
            self.next = start + self.gap
        time.sleep(max(0.0, start - now))


_PACER = Pacer(0.36)          # NCBI allows three requests a second without a key


def _get(endpoint: str, q: dict) -> bytes:
    q = {**q, "tool": TOOL}
    if os.environ.get("NCBI_EMAIL"):
        q["email"] = os.environ["NCBI_EMAIL"]
    body = urllib.parse.urlencode(q).encode()
    last = None
    for attempt in range(6):
        _PACER.wait()
        try:
            with urllib.request.urlopen(urllib.request.Request(f"{EUTILS}/{endpoint}.fcgi", data=body), timeout=60) as r:
                return r.read()
        except Exception as e:  # rate limits, timeouts and resets are all retried
            last = e
            time.sleep(2 ** attempt)
    raise SystemExit(f"NCBI did not answer {endpoint} after 6 attempts ({last}); what was answered before this is kept, re-run to go on")


def esearch(term: str) -> dict:
    f = lit_dir() / "pubmed" / (hashlib.sha256(f"{term}|{RECORDS}|relevance".encode()).hexdigest()[:24] + ".json")
    if f.is_file():
        return json.loads(f.read_text())
    out = parse_esearch(json.loads(_get("esearch", {"db": "pubmed", "term": term, "retmode": "json", "retmax": RECORDS, "sort": "relevance"})))
    out["term"] = term
    if out["error"] is None:
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps(out) + "\n")
    return out


def search_many(terms: list[str], workers: int = 3, every: int = 500) -> dict[str, dict]:
    terms = sorted(set(terms))
    out: dict[str, dict] = {}
    with ThreadPoolExecutor(workers) as ex:
        for k, (t, r) in enumerate(zip(terms, ex.map(esearch, terms)), 1):
            out[t] = r
            if k % every == 0:
                print(f"  searched {k:,} / {len(terms):,}", flush=True)
    return out


def lift(both: int, gene: int, disease: int, total: int) -> float:
    """Papers that mention both, over the number expected if the two were mentioned independently."""
    if not both or not gene or not disease:
        return 0.0
    return both * total / (gene * disease)


# ---- records -----------------------------------------------------------------------------------------------

def parse_efetch(xml: bytes) -> dict[str, dict]:
    out = {}
    for art in ET.fromstring(xml).iter("PubmedArticle"):
        pmid = art.findtext("MedlineCitation/PMID")
        title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
        parts = []
        for a in art.iter("AbstractText"):
            text = "".join(a.itertext()).strip()
            parts.append(f"{a.get('Label')}: {text}" if a.get("Label") else text)
        out[str(pmid)] = {"pmid": str(pmid), "title": " ".join(title.split()), "abstract": " ".join(" ".join(parts).split())}
    return out


def records(pmids: list[str]) -> list[dict]:
    """Title and abstract for each id, in the order asked; an id PubMed returns nothing for is left out."""
    d = lit_dir() / "records"
    d.mkdir(parents=True, exist_ok=True)
    need = [p for p in pmids if not (d / f"{p}.json").is_file()]
    if need:
        got = parse_efetch(_get("efetch", {"db": "pubmed", "id": ",".join(need), "retmode": "xml"}))
        for p in need:
            (d / f"{p}.json").write_text(json.dumps(got.get(p, {"pmid": p, "title": "", "abstract": ""})) + "\n")
    out = [json.loads((d / f"{p}.json").read_text()) for p in pmids]
    return [r for r in out if r["title"] or r["abstract"]]


# ---- the reader --------------------------------------------------------------------------------------------

def shown(rec: dict) -> dict:
    """A record as the reader is given it: each field cut at a fixed length."""
    return {"pmid": rec["pmid"], "title": rec["title"][:CUT], "abstract": rec["abstract"][:CUT]}


def prompt(symbol: str, name: str, recs: list[dict]) -> str:
    body = "\n\n".join(f"PMID {r['pmid']}\nTitle: {r['title']}\nAbstract: {r['abstract'] or '(none)'}" for r in map(shown, recs))
    return f"Gene: {symbol}\nDisease: {name}\n\n{body}"


def normalise(s: str) -> str:
    s = s.translate(str.maketrans({"‘": "'", "’": "'", "“": '"', "”": '"', "–": "-", "—": "-", " ": " "}))
    return " ".join(s.lower().split()).strip(" .\"'")


def parse_answer(text: str) -> dict | None:
    for m in re.finditer(r"\{[^{}]*\}", text or "", re.S):
        try:
            a = json.loads(m.group(0))
        except ValueError:
            continue
        if isinstance(a, dict) and "verdict" in a:
            return a
    return None


def check(text: str, recs: list[dict]) -> dict:
    """The reader's answer as code accepts it: a verdict, or None with the reason it was refused."""
    a = parse_answer(text)
    if a is None:
        return {"verdict": None, "why": "no JSON answer"}
    if a.get("verdict") == "mentions_only":
        return {"verdict": "mentions_only"}
    if a.get("verdict") != "asserts":
        return {"verdict": None, "why": f"verdict {str(a.get('verdict'))[:30]!r} is neither of the two"}
    pmid = str(a.get("pmid", "")).strip()
    rec = next((shown(r) for r in recs if r["pmid"] == pmid), None)
    if rec is None:
        return {"verdict": None, "why": f"PMID {pmid[:20]!r} is not one of the records given"}
    q = normalise(str(a.get("quote", "")))
    if len(q) < MIN_PASSAGE:
        return {"verdict": None, "why": f"the passage is under {MIN_PASSAGE} characters"}
    if q not in normalise(rec["title"] + " " + rec["abstract"]):
        return {"verdict": None, "why": f"the passage is not in PMID {pmid}"}
    return {"verdict": "asserts", "pmid": pmid, "quote": str(a["quote"])}


def spent() -> float:
    d = lit_dir() / "calls"
    return sum(json.loads(f.read_text()).get("cost", 0.0) for f in d.glob("*.json")) if d.is_dir() else 0.0


def read(symbol: str, name: str, recs: list[dict], model: str, key: str | None, base_url: str, prices: tuple[float, float] | None = None, post=None) -> dict:
    """One pair read by one model: asked, checked, asked once more if the check refuses, and kept."""
    if post is None:
        from pharos_gap_kg.usability import post
    ids = [r["pmid"] for r in recs]
    f = lit_dir() / "calls" / (hashlib.sha256(f"{BRIEF_SHA}|{model}|{symbol}|{name}|{','.join(ids)}".encode()).hexdigest()[:24] + ".json")
    if f.is_file():
        return {**json.loads(f.read_text()), "cached": True}
    msgs = [{"role": "system", "content": BRIEF}, {"role": "user", "content": prompt(symbol, name, recs)}]
    row = {"model": model, "symbol": symbol, "name": name, "pmids": ids, "brief": BRIEF_SHA, "in": 0, "out": 0, "cost": 0.0, "answers": [], "error": None}
    res = {"verdict": None, "why": "not asked"}
    for attempt in range(2):
        body = {"model": model, "messages": msgs, "max_tokens": 300, "temperature": 0}
        if "openrouter.ai" in base_url:
            body["usage"] = {"include": True}
        r = post(base_url, key, body)
        if "error" in r or not r.get("choices"):
            row["error"] = str(r.get("error"))[:300]
            break
        u = r.get("usage") or {}
        i, o, c = int(u.get("prompt_tokens") or 0), int(u.get("completion_tokens") or 0), float(u.get("cost") or 0.0)
        row["in"] += i
        row["out"] += o
        row["cost"] += c if c or not prices else (i * prices[0] + o * prices[1]) / 1e6
        text = r["choices"][0]["message"].get("content") or ""
        row["answers"].append(text[:1200])
        res = check(text, recs)
        if res["verdict"]:
            break
        msgs += [{"role": "assistant", "content": text or "(no text)"},
                 {"role": "user", "content": f"That answer was refused: {res['why']}. Answer again in the same format. A quote must be copied exactly from the record you cite."}]
    row.update({"verdict": res["verdict"], "pmid": res.get("pmid"), "quote": res.get("quote"), "why_refused": res.get("why"), "attempts": len(row["answers"])})
    if row["error"] is None:
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps(row, indent=1) + "\n")
    return row


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("a disease name loses its punctuation and any bare operator",
       clean_name("Leukemia, Myeloid, Acute (AML)") == "leukemia myeloid acute aml" and clean_name("Deafness, not otherwise specified") == "deafness otherwise specified"
       and clean_name("Crohn's disease") == "crohn's disease" and clean_name('46,XY "sex reversal" [type 3]: a OR b') == "46 xy sex reversal type 3 a b")
    ck("the pair query asks for the symbol in the title or abstract and leaves the disease to PubMed's mapping",
       pair_term("APOE", "Alzheimer disease") == '"APOE"[tiab] AND (alzheimer disease)')
    good = {"esearchresult": {"count": "42", "idlist": ["11", "22"], "errorlist": {"phrasesnotfound": [], "fieldsnotfound": []}}}
    miss = {"esearchresult": {"count": "7", "idlist": [], "errorlist": {"phrasesnotfound": ["zzyzx"]}, "warninglist": {"quotedphrasesnotfound": ['"QQQ1"[tiab]']}}}
    ck("a count and its ids are read", parse_esearch(good) == {"count": 42, "ids": ["11", "22"], "not_found": [], "error": None})
    ck("phrases PubMed could not find are reported, quoted or not", parse_esearch(miss)["not_found"] == ["zzyzx", '"QQQ1"[tiab]'])
    ck("an error response has no count", parse_esearch({"esearchresult": {"ERROR": "bad"}})["count"] is None and parse_esearch({})["error"] is not None)
    ck("the lift is observed over expected, and zero when no paper mentions both", lift(10, 100, 1000, 1_000_000) == 100.0 and lift(0, 5, 5, 100) == 0.0 and lift(3, 0, 5, 100) == 0.0)
    xml = b"""<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID Version="1">11</PMID><Article><ArticleTitle>APOE and <i>late</i> onset.</ArticleTitle>
<Abstract><AbstractText Label="RESULTS">The APOE e4 allele   raised risk.</AbstractText><AbstractText>Second part.</AbstractText></Abstract></Article></MedlineCitation></PubmedArticle>
<PubmedArticle><MedlineCitation><PMID>22</PMID><Article><ArticleTitle>No abstract here</ArticleTitle></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"""
    recs = parse_efetch(xml)
    ck("a record's title keeps its inline markup's text, and a labelled abstract keeps its labels",
       recs["11"] == {"pmid": "11", "title": "APOE and late onset.", "abstract": "RESULTS: The APOE e4 allele raised risk. Second part."} and recs["22"]["abstract"] == "")
    rl = [recs["11"], recs["22"]]
    ck("the reader is shown the gene, the disease and each record", prompt("APOE", "Alzheimer disease", rl).startswith("Gene: APOE\nDisease: Alzheimer disease\n\nPMID 11\nTitle: APOE and late onset.") and "Abstract: (none)" in prompt("APOE", "x", rl))
    ck("an assertion with a passage from the cited record is accepted, whatever its case and spacing",
       check('{"verdict": "asserts", "pmid": "11", "quote": "the APOE e4  allele raised risk."}', rl)["verdict"] == "asserts")
    ck("a passage that is not in the cited record is refused", check('{"verdict": "asserts", "pmid": "11", "quote": "APOE causes the disease in all carriers"}', rl) == {"verdict": None, "why": "the passage is not in PMID 11"})
    ck("a passage from another record, a PMID not given and a short passage are refused",
       check('{"verdict": "asserts", "pmid": "22", "quote": "The APOE e4 allele raised risk."}', rl)["verdict"] is None
       and check('{"verdict": "asserts", "pmid": "99", "quote": "The APOE e4 allele raised risk."}', rl)["why"].startswith("PMID '99'")
       and check('{"verdict": "asserts", "pmid": "11", "quote": "raised risk"}', rl)["why"].startswith("the passage is under"))
    ck("mentions_only needs nothing more, and text around the JSON is ignored", check('Sure.\n{"verdict": "mentions_only"}\n', rl) == {"verdict": "mentions_only"})
    ck("no JSON, or a third verdict, is refused", check("I think it asserts.", rl)["verdict"] is None and check('{"verdict": "supported"}', rl)["verdict"] is None)
    long = {"pmid": "33", "title": "t", "abstract": "a" * CUT + " the tail states the association clearly"}
    ck("a passage beyond the part the reader was shown is refused", check('{"verdict": "asserts", "pmid": "33", "quote": "the tail states the association clearly"}', [long])["verdict"] is None)

    import tempfile
    answers = iter(['{"verdict": "asserts", "pmid": "11", "quote": "APOE is the cause of everything"}', '{"verdict": "asserts", "pmid": "11", "quote": "The APOE e4 allele raised risk."}'])
    fake = lambda base, key, body: {"choices": [{"message": {"content": next(answers)}}], "usage": {"prompt_tokens": 100, "completion_tokens": 10}}
    old = paths.EVAL
    with tempfile.TemporaryDirectory() as d:
        paths.EVAL = Path(d)
        try:
            row = read("APOE", "Alzheimer disease", rl, "m", None, "http://x", prices=(1.0, 2.0), post=fake)
            ck("a refused answer is asked for once more, and the second is checked the same way", row["attempts"] == 2 and row["verdict"] == "asserts" and row["pmid"] == "11")
            ck("tokens are priced when the provider reports no cost", abs(row["cost"] - 2 * (100 * 1.0 + 10 * 2.0) / 1e6) < 1e-12 and abs(spent() - row["cost"]) < 1e-12)
            again = read("APOE", "Alzheimer disease", rl, "m", None, "http://x", post=lambda *a: (_ for _ in ()).throw(AssertionError("asked twice")))
            ck("a pair a model has read is not asked again, and says so", again == {**row, "cached": True})
            bad = read("APOE", "other", rl, "m", None, "http://x", post=lambda *a: {"choices": [{"message": {"content": "no idea"}}], "usage": {}})
            ck("two refused answers leave no verdict", bad["verdict"] is None and bad["attempts"] == 2 and bad["why_refused"] == "no JSON answer")
        finally:
            paths.EVAL = old
    return 0 if ok else 1


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    raise SystemExit(__doc__)
