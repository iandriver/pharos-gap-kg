"""Does full text reach the disease sources that abstracts could not? A check under a daily allowance.

The second pruning round (findings/PRUNING2.md) left three sources that no rule cleans: more than half
of their real pairs have no PubMed abstract that mentions the gene with the disease. This asks the
same question of full text, through Paperclip (paperclip.gxl.ai), which searches PubMed Central and
preprints. It is an optional source: nothing else in the package needs it, and it needs a key
(`PAPERCLIP_API_KEY`) that is never written to this repository.

THE ALLOWANCE is what shapes the test. The free service answers 100 requests a day and 25 a minute.
Every request is written to a log before it is sent, an answer is kept and never asked for twice, and a
run stops when the day's share is used and goes on another day from where it stopped.

TWO FORMS OF THE QUESTION were tried on hand-picked pairs outside the sample (`pilot`), and one was
fixed in docs/PROTOCOL.md before a sample pair was searched:

    doc    some full text holds the gene symbol and every word of the disease name, anywhere in it
    near   some full text holds the symbol within 200 characters, on one line, of the disease's name
           or of another name the ontology gives it as exact

THE FORM FIXED is `near`, as an exhaustive scan. The whole-paper form needs every word of one name
and stops counting at 20.

Run:  python -m pharos_gap_kg.fulltext pilot pharos319
      python -m pharos_gap_kg.fulltext run pharos319 --calls 30    # the next pairs in the fixed order
      python -m pharos_gap_kg.fulltext report pharos319
      python -m pharos_gap_kg.fulltext --selftest
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import litcheck, ontology, paths, prune2, registry

MCP = "https://paperclip.gxl.ai/mcp"
DAY, MINUTE = 100, 25           # the free allowance, as Ian gave it on 9 Oct 2026
PER_DAY, GAP = 90, 3.0          # what this module lets itself use: some of the day is left, and 20 a minute
SOURCES = ("Expression Atlas", "DrugCentral Indication", "JensenLab Knowledge AmyCo")
FIRST = 50                      # of each source: the first 50 planted and the first 50 real pairs, in drawn order
WINDOW = 200
MAX_NAMES, MIN_NAME = 12, 5     # how many names of a disease go into one question, and the shortest one used
PILOT_SEED, PILOT_EACH = 20261010, 3
KEY = prune2.KEY


class Spent(Exception):
    """The day's allowance, ours or the service's, is used."""


def home() -> Path:
    return paths.EVAL / "paperclip"


def today() -> str:
    return time.strftime("%Y-%m-%d", time.gmtime())


def sent_today() -> int:
    f = home() / "calls.jsonl"
    return sum(1 for line in f.read_text().splitlines() if json.loads(line)["at"].startswith(today())) if f.is_file() else 0


def note(what: str) -> None:
    """Write a request down before it is sent: a request that fails still used the allowance."""
    home().mkdir(parents=True, exist_ok=True)
    with (home() / "calls.jsonl").open("a") as f:
        f.write(json.dumps({"at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "what": what}) + "\n")


# ---- the two forms of the question -------------------------------------------------------------------------

def literal(s: str) -> str:
    """Text as a regular expression that matches it exactly, with any run of spaces matching one or more."""
    return r"\s+".join(re.sub(r"([^\w\s])", r"\\\1", w) for w in s.split())


def doc_question(symbol: str, name: str) -> tuple[str, dict]:
    return "search", {"query": f"{symbol} {litcheck.clean_name(name)}", "source": "pmc", "mode": "all", "count_only": True}


def exact_synonyms(lines) -> dict[str, list[str]]:
    """id -> the names an OBO file gives a term as EXACT synonyms."""
    out: dict[str, list[str]] = {}
    cur = None
    for line in lines:
        if line.startswith("["):
            cur = None
        elif line.startswith("id: "):
            cur = out.setdefault(ontology.norm(line[4:].strip()), [])
        elif cur is not None and line.startswith('synonym: "'):
            m = re.match(r'synonym: "(.*)" EXACT\b', line)
            if m:
                cur.append(m.group(1))
    return out


def names_of(disease: str, name: str, terms: dict, synonyms: dict) -> list[str]:
    """The names a disease is searched under: the release's own, the ontology's label and its exact
    synonyms, cleaned, without repeats, none shorter than MIN_NAME characters, at most MAX_NAMES."""
    c = ontology.canonical(disease, terms)
    extra = ([terms[c]["name"]] + sorted(synonyms.get(c, []))) if c else []
    seen, out = set(), []
    for n in [name] + sorted(extra, key=lambda x: (len(x.split()), x.lower())):
        k = litcheck.clean_name(n) if isinstance(n, str) else ""
        if len(k) >= MIN_NAME and k not in seen:
            seen.add(k)
            out.append(k)
    return out[:MAX_NAMES]


def near_question(symbol: str, names: list[str]) -> tuple[str, dict]:
    g, d = literal(symbol), "(" + "|".join(literal(n) for n in names) + ")"
    return "grep", {"pattern": f"{g}.{{0,{WINDOW}}}{d}|{d}.{{0,{WINDOW}}}{g}", "path": "/papers/", "count_only": True, "ignore_case": True, "exhaustive": True}




def parse_count(tool: str, text: str) -> int | None:
    """How many papers the answer reports; None when the answer is not a count (an error, a refusal)."""
    t = (text or "").strip()
    if tool == "search":
        if t.startswith("No papers found"):
            return 0
        m = re.match(r"Found (\d+) papers?", t)
        return int(m.group(1)) if m else None
    if t.startswith("No matches for"):
        return 0
    m = re.match(r"(\d+)\s*$", t.splitlines()[0]) if t else None
    return int(m.group(1)) if m else None


# ---- the service -------------------------------------------------------------------------------------------

class Client:
    """Paperclip's hosted tools over HTTP. Read tools only: search and grep."""

    def __init__(self, per_day: int = PER_DAY) -> None:
        self.key = os.environ.get("PAPERCLIP_API_KEY")
        if not self.key:
            raise SystemExit("set PAPERCLIP_API_KEY (a Paperclip key; it is read from the environment and stored nowhere)")
        self.per_day, self.sid, self.last, self.n = per_day, None, 0.0, 0

    def _post(self, body: dict) -> tuple[int, dict, object]:
        wait = GAP - (time.monotonic() - self.last)
        if wait > 0:
            time.sleep(wait)
        self.last = time.monotonic()
        h = {"X-API-Key": self.key, "Content-Type": "application/json", "Accept": "application/json, text/event-stream"}
        if self.sid:
            h["Mcp-Session-Id"] = self.sid
        try:
            with urllib.request.urlopen(urllib.request.Request(MCP, data=json.dumps(body).encode(), headers=h), timeout=180) as r:
                raw, hd, st = r.read().decode("utf-8", "replace"), dict(r.headers), r.status
        except urllib.error.HTTPError as e:
            return e.code, dict(e.headers), e.read()[:300].decode("utf-8", "replace")
        if "event-stream" in (hd.get("Content-Type") or hd.get("content-type") or ""):
            data = [line[5:].strip() for line in raw.splitlines() if line.startswith("data:")]
            raw = data[-1] if data else raw
        try:
            return st, hd, json.loads(raw)
        except ValueError:
            return st, hd, raw[:300]

    def _open(self) -> None:
        note("initialize")
        st, hd, r = self._post({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                                "params": {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "pharos-gap-kg", "version": "0"}}})
        if st == 429:
            raise Spent("the service refused the session (429)")
        if st != 200:
            raise SystemExit(f"Paperclip did not open a session: HTTP {st}: {r}")
        self.sid = hd.get("Mcp-Session-Id") or hd.get("mcp-session-id")
        self._post({"jsonrpc": "2.0", "method": "notifications/initialized"})

    def ask(self, tool: str, args: dict) -> dict:
        """One read-only tool call, kept once answered. Raises Spent when the day's share is used."""
        if tool not in ("search", "grep"):
            raise SystemExit(f"REFUSED: this module calls Paperclip's search and grep only, not {tool!r}")
        f = home() / "answers" / (hashlib.sha256(json.dumps([tool, args], sort_keys=True).encode()).hexdigest()[:24] + ".json")
        if f.is_file():
            return {**json.loads(f.read_text()), "kept": True}
        if sent_today() + (2 if self.sid is None else 1) > self.per_day:
            raise Spent(f"{sent_today()} requests were sent today and this module allows itself {self.per_day}")
        if self.sid is None:
            self._open()
        note(f"{tool} {json.dumps(args, sort_keys=True)}")
        self.n += 1
        st, _, r = self._post({"jsonrpc": "2.0", "id": 10 + self.n, "method": "tools/call", "params": {"name": tool, "arguments": args}})
        if st == 429:
            raise Spent("the service says its allowance is used (429)")
        res = r.get("result") if isinstance(r, dict) else None
        if st != 200 or not isinstance(res, dict) or "content" not in res:
            return {"tool": tool, "args": args, "text": None, "count": None, "error": f"HTTP {st}: {str(r)[:200]}"}
        text = "\n".join(c.get("text", "") for c in res["content"])
        out = {"tool": tool, "args": args, "text": text[:4000], "count": parse_count(tool, text), "at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "error": None}
        if out["count"] is None:
            out["error"] = "the answer is not a count: " + text[:160].replace("\n", " ")
            return out
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps(out) + "\n")
        return out


# ---- the pilot: hand-picked pairs outside the sample -------------------------------------------------------

def pilot_pairs(release: str, each: int = PILOT_EACH, seed: int = PILOT_SEED) -> pd.DataFrame:
    """A few Expression Atlas pairs that are in neither the sample nor the planted set: real one-source
    pairs, and crossings of them that no source asserts. They fix the form of the question and are
    never scored."""
    links = pd.read_parquet(paths.context(release) / "disease_links.parquet")
    pl, s = prune2.frozen(release)
    name_of, symbol_of = prune2.names(release)
    single = prune2.one_source(links)
    g = single[single.source == SOURCES[0]].reset_index(drop=True)
    used = set(zip(s.target, s.disease)) | set(zip(pl.target, pl.disease))
    rng = np.random.default_rng(seed)
    rows, real = [], []
    for i in rng.permutation(len(g)):
        t, d = g.target[i], g.disease[i]
        if (t, d) not in used and t in symbol_of:
            real.append((t, d))
            rows.append((t, d, False))
        if len(real) == 2 * each:
            break
    asserted = set(zip(links.target, links.disease))
    crossed = [(real[k][0], real[k + each][1]) for k in range(each)]
    crossed = [c for c in crossed if c not in asserted and c not in used]
    out = pd.DataFrame(rows[:each] + [(t, d, True) for t, d in crossed], columns=["target", "disease", "crossed"])
    out["symbol"] = out.target.map(symbol_of)
    out["name"] = [name_of.get((SOURCES[0], d)) for d in out.disease]
    return out


def with_names(pairs: pd.DataFrame) -> pd.DataFrame:
    """`pairs` (disease, name) with `names`, the list each disease is searched under."""
    path = registry.require("efo_obo")
    terms = ontology.read(path)
    with path.open(encoding="utf-8", errors="replace") as f:
        syn = exact_synonyms(f)
    return pairs.assign(names=[names_of(d, n, terms, syn) for d, n in zip(pairs.disease, pairs.name)])


def pilot(release: str) -> pd.DataFrame:
    c = Client()
    p = with_names(pilot_pairs(release))
    p["doc"] = [c.ask(*doc_question(sym, name)).get("count") for sym, name in zip(p.symbol, p.name)]
    p["near"] = [c.ask(*near_question(sym, names)).get("count") for sym, names in zip(p.symbol, p.names)]
    (home() / "pilot.json").write_text(p.to_json(orient="records", indent=1) + "\n")
    return p


# ---- the sample, in a fixed order, a day's share at a time --------------------------------------------------

def ordered(sample: pd.DataFrame, first: int = FIRST) -> pd.DataFrame:
    """The sample pairs this check asks about, in the order it asks: source by source, and within a source
    planted and real pairs in turn by drawn order, so that wherever a day stops the two groups are level."""
    s = sample[sample.source.isin(SOURCES) & (sample.order < first)].copy()
    s["_src"] = s.source.map({x: i for i, x in enumerate(SOURCES)})
    return s.sort_values(["_src", "order", "planted"], ascending=[True, True, False]).drop(columns="_src").reset_index(drop=True)


def run(release: str, calls: int) -> pd.DataFrame:
    """Ask about the next pairs, at most `calls` new requests, and write what is known so far."""
    _, sample = prune2.frozen(release)
    plan = with_names(ordered(sample))
    c = Client()
    rows, new, why = [], 0, None
    for r in plan.itertuples():
        row = {k: getattr(r, k) for k in KEY + ["planted", "order"]}
        if not isinstance(r.symbol, str) or not r.names:
            rows.append({**row, "count": None, "state": "unsearchable"})
            continue
        q = near_question(r.symbol, r.names)
        if why is None and new >= calls:
            why = f"the {calls} requests asked for were sent"
        try:
            a = c.ask(*q) if why is None else Client(per_day=0).ask(*q)      # once stopped, only kept answers are read
        except Spent as e:
            why = why or str(e)
            rows.append({**row, "count": None, "state": "not asked yet"})
            continue
        new += 0 if a.get("kept") else 1
        rows.append({**row, "count": a["count"], "state": "answered" if a["count"] is not None else "error: " + str(a["error"])[:80]})
    out = pd.DataFrame(rows)
    out["count"] = out["count"].astype("Int64")
    W = prune2.work(release)
    out.to_parquet(W / "fulltext.parquet", index=False)
    print(f"  {new} new requests; stopped because {why or 'every pair is answered'}; requests today {sent_today()}")
    return out


def report(release: str) -> dict:
    W = prune2.work(release)
    ft = pd.read_parquet(W / "fulltext.parquet")
    se = pd.read_parquet(W / "search.parquet")[KEY + ["planted", "both", "unsearchable"]]
    x = ft.merge(se, on=KEY + ["planted"])
    done = x[x.state == "answered"].copy()
    done["f"] = (done["count"] == 0).to_numpy()
    done["l0"] = (done.both == 0).to_numpy()
    out = {"release": release, "first": FIRST, "window": WINDOW, "sample_sha256": prune2.graph.file_sha(W / "sample.parquet"),
           "states": {str(k): int(v) for k, v in ft.state.str.replace(r"^error.*", "error", regex=True).value_counts().items()}, "sources": {}}
    for src in SOURCES:
        g = done[done.source == src]
        p, r = g[g.planted], g[~g.planted]
        if not len(g):
            continue
        fp, fr = prune2.share(p.f.to_numpy()), prune2.share(r.f.to_numpy())
        out["sources"][src] = {"planted_withheld": fp, "real_withheld": fr, **prune2.bar(fp, fr), "complete": bool(len(p) == FIRST and len(r) == FIRST),
                               "abstracts_on_the_same_pairs": {"planted_withheld": prune2.share(p.l0.to_numpy()), "real_withheld": prune2.share(r.l0.to_numpy())},
                               "real_pairs_no_abstract_mentions": {"n": int(r.l0.sum()), "found_in_full_text": int((r.l0 & ~r.f).sum())},
                               "planted_pairs_no_abstract_mentions": {"n": int(p.l0.sum()), "found_in_full_text": int((p.l0 & ~p.f).sum())},
                               "real_pairs_an_abstract_mentions": {"n": int((~r.l0).sum()), "not_found_in_full_text": int((~r.l0 & r.f).sum())}}
    (paths.DOCS_DATA / f"pruning2_fulltext_{release}.json").write_text(prune2.dump(out))
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("a search answer is read as its count, none found as zero, and anything else as no count",
       [parse_count("search", t) for t in ("Found 14 papers  [s_1]\n\n  1. A title", "Found 1 paper [s_2]", "No papers found.\nTry broader terms", "Error: search requires a source flag (-s).")] == [14, 1, 0, None])
    ck("a grep answer is read as its count, no match as zero, and an error as no count",
       [parse_count("grep", t) for t in ("100\n", "0\n", "No matches for /A.{0,200}(b)/ across corpus\n[exit 1]\n", "ERR: vsh: grep: X: Cannot read path", "")] == [100, 0, 0, None, None])
    ck("text becomes an exact pattern, whatever punctuation it holds", literal("HLA-A") == r"HLA\-A" and literal("crohn's disease") == r"crohn\'s\s+disease" and re.fullmatch(literal("a.b c"), "a.b   c") and not re.fullmatch(literal("a.b"), "axb"))
    tool, args = doc_question("CENPV", "Lung carcinoma, NOS")
    ck("the whole-paper form asks for the symbol and every word of the name in PubMed Central", tool == "search" and args == {"query": "CENPV lung carcinoma nos", "source": "pmc", "mode": "all", "count_only": True})
    obo = '''[Term]
id: MONDO:1
name: gastric cancer
synonym: "stomach cancer" EXACT []
synonym: "GC" EXACT []
synonym: "cancer of stomach" EXACT [NCIT:1]
synonym: "gastric tumour" RELATED []

[Term]
id: MONDO:2
name: old name
is_obsolete: true
replaced_by: MONDO:1

[Typedef]
id: part_of
synonym: "x" EXACT []
'''.splitlines()
    syn, terms = exact_synonyms(obo), ontology.load(obo)
    ck("exact synonyms are read per term, and related ones and other stanzas are not", syn["MONDO:1"] == ["stomach cancer", "GC", "cancer of stomach"] and syn.get("part_of", []) == ["x"] and "gastric tumour" not in syn["MONDO:1"])
    nm = names_of("MONDO:2", "Gastric Cancer", terms, syn)
    ck("a disease is searched under the release's name, the label and exact synonyms of the term that replaced it, short names left out",
       nm == ["gastric cancer", "stomach cancer", "cancer of stomach"] and names_of("UMLS:9", "Some disease", terms, syn) == ["some disease"])
    tool, args = near_question("PLOD1", nm)
    ck("the nearby form asks for the symbol within 200 characters of any of the names, either order, any case",
       tool == "grep" and args["ignore_case"] and args["count_only"] and args["exhaustive"] and re.search(args["pattern"], "plod1 was high in stomach  cancer", re.I)
       and re.search(args["pattern"], "Gastric cancer cells lacking PLOD1", re.I) and not re.search(args["pattern"], "PLOD1 " + "x" * 300 + " gastric cancer", re.I)
       and not re.search(args["pattern"], "PLOD1 in gastric tumour", re.I))
    sm = pd.DataFrame([(s, pl_, o) for s in (SOURCES[1], SOURCES[0], "CTD") for pl_ in (False, True) for o in (2, 0, 1)], columns=["source", "planted", "order"])
    od = ordered(sm, first=2)
    ck("pairs are asked source by source, planted and real in turn by drawn order, and only the first of each",
       list(zip(od.source, od.planted, od.order)) == [(SOURCES[0], True, 0), (SOURCES[0], False, 0), (SOURCES[0], True, 1), (SOURCES[0], False, 1),
                                                       (SOURCES[1], True, 0), (SOURCES[1], False, 0), (SOURCES[1], True, 1), (SOURCES[1], False, 1)])
    import tempfile
    old = paths.EVAL
    with tempfile.TemporaryDirectory() as d:
        paths.EVAL = Path(d)
        try:
            ck("no request has been sent before the first is written down", sent_today() == 0)
            note("a")
            note("b")
            ck("a request is counted from the moment it is written down", sent_today() == 2)
            os.environ["PAPERCLIP_API_KEY"] = os.environ.get("PAPERCLIP_API_KEY") or "test"
            c = Client(per_day=3)
            try:
                c.ask("rm", {"path": "/x"})
                ck("only search and grep are ever called", False)
            except SystemExit:
                ck("only search and grep are ever called", True)
            try:
                c.ask("search", {"query": "x", "source": "pmc"})
                ck("a call that would pass the day's share is not sent", False)
            except Spent:
                ck("a call that would pass the day's share is not sent", sent_today() == 2)
            f = home() / "answers" / (hashlib.sha256(json.dumps(["grep", {"pattern": "p"}], sort_keys=True).encode()).hexdigest()[:24] + ".json")
            f.parent.mkdir(parents=True)
            f.write_text(json.dumps({"tool": "grep", "args": {"pattern": "p"}, "text": "7\n", "count": 7, "error": None}))
            ck("an answer already kept is returned without a request", Client(per_day=0).ask("grep", {"pattern": "p"})["count"] == 7 and sent_today() == 2)
        finally:
            paths.EVAL = old
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) < 2 or argv[0] not in ("pilot", "run", "report"):
        raise SystemExit(__doc__)
    if argv[0] == "run":
        calls = int(argv[argv.index("--calls") + 1]) if "--calls" in argv else PER_DAY
        out = run(argv[1], calls)
        print(out.groupby(["source", "planted"]).state.apply(lambda v: int((v == "answered").sum())).unstack().to_string())
        return 0
    if argv[0] == "report":
        out = report(argv[1])
        pc = lambda d: f"{d['share']:.1%} [{d['low']:.1%}, {d['high']:.1%}] n={d['n']}"
        print(f"  states: {out['states']}")
        for src, r in out["sources"].items():
            a = r["abstracts_on_the_same_pairs"]
            print(f"  {src:28s} full text: planted {pc(r['planted_withheld'])}  real {pc(r['real_withheld'])}  {'PASS' if r['passes'] else 'no'}{' (provisional)' if r['provisional'] else ''}{'' if r['complete'] else '  [not complete]'}")
            print(f"  {'':28s} abstracts: planted {pc(a['planted_withheld'])}  real {pc(a['real_withheld'])}; real pairs no abstract mentions: {r['real_pairs_no_abstract_mentions']}; planted: {r['planted_pairs_no_abstract_mentions']}")
        return 0
    if argv[0] == "pilot":
        try:
            p = pilot(argv[1])
        except Spent as e:
            raise SystemExit(f"stopped: {e}. Requests sent today: {sent_today()}")
        print(p[["crossed", "symbol", "name", "doc", "near"]].to_string(index=False))
        print(f"  requests sent today by this module: {sent_today()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
