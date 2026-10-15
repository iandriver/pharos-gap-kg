"""Usability with a person: twelve frozen tasks, asked through an MCP client, graded by code.

The definitions are in docs/PROTOCOL.md, "Usability with a person", and the session is run from
docs/PERSON_TEST.md. The tasks are generated with a fixed seed and frozen before the first session.
Expected answers are computed from the bundle's tables without the question layer, and `check` asks
the question layer each task and refuses a set where the two disagree. No disease, target or drug in
a task appears in README.md or docs/ASK.md.

Run:  python -m pharos_gap_kg.persontest freeze pharos319     # write docs/data/person_tasks.json, once
      python -m pharos_gap_kg.persontest check pharos319      # the question layer gives every expected answer
      python -m pharos_gap_kg.persontest sheet                # the participant's task sheet, without answers
      python -m pharos_gap_kg.persontest template             # the file a facilitator fills in during a session
      python -m pharos_gap_kg.persontest grade answers.json   # one session's answers -> docs/data/person_test_<P>.json
      python -m pharos_gap_kg.persontest --selftest
"""
from __future__ import annotations

import hashlib
import json
import random
import re
import statistics
import sys
from pathlib import Path

import pandas as pd

from pharos_gap_kg import paths

SEED = 20261019
FAKE = "QZRV7"
LIMIT_MINUTES, MEDIAN_MINUTES, RIGHT_AT_LEAST = 8.0, 5.0, 9
SIZE_WITHIN, OR_WITHIN = 0.05, 0.01
SUS_AVERAGE = 68
SUS_ITEMS = (
    "I think that I would like to use this system frequently.",
    "I found the system unnecessarily complex.",
    "I thought the system was easy to use.",
    "I think that I would need the support of a technical person to be able to use this system.",
    "I found the various functions in this system were well integrated.",
    "I thought there was too much inconsistency in this system.",
    "I would imagine that most people would learn to use this system very quickly.",
    "I found the system very cumbersome to use.",
    "I felt very confident using the system.",
    "I needed to learn a lot of things before I could get going with this system.",
)
OPEN = ("What was hardest to find, or to trust?", "Was there an answer you did not believe? Which, and why?",
        "What would you ask this graph in your own work?")


def tasks_path() -> Path:
    return paths.DOCS_DATA / "person_tasks.json"


def digest(tasks: list[dict]) -> str:
    return hashlib.sha256(json.dumps(tasks, sort_keys=True).encode()).hexdigest()


def _mentioned(text: str):
    t = text.casefold()
    return lambda name: bool(re.search(r"(?<![\w-])" + re.escape(str(name).casefold()) + r"(?![\w-])", t))


def make_tasks(bundle: Path, seed: int = SEED) -> list[dict]:
    from pharos_gap_kg import ask, gaptables
    G = ask.Graph(bundle)
    rng = random.Random(seed)
    T = gaptables.read(bundle)
    nodes, links, cn = T["nodes"], T["links"], T["context_nodes"]
    t, d = nodes[nodes.kind == "target"], nodes[nodes.kind == "trait"]
    in_docs = _mentioned((paths.ROOT / "README.md").read_text() + (paths.ROOT / "docs" / "ASK.md").read_text())
    plain = lambda r: r["ok"] and not r["note"]
    n_sym = t.groupby("name").id.transform("size")
    mhc = (t.chr == ask.MHC[0]) & (t.start >= ask.MHC[1]) & (t.start <= ask.MHC[2])
    clean_t = {r.id for r in t[(n_sym == 1) & ~mhc].itertuples() if not in_docs(r.name)
               and (lambda x: plain(x) and x["kind"] == "target" and G.t.id[x["i"]] == r.id)(G.resolve(r.name))}
    deg = links.trait.value_counts()
    clean_d = {r.id for r in d[d.disease_like].itertuples() if deg.get(r.id, 0) > 0 and not in_docs(r.name)
               and (lambda x: plain(x) and x["kind"] == "trait" and G.d.id[x["i"]] == r.id)(G.resolve(r.name))}
    cond = cn[cn.kind == "condition"]
    one_name = cond.groupby(cond.name.map(gaptables.norm)).id.transform("size") == 1
    clean_c = {r.id for r in cond[one_name].itertuples() if not in_docs(r.name)
               and (lambda x: plain(x) and x["kind"] == "condition" and x["key"] == r.id)(G.resolve(r.name))}
    R = gaptables.build(T, only=clean_d | clean_c)
    S, summ = R["same_condition"], R["disease_gap_summary"].set_index("condition_id")
    sym, name = dict(zip(t.id, t.name)), dict(zip(d.id, d.name)) | dict(zip(cn.id, cn.name))
    asked_of = S[S.asked_id.isin(clean_d | clean_c)].groupby("read_id").asked_id.apply(list).to_dict()
    used_c, used_t = set(), set()
    out: list[dict] = []

    def add(kind: str, text: str, record: str, gold: dict) -> None:
        out.append({"id": f"T{len(out) + 1}", "kind": kind, "text": text, "record": record, "gold": gold})

    def disease(keep) -> str:
        ok = sorted(c for c in clean_d if c not in used_c and c in summ.index and keep(summ.loc[c]))
        c = rng.choice(ok)
        used_c.add(c)
        return c

    def pair(rows: pd.DataFrame, keep) -> tuple[str, str, dict]:
        for r in rows.sample(frac=1.0, random_state=rng.randrange(2 ** 31)).itertuples():
            if r.target not in clean_t or r.target in used_t:
                continue
            for c in asked_of.get(r.condition, []):
                if c in used_c:
                    continue
                p = gaptables.pair(T, S, r.target, c)
                if keep(p):
                    used_t.add(r.target), used_c.add(c)
                    return r.target, c, p
        raise SystemExit("REFUSED: no pair of the kind wanted")

    c = disease(lambda r: 15 <= r.genes_tied <= 300 and 5 <= r.no_compound < r.genes_tied)
    add("tied_no_compound", f"How many genes does the graph tie to {name[c]}? Of those, how many have no drug acting on them and no compound at their family's potency cutoff?",
        "two numbers: genes tied, and genes with no drug and no compound", {"tied": int(summ.loc[c].genes_tied), "no_compound": int(summ.loc[c].no_compound)})
    c = disease(lambda r: r.loci_not_resolved >= 3)
    add("loci", f"How many loci of {name[c]} carry its GWAS link on more than one gene?", "a number", {"loci": int(summ.loc[c].loci_not_resolved)})
    c = disease(lambda r: pd.notna(r["cases_needed_for_odds_ratio_1.1"]))
    need = int(summ.loc[c]["cases_needed_for_odds_ratio_1.1"])
    add("study_size", f"About how many cases would a study of {name[c]} need to find a common variant with an odds ratio of 1.1? Did the largest study so far have that many?",
        "a number of cases, and yes or no", {"cases": need, "had": bool(int(summ.loc[c].cases) >= need)})

    ok_links = links[links.target.isin(clean_t) & links.trait.isin(clean_d)]
    r = ok_links[ok_links.n_studies >= 2].sample(1, random_state=rng.randrange(2 ** 31)).iloc[0]
    used_t.add(r.target), used_c.add(r.trait)
    add("first_year", f"Is {sym[r.target]} linked to {name[r.trait]} by GWAS? In what year was the link first reported, and by how many studies in all?",
        "yes or no, a year, a number of studies", {"linked": True, "year": int(r.first_year), "studies": int(r.n_studies)})
    later = ok_links[(ok_links.first_year >= 2012) & (ok_links.first_year <= 2020) & ~ok_links.target.isin(used_t) & ~ok_links.trait.isin(used_c)]
    r = later.sample(1, random_state=rng.randrange(2 ** 31)).iloc[0]
    used_t.add(r.target), used_c.add(r.trait)
    add("as_of", f"As of {int(r.first_year) - 1}, was {sym[r.target]} linked to {name[r.trait]} by GWAS? If not, has a link been reported since, and in what year?",
        "yes or no, and a year or none", {"linked_then": False, "year": int(r.first_year)})

    A, X, DT, DI = T["assertions"], T["expression"], T["drug_targets"], T["drug_indications"]
    many = A[~A.withheld].groupby(["target", "condition"]).source.nunique().reset_index()
    tg_, c, p = pair(many[many.source >= 2], lambda p: p["sources"] >= 2)
    add("sources", f"How many disease sources assert {sym[tg_]} with {name[c]}?", "a number", {"sources": p["sources"]})
    lit = A[A.withheld & A.rule.str.contains("literature_lift")]
    tg_, c, p = pair(lit, lambda p: p["withheld_rows"] == 1 and len(p["papers"]) == 1)
    src = A[(A.target == tg_) & A.withheld & A.condition.isin(set(S[S.asked_id == c].read_id) or {c})].source.iloc[0]
    add("withheld", f"The graph withholds a row for {sym[tg_]} and {name[c]}. Which source's row is it, and how many papers mention both?",
        "a source, and a number of papers", {"source": src, "papers": p["papers"][0]})
    path = DT[DT.mechanism].merge(DI, on="drug")
    tg_, c, p = pair(path[["target", "condition"]], lambda p: p["drugs"] >= 1 and p["mechanism"])
    ids = set(S[S.asked_id == c].read_id) or {c}
    drugs = sorted({name.get(g, g) for g in path[(path.target == tg_) & path.condition.isin(ids)].drug})
    add("drug", f"Name one drug indicated for {name[c]} whose mechanism of action is on {sym[tg_]}.", "a drug name", {"drug": drugs})

    tgap = pd.read_parquet(bundle / "target_gaps.parquet")
    dark = tgap[tgap.flips.isin([1, 2]) & tgap.id.isin(clean_t - used_t)].sort_values("id")
    r = dark.iloc[rng.randrange(len(dark))]
    used_t.add(r.id)
    add("dark", f"Is {sym[r.id]} Tdark? If it is, how many criteria must it still clear to leave Tdark?", "yes or no, and a number", {"tdark": True, "criteria": int(r.flips)})
    if G.resolve(FAKE)["ok"] or G.resolve(FAKE)["problem"] != "not in the graph":
        raise SystemExit(f"REFUSED: {FAKE!r} is in the graph, or looks like something in it")
    add("not_in_graph", f"What does the graph hold on the gene {FAKE}?", "what the graph says", {"problem": "not_in_graph"})

    c = disease(lambda r: pd.notna(r.smallest_odds_ratio_it_could_find))
    linked_to_c = set(links[links.trait == c].target)
    tg_ = rng.choice(sorted(clean_t - used_t - linked_to_c))
    used_t.add(tg_)
    add("absence", f"The graph has no GWAS link between {sym[tg_]} and {name[c]}. Does that mean {sym[tg_]} has no association with {name[c]}? "
        f"For the largest study behind {name[c]}'s links, what is the smallest odds ratio a variant carried by 20% of chromosomes needed in order to be found?",
        "yes or no, and an odds ratio", {"no_association": False, "odds_ratio": float(summ.loc[c].smallest_odds_ratio_it_could_find)})
    tg_, c, p = pair(X, lambda p: p["observations"] == 1 and p["sources"] == 0 and p["withheld_rows"] == 0 and p["drugs"] == 0)
    add("observation", f"The graph holds an expression observation for {sym[tg_]} in {name[c]}. Was its expression higher or lower? Does the graph say {sym[tg_]} is involved in {name[c]}?",
        "higher or lower, and yes or no", {"direction": p["directions"][0], "involved": False})
    if len(out) != 12:
        raise SystemExit(f"REFUSED: {len(out)} tasks, 12 are wanted")
    return out


def via_ask(task: dict, G) -> dict:
    """What the question layer gives for one task, for the check against the expected answer."""
    k, x = task["kind"], task["text"]
    if k in ("tied_no_compound", "loci", "study_size"):
        c = re.match(r"^(?:How many genes does the graph tie to|How many loci of|About how many cases would a study of) (.+?)(?:\?| carry| need)", x).group(1)
        o = G.gaps(disease=c)
        if k == "tied_no_compound":
            return {"tied": o["holds"]["tied_genes"], "no_compound": o["counts"]["no_compound"]}
        if k == "loci":
            return {"loci": o["counts"]["locus_not_resolved"]}
        nd = next(g for g in o["gaps"] if g["kind"] == "no_genetic_support")["what_an_absent_link_could_have_seen"]["the_study_that_would_close_it"]
        c11 = next(v for v in nd["common_variant"] if v["odds_ratio"] == 1.1)
        return {"cases": c11["cases"], "had": c11["the_largest_study_had_them"]}
    if k == "first_year":
        a, b = re.match(r"^Is (\S+) linked to (.+?) by GWAS\?", x).groups()
        w = G.why(a, b, limit=100_000)
        return {"linked": w["found"], "year": min(s_["year"] for s_ in w["studies"]) if w["found"] else None,
                "studies": len({s_["study"] for s_ in w["studies"]}) if w["found"] else None}
    if k == "as_of":
        y, a, b = re.match(r"^As of (\d+), was (\S+) linked to (.+?) by GWAS\?", x).groups()
        w = G.why(a, b, as_of=int(y))
        return {"linked_then": w["found"], "year": w.get("reported_later")}
    if k in ("sources", "withheld", "drug", "observation"):
        a, b = (re.match(r"^How many disease sources assert (\S+) with (.+)\?$", x) or re.match(r"^The graph withholds a row for (\S+) and (.+?)\. Which", x)
                or re.match(r"^Name one drug indicated for (.+?) whose mechanism of action is on (\S+)\.$", x) or re.match(r"^The graph holds an expression observation for (\S+) in (.+?)\. Was", x)).groups()
        if k == "drug":
            a, b = b, a
        o = G.why(a, b)["other_sources"]
        if k == "sources":
            return {"sources": o["sources_asserting"]}
        if k == "withheld":
            w = o["withheld"][0]
            return {"source": w["source"], "papers": w.get("papers_mentioning_both")}
        if k == "drug":
            return {"drug": sorted({p["drug"] for p in o["through_a_drug"] if p["target_is_the_drugs_mechanism"]})}
        return {"direction": o["expression_observations"][0]["direction"], "involved": False}
    if k == "dark":
        a = re.match(r"^Is (\S+) Tdark\?", x).group(1)
        g = G.gaps(a)
        d = next((v for v in g["gaps"] if v["kind"] == "dark"), None)
        return {"tdark": g["target"]["pharos_level"] == "Tdark", "criteria": d["criteria_to_clear"] if d else None}
    if k == "not_in_graph":
        k_ = G.known(FAKE)
        return {"problem": "not_in_graph" if k_.get("problem") == "not in the graph" else "found"}
    if k == "absence":
        a, b = re.match(r"^The graph has no GWAS link between (\S+) and (.+?)\. Does", x).groups()
        w = G.why(a, b)
        return {"no_association": False if not w["found"] else None, "odds_ratio": w["what_an_absent_link_could_have_seen"]["largest"]["smallest_odds_ratio_it_could_find"]}
    raise SystemExit(f"unknown kind {k}")


def check(bundle: Path, tasks: list[dict]) -> list[str]:
    from pharos_gap_kg import ask
    G = ask.Graph(bundle)
    bad = []
    for t in tasks:
        got = via_ask(t, G)
        want = t["gold"]
        if t["kind"] == "drug":                         # the answer lists the first drugs; any one of the set is right
            ok = bool(got["drug"]) and set(got["drug"]) <= set(want["drug"])
        else:
            ok = all(got.get(k) == v if not isinstance(v, float) else abs(got.get(k) - v) < 1e-9 for k, v in want.items())
        if not ok:
            bad.append(f"{t['id']} {t['kind']}: expected {want}, the question layer gives {got}")
    return bad


# ---- grading one session ------------------------------------------------------------------------------------

def _yes(v) -> bool | None:
    if isinstance(v, bool):
        return v
    s = str(v).strip().casefold()
    return True if s in ("yes", "y", "true") else False if s in ("no", "n", "false") else None


def _num(v) -> float | None:
    try:
        return float(str(v).replace(",", "").strip())
    except (TypeError, ValueError):
        return None


def right(task: dict, answer: dict | None) -> tuple[bool, bool]:
    """(right, invented) for one recorded answer."""
    g, k = task["gold"], task["kind"]
    if not answer:
        return False, False
    if k == "not_in_graph":
        said = re.sub(r"[_\s]+", " ", str(answer.get("problem", answer.get("answer", ""))).casefold()).strip()
        ok = bool(re.search(r"not (in|found in) (the )?graph|^nothing|no such", said))
        return ok, bool(said) and not ok
    if k == "drug":
        return str(answer.get("drug", "")).strip().casefold() in {x.casefold() for x in g["drug"]}, False
    for f, v in g.items():
        a = answer.get(f)
        if isinstance(v, bool):
            if _yes(a) is not v:
                return False, False
        elif f == "cases":
            n = _num(a)
            if n is None or abs(n - v) > SIZE_WITHIN * v:
                return False, False
        elif f == "odds_ratio":
            n = _num(a)
            if n is None or abs(n - v) > OR_WITHIN:
                return False, False
        elif isinstance(v, int):
            if _num(a) != v:
                return False, False
        elif str(a).strip().casefold() != str(v).casefold():
            return False, False
    return True, False


def sus(items: list[int]) -> float:
    if len(items) != 10 or not all(1 <= int(i) <= 5 for i in items):
        raise SystemExit("REFUSED: the System Usability Scale has ten items, each scored 1 to 5")
    return 2.5 * sum((int(v) - 1) if i % 2 == 0 else (5 - int(v)) for i, v in enumerate(items))


def grade(session: dict, tasks: list[dict], frozen_digest: str) -> dict:
    if not re.fullmatch(r"P\d+", str(session.get("participant", ""))):
        raise SystemExit("REFUSED: a participant is recorded by a code such as P1, and by nothing that names them")
    rows, times = [], []
    for t in tasks:
        a = session["tasks"].get(t["id"], {})
        m = _num(a.get("minutes"))
        ok, invented = right(t, a.get("answer")) if m is not None and m <= LIMIT_MINUTES else (False, False)
        rows.append({"id": t["id"], "kind": t["kind"], "right": ok, "invented": invented, "minutes": m})
        if m is not None:
            times.append(min(m, LIMIT_MINUTES))
    n_right = sum(r["right"] for r in rows)
    med = statistics.median(times) if times else None
    score = sus(session["sus"]) if session.get("sus") else None
    usable = n_right >= RIGHT_AT_LEAST and not any(r["invented"] for r in rows) and med is not None and med <= MEDIAN_MINUTES
    return {"participant": session["participant"], "field": session.get("field"), "years_in_research": session.get("years_in_research"),
            "client": session.get("client"), "model": session.get("model"), "date": session.get("date"), "tasks_sha256": frozen_digest,
            "right": n_right, "of": len(tasks), "invented": sum(r["invented"] for r in rows), "median_minutes": med,
            "usable_by_this_researcher": bool(usable), "sus": score, "sus_against_average": None if score is None else score - SUS_AVERAGE,
            "tasks": rows, "open_answers": session.get("open", [])}


def sheet(tasks: list[dict]) -> str:
    """The participant's sheet: the tasks and what to answer, with no expected answer in it."""
    lines = ["# Tasks", "", f"Take them in order. Each has at most {LIMIT_MINUTES:g} minutes; say your final answer aloud when you have it.", ""]
    for t in tasks:
        lines += [f"**{t['id']}.** {t['text']}", "", f"*Answer:* {t['record']}", ""]
    return "\n".join(lines)


def template(tasks: list[dict]) -> dict:
    """The file a facilitator fills in during a session: one answer field per expected field, and nothing that names the participant."""
    return {"participant": "P1", "field": "", "years_in_research": "", "client": "", "model": "", "date": "",
            "tasks": {t["id"]: {"answer": {k: None for k in t["gold"]}, "minutes": None} for t in tasks},
            "sus": [None] * 10, "open": ["", "", ""]}


def _selftest() -> int:
    ok = True

    def ck(label: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {label}")

    tasks = [{"id": "T1", "kind": "study_size", "gold": {"cases": 14316, "had": True}}, {"id": "T2", "kind": "drug", "gold": {"drug": ["donepezil", "galantamine"]}},
             {"id": "T3", "kind": "not_in_graph", "gold": {"problem": "not_in_graph"}}, {"id": "T4", "kind": "absence", "gold": {"no_association": False, "odds_ratio": 1.084}},
             {"id": "T5", "kind": "withheld", "gold": {"source": "CTD", "papers": 0}}]
    ck("a study size said to two figures is right within 5%, and a yes is read as yes", right(tasks[0], {"cases": "14,000", "had": "yes"}) == (True, False)
       and right(tasks[0], {"cases": 13000, "had": "yes"})[0] is False)
    ck("any qualifying drug is right, in any case", right(tasks[1], {"drug": "Galantamine"}) == (True, False) and right(tasks[1], {"drug": "aspirin"})[0] is False)
    ck("a name not in the graph: saying so is right, anything else is invented", right(tasks[2], {"problem": "not in the graph"}) == (True, False)
       and right(tasks[2], {"problem": "a kinase on chromosome 3"}) == (False, True))
    ck("an odds ratio is right within 0.01, and a no for no association", right(tasks[3], {"no_association": "no", "odds_ratio": "1.08"}) == (True, False)
       and right(tasks[3], {"no_association": "yes", "odds_ratio": 1.084})[0] is False)
    ck("a count of zero is a count, not a missing answer", right(tasks[4], {"source": "ctd", "papers": 0}) == (True, False) and right(tasks[4], {"source": "CTD"})[0] is False)
    ck("SUS: all agree on the positive items and disagree on the negative ones is 100; the reverse is 0", sus([5, 1] * 5) == 100 and sus([1, 5] * 5) == 0 and sus([3] * 10) == 50)
    s = {"participant": "P1", "sus": [4, 2] * 5, "tasks": {"T1": {"answer": {"cases": 14000, "had": "yes"}, "minutes": 3}, "T2": {"answer": {"drug": "donepezil"}, "minutes": 9},
                                                           "T3": {"answer": {"problem": "not in the graph"}, "minutes": 1}}}
    g = grade(s, tasks, "x")
    ck("a task past the time limit is wrong whatever was said, and an unanswered one is wrong", [r["right"] for r in g["tasks"]] == [True, False, True, False, False]
       and g["right"] == 2 and g["median_minutes"] == 3 and g["sus"] == 75)
    try:
        grade({**s, "participant": "Jane"}, tasks, "x")
        ck("a participant is recorded by a code, never a name", False)
    except SystemExit:
        ck("a participant is recorded by a code, never a name", True)
    ck("the sheet holds the tasks and no expected answer", "14316" not in sheet([{**t, "text": "q", "record": "r"} for t in tasks]))
    ck("a name in the docs is found as a whole word only", _mentioned("APOE and Parkinson disease")("APOE") and not _mentioned("APOE4")("APOE"))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) >= 2 and argv[0] in ("freeze", "check"):
        b = paths.bundle(argv[1])
        if argv[0] == "freeze":
            if tasks_path().is_file():
                raise SystemExit(f"REFUSED: {tasks_path().name} is frozen and is not written again")
            ts = make_tasks(b)
            bad = check(b, ts)
            if bad:
                raise SystemExit("REFUSED: the question layer disagrees with the expected answers:\n  " + "\n  ".join(bad))
            tasks_path().write_text(json.dumps({"seed": SEED, "sha256": digest(ts), "tasks": ts}, indent=1) + "\n")
            print(f"  {len(ts)} tasks frozen -> {tasks_path()}  sha256 {digest(ts)}")
            return 0
        f = json.loads(tasks_path().read_text())
        bad = check(b, f["tasks"])
        print("  the question layer gives every expected answer" if not bad else "\n".join(bad))
        return 0 if not bad else 1
    if argv and argv[0] == "template":
        ts = json.loads(tasks_path().read_text())["tasks"]
        print(json.dumps(template(ts), indent=1))
        return 0
    if argv and argv[0] == "sheet":
        print(sheet(json.loads(tasks_path().read_text())["tasks"]))
        return 0
    if len(argv) == 2 and argv[0] == "grade":
        f = json.loads(tasks_path().read_text())
        if digest(f["tasks"]) != f["sha256"]:
            raise SystemExit("REFUSED: the frozen tasks do not match their digest")
        g = grade(json.loads(Path(argv[1]).read_text()), f["tasks"], f["sha256"])
        out = paths.DOCS_DATA / f"person_test_{g['participant']}.json"
        if out.is_file():
            raise SystemExit(f"REFUSED: {out.name} exists; a session is graded once")
        out.write_text(json.dumps(g, indent=1) + "\n")
        print(f"  {g['participant']}: {g['right']} of {g['of']} right, {g['invented']} invented, median {g['median_minutes']} min; usable: {g['usable_by_this_researcher']}; SUS {g['sus']}")
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
