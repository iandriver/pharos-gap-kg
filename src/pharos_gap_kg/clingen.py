"""Kept and withheld rows against ClinGen's expert curations of gene-disease validity.

The pruning rule withholds a disease source's row when the literature does not mention the gene and the
disease together more often than chance (and Monarch's when no other source corroborates it). Neither
the rule nor the sources read ClinGen. This puts both sides of the rule beside ClinGen's curations
(docs/PROTOCOL.md, "Kept and withheld rows against ClinGen", defined on 26 Oct 2026 before any row was
matched to a classification):

  in scope    a row whose target's symbol is a gene ClinGen has curated and whose condition is a MONDO id
  judged      an in-scope row whose disease is one ClinGen curated for that gene: the same MONDO id, or, as
              the secondary reading, the same live term in EFO or one within two is_a steps above or below
  supported   judged, and a curation of the pair is Definitive, Strong or Moderate
  limited     judged, not supported, and a curation is Limited
  against     judged, and every curation of the pair is Disputed, Refuted or No Known Disease Relationship

ClinGen curates a few thousand mostly Mendelian genes, so this measures agreement with expert curation
where it exists, not the precision of the whole graph, and a pair ClinGen has not curated is not
evidence against the row.

Run:  python -m pharos_gap_kg.clingen fetch clingen_gv_20261026 [--from FILE]   # the pinned file, checked
      python -m pharos_gap_kg.clingen scope pharos319     # rows in scope by source; reads no classification
      python -m pharos_gap_kg.clingen audit pharos319     # writes docs/data/clingen_audit_<release>.json
      python -m pharos_gap_kg.clingen --selftest
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import sys
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import ontology, paths, registry

KEY = "clingen_gv_20261026"
SUPPORTED = {"Definitive", "Strong", "Moderate"}
LIMITED = {"Limited"}
AGAINST = {"Disputed", "Refuted", "No Known Disease Relationship"}
RANK = {"supported": 2, "limited": 1, "against": 0}          # several curations match a row: the best class counts
PRUNED = ["CTD", "JensenLab Knowledge MedlinePlus", "JensenLab Knowledge UniProtKB-KW", "Monarch", "UniProt Disease"]
STEPS = 2                 # the secondary reading: two is_a steps above or below, never across to a sibling
MIN_JUDGED = 10           # a pruned source gets a verdict of its own only with this many judged withheld rows
SEED, DRAWS = 20261026, 2000
UA = {"User-Agent": "pharos-gap-kg (research; github.com/iandriver)"}


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def fetch(s: registry.Source, src: Path | None = None) -> Path:
    """The pinned file, from ClinGen or from a copy already downloaded; refused unless its size and digest are the registry's."""
    side = s.path.with_name(s.filename + ".lineage.json")
    if s.path.is_file() and side.is_file() and sha(s.path.read_bytes()) == s.sha256:
        print(f"{s.filename}: already fetched ({s.sha256[:12]})")
        return s.path
    modified = ""
    if src is None:
        with urllib.request.urlopen(urllib.request.Request(s.url, headers=UA), timeout=120) as r:
            data, modified = r.read(), r.headers.get("Last-Modified", "")
    else:
        data = Path(src).read_bytes()
    if len(data) != s.bytes or sha(data) != s.sha256:
        raise SystemExit(f"REFUSED {s.filename}: {len(data):,} bytes, sha256 {sha(data)[:16]}; the registry pins {s.bytes:,} bytes, "
                         f"{s.sha256[:16]}. ClinGen remakes this file as curations are added; a newer one is a new registry entry.")
    s.path.parent.mkdir(parents=True, exist_ok=True)
    tmp = s.path.with_name(s.filename + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, s.path)
    side.write_text(json.dumps({"key": s.key, "url": s.url, "bytes": len(data), "sha256": s.sha256, "server_last_modified": modified,
                                "published": s.published, "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                                "from": str(src) if src else s.url}, indent=1) + "\n")
    print(f"{s.filename}: {len(data):,} bytes, sha256 {s.sha256}")
    return s.path


def read(text: str) -> pd.DataFrame:
    """ClinGen's summary CSV: a few lines about the file, the header, rows of '+' between, then one row per curation."""
    rows = list(csv.reader(io.StringIO(text)))
    i = next(k for k, r in enumerate(rows) if r and r[0] == "GENE SYMBOL")
    body = [r for r in rows[i + 1:] if r and not r[0].startswith("+")]
    df = pd.DataFrame(body, columns=rows[i])
    out = pd.DataFrame({"symbol": df["GENE SYMBOL"], "mondo": df["DISEASE ID (MONDO)"], "label": df["DISEASE LABEL"],
                        "classification": df["CLASSIFICATION"], "moi": df["MOI"], "date": df["CLASSIFICATION DATE"].str[:10],
                        "report": df["ONLINE REPORT"]})
    unknown = set(out.classification) - SUPPORTED - LIMITED - AGAINST
    if unknown:
        raise SystemExit(f"REFUSED: classifications the protocol does not define: {sorted(unknown)}")
    return out


def created(text: str) -> str:
    line = next((l for l in text.splitlines()[:5] if "FILE CREATED" in l), "")
    return line.split("FILE CREATED:")[-1].split('"')[0].strip()


def pair_classes(cg: pd.DataFrame) -> pd.DataFrame:
    def cls(s) -> str:
        s = set(s)
        return "supported" if s & SUPPORTED else "limited" if s & LIMITED else "against"
    return cg.groupby(["symbol", "mondo"]).agg(cls=("classification", cls), label=("label", "first"), report=("report", "first")).reset_index()


def near(term: str, terms: dict, children: dict, steps: int = STEPS) -> set[str]:
    """The live term and every term within `steps` is_a steps above it or below it."""
    c = ontology.canonical(term, terms)
    if c is None:
        return set()
    out = {c}
    for nxt in (lambda t: terms[t]["parents"] if t in terms else (), lambda t: children.get(t, ())):
        frontier = {c}
        for _ in range(steps):
            frontier = {n for t in frontier for n in nxt(t)} - out
            out |= frontier
    return out


def children_of(terms: dict) -> dict[str, set[str]]:
    ch: dict[str, set[str]] = defaultdict(set)
    for t, d in terms.items():
        for p in d["parents"]:
            ch[p].add(t)
    return ch


def judge(rows: pd.DataFrame, pairs: pd.DataFrame, terms: dict | None) -> pd.DataFrame:
    """For each row: in scope or not, its class by the same MONDO id, and its class by the secondary reading."""
    genes = set(pairs.symbol)
    rows = rows.assign(in_scope=rows.symbol.isin(genes) & rows.condition.str.startswith("MONDO:"))
    exact = {(s, m): c for s, m, c in zip(pairs.symbol, pairs.mondo, pairs.cls)}
    rows["exact"] = [exact.get((s, c)) if i else None for s, c, i in zip(rows.symbol, rows.condition, rows.in_scope)]
    if terms is None:
        rows["near"] = rows["exact"]
        return rows
    ch, cache = children_of(terms), {}
    by_gene: dict[str, list] = defaultdict(list)
    for s, m, c in zip(pairs.symbol, pairs.mondo, pairs.cls):
        by_gene[s].append((m, c))

    def best(s: str, cond: str, ex: str | None) -> str | None:
        cc, found = ontology.canonical(cond, terms), [ex] if isinstance(ex, str) else []
        if cc is not None:
            for m, c in by_gene.get(s, ()):
                if m not in cache:
                    cache[m] = near(m, terms, ch)
                if cc in cache[m]:
                    found.append(c)
        return max(found, key=RANK.get) if found else None
    rows["near"] = [best(s, c, e) if i else None for s, c, e, i in zip(rows.symbol, rows.condition, rows.exact, rows.in_scope)]
    return rows


def wilson(k: int, n: int, z: float = 1.959964) -> list:
    if n == 0:
        return [None, None, None]
    p = k / n
    d = 1 + z * z / n
    mid, half = (p + z * z / (2 * n)) / d, z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(p, 4), round(mid - half, 4), round(mid + half, 4)]


def newcombe(k1: int, n1: int, k2: int, n2: int) -> list:
    """p1 - p2 with Newcombe's hybrid score interval (method 10)."""
    if not n1 or not n2:
        return [None, None, None]
    (p1, l1, u1), (p2, l2, u2) = [[k / n] + wilson(k, n)[1:] for k, n in ((k1, n1), (k2, n2))]
    d = p1 - p2
    return [round(d, 4), round(d - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), 4), round(d + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2), 4)]


def counts(c: pd.Series) -> dict:
    """Counts and rates for one group of rows, given each row's class (None for a row ClinGen did not judge)."""
    n_scope = int(len(c))
    k = {x: int((c == x).sum()) for x in RANK}
    judged = sum(k.values())
    return {"in_scope": n_scope, "judged": judged, **k, "agreement": wilson(k["supported"], judged),
            "against_rate": wilson(k["against"], judged), "confirmation": wilson(k["supported"], n_scope)}


def pooled(cells: dict[str, tuple], rng: np.random.Generator) -> list:
    """sum_s w_s (p_kept,s - p_withheld,s) / sum_s w_s, w_s the withheld denominator; 95% interval from a parametric bootstrap of each cell."""
    use = {s: c for s, c in cells.items() if c[1] and c[3]}
    if not use:
        return [None, None, None]
    w = np.array([c[3] for c in use.values()], float)
    point = sum(wi * (c[0] / c[1] - c[2] / c[3]) for wi, c in zip(w, use.values())) / w.sum()
    draws = np.zeros(DRAWS)
    for wi, (k1, n1, k2, n2) in zip(w, use.values()):
        draws += wi * (rng.binomial(n1, k1 / n1, DRAWS) / n1 - rng.binomial(n2, k2 / n2, DRAWS) / n2)
    draws /= w.sum()
    return [round(point, 4), round(float(np.percentile(draws, 2.5)), 4), round(float(np.percentile(draws, 97.5)), 4)]


def compare(J: pd.DataFrame, col: str, rng: np.random.Generator) -> dict:
    """Kept against withheld for each pruned source, and pooled over them, with and without Monarch."""
    out, cells_a, cells_c = {}, {}, {}
    for s in PRUNED:
        g = J[(J.source == s) & J.in_scope]
        kc, wc = counts(g[~g.withheld][col]), counts(g[g.withheld][col])
        cells_a[s] = (kc["supported"], kc["judged"], wc["supported"], wc["judged"])
        cells_c[s] = (kc["supported"], kc["in_scope"], wc["supported"], wc["in_scope"])
        out[s] = {"agreement_kept_minus_withheld": newcombe(*cells_a[s]), "confirmation_kept_minus_withheld": newcombe(*cells_c[s]),
                  "judged_withheld": wc["judged"], "verdict_of_its_own": wc["judged"] >= MIN_JUDGED}
    no_m = lambda d: {s: c for s, c in d.items() if s != "Monarch"}
    out["pooled"] = {"agreement_kept_minus_withheld": pooled(cells_a, rng), "confirmation_kept_minus_withheld": pooled(cells_c, rng),
                     "agreement_without_monarch": pooled(no_m(cells_a), rng), "confirmation_without_monarch": pooled(no_m(cells_c), rng)}
    return out


def rows_of(bundle: Path) -> pd.DataFrame:
    A = pd.read_parquet(bundle / "assertions.parquet")
    N = pd.read_parquet(bundle / "nodes.parquet", columns=["id", "name"])
    return A.assign(symbol=A.target.map(dict(zip(N.id, N.name))))


def examples(J: pd.DataFrame, pairs: pd.DataFrame, mask, n: int = 20) -> list:
    lab = {(s, m): (l, r) for s, m, l, r in zip(pairs.symbol, pairs.mondo, pairs.label, pairs.report)}
    out = []
    for r in J[mask].sort_values(["source", "symbol", "condition"]).head(n).itertuples(index=False):
        l, u = lab.get((r.symbol, r.condition), ("", ""))
        out.append({"source": r.source, "symbol": r.symbol, "condition": r.condition, "clingen_label": l, "clingen_report": u,
                    "rule": r.rule if r.withheld else ""})
    return out


def audit(release: str) -> dict:
    text = registry.require(KEY).read_text(encoding="utf-8")
    cg = read(text)
    pairs = pair_classes(cg)
    terms = ontology.read(registry.require("efo_obo"))
    J = judge(rows_of(paths.bundle(release)), pairs, terms)
    rng = np.random.default_rng(SEED)
    by_source = []
    for (s, w), g in J.groupby(["source", "withheld"]):
        sc = g[g.in_scope]
        by_source.append({"source": s, "side": "withheld" if w else "kept", "rows": int(len(g)),
                          "exact": counts(sc.exact), "near": counts(sc.near)})
    K = J[~J.withheld & J.in_scope]
    return {
        "release": release, "protocol": "docs/PROTOCOL.md, Kept and withheld rows against ClinGen (26 Oct 2026)",
        "clingen": {"key": KEY, "file": registry.SOURCES[KEY].filename, "sha256": registry.SOURCES[KEY].sha256, "file_created": created(text),
                    "curations": int(len(cg)), "genes": int(cg.symbol.nunique()), "diseases": int(cg.mondo.nunique()), "pairs": int(len(pairs)),
                    "pairs_by_class": {c: int((pairs.cls == c).sum()) for c in RANK}},
        "ontology_for_the_secondary_reading": registry.SOURCES["efo_obo"].filename, "steps": STEPS,
        "rows": int(len(J)), "in_scope": int(J.in_scope.sum()), "in_scope_withheld": int((J.in_scope & J.withheld).sum()),
        "clingen_genes_among_targets": int(len(set(pairs.symbol) & set(J.symbol.dropna()))),
        "kept_rows_of_all_sources": {"exact": counts(K.exact), "near": counts(K.near)},
        "withheld_rows_of_all_sources": {"exact": counts(J[J.withheld & J.in_scope].exact), "near": counts(J[J.withheld & J.in_scope].near)},
        "by_source": by_source,
        "pruned_sources": {"exact": compare(J, "exact", rng), "near": compare(J, "near", rng), "min_judged_withheld_for_a_verdict": MIN_JUDGED,
                           "seed": SEED, "draws": DRAWS},
        "withheld_rows_clingen_supports": {"by_source": {s: int(((J.source == s) & J.withheld & (J.exact == "supported")).sum()) for s in PRUNED},
                                           "examples": examples(J, pairs, J.withheld & (J.exact == "supported"))},
        "kept_rows_clingen_is_against": {"by_source": {s: int(((J.source == s) & ~J.withheld & (J.exact == "against")).sum()) for s in sorted(J.source.unique())},
                                         "examples": examples(J, pairs, ~J.withheld & (J.exact == "against"))},
    }


def scope(release: str) -> pd.DataFrame:
    """Rows in scope by source and side. Reads the curated genes and nothing about their classifications."""
    genes = set(read(registry.require(KEY).read_text(encoding="utf-8")).symbol)
    R = rows_of(paths.bundle(release))
    R["in_scope"] = R.symbol.isin(genes) & R.condition.str.startswith("MONDO:")
    return R.groupby(["source", "withheld"]).in_scope.agg(["size", "sum"])


SAMPLE = '''"CLINGEN GENE DISEASE VALIDITY CURATIONS","","","","","","","","",""
"FILE CREATED: 2026-10-26","","","","","","","","",""
"WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity","","","","","","","","",""
"+++++","+++","+++","+++","+++","+++","+++","+++","+++","+++"
"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"
"+++++","+++","+++","+++","+++","+++","+++","+++","+++","+++"
"AAA","HGNC:1","disease one","MONDO:0000001","AD","SOP10","Definitive","u1","2024-03-14T16:00:00.000Z","p"
"AAA","HGNC:1","disease one","MONDO:0000001","AR","SOP10","Disputed","u2","2024-03-14T16:00:00.000Z","p"
"BBB","HGNC:2","disease two","MONDO:0000002","AD","SOP10","Disputed","u3","2024-03-14T16:00:00.000Z","p"
"BBB","HGNC:2","disease two","MONDO:0000002","AD","SOP10","Refuted","u4","2024-03-14T16:00:00.000Z","p"
"CCC","HGNC:3","disease three","MONDO:0000003","AD","SOP10","Limited","u5","2024-03-14T16:00:00.000Z","p"
"CCC","HGNC:3","disease three","MONDO:0000003","AD","SOP10","No Known Disease Relationship","u6","2024-03-14T16:00:00.000Z","p"
"DDD","HGNC:4","leaf","MONDO:0000013","AD","SOP10","Strong","u7","2024-03-14T16:00:00.000Z","p"
'''


def _selftest() -> int:
    import tempfile
    ok = True

    def ck(label: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {label}")

    cg = read(SAMPLE)
    ck("the header lines and the '+' rows are skipped", len(cg) == 7 and created(SAMPLE) == "2026-10-26" and cg.date.iloc[0] == "2024-03-14")
    P = pair_classes(cg).set_index(["symbol", "mondo"]).cls
    ck("a pair is supported by any strong curation, against only if every one is", P[("AAA", "MONDO:0000001")] == "supported"
       and P[("BBB", "MONDO:0000002")] == "against" and P[("CCC", "MONDO:0000003")] == "limited")
    try:
        read(SAMPLE.replace('"Strong"', '"Probable"'))
        ck("a classification the protocol does not define is refused", False)
    except SystemExit as e:
        ck("a classification the protocol does not define is refused", "Probable" in str(e))
    # MONDO:10 > 11 > 12 > 13 (13 is ClinGen's), MONDO:14 a sibling of 13 under 12, MONDO:15 below 13, MONDO:16 below 15, MONDO:9 obsolete for 12
    T = lambda parents, obs=False, rb="": {"parents": set(parents), "obsolete": obs, "replaced_by": rb}
    terms = {"MONDO:0000010": T([]), "MONDO:0000011": T(["MONDO:0000010"]), "MONDO:0000012": T(["MONDO:0000011"]),
             "MONDO:0000013": T(["MONDO:0000012"]), "MONDO:0000014": T(["MONDO:0000012"]), "MONDO:0000015": T(["MONDO:0000013"]),
             "MONDO:0000016": T(["MONDO:0000015"]), "MONDO:0000017": T(["MONDO:0000016"]), "MONDO:0000009": T([], True, "MONDO:0000012")}
    ck("two steps above and below, no siblings, nothing three steps away",
       near("MONDO:0000013", terms, children_of(terms)) == {f"MONDO:00000{i}" for i in (11, 12, 13, 15, 16)})
    rows = pd.DataFrame([
        ("AAA", "MONDO:0000001", "S", False), ("BBB", "MONDO:0000002", "S", True), ("CCC", "MONDO:0000003", "S", False),
        ("DDD", "MONDO:0000013", "S", False), ("DDD", "MONDO:0000011", "S", False), ("DDD", "MONDO:0000009", "S", True),
        ("DDD", "MONDO:0000014", "S", False), ("DDD", "MONDO:0000010", "S", False), ("DDD", "MESH:D000001", "S", False),
        ("EEE", "MONDO:0000001", "S", False)], columns=["symbol", "condition", "source", "withheld"])
    J = judge(rows, pair_classes(cg), terms)
    ck("in scope: a curated gene and a MONDO condition", J.in_scope.tolist() == [True] * 8 + [False, False])
    cls = lambda col: [x if isinstance(x, str) else None for x in J[col]]
    ck("the same MONDO id, read exactly", cls("exact")[:5] == ["supported", "against", "limited", "supported", None])
    ck("the secondary reading reaches two steps up and an obsolete id's replacement, not a sibling or three steps",
       cls("near")[3:8] == ["supported", "supported", "supported", None, None])
    c = counts(J[J.in_scope].exact)
    ck("counts and rates over in-scope rows", (c["in_scope"], c["judged"], c["supported"], c["limited"], c["against"]) == (8, 4, 2, 1, 1)
       and c["agreement"][0] == 0.5 and c["confirmation"][0] == 0.25)
    ck("Wilson's interval, 0 of 10", wilson(0, 10) == [0.0, 0.0, 0.2775])
    ck("Newcombe's hybrid interval, his example 56/70 against 48/80", newcombe(56, 70, 48, 80) == [0.2, 0.0524, 0.3339])
    a = pooled({"x": (80, 100, 40, 50), "y": (30, 60, 5, 10)}, np.random.default_rng(1))
    b = pooled({"x": (80, 100, 40, 50), "y": (30, 60, 5, 10)}, np.random.default_rng(1))
    ck("the pooled difference is weighted by the withheld rows, and the same seed gives the same interval",
       a == b and a[0] == round((50 * 0.0 + 10 * 0.0) / 60, 4) and a[1] < 0 < a[2])
    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        old = paths.RAW
        paths.RAW = t / "raw"
        try:
            data = SAMPLE.encode()
            s = registry.Source("x", "https://example.org/download", len(data), "2026-10-26", "sample", sha256=sha(data), file="sample.csv", command="x")
            (t / "copy.csv").write_bytes(data)
            fetch(s, t / "copy.csv")
            side = json.loads((t / "raw" / "sample.csv.lineage.json").read_text())
            ck("a copy with the pinned digest is kept, with its lineage", (t / "raw" / "sample.csv").read_bytes() == data and side["sha256"] == sha(data))
            (t / "other.csv").write_bytes(data + b"\n")
            (t / "raw" / "sample.csv").unlink()
            try:
                fetch(s, t / "other.csv")
                ck("a file that is not the pinned one is refused and not kept", False)
            except SystemExit as e:
                ck("a file that is not the pinned one is refused and not kept", str(e).startswith("REFUSED") and not (t / "raw" / "sample.csv").exists())
        finally:
            paths.RAW = old
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) >= 2 and argv[0] == "fetch":
        fetch(registry.SOURCES[argv[1]], Path(argv[argv.index("--from") + 1]) if "--from" in argv else None)
        return 0
    if len(argv) == 2 and argv[0] == "scope":
        print(scope(argv[1]).unstack().to_string())
        return 0
    if len(argv) == 2 and argv[0] == "audit":
        out = audit(argv[1])
        dest = paths.DOCS_DATA / f"clingen_audit_{argv[1]}.json"
        tmp = dest.with_name(dest.name + ".tmp")
        tmp.write_text(json.dumps(out, indent=1) + "\n")
        os.replace(tmp, dest)
        print(f"  {dest}")
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
