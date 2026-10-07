"""The evidence each disease source carries for its own pairs, as one number per (target, disease, field).

A source that grades its rows does so in its own column: DisGeNET a score, Expression Atlas a p-value
and a fold change, JensenLab text mining a z-score. This reads those columns into one long table. A
pair a source asserts in several rows takes the strongest value. Which fields are read is fixed in
docs/PROTOCOL.md; four sources carry none that varies and have no rows here.

Run:  python -m pharos_gap_kg.evidence pharos319
      python -m pharos_gap_kg.evidence --selftest
"""
from __future__ import annotations

import json
import re
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import ctxtables, graph, paths

PUBMED_IDS = re.compile(r"(\d+)\s+PubMed IDs")
# fields a paper-count label is not independent of (docs/PROTOCOL.md)
NOT_INDEPENDENT = {("JensenLab Text Mining", "zscore"), ("DisGeNET", "score"), ("DisGeNET", "pubmed_ids")}
NO_FIELD = ("UniProt Disease", "JensenLab Knowledge AmyCo", "JensenLab Knowledge MedlinePlus", "JensenLab Knowledge UniProtKB-KW")


def _count(s: pd.Series, sep: str) -> pd.Series:
    """How many distinct non-empty parts each delimited string has."""
    return s.fillna("").astype(str).map(lambda x: len({p for p in x.split(sep) if p.strip()})).astype(float)


def fields(rows: pd.DataFrame) -> pd.DataFrame:
    """source, protein_id, did, mondoid, field, value: one row per source row and field, higher meaning stronger."""
    out = []

    def add(src: str, field: str, value: pd.Series) -> None:
        g = rows.loc[value.index, ["protein_id", "did", "mondoid"]].assign(source=src, field=field, value=value.astype(float))
        out.append(g[np.isfinite(g.value)])

    d = rows[rows.dtype == "DisGeNET"]
    add("DisGeNET", "score", d.score)
    add("DisGeNET", "pubmed_ids", d.evidence.fillna("").str.extract(PUBMED_IDS)[0].fillna(0).astype(float))
    add("DisGeNET", "upstream_sources", _count(d.source, ";"))
    e = rows[rows.dtype == "Expression Atlas"]
    p = pd.to_numeric(e.pvalue, errors="coerce")
    add("Expression Atlas", "minus_log10_p", -np.log10(p.clip(lower=1e-300)))
    add("Expression Atlas", "abs_log2_fold_change", e.log2foldchange.abs())
    t = rows[rows.dtype == "JensenLab Text Mining"]
    add("JensenLab Text Mining", "zscore", t.zscore)
    m = rows[rows.dtype == "Monarch"]
    add("Monarch", "s2o", m.S2O)
    r = rows[rows.dtype == "eRAM"]
    add("eRAM", "upstream_sources", _count(r.source, "|"))
    c = rows[rows.dtype == "CTD"]
    add("CTD", "therapeutic", c.evidence.fillna("").str.contains("therapeutic").astype(float))
    g = rows[rows.dtype == "DrugCentral Indication"]
    n = g.groupby(["protein_id", "did"], dropna=False).drug_name.transform("nunique")
    add("DrugCentral Indication", "drugs", n)
    return pd.concat(out, ignore_index=True)


def per_pair(long: pd.DataFrame, target_of: dict) -> pd.DataFrame:
    """target, disease, source, field, value: the strongest value a source gives a pair, ids as in disease_links."""
    x = long.copy()
    x["disease"] = [ctxtables.disease_id(m, d) for m, d in zip(x.mondoid, x.did)]
    x["target"] = x.protein_id.map(target_of)
    x = x.dropna(subset=["disease", "target"])
    return x.groupby(["target", "disease", "source", "field"], as_index=False).value.max()


def build(release: str) -> dict:
    d = paths.tables(release)
    cols = ["dtype", "protein_id", "did", "mondoid", "evidence", "zscore", "drug_name", "log2foldchange", "pvalue", "score", "source", "S2O"]
    rows = pq.read_table(d / "disease.parquet", columns=cols).to_pandas()
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot"]).to_pandas()
    ev = per_pair(fields(rows), dict(zip(prot.id, "UniProtKB:" + prot.uniprot)))
    C = paths.context(release)
    links = pd.read_parquet(C / "disease_links.parquet")
    have = set(zip(links.target, links.disease, links.source))
    stray = [k for k in zip(ev.target, ev.disease, ev.source) if k not in have]
    if stray:
        raise SystemExit(f"REFUSED: {len(stray)} evidence rows are for pairs disease_links does not hold, e.g. {stray[0]}; rebuild the context tables")
    ev.to_parquet(C / "disease_evidence.parquet", index=False)
    per = {}
    for (s, f), g in ev.groupby(["source", "field"]):
        n = int((links.source == s).sum())
        per[f"{s} | {f}"] = {"pairs_with_value": len(g), "pairs_in_source": n, "distinct_values": int(g.value.nunique()),
                             "min": float(g.value.min()), "median": float(g.value.median()), "max": float(g.value.max())}
        if len(g) != n:
            raise SystemExit(f"REFUSED: {s} has {n} pairs and {len(g)} have a value for {f}; a field in the protocol covers every pair of its source")
    man = {"release": release, "disease_links_sha256": json.loads((C / "manifest.json").read_text())["tables"]["disease_links"]["sha256"],
           "rows": len(ev), "sha256": graph.file_sha(C / "disease_evidence.parquet"), "fields": per, "sources_with_no_field": list(NO_FIELD)}
    (paths.DOCS_DATA / f"evidence_fields_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    na = float("nan")
    rows = pd.DataFrame([
        ("DisGeNET", 1, "UMLS:C1", "MONDO:1", "3 PubMed IDs; 1 SNPs", na, None, na, None, 0.4, "CTD_human;UNIPROT", na),
        ("DisGeNET", 1, "UMLS:C1", "MONDO:1", "0 SNPs", na, None, na, None, 0.1, "CTD_human", na),
        ("Expression Atlas", 2, "DOID:2", "MONDO:2", None, na, None, -2.0, "0.01", na, None, na),
        ("Expression Atlas", 2, "DOID:2", "MONDO:2", None, na, None, 1.5, "0.0001", na, None, na),
        ("JensenLab Text Mining", 3, "DOID:3", None, None, 6.5, None, na, None, na, None, na),
        ("Monarch", 3, "MIM:100", None, "ECO:1", na, None, na, None, na, None, 88.0),
        ("eRAM", 4, "DOID:4", "MONDO:4", None, na, None, na, None, na, "CLINVAR|GHR|GHR", na),
        ("CTD", 4, "MESH:D4", "MONDO:4", "marker/mechanism|therapeutic", na, None, na, None, na, None, na),
        ("CTD", 5, "MESH:D4", "MONDO:4", "marker/mechanism", na, None, na, None, na, None, na),
        ("DrugCentral Indication", 5, "DOID:5", None, None, na, "a", na, None, na, None, na),
        ("DrugCentral Indication", 5, "DOID:5", None, None, na, "b", na, None, na, None, na),
        ("UniProt Disease", 5, "MIM:5", None, "7", na, None, na, None, na, None, na),
    ], columns=["dtype", "protein_id", "did", "mondoid", "evidence", "zscore", "drug_name", "log2foldchange", "pvalue", "score", "source", "S2O"])
    T = {i: f"UniProtKB:P{i}" for i in range(1, 6)}
    ev = per_pair(fields(rows), T)
    v = {(s, f, t): x for s, f, t, x in zip(ev.source, ev.field, ev.target, ev.value)}
    ck("a pair asserted in two rows takes the strongest score, paper count and source count",
       v[("DisGeNET", "score", "UniProtKB:P1")] == 0.4 and v[("DisGeNET", "pubmed_ids", "UniProtKB:P1")] == 3 and v[("DisGeNET", "upstream_sources", "UniProtKB:P1")] == 2)
    ck("a p-value is read as minus log10 and a fold change by its size",
       abs(v[("Expression Atlas", "minus_log10_p", "UniProtKB:P2")] - 4.0) < 1e-9 and v[("Expression Atlas", "abs_log2_fold_change", "UniProtKB:P2")] == 2.0)
    ck("z-score and S2O are read as they stand", v[("JensenLab Text Mining", "zscore", "UniProtKB:P3")] == 6.5 and v[("Monarch", "s2o", "UniProtKB:P3")] == 88.0)
    ck("upstream sources are counted once each", v[("eRAM", "upstream_sources", "UniProtKB:P4")] == 2)
    ck("therapeutic is 1 where the evidence names it and 0 where it does not", v[("CTD", "therapeutic", "UniProtKB:P4")] == 1 and v[("CTD", "therapeutic", "UniProtKB:P5")] == 0)
    ck("drugs are counted per pair", v[("DrugCentral Indication", "drugs", "UniProtKB:P5")] == 2)
    ck("a source with no field has no rows, and ids are written as in disease_links",
       "UniProt Disease" not in set(ev.source) and set(ev[ev.source == "Monarch"].disease) == {"OMIM:100"})
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 1:
        raise SystemExit(__doc__)
    man = build(argv[0])
    for k, r in man["fields"].items():
        print(f"  {k:46s} {r['pairs_with_value']:>8,} pairs  {r['distinct_values']:>7,} values  min {r['min']:.3g}  median {r['median']:.3g}  max {r['max']:.3g}")
    print(f"  no field: {', '.join(man['sources_with_no_field'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
