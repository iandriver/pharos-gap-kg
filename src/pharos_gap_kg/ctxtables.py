"""The rest of a Pharos release as link tables: disease sources, interactions, pathways, GO, compounds.

These are what the database holds besides the dated GWAS links. Each is written as one small table of
CURIEs with the source beside every row, and a manifest that counts what was left out and why.

LEFT OUT, each counted in the manifest:
  - the disease source "JensenLab Experiment TIGA", which is the GWAS data again (findings/CONTEXT.md);
  - disease rows with no identifier for the disease (a name only);
  - interaction rows that pair a protein with itself (the release marks them "mock");
  - STRING interactions below a combined score of 700, unless BioPlex or Reactome also report the pair;
  - compounds active on one target only, which tie nothing together.

Run:  python -m pharos_gap_kg.ctxtables pharos319
      python -m pharos_gap_kg.ctxtables --selftest
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import graph, paths, sqldump

GWAS_AGAIN = "JensenLab Experiment TIGA"
STRING_MIN = 700
PATHWAY_SOURCES = ("Reactome", "KEGG", "WikiPathways")
FILES = ["disease_links", "partners", "pathways", "go_terms", "compounds"]


def disease_id(mondo, did) -> str | None:
    """The MONDO id where the release gives one, else the source's own id, with OMIM written one way."""
    if isinstance(mondo, str) and mondo:
        return mondo
    if not isinstance(did, str) or not did:
        return None
    return "OMIM:" + did[4:] if did.startswith("MIM:") else did


def disease_links(rows: pd.DataFrame, target_of: dict) -> tuple[pd.DataFrame, dict]:
    """target, disease, source: one row per (target, disease, source). `rows` has dtype, protein_id, did, mondoid."""
    left = {"gwas_again": int((rows.dtype == GWAS_AGAIN).sum())}
    r = rows[rows.dtype != GWAS_AGAIN].copy()
    r["disease"] = [disease_id(m, d) for m, d in zip(r.mondoid, r.did)]
    left["no_disease_id"] = int(r.disease.isna().sum())
    r = r.dropna(subset=["disease"])
    r["target"] = r.protein_id.map(target_of)
    out = r.rename(columns={"dtype": "source"})[["target", "disease", "source"]].dropna().drop_duplicates().reset_index(drop=True)
    return out, left


def partners(rows: pd.DataFrame, target_of: dict) -> tuple[pd.DataFrame, dict]:
    """a, b (a < b), sources, string_score: undirected interaction pairs that pass the bar."""
    self_pairs = rows.protein_id == rows.other_id
    r = rows[~self_pairs]
    curated = r.ppitypes.str.contains("BioPlex|Reactome", na=False)
    keep = curated | (r.score >= STRING_MIN)
    left = {"self_pairs": int(self_pairs.sum()), "string_below_bar": int((~keep).sum())}
    r = r[keep]
    a, b = r.protein_id.map(target_of), r.other_id.map(target_of)
    lo, hi = np.where(a < b, a, b), np.where(a < b, b, a)
    out = pd.DataFrame({"a": lo, "b": hi, "sources": r.ppitypes.to_numpy(), "string_score": r.score.to_numpy()}).dropna(subset=["a", "b"])
    out = out.sort_values("string_score", ascending=False).drop_duplicates(["a", "b"]).sort_values(["a", "b"]).reset_index(drop=True)
    return out, left


def memberships(rows: pd.DataFrame, target_col: str, group: pd.Series, target_of: dict, min_size: int = 2) -> pd.DataFrame:
    """target, grp: one row per membership, for groups with at least `min_size` targets."""
    out = pd.DataFrame({"target": rows[target_col].map(target_of), "grp": group.to_numpy()}).dropna().drop_duplicates()
    size = out.groupby("grp").target.transform("size")
    return out[size >= min_size].reset_index(drop=True)


def build(release: str) -> dict:
    sqldump.need(release, ["protein", "t2tc", "disease", "ncats_ppi", "pathway", "goa", "ncats_ligand_activity"])
    d = paths.tables(release)
    rd = lambda t, cols: pq.read_table(d / f"{t}.parquet", columns=cols).to_pandas()
    prot = rd("protein", ["id", "uniprot"])
    by_protein = dict(zip(prot.id, "UniProtKB:" + prot.uniprot))
    t2 = rd("t2tc", ["target_id", "protein_id"])
    by_target = {t: by_protein[p] for t, p in zip(t2.target_id, t2.protein_id) if p in by_protein}
    out = paths.context(release)
    out.mkdir(parents=True, exist_ok=True)
    man = {"release": release, "source_sha256": json.loads((d / "_manifest.json").read_text())["sha256"], "left_out": {}, "tables": {}}

    dl, left = disease_links(rd("disease", ["dtype", "protein_id", "did", "mondoid"]), by_protein)
    man["left_out"]["disease_links"] = left
    pp, left = partners(rd("ncats_ppi", ["ppitypes", "protein_id", "other_id", "score"]), by_protein)
    man["left_out"]["partners"] = left
    pw = rd("pathway", ["protein_id", "pwtype", "id_in_source"])
    pw = pw[pw.pwtype.isin(PATHWAY_SOURCES) & pw.id_in_source.notna()]
    pw = memberships(pw, "protein_id", pw.pwtype + ":" + pw.id_in_source.astype(str), by_protein)
    go = rd("goa", ["protein_id", "go_id"])
    go = memberships(go, "protein_id", go.go_id, by_protein)
    la = rd("ncats_ligand_activity", ["ncats_ligand_id", "target_id"])
    one_target = int((la.drop_duplicates().groupby("ncats_ligand_id").target_id.nunique() < 2).sum())
    cp = memberships(la, "target_id", "ligand:" + la.ncats_ligand_id.astype(str), by_target)
    man["left_out"]["compounds"] = {"active_on_one_target_only": one_target}
    for name, df in zip(FILES, (dl, pp, pw, go, cp)):
        df.to_parquet(out / f"{name}.parquet", index=False)
        man["tables"][name] = {"rows": len(df), "sha256": graph.file_sha(out / f"{name}.parquet")}
    pairs = dl.groupby(["target", "disease"]).source.nunique()
    man["disease_pairs"] = {"pairs": int(len(pairs)), "one_source": int((pairs == 1).sum()), "two_or_more": int((pairs > 1).sum()),
                            "sources": int(dl.source.nunique()), "targets": int(dl.target.nunique()), "diseases": int(dl.disease.nunique())}
    man["groups"] = {"pathways": int(pw.grp.nunique()), "go_terms": int(go.grp.nunique()), "compounds": int(cp.grp.nunique()),
                     "partner_pairs": int(len(pp)), "targets_with_a_partner": int(pd.concat([pp.a, pp.b]).nunique())}
    (out / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"context_tables_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    T = {1: "UniProtKB:A", 2: "UniProtKB:B", 3: "UniProtKB:C"}
    ck("a disease takes its MONDO id, else the source's id, and MIM is written OMIM",
       [disease_id("MONDO:1", "DOID:9"), disease_id(None, "MIM:219500"), disease_id(float("nan"), "UMLS:C1"), disease_id(None, None)] == ["MONDO:1", "OMIM:219500", "UMLS:C1", None])
    rows = pd.DataFrame({"dtype": ["DisGeNET", "DisGeNET", "CTD", GWAS_AGAIN, "DrugCentral Indication", "CTD"],
                         "protein_id": [1, 1, 1, 2, 2, 9], "did": ["UMLS:C1", "UMLS:C1", "MESH:D1", "DOID:5", None, "MESH:D2"],
                         "mondoid": ["MONDO:1", "MONDO:1", "MONDO:1", "MONDO:5", None, None]})
    dl, left = disease_links(rows, T)
    ck("a source asserting a pair twice is one row, and two sources are two", dl.values.tolist() == [["UniProtKB:A", "MONDO:1", "DisGeNET"], ["UniProtKB:A", "MONDO:1", "CTD"]])
    ck("the GWAS source, a row with no disease id and an unknown protein are left out and the first two counted", left == {"gwas_again": 1, "no_disease_id": 1})
    pr = pd.DataFrame({"ppitypes": ["mock", "STRINGDB", "STRINGDB", "BioPlex", "BioPlex,STRINGDB", "STRINGDB"],
                       "protein_id": [1, 1, 2, 1, 3, 3], "other_id": [1, 2, 1, 3, 1, 2], "score": [None, 900.0, 900.0, None, 150.0, 400.0]})
    pp, left = partners(pr, T)
    ck("a pair is undirected and listed once", pp[["a", "b"]].values.tolist() == [["UniProtKB:A", "UniProtKB:B"], ["UniProtKB:A", "UniProtKB:C"]])
    ck("a self pair and a weak STRING-only pair are left out and counted; a curated pair needs no score", left == {"self_pairs": 1, "string_below_bar": 1})
    m = memberships(pd.DataFrame({"p": [1, 2, 3, 1, 1]}), "p", pd.Series(["g1", "g1", "g2", "g2", "g2"]), T)
    ck("a group needs two distinct targets, and a repeated membership is one", sorted(map(tuple, m.values.tolist())) == [("UniProtKB:A", "g1"), ("UniProtKB:A", "g2"), ("UniProtKB:B", "g1"), ("UniProtKB:C", "g2")])
    ck("a group of one is dropped", memberships(pd.DataFrame({"p": [1]}), "p", pd.Series(["solo"]), T).empty)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 1:
        raise SystemExit(__doc__)
    man = build(argv[0])
    for k, v in man["tables"].items():
        print(f"  {k:14s} {v['rows']:>10,} rows")
    print(f"  disease pairs: {man['disease_pairs']}")
    print(f"  groups: {man['groups']}")
    print(f"  left out: {man['left_out']}")
    print(f"  -> {paths.context(man['release'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
