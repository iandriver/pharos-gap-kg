"""Studies added after the release: their links, kept beside the bundle and apart from it.

A study published after the release is added by hand, one protocol entry each. Its links go to
data/later/<release>/<key>.parquet with a manifest of digests. The bundle, its manifest and the release
archive are unchanged, and only the answers that say what was reported after the release read the
table (docs/PROTOCOL.md, "A study added after the release" and its amendment of 29 Oct).

The one study so far: Uffelmann et al., Nature, 28 Oct 2026 (doi:10.1038/s41586-026-11206-y), an
Alzheimer disease GWAS. From table SD2 of its Supplementary Data, each predicted effector gene that
names a graph gene, by the matching rule of `adstudy`, linked to Alzheimer disease (MONDO:0004975) with
its locus. The file is licensed CC BY-NC-ND 4.0, so the table is built from the local file, is not
committed, and is not in the release archive; docs/data/later_<release>.json records its counts and digest.

Run:  python -I -m pharos_gap_kg.later build pharos319      # reads the registered file in memory
      python -m pharos_gap_kg.later --selftest
"""
from __future__ import annotations

import json
import sys

import pandas as pd

from pharos_gap_kg import adstudy, paths, registry

STUDIES = {
    "uffelmann2026_supp": {"citation": "Uffelmann et al., Nature (2026)", "doi": "10.1038/s41586-026-11206-y", "published": "2026-10-28", "added": "2026-10-29",
                           "year": 2026, "trait": adstudy.TRAIT, "licence": "CC BY-NC-ND 4.0", "table": "SD2", "loci": adstudy.LOCI, "new_loci": adstudy.NEW_LOCI,
                           "what": "a genome-wide association study of Alzheimer disease: each locus with the gene or genes the study predicts act there"},
}
NAMED_BY = ("FLAMES", "Exome", "Coloc_MR", "Literature")
COLUMNS = ["target", "symbol", "study_symbol", "trait", "year", "study", "locus", "rsid", "chromosome", "start", "end", "p_value", "new_locus", "also_named_by"]


def later_dir(release: str):
    return paths.DATA / "later" / release


def loci(rows: list[list[str]]) -> pd.DataFrame:
    """SD2's loci as `adstudy.study_loci` reads them (which refuses unless 127 loci, 48 new), with each
    locus's place and the genes of the columns that may also name an effector gene."""
    L = adstudy.study_loci(rows)
    head = next(i for i, r in enumerate(rows) if len(r) > 1 and r[1] == "Predicted effector gene(s)")
    cols = rows[head]
    t = pd.DataFrame([(r + [""] * len(cols))[:len(cols)] for r in rows[head + 1:] if len(r) > 2 and r[2].strip()], columns=cols)
    place = pd.DataFrame({"locus": t["locus"].astype(int), "rsid": t["rsid"], "chromosome": t["chromosome"].astype(str),
                          "start": t["start"].astype(int), "end": t["end"].astype(int), "p_value": pd.to_numeric(t["p_value"], errors="raise"),
                          **{f"by_{c}": t[c].map(adstudy.genes) for c in NAMED_BY}})
    return L.merge(place, on="locus", validate="one_to_one")


def links(L: pd.DataFrame, m: dict[str, set[str]], ids: dict[str, list[str]], key: str) -> pd.DataFrame:
    """One row per graph target and locus: each predicted effector gene that names a graph gene.
    m: study symbol -> Pharos symbols (`adstudy.symbol_map`); ids: Pharos symbol -> target ids."""
    s = STUDIES[key]
    to = lambda gs: set().union(*(m.get(g, set()) for g in gs)) if gs else set()
    out = []
    for r in L.itertuples():
        by = {c: to(getattr(r, f"by_{c}")) for c in NAMED_BY}
        for g in r.effector:
            for p in sorted(m.get(g, ())):
                for t in ids.get(p, ()):
                    out.append((t, p, g, s["trait"], s["year"], key, int(r.locus), r.rsid, r.chromosome, int(r.start), int(r.end), float(r.p_value), bool(r.new),
                                [c for c in NAMED_BY if p in by[c]]))
    return pd.DataFrame(out, columns=COLUMNS).drop_duplicates(["target", "locus"]).sort_values(["locus", "symbol", "target"]).reset_index(drop=True)


def build(release: str, key: str = adstudy.KEY) -> dict:
    from pharos_gap_kg import ask, graph
    src = registry.require(key)
    blob = src.read_bytes()
    if adstudy.sha(blob) != registry.SOURCES[key].sha256:
        raise SystemExit(f"REFUSED: {src.name} is not the file the registry pins")
    L = loci(adstudy.sheet(blob, STUDIES[key]["table"]))
    G = ask.Graph(paths.bundle(release))
    m = adstudy.symbol_map(G, release)
    ids: dict[str, list[str]] = {}
    for t, p in zip(G.t.id, G.t["name"]):
        ids.setdefault(p, []).append(t)
    T = links(L, m, ids, key)
    eff = sorted({g for gs in L.effector for g in gs})
    d = later_dir(release)
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{key}.parquet"
    T.to_parquet(f, index=False)
    mf = d / "manifest.json"
    man = json.loads(mf.read_text()) if mf.is_file() else {"release": release, "files": {}, "studies": {}}
    if man["release"] != release:
        raise SystemExit(f"REFUSED: {mf} is for {man['release']}")
    man["files"][f.name] = graph.file_sha(f)
    man["studies"][key] = STUDIES[key]
    mf.write_text(json.dumps(man, indent=1, sort_keys=True) + "\n")
    rec = {"release": release, "study": key, **STUDIES[key], "source_sha256": registry.SOURCES[key].sha256, "kept_in": f"data/later/{release}/{f.name}, not committed",
           "loci_read": len(L), "new_loci_read": int(L.new.sum()), "effector_symbols": len(eff), "effector_symbols_naming_a_graph_gene": sum(1 for g in eff if m.get(g)),
           "rows": len(T), "targets": int(T.target.nunique()), "pharos_symbols": int(T.symbol.nunique()), "loci_with_a_row": int(T.locus.nunique()),
           "pharos_symbols_at_new_loci": int(T[T.new_locus].symbol.nunique()),
           "rows_also_named_by": {c: int(T.also_named_by.map(lambda x, c=c: c in x).sum()) for c in NAMED_BY}, "file_sha256": man["files"][f.name]}
    (paths.DOCS_DATA / f"later_{release}.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    L = pd.DataFrame({"locus": [1, 2], "new": [True, False], "effector": [["AAA", "BBB"], ["CCC"]], "rsid": ["rs1", "rs2"], "chromosome": ["1", "2"],
                      "start": [10, 20], "end": [15, 25], "p_value": [1e-9, 2e-8], "by_FLAMES": [["BBB"], []], "by_Exome": [[], ["C-OLD"]], "by_Coloc_MR": [["AAA", "BBB"], []],
                      "by_Literature": [[], []]})
    m = {"AAA": {"AAA"}, "BBB": {"BBB"}, "CCC": {"CCC"}, "C-OLD": {"CCC"}}
    T = links(L, m, {"AAA": ["t:A"], "BBB": ["t:B", "t:B2"], "CCC": ["t:C"]}, adstudy.KEY)
    ck("one row per target and locus, two entries of one symbol both kept", T[["target", "locus"]].values.tolist() == [["t:A", 1], ["t:B", 1], ["t:B2", 1], ["t:C", 2]])
    ck("the columns that also name the gene are read through the same symbol map", T.also_named_by.tolist() == [["Coloc_MR"], ["FLAMES", "Coloc_MR"], ["FLAMES", "Coloc_MR"], ["Exome"]])
    ck("the trait, year, study and whether the locus is new travel on each row",
       set(T.trait) == {"MONDO:0004975"} and set(T.year) == {2026} and set(T.study) == {adstudy.KEY} and T.new_locus.tolist() == [True, True, True, False])
    ck("a study symbol that names no graph gene adds no row", links(L.assign(effector=[["ZZZ"], []]), m, {}, adstudy.KEY).empty)
    ck("the table is kept beside the bundle, not in it", later_dir("x").parent.parent == paths.DATA and later_dir("x") != paths.bundle("x") and later_dir("x").parent != paths.bundle("x"))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 2 or argv[0] != "build":
        raise SystemExit(__doc__)
    r = build(argv[1])
    print(f"{r['release']}: {r['study']}: {r['loci_read']} loci ({r['new_loci_read']} new); {r['effector_symbols_naming_a_graph_gene']} of {r['effector_symbols']} effector "
          f"symbols name a graph gene; {r['rows']} rows, {r['pharos_symbols']} Pharos symbols at {r['loci_with_a_row']} loci -> data/later/{r['release']}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
