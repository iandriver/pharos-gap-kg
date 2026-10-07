"""Development levels at two dates, joined by UniProt accession: who moved, and who gained a ligand.

Pharos grades every human protein target: Tclin (an approved drug with a known mechanism), Tchem (a
small molecule above a potency cutoff), Tbio (biological evidence only), Tdark (little of anything).
The archived database gives each target's level as of its own contents; the Pharos 4.0.0 table gives
the level in May 2026. A target with no qualifying ligand at the first date and one at the second is
a dated, public label for "this target-ligand link did not exist and now does".

THE RELEASE NOTES ARE A SECOND SOURCE, so they are read and compared here on every run. The notes
count upward moves only. This module counts both directions, because a level that can fall is a level
some of whose rises are the pipeline changing and not the knowledge.

Run:  python -m pharos_gap_kg.levels pharos319              # against Pharos 4.0.0 and its release notes
      python -m pharos_gap_kg.levels pharos316 pharos319    # two database releases against each other
      python -m pharos_gap_kg.levels --selftest
"""
from __future__ import annotations

import json
import re
import sys

import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import paths, registry

ORDER = ["Tdark", "Tbio", "Tchem", "Tclin"]
RANK = {k: i for i, k in enumerate(ORDER)}
HAS_LIGAND = {"Tchem", "Tclin"}
NOTE_ROW = re.compile(r"(T(?:dark|bio|chem|clin))\s*-(?:&gt;|>)\s*(T(?:dark|bio|chem|clin))\s+([\d,]+)")


def old_levels(release: str) -> pd.DataFrame:
    """uniprot, sym, tdl, fam for every target of an extracted database release."""
    d = paths.tables(release)
    need = ["target", "t2tc", "protein"]
    missing = [t for t in need if not (d / f"{t}.parquet").is_file()]
    if missing:
        raise SystemExit(f"{release} has no extracted {missing}; run:  "
                         f"python -m pharos_gap_kg.sqldump extract {release} --tables {','.join(need)}")
    rd = lambda t, cols: pq.read_table(d / f"{t}.parquet", columns=cols).to_pandas()
    t = rd("target", ["id", "tdl", "fam"]).rename(columns={"id": "target_id"})
    link = rd("t2tc", ["target_id", "protein_id"])
    p = rd("protein", ["id", "uniprot", "sym"]).rename(columns={"id": "protein_id"})
    out = t.merge(link, on="target_id").merge(p, on="protein_id")
    if len(out) != len(t) or out.uniprot.duplicated().any():
        raise SystemExit(f"REFUSED {release}: targets and proteins are not one to one "
                         f"({len(t):,} targets, {len(out):,} joined, {int(out.uniprot.duplicated().sum())} repeated accessions)")
    return out[["uniprot", "sym", "tdl", "fam", "target_id"]]


def new_levels() -> pd.DataFrame:
    d = pd.read_csv(registry.require("pharos400_canonicals"), dtype=str, keep_default_na=False,
                    usecols=["uniprot_id", "tdl", "symbol", "idg_family", "tdl_ligand_count", "tdl_drug_count"])
    if d.uniprot_id.duplicated().any():
        raise SystemExit("REFUSED: the Pharos 4.0.0 canonical table repeats an accession")
    return d.rename(columns={"uniprot_id": "uniprot"})


def transitions(old: pd.DataFrame, new: pd.DataFrame) -> dict:
    m = old.merge(new, on="uniprot", how="outer", suffixes=("_old", "_new"), indicator=True)
    both = m[m._merge == "both"]
    bad = sorted((set(both.tdl_old) | set(both.tdl_new)) - set(ORDER))
    if bad:
        raise SystemExit(f"REFUSED: unknown development levels {bad}")
    matrix = {a: {b: int(((both.tdl_old == a) & (both.tdl_new == b)).sum()) for b in ORDER} for a in ORDER}
    up = both.tdl_new.map(RANK) > both.tdl_old.map(RANK)
    down = both.tdl_new.map(RANK) < both.tdl_old.map(RANK)
    cand = both[~both.tdl_old.isin(HAS_LIGAND)]
    gained = cand.tdl_new.isin(HAS_LIGAND)
    had = both[both.tdl_old.isin(HAS_LIGAND)]
    lost = ~had.tdl_new.isin(HAS_LIGAND)
    fam = cand.assign(gained=gained, fam=cand.fam.fillna("none")).groupby("fam").gained.agg(["size", "sum"])
    return {
        "shared": len(both), "only_old": int((m._merge == "left_only").sum()), "only_new": int((m._merge == "right_only").sum()),
        "matrix": matrix, "up": int(up.sum()), "down": int(down.sum()), "unchanged": int((~up & ~down).sum()),
        "ligand_candidates": len(cand), "ligand_gained": int(gained.sum()),
        "ligand_holders": len(had), "ligand_lost": int(lost.sum()),
        "gained_by_family": {k: {"candidates": int(r["size"]), "gained": int(r["sum"])} for k, r in fam.iterrows()},
    }


def level_against_ligands(release: str) -> dict:
    """Does a release's level agree with its own ligand table? Targets graded as holding a ligand
    with no activity row, and targets with activity rows graded as holding none."""
    d = paths.tables(release)
    p = d / "ncats_ligand_activity.parquet"
    if not p.is_file():
        raise SystemExit(f"{release} has no extracted ncats_ligand_activity; run:  "
                         f"python -m pharos_gap_kg.sqldump extract {release} --tables ncats_ligand_activity")
    with_rows = set(pq.read_table(p, columns=["target_id"]).to_pandas().target_id)
    t = pq.read_table(d / "target.parquet", columns=["id", "tdl"]).to_pandas()
    graded = t.tdl.isin(HAS_LIGAND)
    has = t.id.isin(with_rows)
    return {"targets_with_activity_rows": int(has.sum()), "graded_as_holding": int(graded.sum()),
            "rows_but_graded_without": int((has & ~graded).sum()),
            "rows_but_graded_without_by_level": {k: int(v) for k, v in t.tdl[has & ~graded].value_counts().items()},
            "graded_holding_but_no_rows": int((graded & ~has).sum()),
            "graded_holding_but_no_rows_by_level": {k: int(v) for k, v in t.tdl[graded & ~has].value_counts().items()}}


def ligand_label(release: str, old: pd.DataFrame, new: pd.DataFrame) -> dict:
    """The label "gained a qualifying ligand by Pharos 4.0.0", with candidates defined by the release's
    ligand TABLE and not its grade: a target that already has activity rows did not gain a ligand,
    whatever it was graded."""
    with_rows = set(pq.read_table(paths.tables(release) / "ncats_ligand_activity.parquet", columns=["target_id"]).to_pandas().target_id)
    both = old.merge(new, on="uniprot", suffixes=("_old", "_new"))
    by_grade = both[~both.tdl_old.isin(HAS_LIGAND)]
    had_rows = by_grade.target_id.isin(with_rows)
    gained = by_grade.tdl_new.isin(HAS_LIGAND)
    clean = by_grade[~had_rows]
    return {"candidates_by_grade": len(by_grade), "gained_by_grade": int(gained.sum()),
            "candidates_with_activity_rows": int(had_rows.sum()), "of_those_counted_as_gained": int((had_rows & gained).sum()),
            "candidates": len(clean), "gained": int(clean.tdl_new.isin(HAS_LIGAND).sum()),
            "gained_from": {k: int(v) for k, v in clean.tdl_old[clean.tdl_new.isin(HAS_LIGAND)].value_counts().items()}}


def notes_transitions() -> dict[str, int]:
    """The upward moves the Pharos 4.0.0 release notes report, as {'Tbio>Tchem': 1320, ...}."""
    text = re.sub(r"<[^>]*>", " ", registry.require("pharos400_readme").read_text(errors="replace"))
    return {f"{a}>{b}": int(n.replace(",", "")) for a, b, n in NOTE_ROW.findall(text)}


def against_notes(r: dict, notes: dict[str, int]) -> list[dict]:
    return [{"move": k, "notes": n, "measured": r["matrix"][k.split(">")[0]][k.split(">")[1]],
             "difference": r["matrix"][k.split(">")[0]][k.split(">")[1]] - n} for k, n in sorted(notes.items())]


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    old = pd.DataFrame({"uniprot": list("ABCDEFG"), "sym": list("abcdefg"),
                        "tdl": ["Tdark", "Tdark", "Tbio", "Tbio", "Tchem", "Tclin", "Tbio"],
                        "fam": ["Kinase", None, "Kinase", None, None, None, None]})
    new = pd.DataFrame({"uniprot": list("ABCDEFH"), "tdl": ["Tbio", "Tchem", "Tchem", "Tdark", "Tbio", "Tclin", "Tdark"]})
    r = transitions(old, new)
    ck("accessions are matched, and the unmatched counted on each side", (r["shared"], r["only_old"], r["only_new"]) == (6, 1, 1))
    ck("rises and falls are both counted", (r["up"], r["down"], r["unchanged"]) == (3, 2, 1))
    ck("the matrix holds each move once", r["matrix"]["Tdark"]["Tchem"] == 1 and r["matrix"]["Tchem"]["Tbio"] == 1 and sum(sum(v.values()) for v in r["matrix"].values()) == 6)
    ck("a ligand is gained only from a level without one", (r["ligand_candidates"], r["ligand_gained"]) == (4, 2))
    ck("a ligand lost is counted against those that held one", (r["ligand_holders"], r["ligand_lost"]) == (2, 1))
    ck("family counts sum to the candidates", sum(v["candidates"] for v in r["gained_by_family"].values()) == 4 and r["gained_by_family"]["Kinase"] == {"candidates": 2, "gained": 1})
    cmp_ = against_notes(r, {"Tdark>Tchem": 1, "Tbio>Tchem": 3})
    ck("a move is compared with the notes' count", [c["difference"] for c in cmp_] == [-2, 0])
    ck("the notes' arrow is read with or without HTML escaping", NOTE_ROW.findall("Tbio -&gt; Tchem 1,320 Tdark -> Tbio 743") == [("Tbio", "Tchem", "1,320"), ("Tdark", "Tbio", "743")])
    r2 = transitions(old.assign(target_id=range(7)), new)
    ck("an extra column on the old side changes nothing", r2["matrix"] == r["matrix"])
    try:
        transitions(old.assign(tdl="Tnew"), new)
        ck("an unknown level is refused", False)
    except SystemExit:
        ck("an unknown level is refused", True)
    return 0 if ok else 1


def _between(a: str, b: str) -> int:
    """Two database releases against each other: the same counts, with no release notes to compare."""
    r = transitions(old_levels(a), old_levels(b)[["uniprot", "tdl"]])
    print(f"{a} -> {b}: {r['shared']:,} shared accessions, {r['only_old']:,} only in {a}, {r['only_new']:,} only in {b}")
    print("  from \\ to " + "".join(f"{c:>8s}" for c in ORDER))
    for x in ORDER:
        print(f"  {x:9s} " + "".join(f"{r['matrix'][x][y]:>8,}" for y in ORDER))
    print(f"  moved up {r['up']:,}, moved down {r['down']:,}, unchanged {r['unchanged']:,}")
    print(f"  no qualifying ligand in {a}: {r['ligand_candidates']:,}; with one in {b}: {r['ligand_gained']:,}")
    print(f"  held one in {a}: {r['ligand_holders']:,}; without one in {b}: {r['ligand_lost']:,}")
    sha = lambda k: json.loads((paths.tables(k) / "_manifest.json").read_text())["sha256"]
    out = paths.DOCS_DATA / f"levels_{a}_{b}.json"
    out.write_text(json.dumps({"old": a, "old_sha256": sha(a), "new": b, "new_sha256": sha(b), **r}, indent=1) + "\n")
    print(f"  -> {out}")
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) == 2:
        return _between(argv[0], argv[1])
    if len(argv) != 1:
        raise SystemExit(__doc__)
    release = argv[0]
    old, new = old_levels(release), new_levels()
    r = transitions(old, new)
    label = ligand_label(release, old, new)
    notes = notes_transitions()
    cmp_ = against_notes(r, notes)
    print(f"{release} -> Pharos 4.0.0: {r['shared']:,} shared accessions, {r['only_old']:,} only in {release}, {r['only_new']:,} only in 4.0.0")
    print("  from \\ to " + "".join(f"{c:>8s}" for c in ORDER))
    for a in ORDER:
        print(f"  {a:9s} " + "".join(f"{r['matrix'][a][b]:>8,}" for b in ORDER))
    print(f"  moved up {r['up']:,}, moved down {r['down']:,}, unchanged {r['unchanged']:,}")
    print(f"  no qualifying ligand in {release}: {r['ligand_candidates']:,}; with one in 4.0.0: {r['ligand_gained']:,} "
          f"({r['ligand_gained'] / r['ligand_candidates']:.1%})")
    print(f"  held one in {release}: {r['ligand_holders']:,}; without one in 4.0.0: {r['ligand_lost']:,} "
          f"({r['ligand_lost'] / r['ligand_holders']:.1%})")
    print(f"  release notes, upward moves ({sum(notes.values()):,} in all):")
    for c in cmp_:
        print(f"    {c['move']:12s} notes {c['notes']:>6,}  measured {c['measured']:>6,}  difference {c['difference']:+d}")
    own = level_against_ligands(release)
    print(f"  {release} against its own ligand table: {own['targets_with_activity_rows']:,} targets have activity rows, "
          f"{own['graded_as_holding']:,} are graded Tchem or Tclin; {own['rows_but_graded_without']:,} have rows and a lower grade, "
          f"{own['graded_holding_but_no_rows']:,} have the grade and no rows")
    print(f"  label, candidates taken from the ligand table: {label['candidates']:,} targets with no activity row, "
          f"{label['gained']:,} with a ligand in 4.0.0 ({label['gained'] / label['candidates']:.1%}); "
          f"{label['of_those_counted_as_gained']} of the {label['candidates_with_activity_rows']} already holding rows would have counted as gains")
    man = json.loads((paths.tables(release) / "_manifest.json").read_text())
    side = lambda k: json.loads(registry.SOURCES[k].path.with_name(registry.SOURCES[k].filename + ".lineage.json").read_text())["sha256"]
    out = paths.DOCS_DATA / f"levels_{release}_pharos400.json"
    out.write_text(json.dumps({"old": release, "old_sha256": man["sha256"], "new": "pharos400_canonicals",
                               "new_sha256": side("pharos400_canonicals"), "notes_sha256": side("pharos400_readme"),
                               **r, "release_notes": cmp_, "level_against_own_ligand_table": own, "ligand_label": label}, indent=1) + "\n")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
