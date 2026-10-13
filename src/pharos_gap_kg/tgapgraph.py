"""Target-level gaps as KGX: each gap an edge that says what the target lacks and what would close it.

`tgaps` reads, from a release's own evidence, what each target lacks for the next development level.
Here each of those is written as a gap edge in the form `kgx-gaps` defines: from the target to the
kind of gap, with the reason nothing is asserted, the measurement that closes it and the observation
that would retire it unresolved. Nothing is asserted, so the export has gaps and no associations.

    below_tbio_evidence      a Tdark target: how many criteria must clear, the nearest route and its
                             cost, and the share of its band that was promoted in the back-test
    no_qualifying_compound   a target with no compound at its family's potency cutoff

Two fields `kgx-gaps` defines are left empty because they do not apply: `detection_floor` (no effect
size was measured) and `n_required` (no sample size would decide it).

Run:  python -m pharos_gap_kg.tgapgraph pharos319
      python -m pharos_gap_kg.tgapgraph --selftest
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
from kgx_gaps import Mapping, conformance
from kgx_gaps.export import Exporter

from pharos_gap_kg import graph, paths

PREFIX = graph.PREFIX
SOURCE = "infores:pharos"
DARK, LIGAND = "below_tbio_evidence", "no_qualifying_compound"
GAP_NAMES = {DARK: "evidence short of the Tbio level", LIGAND: "no compound at the potency cutoff"}
KILL = "the target is withdrawn from, or merged in, a later release"
OTHER_ROUTE = "an experimental Gene Ontology function or process annotation"
EXTRA = [f"{PREFIX}:{c}" for c in ("level", "criteria_to_clear", "route_cost", "potency_cutoff", "band", "closed_in_backtest", "backtest")]


def release_name(release: str) -> str:
    m = re.fullmatch(r"pharos(\d)(\d+)", release)
    return f"Pharos {m.group(1)}.{m.group(2)}" if m else release


def mapping() -> Mapping:
    m = Mapping(prefix=PREFIX, knowledge_source=SOURCE, gap_types={DARK, LIGAND}, gap_names=dict(GAP_NAMES))
    problems = m.validate()
    if problems:
        raise SystemExit("REFUSED mapping: " + "; ".join(problems))
    return m


def rows(dark: pd.DataFrame, lig: pd.DataFrame, summary: dict, spans: str) -> list[dict]:
    """One row per gap. `dark`: uniprot, flips, cost, route, proposal, band. `lig`: uniprot, tdl, cutoff."""
    bands, lg = summary["dark_gap"]["bands"], summary["ligand_gap"]
    out = []
    for r in dark.sort_values("uniprot").itertuples(index=False):
        n = int(r.flips) + 1
        out.append({"type": DARK, "subject": f"UniProtKB:{r.uniprot}", "kill_condition": KILL, "level": "Tdark",
                    "reason": (f"Tdark: it meets {n} of the three dark criteria (PubMed score under 5, at most 3 GeneRIFs, at most 50 antibodies) "
                               f"and has no experimental Gene Ontology function or process leaf annotation"),
                    "proposal": f"{r.proposal}; or {OTHER_ROUTE}",
                    "criteria_to_clear": int(r.flips), "route_cost": round(float(r.cost), 6), "band": r.band,
                    "closed_in_backtest": round(float(bands[r.band]["share"]), 4), "backtest": spans})
    rate = round(lg["gained"] / lg["candidates"], 4)
    for r in lig.sort_values("uniprot").itertuples(index=False):
        out.append({"type": LIGAND, "subject": f"UniProtKB:{r.uniprot}", "kill_condition": KILL, "level": r.tdl,
                    "reason": "no compound activity in this release passes the potency cutoff for the target's family",
                    "proposal": f"a compound active against this target at {r.cutoff} or better",
                    "potency_cutoff": r.cutoff, "closed_in_backtest": rate, "backtest": spans})
    return out


def export(gap_rows: list[dict], names: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    def resolver(label):
        return label, "biolink:Protein", True, names.get(label, "")

    def extras(row, kind):
        return {f"{PREFIX}:{c.split(':', 1)[1]}": row.get(c.split(":", 1)[1]) for c in EXTRA} if kind == "gap" else {}

    ex = Exporter(mapping=mapping(), resolver=resolver, extra_gap_columns=EXTRA, extras=extras)
    for row in gap_rows:
        if ex.add(row) != "gap":
            raise SystemExit(f"REFUSED: a target-level gap for {row['subject']} was written as an association")
    n, e, g = ex.frames()
    return n, e, g, ex.merged


def load_rows(release: str) -> tuple[list[dict], dict, dict, Path]:
    """One row per target-level gap, the target names, the summary the rows were read with, and its file."""
    C, summary_file = paths.context(release), paths.DOCS_DATA / f"target_gaps_{release}.json"
    if not summary_file.is_file():
        raise SystemExit(f"no target-level gaps; run:  python -m pharos_gap_kg.tgaps {release}")
    summary = json.loads(summary_file.read_text())
    for f, want in summary["files"].items():
        if graph.file_sha(C / f) != want:
            raise SystemExit(f"REFUSED: {C / f} is not the table {summary_file.name} describes; run the target-level gaps again")
    dark, lig = pd.read_parquet(C / "target_gaps_dark.parquet"), pd.read_parquet(C / "target_gaps_ligand.parquet")
    names = {f"UniProtKB:{u}": (s if isinstance(s, str) and s else u) for u, s in zip(lig.uniprot, lig.sym)}
    names.update({f"UniProtKB:{u}": (s if isinstance(s, str) and s else u) for u, s in zip(dark.uniprot, dark.sym)})
    spans = f"{release_name(release)} to Pharos 4.0.0 (May 2026)"
    return rows(dark, lig, summary, spans), names, summary, summary_file


def build(release: str) -> dict:
    gap_rows, names, summary, summary_file = load_rows(release)
    nodes, edges, gaps, merged = export(gap_rows, names)
    rep = conformance.check(nodes, edges, gaps, mapping(), n_input_rows=len(gap_rows), n_merged=merged)
    if not rep.passed or rep.n_skipped:
        raise SystemExit(f"REFUSED: the target-gap graph does not pass conformance\n{rep}\nNothing was written.")
    out = paths.GRAPH / f"{release}_target_gaps"
    files = graph.write(out, nodes, edges, gaps)
    man = {"release": release, "summary_sha256": graph.file_sha(summary_file), "tables_sha256": summary["files"],
           "rows": len(gap_rows), "edges": len(edges), "gaps": len(gaps), "merged": merged, "nodes": len(nodes),
           "targets": int(gaps.subject.nunique()), "targets_with_both_gaps": int((gaps.groupby("subject").size() == 2).sum()),
           "gaps_by_type": {k: int(v) for k, v in gaps[f"{PREFIX}:gap_type"].value_counts().items()},
           "gap_predicate": mapping().gap_predicate, "knowledge_source": SOURCE,
           "conformance": "10 of 10 checks passed", "files": {k: graph.file_sha(p) for k, p in files.items()}}
    (out / "graph.json").write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"graph_{release}_target_gaps.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    dark = pd.DataFrame({"uniprot": ["P2", "P1"], "sym": ["BBB", "AAA"], "flips": [2, 1], "cost": [1.2345678, 0.2], "route": ["pubmed_score|generifs", "generifs"],
                         "proposal": ["PubMed score up by 3.0 (now 2.0, clears at 5); and 4 more GeneRIFs (now 0, clears at 4)", "1 more GeneRIF (now 3, clears at 4)"],
                         "band": ["2 flips, nearest third", "1 flip, nearest third"]})
    lig = pd.DataFrame({"uniprot": ["P1", "P2", "P3"], "sym": ["AAA", "BBB", None], "tdl": ["Tdark", "Tdark", "Tbio"], "cutoff": ["30 nM", "1 uM", "1 uM"]})
    summary = {"dark_gap": {"bands": {"1 flip, nearest third": {"share": 0.5361}, "2 flips, nearest third": {"share": 0.11412}}}, "ligand_gap": {"gained": 20, "candidates": 200}}
    R = rows(dark, lig, summary, "Tiny 1 to Tiny 2")
    ck("every gap is one row: a dark gap for each dark target and a ligand gap for each target without a compound", [r["type"] for r in R] == [DARK, DARK, LIGAND, LIGAND, LIGAND])
    n, e, g, merged = export(R, {"UniProtKB:P1": "AAA", "UniProtKB:P2": "BBB"})
    rep = conformance.check(n, e, g, mapping(), n_input_rows=len(R), n_merged=merged)
    ck("the export is gaps only and passes conformance with nothing skipped", len(e) == 0 and len(g) == 5 and merged == 0 and rep.passed and not rep.n_skipped)
    ck("a gap runs from the target to the kind of gap, under the gap predicate and Pharos as its source",
       set(g.predicate) == {f"{PREFIX}:evidence_missing_for"} and set(g.object) == {f"{PREFIX}:GAP:{DARK}", f"{PREFIX}:GAP:{LIGAND}"} and set(g.primary_knowledge_source) == {SOURCE})
    a = g[(g.subject == "UniProtKB:P1") & (g[f"{PREFIX}:gap_type"] == DARK)].iloc[0]
    ck("a dark gap says how many criteria it meets, what closes it either way, and what would retire it",
       "meets 2 of the three" in a[f"{PREFIX}:gap_reason"] and a[f"{PREFIX}:proposal"] == f"1 more GeneRIF (now 3, clears at 4); or {OTHER_ROUTE}" and a[f"{PREFIX}:kill_condition"] == KILL)
    ck("it carries its flips, cost, band and the share of the band that closed, and the releases that share spans",
       (a[f"{PREFIX}:criteria_to_clear"], a[f"{PREFIX}:route_cost"], a[f"{PREFIX}:band"], a[f"{PREFIX}:closed_in_backtest"], a[f"{PREFIX}:backtest"]) == (1, 0.2, "1 flip, nearest third", 0.5361, "Tiny 1 to Tiny 2"))
    b = g[(g.subject == "UniProtKB:P3")].iloc[0]
    ck("a ligand gap names the cutoff it lacks and has no flips or cost", b[f"{PREFIX}:potency_cutoff"] == "1 uM" and b[f"{PREFIX}:criteria_to_clear"] == "" and b[f"{PREFIX}:route_cost"] == ""
       and b[f"{PREFIX}:closed_in_backtest"] == 0.1 and "1 uM or better" in b[f"{PREFIX}:proposal"])
    ck("the two fields that do not apply are left empty", set(g[f"{PREFIX}:detection_floor"]) == {""} and set(g[f"{PREFIX}:n_required"]) == {""} and set(g[f"{PREFIX}:withheld_from"]) == {""})
    ck("one target's two gaps have two ids, and an id does not move when the gap's numbers do",
       g.id.is_unique and g.id.str.startswith(f"{PREFIX}:g").all()
       and export(rows(dark.assign(cost=[9.0, 9.0]), lig, summary, "other"), {})[2].id.tolist() == g.id.tolist())
    ck("the gap nodes are named, and a target is named by its symbol or its accession", n.set_index("id").name[f"{PREFIX}:GAP:{DARK}"] == GAP_NAMES[DARK]
       and n.set_index("id").name["UniProtKB:P1"] == "AAA" and n.set_index("id").name["UniProtKB:P3"] == "UniProtKB:P3")
    ck("a release is named as a person would write it", release_name("pharos319") == "Pharos 3.19" and release_name("tiny") == "tiny")
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 1:
        raise SystemExit(__doc__)
    man = build(argv[0])
    print(f"{man['release']} target-level gaps: {man['rows']:,} rows -> {man['edges']:,} edges + {man['gaps']:,} gaps, {man['nodes']:,} nodes; {man['conformance']}")
    print(f"  by type: {man['gaps_by_type']}; {man['targets']:,} targets, {man['targets_with_both_gaps']:,} with both")
    return 0


if __name__ == "__main__":
    sys.exit(main())
