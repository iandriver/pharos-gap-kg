"""Everything the question-answering front end reads, built once into one directory with a manifest.

The front end answers from the bundle and from nothing else, so what a user is told can be traced to
these files and their digests:

    nodes.parquet        every Pharos target and every GWAS trait, with what names and grades them
    links.parquet        one row per target-trait link: its edge id, first and last year, study counts
    evidence.parquet     one row per (link, study): the study, its paper and the paper's year
    lineage.parquet      pairs of traits that are one term, or broader and narrower
    neighbours.parquet   pairs of targets within 500,000 bases, and the bases between them
    backtest.json        how the ranking the front end offers has held on the frozen split, when measured
    target_gaps.*        what each target lacks for the next development level, and how often targets in
                         each band were promoted in the back-test (where `tgaps` has been run)
    manifest.json        the inputs' digests, each file's digest, and the counts

The links are every link in the release, whatever its year. "As of 2018" is a filter the front end
applies; the bundle does not choose a date. Edge ids come from the same exporter as the KGX graph,
and where a dated graph has been built its edges must be among these under the same ids.

Run:  python -m pharos_gap_kg.bundle pharos319
      python -m pharos_gap_kg.bundle --selftest
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from kgx_gaps import conformance

from pharos_gap_kg import context, graph, ontology, paths, position, pubyear, registry

FILES = ["nodes.parquet", "links.parquet", "evidence.parquet", "lineage.parquet", "neighbours.parquet"]
GAP_FILES = ["target_gaps.parquet", "target_gaps.json"]


def gap_files(release: str, out: Path) -> list[str]:
    """Write what the `gaps` question reads, where the target-level gaps have been computed
    (`python -m pharos_gap_kg.tgaps <release>`), and return the files written. One table of what each
    target lacks, and the rates at which targets in each band were promoted in the back-test."""
    C, summary = paths.context(release), paths.DOCS_DATA / f"target_gaps_{release}.json"
    dark, lig = C / "target_gaps_dark.parquet", C / "target_gaps_ligand.parquet"
    if not (dark.is_file() and lig.is_file() and summary.is_file()):
        return []
    s = json.loads(summary.read_text())
    for f, want in s["files"].items():
        if graph.file_sha(C / f) != want:
            raise SystemExit(f"REFUSED: {C / f} is not the table {summary.name} describes; run the target-level gaps again")
    d = pd.read_parquet(dark)[["uniprot", "pubmed_score", "generifs", "antibodies", "flips", "cost", "route", "proposal", "band"]]
    g = pd.read_parquet(lig)[["uniprot", "tdl", "fam", "cutoff", "partners_with_ligand"]].assign(lacks_ligand=True)
    t = g.merge(d, on="uniprot", how="outer")
    t.insert(0, "id", "UniProtKB:" + t.uniprot)
    t = t.drop(columns="uniprot").rename(columns={"tdl": "level", "fam": "family"}).sort_values("id").reset_index(drop=True)
    t["flips"] = t.flips.astype("Int64")
    t.to_parquet(out / GAP_FILES[0], index=False)
    dg, lg = s["dark_gap"], s["ligand_gap"]
    m = re.fullmatch(r"pharos(\d)(\d+)", release)
    rates = {"earlier": f"Pharos {m.group(1)}.{m.group(2)}" if m else release, "later": "Pharos 4.0.0 (May 2026)", "summary_sha256": graph.file_sha(summary),
             "dark": {"targets": dg["with_a_gap"], "in_both_releases": dg["in_both_releases"], "promoted": dg["promoted"],
                      "bands": {b: {"promoted": v, "fell_in_the_mirror": dg["mirror"]["bands"][b]} for b, v in dg["bands"].items()},
                      "route_cost_auc": dg["scores"]["route cost"]["auc"], "best_single_count": dg["best_single_count"],
                      "best_single_count_auc": dg["scores"][dg["best_single_count"]]["auc"],
                      "named_criterion_cleared": dg["proposal"]["named_criterion_cleared"], "by_naming_at_random": dg["proposal"]["by_naming_at_random"]},
             "ligand": {"candidates": lg["candidates"], "gained": lg["gained"]}}
    (out / GAP_FILES[1]).write_text(json.dumps(rates, indent=1) + "\n")
    return list(GAP_FILES)


def add_gaps(release: str) -> dict:
    """Add the target-level gaps to a bundle that is already built, and record them in its manifest."""
    out = paths.bundle(release)
    mf = out / "manifest.json"
    if not mf.is_file():
        raise SystemExit(f"no bundle at {out}; build one with:  python -m pharos_gap_kg.bundle {release}")
    files = gap_files(release, out)
    if not files:
        raise SystemExit(f"no target-level gaps to add; run:  python -m pharos_gap_kg.tgaps {release}")
    man = json.loads(mf.read_text())
    man["files"].update({f: graph.file_sha(out / f) for f in files})
    man["target_gaps"] = int(len(pd.read_parquet(out / GAP_FILES[0])))
    mf.write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"bundle_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def lineage_table(traits: list[str], terms: dict) -> pd.DataFrame:
    """a, b, relation: for each trait a, each other trait b in its lineage, and what b is to a."""
    lin = ontology.lineage(traits, terms)
    canon = {t: ontology.canonical(t, terms) for t in traits}
    anc = {t: (ontology.ancestors(c, terms) if c else set()) for t, c in canon.items()}
    rows = []
    for a in traits:
        for b in sorted(lin[a]):
            rel = "same term" if canon[a] == canon[b] else ("broader" if canon[b] in anc[a] else "narrower")
            rows.append((a, b, rel))
    return pd.DataFrame(rows, columns=["a", "b", "relation"])


def node_table(release: str, known: dict, terms: dict, where: pd.DataFrame | None = None) -> pd.DataFrame:
    d = paths.tables(release)
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot", "sym", "description", "geneid"]).to_pandas()
    tgt = pq.read_table(d / "target.parquet", columns=["id", "tdl", "fam"]).to_pandas().rename(columns={"id": "target_id"})
    link = pq.read_table(d / "t2tc.parquet", columns=["target_id", "protein_id"]).to_pandas()
    p = prot.merge(link, left_on="id", right_on="protein_id").merge(tgt, on="target_id")
    p["node"] = "UniProtKB:" + p.uniprot
    p = p.set_index("node")
    rows = []
    at = {} if where is None else where.sort_values(["target", "start"]).drop_duplicates("target").set_index("target")[["chr", "start", "end"]].to_dict("index")
    traits = [c for c, v in known.items() if v[1] != "biolink:Protein"]
    canon, dis = context.kinds(traits, terms)
    live = dict(zip(traits, canon))
    disease = dict(zip(traits, dis.tolist()))
    for c, (curie, cat, _, name) in known.items():
        if cat == "biolink:Protein":
            r = p.loc[c]
            rows.append({"id": c, "kind": "target", "category": cat, "name": name, "full_name": r.description or "",
                         "level": r.tdl or "", "family": r.fam if isinstance(r.fam, str) else "",
                         "xref": f"NCBIGene:{int(r.geneid)}" if pd.notna(r.geneid) else "", "live_term": "", "disease_like": False,
                         "chr": at.get(c, {}).get("chr", ""), "start": int(at.get(c, {}).get("start", 0)), "end": int(at.get(c, {}).get("end", 0))})
        else:
            rows.append({"id": c, "kind": "trait", "category": cat, "name": name, "full_name": "", "level": "", "family": "",
                         "xref": "", "live_term": live[c] or "", "disease_like": bool(disease[c]), "chr": "", "start": 0, "end": 0})
    return pd.DataFrame(rows)


def check_against_dated_graphs(release: str, links: pd.DataFrame) -> list[str]:
    """Every dated graph already built for this release must have its edges here, under the same ids."""
    said = []
    have = links.set_index("edge_id")[["target", "trait"]]
    for g in sorted(paths.GRAPH.glob(f"{release}_asof*")):
        e = pd.read_csv(g / "edges.tsv", sep="\t", dtype=str, keep_default_na=False, usecols=["id", "subject", "object"])
        j = e.join(have, on="id")
        bad = j[(j.target != j.subject) | (j.trait != j.object)]
        if len(bad):
            raise SystemExit(f"REFUSED: {len(bad):,} edges of {g.name} are missing from the bundle or sit under another id")
        said.append(f"{g.name}: {len(e):,} edges found under the same ids")
    return said


def backtest() -> dict | None:
    """What the front end may say about how well its ranking has held: the headline scores of the
    frozen backtest and the calibration by rank. Read from the committed summaries, newest cutoff."""
    cal = sorted(paths.DOCS_DATA.glob("calibration_gene_trait_*.json"))
    if not cal:
        return None
    c = json.loads(cal[-1].read_text())
    b = paths.DOCS_DATA / f"baselines_gene_trait_{c['cutoff']}.json"
    k = paths.DOCS_DATA / f"controls_gene_trait_{c['cutoff']}.json"
    if not (b.is_file() and k.is_file()):
        return None
    base, ctl = json.loads(b.read_text()), json.loads(k.read_text())
    if base["split"] != c["split"] or ctl["split"] != c["split"]:
        raise SystemExit("REFUSED: the baseline, control and calibration summaries are not from one split")
    pair = ctl["paired"]["B3-B1"]
    head = {"cutoff": c["cutoff"], "data_ends": c["data_ends"], "measure": "mean average precision over traits",
            "chance": base["results"]["B0"]["primary"][0], "degree": base["results"]["B1"]["primary"][0],
            "neighbours": base["results"]["B3"]["primary"][0], "positives": base["positives"], "traits_scored": base["traits_scored"],
            "share_of_later_links_explained_by_wording_or_position": ctl["shares"]["either"] / ctl["shares"]["positives"],
            "neighbours_minus_degree_all": [pair["full"]["difference"], pair["full"]["low"], pair["full"]["high"]],
            "neighbours_minus_degree_unexplained_lower_bound": [pair["low"]["difference"], pair["low"]["low"], pair["low"]["high"]],
            "neighbours_minus_degree_unexplained_neutral": [pair["neutral"]["difference"], pair["neutral"]["low"], pair["neutral"]["high"]],
            "split": c["split"]}
    return {"headline": head, "calibration": {k_: c[k_] for k_ in ("cutoff", "data_ends", "base_rate", "by_trait", "by_target")}}


def build(release: str) -> dict:
    years = pubyear.load().set_index("pmid").year
    rows = graph.evidence_rows(release, years)
    latest = int(rows.year.max())
    ev = graph.evidence(release, years, latest)
    known, xref = graph.nodes_of(release)
    ex = graph.export(ev, known, xref)
    n, e, g = graph.frames(ex)
    rep = conformance.check(n, e, g, ex.mapping, n_input_rows=len(ev), n_merged=ex.merged)
    if not rep.passed or rep.n_skipped:
        raise SystemExit(f"REFUSED: the full graph does not pass conformance\n{rep}")
    grp = rows.groupby(["subject", "object"], sort=True)
    stats = pd.DataFrame({"first_year": grp.year.min(), "last_year": grp.year.max(),
                          "n_studies": grp.study_acc.nunique(), "n_papers": grp.pubmedid.nunique()}).reset_index()
    links = e[["id", "subject", "object"]].merge(stats, on=["subject", "object"], how="left")
    links = links.rename(columns={"id": "edge_id", "subject": "target", "object": "trait"})
    if links.first_year.isna().any() or len(links) != len(stats):
        raise SystemExit("REFUSED: the exported edges and the dated links are not the same set")
    for c in ("first_year", "last_year", "n_studies", "n_papers"):
        links[c] = links[c].astype(int)
    terms = ontology.read(registry.require("efo_obo"))
    genes = position.genes(registry.require("ensembl_gtf_113"))
    where = position.target_genes(release, genes)
    nodes = node_table(release, known, terms, where)
    traits = nodes[nodes.kind == "trait"].id.tolist()
    targets = nodes[nodes.kind == "target"].id.tolist()
    lin = lineage_table(traits, terms)
    nb = position.neighbour_table(where, targets)
    evidence = rows.rename(columns={"subject": "target", "object": "trait", "study_acc": "study", "pubmedid": "pmid"})
    out = paths.bundle(release)
    out.mkdir(parents=True, exist_ok=True)
    for name, df in (("nodes", nodes), ("links", links), ("evidence", evidence), ("lineage", lin), ("neighbours", nb)):
        df.to_parquet(out / f"{name}.parquet", index=False)
    files = list(FILES)
    bt = backtest()
    if bt is not None:
        (out / "backtest.json").write_text(json.dumps(bt, indent=1) + "\n")
        files.append("backtest.json")
    elif (out / "backtest.json").exists():
        (out / "backtest.json").unlink()
    files += gap_files(release, out)
    side = lambda k: json.loads(registry.SOURCES[k].path.with_name(registry.SOURCES[k].filename + ".lineage.json").read_text())["sha256"]
    src = registry.SOURCES[release]
    man = {
        "release": release, "source_file": src.filename, "source_url": src.url, "source_published": src.published,
        "source_sha256": json.loads((paths.tables(release) / "_manifest.json").read_text())["sha256"],
        "source_tables": graph.TABLES, "ontology_sha256": side("efo_obo"), "coordinates_sha256": side("ensembl_gtf_113"),
        "window": position.WINDOW, "first_year": int(links.first_year.min()), "latest_year": latest,
        "targets": len(targets), "targets_with_a_link": int(links.target.nunique()), "traits": len(traits),
        "links": len(links), "studies": int(evidence.study.nunique()), "papers": int(evidence.pmid.nunique()),
        "lineage_pairs": int(len(lin) // 2), "neighbour_pairs": int(len(nb) // 2),
        "conformance": "10 of 10 checks passed" if rep.passed else "failed",
        "dated_graphs": check_against_dated_graphs(release, links),
        "backtest": (f"split frozen at {bt['headline']['cutoff']}" if bt else "none"),
        "files": {f: graph.file_sha(out / f) for f in files},
    }
    if GAP_FILES[0] in files:
        man["target_gaps"] = int(len(pd.read_parquet(out / GAP_FILES[0])))
    (out / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"bundle_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    terms = ontology.load("""[Term]
id: MONDO:1
name: disease

[Term]
id: MONDO:2
name: heart disease
is_a: MONDO:1

[Term]
id: efo:EFO_0000009
name: obsolete_heart disease
is_obsolete: true
replaced_by: http://purl.obolibrary.org/obo/MONDO_2

[Term]
id: MONDO:4
name: lung disease
is_a: MONDO:1
""".splitlines())
    L = lineage_table(["MONDO:1", "MONDO:2", "EFO:0000009", "MONDO:4"], terms)
    rel = {(r.a, r.b): r.relation for r in L.itertuples()}
    ck("a parent is broader and a child narrower, seen from each side", rel[("MONDO:2", "MONDO:1")] == "broader" and rel[("MONDO:1", "MONDO:2")] == "narrower")
    ck("a retired id and its replacement are the same term", rel[("EFO:0000009", "MONDO:2")] == "same term" and rel[("MONDO:2", "EFO:0000009")] == "same term")
    ck("a retired id inherits its replacement's parents", rel[("EFO:0000009", "MONDO:1")] == "broader")
    ck("siblings are not listed, and every pair appears from both sides", ("MONDO:2", "MONDO:4") not in rel and all((b, a) in rel for a, b in rel))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) == 2 and argv[1] == "--gaps-only":
        man = add_gaps(argv[0])
        print(f"{man['release']}: {man['target_gaps']:,} target-level gaps added to {paths.bundle(man['release'])}")
        return 0
    if len(argv) != 1:
        raise SystemExit(__doc__)
    man = build(argv[0])
    print(f"{man['release']}: {man['targets']:,} targets ({man['targets_with_a_link']:,} with a link), {man['traits']:,} traits, "
          f"{man['links']:,} links from {man['studies']:,} studies in {man['papers']:,} papers, {man['first_year']} to {man['latest_year']}")
    print(f"  lineage pairs {man['lineage_pairs']:,}; neighbour pairs {man['neighbour_pairs']:,}; {man['conformance']}")
    for s in man["dated_graphs"]:
        print(f"  {s}")
    print(f"  -> {paths.bundle(man['release'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
