"""The graph as it could have been known at a cutoff year, as KGX nodes, edges and gaps.

WHAT IS IN IT. Every Pharos target is a node, linked or not: a target nothing is known about is the
case the graph exists for. Every GWAS trait with a link reported by the cutoff is a node. An edge is
one target-trait link first reported in or before the cutoff, and it carries only what was known
then: the studies and papers dated up to the cutoff. A study published later is not on the edge, even
though the database lists it, because a model shown the 2021 replication of a 2017 link has been
shown the future.

THE RULES ARE NOT RESTATED HERE. Rows go through the `kgx-gaps` exporter, which derives each id from
what the edge claims, and the result must pass its conformance suite or nothing is written. The edge
set must also equal the frozen training part of the split, link for link.

LINEAGE. `lineage` answers "why does A connect to B" for one edge: the claim, the studies behind it
with their papers and years, the tables and the published file they were read from with its digest,
and the rule that admitted it.

Run:  python -m pharos_gap_kg.graph build pharos319 --cutoff 2018
      python -m pharos_gap_kg.graph lineage pharos319 --cutoff 2018 --target APOE --trait MONDO:0004975
      python -m pharos_gap_kg.graph --selftest
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from kgx_gaps import Mapping, Rule, conformance
from kgx_gaps.export import Exporter

from pharos_gap_kg import paths, pubyear, registry, split

PREFIX = "pharosgapkg"
SOURCE = "infores:tiga"          # Target Illumination GWAS Analytics, which aggregated the GWAS Catalog
AGGREGATOR = "infores:pharos"
TABLES = ["tiga", "tiga_provenance", "protein"]

# What a trait id's prefix says the node is. An id under a prefix not listed here is refused: a
# category guessed from an unknown vocabulary would be a claim nobody checked.
TRAIT_CATEGORY = {
    "MONDO": "biolink:Disease", "Orphanet": "biolink:Disease",
    "HP": "biolink:PhenotypicFeature",
    "EFO": "biolink:DiseaseOrPhenotypicFeature", "NCIT": "biolink:DiseaseOrPhenotypicFeature",
    "GO": "biolink:BiologicalProcess", "PATO": "biolink:PhenotypicQuality",
}
CONDITION = {"biolink:Disease", "biolink:PhenotypicFeature", "biolink:DiseaseOrPhenotypicFeature"}
EXTRA_EDGE = ["aggregator_knowledge_source", "publications", f"{PREFIX}:studies",
              f"{PREFIX}:first_year", f"{PREFIX}:n_studies"]


def repair(name: str) -> str:
    """Undo text that was UTF-8, read as Latin-1 and saved again, once or more ("Ã¢Â\x80Â\x99" for an
    apostrophe). Each round is applied only if it decodes cleanly, so text that was never damaged is
    returned as it came. A no-break space becomes a space."""
    s = name
    for _ in range(3):
        try:
            r = s.encode("latin-1").decode("utf-8")
        except (UnicodeEncodeError, UnicodeDecodeError):
            break
        if r == s:
            break
        s = r
    return " ".join(s.replace("\xa0", " ").split())


def trait_category(curie: str) -> str:
    prefix = curie.split(":", 1)[0]
    if prefix not in TRAIT_CATEGORY:
        raise SystemExit(f"REFUSED: trait {curie} is under the prefix {prefix!r}, which has no declared category")
    return TRAIT_CATEGORY[prefix]


def mapping() -> Mapping:
    m = Mapping(prefix=PREFIX, knowledge_source=SOURCE)
    common = dict(predicate="biolink:genetically_associated_with", knowledge_level="statistical_association",
                  agent_type="data_analysis_pipeline")
    # Biolink 4.4.3 has no gene-to-disease-or-phenotype class, and most traits are DiseaseOrPhenotypicFeature;
    # the two types are kept apart because the type is part of an edge's id (docs/PROTOCOL.md, 19 Oct).
    m.rules["gwas_condition"] = Rule(category="biolink:Association", **common)
    m.rules["gwas_other"] = Rule(category="biolink:Association", **common)   # a process or a quality, not a condition
    problems = m.validate()
    if problems:
        raise SystemExit("REFUSED mapping: " + "; ".join(problems))
    return m


def evidence_rows(release: str, years: pd.Series) -> pd.DataFrame:
    """One row per (target, trait, study): subject, object, study_acc, pubmedid, year. Every year."""
    d = paths.tables(release)
    t = pq.read_table(d / "tiga.parquet", columns=["protein_id", "ensg", "efoid"]).to_pandas().drop_duplicates()
    p = pq.read_table(d / "tiga_provenance.parquet", columns=["ensg", "efoid", "study_acc", "pubmedid"]).to_pandas().drop_duplicates()
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot"]).to_pandas().rename(columns={"id": "protein_id"})
    p["year"] = p.pubmedid.map(years)
    if p.year.isna().any():
        raise SystemExit(f"REFUSED: {int(p.year.isna().sum()):,} cited studies have no publication year; "
                         f"run:  python -m pharos_gap_kg.pubyear {release}")
    m = t.merge(p, on=["ensg", "efoid"]).merge(prot, on="protein_id")
    m["subject"] = "UniProtKB:" + m.uniprot
    m["object"] = m.efoid.map(split.curie)
    m["year"] = m.year.astype(int)
    return m.drop_duplicates(["subject", "object", "study_acc", "pubmedid"])[["subject", "object", "study_acc", "pubmedid", "year"]]


def evidence(release: str, years: pd.Series, cutoff: int) -> pd.DataFrame:
    """One row per link known by the cutoff, with the studies and papers dated up to it."""
    m = evidence_rows(release, years)
    m = m[m.year <= cutoff]
    g = m.sort_values(["year", "study_acc"]).groupby(["subject", "object"], sort=True)
    out = pd.DataFrame({
        "first_year": g.year.min().astype(int),
        "n_studies": g.study_acc.nunique(),
        "studies": g.study_acc.agg(lambda s: "|".join(dict.fromkeys(s))),
        "publications": g.pubmedid.agg(lambda s: "|".join(f"PMID:{i}" for i in dict.fromkeys(s))),
    }).reset_index()
    out["type"] = ["gwas_condition" if trait_category(o) in CONDITION else "gwas_other" for o in out.object]
    return out


def nodes_of(release: str) -> tuple[dict, dict]:
    """label -> (curie, category, grounded, name) for every target and trait, and target xrefs."""
    d = paths.tables(release)
    prot = pq.read_table(d / "protein.parquet", columns=["uniprot", "sym", "geneid", "description"]).to_pandas()
    tg = {}
    xref = {}
    for r in prot.itertuples(index=False):
        c = f"UniProtKB:{r.uniprot}"
        tg[c] = (c, "biolink:Protein", True, r.sym if isinstance(r.sym, str) and r.sym else r.uniprot)
        xref[c] = f"NCBIGene:{int(r.geneid)}" if pd.notna(r.geneid) else ""
    t = pq.read_table(d / "tiga.parquet", columns=["efoid", "trait"]).to_pandas().sort_values("trait").drop_duplicates("efoid")
    for r in t.itertuples(index=False):
        c = split.curie(r.efoid)
        tg[c] = (c, trait_category(c), True, repair(r.trait))
    return tg, xref


def export(ev: pd.DataFrame, known: dict, xref: dict, all_targets: bool = True) -> Exporter:
    def resolver(label):
        if label not in known:
            raise SystemExit(f"REFUSED: {label} is on an edge and is not a Pharos target or a GWAS trait")
        return known[label]

    def extras(row, kind):
        if kind != "edge":
            return {}
        return {"aggregator_knowledge_source": AGGREGATOR, "publications": row["publications"],
                f"{PREFIX}:studies": row["studies"], f"{PREFIX}:first_year": row["first_year"],
                f"{PREFIX}:n_studies": row["n_studies"]}

    ex = Exporter(mapping=mapping(), resolver=resolver, extra_edge_columns=EXTRA_EDGE, extras=extras)
    if all_targets:
        for c, (curie, cat, grounded, name) in known.items():
            if cat == "biolink:Protein":
                ex.node(curie, cat, name, grounded, xref.get(c, ""))
    for row in ev.to_dict("records"):
        ex.add(row)
    return ex


INT_COLUMNS = [f"{PREFIX}:first_year", f"{PREFIX}:n_studies"]


def frames(ex: Exporter) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """The exporter's frames, with whole-number columns written as whole numbers. The exporter turns
    any value that parses as a number into a float, so a year arrives as 2018.0."""
    n, e, g = ex.frames()
    for c in INT_COLUMNS:
        e[c] = [("" if v == "" else str(int(float(v)))) for v in e[c]]
    return n, e, g


def write(out: Path, n: pd.DataFrame, e: pd.DataFrame, g: pd.DataFrame) -> dict[str, Path]:
    out.mkdir(parents=True, exist_ok=True)
    files = {"nodes": out / "nodes.tsv", "edges": out / "edges.tsv", "gaps": out / "gaps.tsv"}
    for k, df in (("nodes", n), ("edges", e), ("gaps", g)):
        df.to_csv(files[k], sep="\t", index=False, na_rep="")
    return files


def file_sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        while buf := f.read(1 << 20):
            h.update(buf)
    return h.hexdigest()


def graph_dir(release: str, cutoff: int) -> Path:
    return paths.GRAPH / f"{release}_asof{cutoff}"


def build(release: str, cutoff: int) -> dict:
    years = pubyear.load().set_index("pmid").year
    ev = evidence(release, years, cutoff)
    train = split.load(cutoff)["train"]
    mine = hashlib.sha256(("\n".join(sorted(ev.subject + "\t" + ev.object)) + "\n").encode()).hexdigest()
    if mine != split.digest(train):
        raise SystemExit(f"REFUSED: the {len(ev):,} links this graph would assert are not the {len(train):,} "
                         f"frozen training links. Nothing was written.")
    known, xref = nodes_of(release)
    raw = pq.read_table(paths.tables(release) / "tiga.parquet", columns=["trait"]).to_pandas().trait.drop_duplicates()
    names_repaired = int((raw.map(repair) != raw).sum())
    ex = export(ev, known, xref)
    n, e, g = frames(ex)
    rep = conformance.check(n, e, g, ex.mapping, n_input_rows=len(ev), n_merged=ex.merged)
    if not rep.passed or rep.n_skipped:
        raise SystemExit(f"REFUSED: the graph does not pass conformance\n{rep}\nNothing was written.")
    out = graph_dir(release, cutoff)
    files = write(out, n, e, g)
    man = {
        "release": release, "cutoff": cutoff,
        "source_sha256": json.loads((paths.tables(release) / "_manifest.json").read_text())["sha256"],
        "source_tables": TABLES, "split_train_sha256": split.digest(train),
        "nodes": len(n), "edges": len(e), "gaps": len(g),
        "node_categories": {k: int(v) for k, v in n.category.value_counts().items()},
        "edge_categories": {k: int(v) for k, v in e.category.value_counts().items()},
        "predicates": {k: int(v) for k, v in e.predicate.value_counts().items()},
        "names_repaired": names_repaired,
        "targets_with_no_edge": int((~n[n.category == "biolink:Protein"].id.isin(set(e.subject))).sum()),
        "edges_by_first_year": {int(k): int(v) for k, v in e[f"{PREFIX}:first_year"].astype(int).value_counts().sort_index().items()},
        "conformance": [str(c) for c in rep.checks] if hasattr(rep, "checks") else str(rep).splitlines(),
        "files": {k: file_sha(p) for k, p in files.items()},
    }
    (out / "graph.json").write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"graph_{release}_asof{cutoff}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def lineage(release: str, cutoff: int, edge_id: str | None = None, target: str | None = None, trait: str | None = None) -> dict:
    out = graph_dir(release, cutoff)
    if not (out / "edges.tsv").is_file():
        raise SystemExit(f"no graph yet; build it with:  python -m pharos_gap_kg.graph build {release} --cutoff {cutoff}")
    e = pd.read_csv(out / "edges.tsv", sep="\t", dtype=str, keep_default_na=False)
    n = pd.read_csv(out / "nodes.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("id")
    if edge_id:
        hit = e[e.id == edge_id]
    else:
        subj = set(n.index[(n.name == target) | (n.index == target)])
        hit = e[e.subject.isin(subj) & (e.object == trait)]
    if hit.empty:
        asked = edge_id or f"{target} and {trait}"
        return {"found": False, "asked": asked,
                "answer": f"no edge for {asked} in the graph as of {cutoff}. That is not evidence of no association: "
                          f"it means no study dated {cutoff} or earlier in this source reported one."}
    r = hit.iloc[0]
    years = pubyear.load().set_index("pmid").year
    pmids = [int(x.split(":")[1]) for x in r.publications.split("|") if x]
    man = json.loads((out / "graph.json").read_text())
    src = registry.SOURCES[release]
    return {
        "found": True, "edge": r.id,
        "claim": {"subject": r.subject, "subject_name": n.name.get(r.subject, ""), "predicate": r.predicate,
                  "object": r.object, "object_name": n.name.get(r.object, ""), "category": r.category},
        "why": {"first_reported": int(r[f"{PREFIX}:first_year"]), "studies_by_cutoff": r[f"{PREFIX}:studies"].split("|"),
                "papers": [{"pmid": i, "year": int(years[i])} for i in pmids]},
        "knowledge": {"primary_source": r.primary_knowledge_source, "aggregator": r.aggregator_knowledge_source,
                      "level": r.knowledge_level, "agent": r.agent_type},
        "read_from": {"tables": man["source_tables"], "file": src.filename, "url": src.url,
                      "published": src.published, "sha256": man["source_sha256"]},
        "admitted_by": f"rule {r[f'{PREFIX}:source_edge_type']!r}: a GWAS association TIGA lists for this target and trait, "
                       f"with at least one study dated {cutoff} or earlier",
        "not_on_this_edge": f"studies dated after {cutoff}, which the database lists and the graph as of {cutoff} leaves out",
    }


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    known = {"UniProtKB:P1": ("UniProtKB:P1", "biolink:Protein", True, "AAA"),
             "UniProtKB:P2": ("UniProtKB:P2", "biolink:Protein", True, "BBB"),
             "UniProtKB:P3": ("UniProtKB:P3", "biolink:Protein", True, "CCC"),
             "MONDO:1": ("MONDO:1", "biolink:Disease", True, "a disease"),
             "GO:2": ("GO:2", "biolink:BiologicalProcess", True, "a process")}
    ev = pd.DataFrame({"subject": ["UniProtKB:P1", "UniProtKB:P2", "UniProtKB:P1"], "object": ["MONDO:1", "MONDO:1", "GO:2"],
                       "type": ["gwas_condition", "gwas_condition", "gwas_other"], "first_year": [2010, 2015, 2012],
                       "n_studies": [2, 1, 1], "studies": ["GCST1|GCST2", "GCST3", "GCST4"],
                       "publications": ["PMID:1|PMID:2", "PMID:3", "PMID:4"]})
    ex = export(ev, known, {"UniProtKB:P1": "NCBIGene:1"})
    n, e, g = frames(ex)
    rep = conformance.check(n, e, g, ex.mapping, n_input_rows=len(ev), n_merged=ex.merged)
    ck("a small graph passes every conformance check with none skipped", rep.passed and not rep.n_skipped)
    ck("a target with no link is still a node", "UniProtKB:P3" in set(n.id) and "UniProtKB:P3" not in set(e.subject))
    ck("a node is named by its symbol and keeps its cross-reference", n.set_index("id").loc["UniProtKB:P1", ["name", "xref"]].tolist() == ["AAA", "NCBIGene:1"])
    ck("a link is a Biolink 4.4.3 Association, and a link to a process keeps a type of its own",
       set(e.category) == {"biolink:Association"} and e.set_index("object").loc["GO:2", f"{PREFIX}:source_edge_type"] == "gwas_other"
       and (e[e.object == "MONDO:1"][f"{PREFIX}:source_edge_type"] == "gwas_condition").all())
    one = e[(e.subject == "UniProtKB:P1") & (e.object == "MONDO:1")].iloc[0]
    ck("an edge carries its studies, papers, first year and aggregator",
       [str(one[c]) for c in ("publications", f"{PREFIX}:studies", f"{PREFIX}:first_year", f"{PREFIX}:n_studies", "aggregator_knowledge_source")]
       == ["PMID:1|PMID:2", "GCST1|GCST2", "2010", "2", AGGREGATOR])
    ck("the source is the primary knowledge source on every edge", set(e.primary_knowledge_source) == {SOURCE})
    ids = set(e.id)
    e2 = frames(export(ev.iloc[::-1], known, {}))[1]
    ck("an edge id does not depend on row order", set(e2.id) == ids and len(ids) == 3)
    e3 = frames(export(ev.assign(n_studies=[9, 9, 9], publications=["PMID:9"] * 3), known, {}))[1]
    ck("an edge id does not depend on the evidence counted for it", set(e3.id) == ids)
    ck("text mis-decoded twice is repaired, and clean text is left alone",
       repair("L\u00c3\u0083\u00c2\u00b6fgren\u00c3\u00a2\u00c2\u0080\u00c2\u0099s") == "L\u00f6fgren\u2019s" and repair("Sj\u00f6gren syndrome") == "Sj\u00f6gren syndrome"
       and repair("x\u00c3\u0082\u00c2\u00a0y") == "x y" and repair("Alzheimer\u00c3\u00a2\u00c2\u0080\u00c2\u0099s") == "Alzheimer\u2019s")
    try:
        trait_category("XYZ:1")
        ck("a trait under an undeclared prefix is refused", False)
    except SystemExit:
        ck("a trait under an undeclared prefix is refused", True)
    try:
        export(pd.concat([ev, ev.iloc[[0]].assign(object="MONDO:404")]), known, {})
        ck("an edge to a node nobody declared is refused", False)
    except SystemExit:
        ck("an edge to a node nobody declared is refused", True)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    l = sub.add_parser("lineage")
    for p in (b, l):
        p.add_argument("release")
        p.add_argument("--cutoff", type=int, required=True)
    l.add_argument("--edge")
    l.add_argument("--target")
    l.add_argument("--trait")
    a = ap.parse_args(argv)
    if a.cmd == "build":
        man = build(a.release, a.cutoff)
        print(f"{a.release} as of {a.cutoff}: {man['nodes']:,} nodes, {man['edges']:,} edges, {man['gaps']:,} gaps")
        print(f"  nodes: {man['node_categories']}")
        print(f"  edges: {man['edge_categories']}")
        print(f"  targets with no edge: {man['targets_with_no_edge']:,}")
        for c in man["conformance"]:
            print(f"  {c}")
        print(f"  -> {graph_dir(a.release, a.cutoff)}")
    else:
        if not a.edge and not (a.target and a.trait):
            raise SystemExit("lineage needs --edge ID, or --target SYMBOL and --trait CURIE")
        print(json.dumps(lineage(a.release, a.cutoff, a.edge, a.target, a.trait), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
