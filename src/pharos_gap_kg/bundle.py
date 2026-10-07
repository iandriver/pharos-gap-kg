"""Everything the question-answering front end reads, built once into one directory with a manifest.

The front end answers from the bundle and from nothing else, so what a user is told can be traced to
six files and their digests:

    nodes.parquet        every Pharos target and every GWAS trait, with what names and grades them
    links.parquet        one row per target-trait link: its edge id, first and last year, study counts
    evidence.parquet     one row per (link, study): the study, its paper and the paper's year
    lineage.parquet      pairs of traits that are one term, or broader and narrower
    neighbours.parquet   pairs of targets within 500,000 bases, and the bases between them
    manifest.json        the inputs' digests, each file's digest, and the counts

The links are every link in the release, whatever its year. "As of 2018" is a filter the front end
applies; the bundle does not choose a date. Edge ids come from the same exporter as the KGX graph,
and where a dated graph has been built its edges must be among these under the same ids.

Run:  python -m pharos_gap_kg.bundle pharos319
      python -m pharos_gap_kg.bundle --selftest
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from kgx_gaps import conformance

from pharos_gap_kg import context, graph, ontology, paths, position, pubyear, registry

FILES = ["nodes.parquet", "links.parquet", "evidence.parquet", "lineage.parquet", "neighbours.parquet"]


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


def node_table(release: str, known: dict, terms: dict) -> pd.DataFrame:
    d = paths.tables(release)
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot", "sym", "description", "geneid"]).to_pandas()
    tgt = pq.read_table(d / "target.parquet", columns=["id", "tdl", "fam"]).to_pandas().rename(columns={"id": "target_id"})
    link = pq.read_table(d / "t2tc.parquet", columns=["target_id", "protein_id"]).to_pandas()
    p = prot.merge(link, left_on="id", right_on="protein_id").merge(tgt, on="target_id")
    p["node"] = "UniProtKB:" + p.uniprot
    p = p.set_index("node")
    rows = []
    traits = [c for c, v in known.items() if v[1] != "biolink:Protein"]
    canon, dis = context.kinds(traits, terms)
    live = dict(zip(traits, canon))
    disease = dict(zip(traits, dis.tolist()))
    for c, (curie, cat, _, name) in known.items():
        if cat == "biolink:Protein":
            r = p.loc[c]
            rows.append({"id": c, "kind": "target", "category": cat, "name": name, "full_name": r.description or "",
                         "level": r.tdl or "", "family": r.fam if isinstance(r.fam, str) else "",
                         "xref": f"NCBIGene:{int(r.geneid)}" if pd.notna(r.geneid) else "", "live_term": "", "disease_like": False})
        else:
            rows.append({"id": c, "kind": "trait", "category": cat, "name": name, "full_name": "", "level": "", "family": "",
                         "xref": "", "live_term": live[c] or "", "disease_like": bool(disease[c])})
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
    nodes = node_table(release, known, terms)
    traits = nodes[nodes.kind == "trait"].id.tolist()
    targets = nodes[nodes.kind == "target"].id.tolist()
    lin = lineage_table(traits, terms)
    genes = position.genes(registry.require("ensembl_gtf_113"))
    nb = position.neighbour_table(position.target_genes(release, genes), targets)
    evidence = rows.rename(columns={"subject": "target", "object": "trait", "study_acc": "study", "pubmedid": "pmid"})
    out = paths.bundle(release)
    out.mkdir(parents=True, exist_ok=True)
    for name, df in (("nodes", nodes), ("links", links), ("evidence", evidence), ("lineage", lin), ("neighbours", nb)):
        df.to_parquet(out / f"{name}.parquet", index=False)
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
        "files": {f: graph.file_sha(out / f) for f in FILES},
    }
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
