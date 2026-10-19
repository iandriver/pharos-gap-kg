"""The eleven disease sources as KGX: what each asserts, and what the pruning rule withheld, as gaps.

Every row of the disease-source table becomes an association or a gap, and none is dropped. A row is
a gap when (1) a pruning rule withheld it (findings/PRUNING.md, findings/PRUNING2.md), or (2) the
source gives a gene identifier where a disease belongs. A gap names the rule or rules that made it one
and what would restore it. A pair withheld by the literature rule also carries what was counted: the
papers that mention both, the lift, and the cut for its source.

Each source is exported under its own knowledge source, so one (target, disease) pair asserted by
three sources is three edges with three ids, and a reader can see who said what.

THE SUBJECT IS THE GENE. Biolink's `gene_associated_with_condition` and `target_for` take a gene, and a
Pharos row names a protein, so the subject is the gene (`NCBIGene`) of that protein and the protein is
kept on the edge as `pharosgapkg:gene_product`, part of the edge's identity: two proteins of one gene
keep two edges. A row whose protein has no gene in the release stays on the protein under
`biolink:associated_with` (docs/PROTOCOL.md, 19 Oct).

Run:  python -m pharos_gap_kg.ctxgraph pharos319
      python -m pharos_gap_kg.ctxgraph --selftest
"""
from __future__ import annotations

import json
import sys

import pandas as pd
import pyarrow.parquet as pq
from kgx_gaps import Mapping, Rule, conformance
from kgx_gaps.export import Exporter

from pharos_gap_kg import ctxtables, graph, paths

PREFIX = graph.PREFIX
ASSOC = "biolink:GeneToDiseaseAssociation"
# Where a source's predicate needs a class of its own: an observation of expression, a drug's target for its indication.
CATEGORY = {"Expression Atlas": "biolink:CorrelatedGeneToDiseaseAssociation", "DrugCentral Indication": "biolink:DruggableGeneToDiseaseAssociation"}
NO_GENE = "biolink:associated_with"
GENE_ONLY = {"biolink:gene_associated_with_condition", "biolink:target_for"}     # predicates whose domain is a gene
GENE_PRODUCT = f"{PREFIX}:gene_product"
# source -> (knowledge source, predicate, knowledge level, agent type). "not_provided" where a source
# mixes curated and mined rows and the release does not say which a row is.
SOURCES = {
    "DisGeNET": ("infores:disgenet", "biolink:gene_associated_with_condition", "not_provided", "not_provided"),
    "CTD": ("infores:ctd", "biolink:gene_associated_with_condition", "knowledge_assertion", "manual_agent"),
    "eRAM": ("infores:eram", "biolink:gene_associated_with_condition", "not_provided", "not_provided"),
    "Monarch": ("infores:monarchinitiative", "biolink:gene_associated_with_condition", "not_provided", "not_provided"),
    "UniProt Disease": ("infores:uniprot", "biolink:gene_associated_with_condition", "knowledge_assertion", "manual_agent"),
    "JensenLab Text Mining": ("infores:diseases", "biolink:gene_associated_with_condition", "statistical_association", "text_mining_agent"),
    "JensenLab Knowledge MedlinePlus": ("infores:diseases", "biolink:gene_associated_with_condition", "knowledge_assertion", "manual_agent"),
    "JensenLab Knowledge UniProtKB-KW": ("infores:diseases", "biolink:gene_associated_with_condition", "knowledge_assertion", "manual_agent"),
    "JensenLab Knowledge AmyCo": ("infores:diseases", "biolink:gene_associated_with_condition", "knowledge_assertion", "manual_agent"),
    "Expression Atlas": ("infores:expression-atlas", "biolink:correlated_with", "statistical_association", "data_analysis_pipeline"),
    "DrugCentral Indication": ("infores:drugcentral", "biolink:target_for", "knowledge_assertion", "manual_agent"),
}
NOT_A_DISEASE = ("NCBIGene:",)
RULES = {"uncorroborated": ("one source, and nothing else in the graph corroborates it",
                            "a second source, the same target with a broader or narrower disease, a GWAS link, or an interaction partner with this disease"),
         "literature_lift": ("PubMed mentions this gene with this disease no more often than the cut for its source",
                             "a lift above the cut for its source, or a second source"),
         "not_a_disease": ("the source gives a gene identifier where a disease belongs", "a disease identifier from the source")}
# A gap edge runs from the target to the kind of gap, so by itself it does not say WHICH disease the
# withheld association was about. That is carried here, and it is part of what the gap is: two
# diseases withheld for one target are two gaps.
WITHHELD_OBJECT = f"{PREFIX}:withheld_object"
PAPERS, LIFT, CUT = f"{PREFIX}:papers_mentioning_both", f"{PREFIX}:lift", f"{PREFIX}:lift_cut"
EXTRA_GAP = ["aggregator_knowledge_source", WITHHELD_OBJECT, GENE_PRODUCT, f"{PREFIX}:rule", f"{PREFIX}:restored_by", PAPERS, LIFT, CUT]


def rules_of(row: dict) -> list[str]:
    """The rules that make this row a gap, none if it is an association. A withheld row that names no
    rule is round one's, the in-graph rule."""
    if str(row["object"]).startswith(NOT_A_DISEASE):
        return ["not_a_disease"]
    if not row.get("withheld"):
        return []
    named = row.get("rules")
    return named.split("|") if isinstance(named, str) and named else ["uncorroborated"]


def why_withheld(row: dict) -> str:
    return "; ".join(RULES[r][0] for r in rules_of(row))


def figure(v, integer: bool = False) -> str:
    """A counted figure as text, empty where none was counted."""
    if v is None or pd.isna(v):
        return ""
    return str(int(v)) if integer else format(float(v), ".4g")


def gap_extras(row: dict) -> dict:
    """What a withheld row carries beside its reason: which disease, which rule, what restores it, what was counted."""
    rules = rules_of(row)
    counted = "literature_lift" in rules
    return {"aggregator_knowledge_source": graph.AGGREGATOR, WITHHELD_OBJECT: row["object"], GENE_PRODUCT: row.get("gene_product") or "",
            f"{PREFIX}:rule": "|".join(rules),
            f"{PREFIX}:restored_by": "; or ".join(RULES[r][1] for r in rules),
            PAPERS: figure(row.get("papers"), integer=True) if counted else "", LIFT: figure(row.get("lift")) if counted else "", CUT: figure(row.get("cut")) if counted else ""}


def predicate_for(pred: str):
    """The source's predicate on a gene; on a protein the release gives no gene, `associated_with` where
    the source's predicate takes only a gene."""
    if pred not in GENE_ONLY:
        return pred
    return lambda row: pred if str(row["subject"]).startswith("NCBIGene:") else NO_GENE


def mapping(ks: str) -> Mapping:
    m = Mapping(prefix=PREFIX, knowledge_source=ks, gap_identity_columns=[WITHHELD_OBJECT, GENE_PRODUCT], edge_identity_columns=[GENE_PRODUCT])
    for s, (k, pred, level, agent) in SOURCES.items():
        if k == ks:
            m.rules[s] = Rule(predicate=predicate_for(pred), category=CATEGORY.get(s, ASSOC), knowledge_level=level, agent_type=agent, precondition=why_withheld)
            if pred in GENE_ONLY:
                m.declare_predicates(s, pred, NO_GENE)
    problems = m.validate()
    if problems:
        raise SystemExit("REFUSED mapping: " + "; ".join(problems))
    return m


def node_category(label: str) -> str:
    if label.startswith("NCBIGene:"):
        return "biolink:Gene"
    if label.startswith("UniProtKB:"):
        return "biolink:Protein"
    return "biolink:Disease"


def genes(release: str) -> tuple[dict, dict]:
    """Protein -> its gene, where the release gives one; and each gene's name, that of its first protein by accession."""
    prot = pq.read_table(paths.tables(release) / "protein.parquet", columns=["uniprot", "geneid", "sym"]).to_pandas().dropna(subset=["geneid"]).sort_values("uniprot")
    gene_of = {f"UniProtKB:{u}": f"NCBIGene:{int(g)}" for u, g in zip(prot.uniprot, prot.geneid)}
    first = prot.drop_duplicates("geneid")
    return gene_of, {f"NCBIGene:{int(g)}": (s if isinstance(s, str) and s else u) for g, s, u in zip(first.geneid, first.sym, first.uniprot)}


def export(rows: pd.DataFrame, names: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    """rows: subject (the gene, or the protein where it has none), gene_product (the protein), object, type
    (the source), withheld. Returns nodes, edges, gaps and rows merged as repeats."""
    unknown = sorted(set(rows.type) - set(SOURCES))
    if unknown:
        raise SystemExit(f"REFUSED: no mapping is declared for the source(s) {unknown}")

    def resolver(label):
        return label, node_category(label), True, names.get(label, "")

    def extras(row, kind):
        return {"aggregator_knowledge_source": graph.AGGREGATOR, GENE_PRODUCT: row.get("gene_product") or ""} if kind == "edge" else gap_extras(row)

    N, E, G, merged = [], [], [], 0
    for ks in sorted({v[0] for v in SOURCES.values()}):
        mine = rows[rows.type.map(lambda s: SOURCES[s][0]) == ks]
        if mine.empty:
            continue
        ex = Exporter(mapping=mapping(ks), resolver=resolver, extra_edge_columns=["aggregator_knowledge_source", GENE_PRODUCT],
                      extra_gap_columns=EXTRA_GAP, extras=extras)
        for row in mine.to_dict("records"):
            ex.add(row)
        n, e, g = ex.frames()
        N.append(n), E.append(e), G.append(g)
        merged += ex.merged
    nodes = pd.concat(N, ignore_index=True).drop_duplicates("id").reset_index(drop=True)
    return nodes, fallback_class(pd.concat(E, ignore_index=True)), pd.concat(G, ignore_index=True), merged


def fallback_class(edges: pd.DataFrame) -> pd.DataFrame:
    """An edge that fell back to `associated_with` on a protein is a plain gene-to-disease association,
    whatever class its source's predicate takes. The class is not part of an edge's id."""
    if len(edges):
        edges.loc[edges.predicate == NO_GENE, "category"] = ASSOC
    return edges


def load_rows(release: str) -> tuple[pd.DataFrame, dict]:
    """Every row of the disease-source table as subject, object, type and whether a rule withheld it; and node names."""
    C = paths.context(release)
    links = pd.read_parquet(C / "disease_links.parquet")
    wp = C / "disease_links_withheld.parquet"
    if not wp.is_file():
        raise SystemExit(f"nothing says which pairs are withheld; run:  python -m pharos_gap_kg.planted {release}  (and applycut write, for the literature rule)")
    withheld = pd.read_parquet(wp)
    for col in ("rules", "papers", "lift", "cut"):       # round one's file names no rule and counted nothing
        if col not in withheld:
            withheld[col] = None
    withheld = withheld[["target", "disease", "source", "rules", "papers", "lift", "cut"]].assign(withheld=True)
    rows = links.merge(withheld, on=["target", "disease", "source"], how="left")
    rows["withheld"] = rows.withheld.fillna(False).astype(bool)
    rows = rows.rename(columns={"target": "gene_product", "disease": "object", "source": "type"})
    gene_of, gene_names = genes(release)
    rows["subject"] = rows.gene_product.map(gene_of).fillna(rows.gene_product)
    d = paths.tables(release)
    prot = pq.read_table(d / "protein.parquet", columns=["uniprot", "sym"]).to_pandas()
    names = {"UniProtKB:" + u: (s if isinstance(s, str) and s else u) for u, s in zip(prot.uniprot, prot.sym)}
    names.update(gene_names)
    dz = pq.read_table(d / "disease.parquet", columns=["did", "mondoid", "name"]).to_pandas().dropna(subset=["name"])
    dz["id"] = [ctxtables.disease_id(m, x) for m, x in zip(dz.mondoid, dz.did)]
    best = dz.dropna(subset=["id"]).groupby(["id", "name"]).size().reset_index(name="n").sort_values(["id", "n", "name"], ascending=[True, False, True])
    names.update({i: graph.repair(n) for i, n in zip(best.drop_duplicates("id").id, best.drop_duplicates("id").name)})
    return rows, names


def build(release: str) -> dict:
    rows, names = load_rows(release)
    d = paths.tables(release)
    nodes, edges, gaps, merged = export(rows, names)
    rep = conformance.check(nodes, edges, gaps, mapping("infores:disgenet"), n_input_rows=len(rows), n_merged=merged)
    if not rep.passed or rep.n_skipped:
        raise SystemExit(f"REFUSED: the context graph does not pass conformance\n{rep}\nNothing was written.")
    out = paths.GRAPH / f"{release}_context_disease"
    files = graph.write(out, nodes, edges, gaps)
    man = {"release": release, "source_sha256": json.loads((d / "_manifest.json").read_text())["sha256"], "source_tables": ["disease", "protein"],
           "rows": len(rows), "edges": len(edges), "gaps": len(gaps), "merged": merged, "nodes": len(nodes),
           "rows_on_a_protein_with_no_gene": int((~rows.subject.str.startswith("NCBIGene:")).sum()),
           "edges_by_source": {k: int(v) for k, v in edges[f"{PREFIX}:source_edge_type"].value_counts().items()},
           "gaps_by_rule": {k: int(v) for k, v in gaps[f"{PREFIX}:rule"].value_counts().items()},
           "gaps_by_source": {k: int(v) for k, v in gaps[f"{PREFIX}:source_edge_type"].value_counts().items()},
           "predicates": {k: int(v) for k, v in edges.predicate.value_counts().items()},
           "conformance": "10 of 10 checks passed", "files": {k: graph.file_sha(p) for k, p in files.items()}}
    (out / "graph.json").write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"graph_{release}_context_disease.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    rows = pd.DataFrame([("UniProtKB:A", "MONDO:1", "DisGeNET", False), ("UniProtKB:A", "MONDO:1", "CTD", False), ("UniProtKB:B", "MONDO:2", "Monarch", True),
                         ("UniProtKB:B", "MONDO:3", "Monarch", True), ("UniProtKB:B", "NCBIGene:8", "Monarch", False),
                         ("UniProtKB:B", "NCBIGene:7", "CTD", False), ("UniProtKB:C", "MONDO:1", "Expression Atlas", False),
                         ("UniProtKB:C", "MONDO:2", "JensenLab Text Mining", False), ("UniProtKB:C", "MONDO:2", "JensenLab Knowledge AmyCo", False)],
                        columns=["subject", "object", "type", "withheld"])
    gene = {"UniProtKB:A": "NCBIGene:1", "UniProtKB:B": "NCBIGene:2", "UniProtKB:C": "NCBIGene:3", "UniProtKB:D": "NCBIGene:4", "UniProtKB:A2": "NCBIGene:1"}

    def on_genes(df):
        return df.assign(gene_product=df.subject, subject=df.subject.map(gene).fillna(df.subject))
    rows = on_genes(rows)
    counted = pd.DataFrame([("UniProtKB:D", "MONDO:4", "CTD", True, "literature_lift", 3, 1.23456, 2.626), ("UniProtKB:D", "MONDO:5", "Monarch", True, "literature_lift|uncorroborated", 0, 0.0, 0.0),
                            ("UniProtKB:D", "MONDO:6", "Monarch", True, "uncorroborated", None, None, None)], columns=["subject", "object", "type", "withheld", "rules", "papers", "lift", "cut"])
    _, e2, g2, _ = export(pd.concat([rows, on_genes(counted)], ignore_index=True), {})
    g2 = g2[g2.subject == "NCBIGene:4"].set_index(WITHHELD_OBJECT)
    ck("a pair withheld by the literature rule carries the papers, the lift and the cut it was judged on",
       (g2.at["MONDO:4", PAPERS], g2.at["MONDO:4", LIFT], g2.at["MONDO:4", CUT]) == ("3", "1.235", "2.626") and g2.at["MONDO:5", PAPERS] == "0" and g2.at["MONDO:5", LIFT] == "0")
    ck("a pair two rules withhold names both, with both reasons and both ways back",
       g2.at["MONDO:5", f"{PREFIX}:rule"] == "literature_lift|uncorroborated" and g2.at["MONDO:5", f"{PREFIX}:gap_reason"] == RULES["literature_lift"][0] + "; " + RULES["uncorroborated"][0]
       and g2.at["MONDO:5", f"{PREFIX}:restored_by"] == RULES["literature_lift"][1] + "; or " + RULES["uncorroborated"][1])
    ck("a pair the in-graph rule alone withholds carries no counts", (g2.at["MONDO:6", PAPERS], g2.at["MONDO:6", LIFT], g2.at["MONDO:6", CUT]) == ("", "", ""))
    n, e, g, merged = export(rows, {"NCBIGene:1": "AAA", "MONDO:1": "a disease"})
    rep = conformance.check(n, e, g, mapping("infores:ctd"), n_input_rows=len(rows), n_merged=merged)
    ck("every row is an association or a gap, and the result passes conformance with nothing skipped", len(e) + len(g) == len(rows) and rep.passed and not rep.n_skipped)
    ck("one pair from two sources is two edges under two knowledge sources and two ids",
       sorted(e[(e.subject == "NCBIGene:1")].primary_knowledge_source) == ["infores:ctd", "infores:disgenet"] and e.id.is_unique)
    ck("two sources that share a knowledge source still get two ids for one pair", len(e[(e.subject == "NCBIGene:3") & (e.object == "MONDO:2")]) == 2)
    ck("the subject is the gene, a Biolink Gene, and the edge names the protein the row was for",
       set(e[e.subject == "NCBIGene:1"][GENE_PRODUCT]) == {"UniProtKB:A"} and n.set_index("id").category["NCBIGene:1"] == "biolink:Gene")
    two = on_genes(pd.DataFrame([("UniProtKB:A", "MONDO:1", "CTD", False), ("UniProtKB:A2", "MONDO:1", "CTD", False), ("UniProtKB:Z", "MONDO:1", "CTD", False)],
                                columns=["subject", "object", "type", "withheld"]))
    _, e3, _, m3 = export(two, {})
    ck("two proteins of one gene keep two edges, and a protein with no gene stays the subject under associated_with",
       len(e3) == 3 and m3 == 0 and e3.id.is_unique and len(e3[e3.subject == "NCBIGene:1"]) == 2
       and e3.set_index("subject").predicate["UniProtKB:Z"] == NO_GENE and set(e3[e3.subject == "NCBIGene:1"].predicate) == {"biolink:gene_associated_with_condition"})
    ck("a withheld row is a gap naming the rule and what restores it", set(g[f"{PREFIX}:rule"]) == {"uncorroborated", "not_a_disease"}
       and g[f"{PREFIX}:restored_by"].str.len().gt(10).all() and g[f"{PREFIX}:gap_reason"].str.len().gt(10).all())
    mine = g[(g.subject == "NCBIGene:2") & (g[f"{PREFIX}:source_edge_type"] == "Monarch")]
    ck("two diseases withheld for one target from one source are two gaps, each naming its disease",
       len(mine) == 3 and mine.id.is_unique and sorted(mine[WITHHELD_OBJECT]) == ["MONDO:2", "MONDO:3", "NCBIGene:8"])
    ck("a gap names the predicate it was withheld from", set(g[f"{PREFIX}:withheld_from"]) == {"biolink:gene_associated_with_condition"})
    x = e[e[f"{PREFIX}:source_edge_type"] == "Expression Atlas"]
    ck("Expression Atlas rows are a correlation, not a gene-condition association, in the class for one",
       x.predicate.tolist() == ["biolink:correlated_with"] and x.category.tolist() == ["biolink:CorrelatedGeneToDiseaseAssociation"])
    ck("a node is named where a name is known", n.set_index("id").name["NCBIGene:1"] == "AAA" and n.set_index("id").name["MONDO:1"] == "a disease")
    try:
        export(rows.assign(type="Unknown Source"), {})
        ck("a source with no declared mapping is refused", False)
    except SystemExit:
        ck("a source with no declared mapping is refused", True)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 1:
        raise SystemExit(__doc__)
    man = build(argv[0])
    print(f"{man['release']} disease sources: {man['rows']:,} rows -> {man['edges']:,} edges + {man['gaps']:,} gaps (+{man['merged']} repeats), {man['nodes']:,} nodes; {man['conformance']}")
    print(f"  gaps by rule: {man['gaps_by_rule']}; by source: {man['gaps_by_source']}")
    print(f"  predicates: {man['predicates']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
