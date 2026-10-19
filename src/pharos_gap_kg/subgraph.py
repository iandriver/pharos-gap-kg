"""The subgraph behind an answer, as KGX.

WHAT IT HOLDS. An answer from `ask` names nodes (targets, traits, conditions, drugs) and cites edges and
gaps by id. Its subgraph is every edge and gap the answer cites, every edge and gap of the graph
between two nodes the answer names, and the nodes at their ends. A dated GWAS link enters only if it
was first reported by the year the answer is as of; the other sources' edges are undated and enter as
of the release, as the answer reads them.

GENES. A disease source's assertion runs from the gene of the protein it names, and a gene has the
protein as its product. A subgraph holds the gene of every target it names and of every assertion it
holds, joined to those targets and to the proteins the assertions name.

WHAT IT IS BUILT FROM. The bundle, so it works from the release archive with no graph files. Each row's
predicate, category, knowledge level, agent type and knowledge source come from the same mapping
objects that write the graph's own export, and the ids are the export's: a check compares every row
with the export where its files are present (`--check`).

THE FORMAT. KGX JSON: {"nodes": [...], "edges": [...]} with KGX's column names, a category and the
publications as lists, and empty fields left out. Gaps, kgx-gaps' extension, are under "gaps", which a
KGX reader ignores. "graph" says what was asked, the release and its digest, and how many rows each
part holds.

Run:  python -m pharos_gap_kg.subgraph why APOE "Alzheimer disease" --out apoe_ad.json
      python -m pharos_gap_kg.subgraph --check pharos319     # every row of a set of answers, and every row of the bundle, against the export
      python -m pharos_gap_kg.subgraph --tsv                 # the set's subgraphs as KGX TSV, for the reference validator (standards.py --dir)
      python -m pharos_gap_kg.subgraph --selftest
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import ask, graph, paths

PREFIX = graph.PREFIX
CURIE = re.compile(r"^[A-Za-z][\w.]*:\S+$")
GAP_NODE = f"{PREFIX}:GAP:"
GAP_CATEGORY = f"{PREFIX}:KnowledgeGap"
GAP_PREDICATE = f"{PREFIX}:evidence_missing_for"
GAP_ASSOC = f"{PREFIX}:KnowledgeGapAssociation"
LIMIT = 500
# The columns compared with the export, for edges and for gaps.
EMPTY = {"assertions": ["target", "gene", "condition", "source", "edge", "withheld", "rule", "papers", "lift", "cut", "evidence"],
         "expression": ["target", "condition", "label", "log2fc", "p", "edge"],
         "drug_targets": ["drug", "target", "action", "mechanism", "activity", "source", "edge"],
         "drug_indications": ["drug", "condition", "edge"], "context_nodes": ["id", "kind", "category", "name", "grounded", "xref"],
         "target_gaps": ["id", "level", "cutoff", "proposal", "dark_gap_id", "ligand_gap_id"]}
CORE = ["id", "subject", "predicate", "object", "category", "primary_knowledge_source", "knowledge_level", "agent_type", "aggregator_knowledge_source"]


class Tables:
    """The bundle's tables that hold edges, and the mapping objects that say what each edge is."""

    def __init__(self, bundle_dir: Path):
        from pharos_gap_kg import ctxgraph, kg, tgapgraph  # noqa: F401
        def rd(n: str) -> pd.DataFrame:
            # A bundle built without the other sources holds only the dated links; their tables are then empty here.
            f = bundle_dir / f"{n}.parquet"
            return pd.read_parquet(f) if f.exists() else pd.DataFrame(columns=EMPTY[n])
        self.dir = bundle_dir
        self.links, self.assertions, self.expression = rd("links"), rd("assertions"), rd("expression")
        self.drug_targets, self.drug_indications, self.evidence = rd("drug_targets"), rd("drug_indications"), rd("evidence")
        self.target_gaps = rd("target_gaps")
        nodes, ctx = rd("nodes"), rd("context_nodes")
        self.nodes = pd.concat([nodes[["id", "category", "name"]], ctx[["id", "category", "name"]]], ignore_index=True).drop_duplicates("id").set_index("id")
        self.manifest = json.loads((bundle_dir / "manifest.json").read_text())
        g = graph.mapping()
        self.link_rule = {k: g.rules[k] for k in ("gwas_condition", "gwas_other")}
        self.link_source = g.knowledge_source
        self.sources = ctxgraph.SOURCES
        self.assoc = ctxgraph.ASSOC
        self.category, self.no_gene, self.gene_only, self.GENE_PRODUCT = ctxgraph.CATEGORY, ctxgraph.NO_GENE, ctxgraph.GENE_ONLY, ctxgraph.GENE_PRODUCT
        if "gene" not in self.assertions:          # a bundle built before the assertions moved to genes
            self.assertions["gene"] = self.assertions.target
        t = nodes[nodes.kind == "target"].sort_values("id")
        t = t[t.xref.astype(str).str.startswith("NCBIGene:")]
        self.gene_of = dict(zip(t.id, t.xref))
        self.gene_names = dict(zip(t.drop_duplicates("xref").xref, t.drop_duplicates("xref").name))
        self.gene_frames, self.gene_rule = kg.gene_frames, kg.gene_mapping().rules[kg.GENE_PRODUCT_TYPE]
        self.rules = ctxgraph.RULES
        em = kg.expression_mapping()
        self.expr_rule, self.expr_source = em.rules[kg.EXPRESSION], em.knowledge_source
        dm = kg.drug_mapping()
        self.drug_rules, self.drug_source = dm.rules, dm.knowledge_source
        self.DT, self.DI, self.DIRECTION = kg.DRUG_TARGET, kg.DRUG_INDICATION, kg.DIRECTION
        tm = tgapgraph.mapping()
        self.tgap_source, self.gap_names = tm.knowledge_source, {**tm.gap_names, "precondition_unmet": "precondition unmet"}
        self.DARK, self.LIGAND = tgapgraph.DARK, tgapgraph.LIGAND


@lru_cache(maxsize=2)
def tables(bundle_dir: str) -> Tables:
    return Tables(Path(bundle_dir))


def ids_in(answer) -> set[str]:
    """Every string in the answer that is written as a compact identifier."""
    out: set[str] = set()

    def walk(x):
        if isinstance(x, dict):
            for v in x.values():
                walk(v)
        elif isinstance(x, (list, tuple)):
            for v in x:
                walk(v)
        elif isinstance(x, str) and CURIE.match(x):
            out.add(x)
    walk(answer)
    return out


def _edge(**kw) -> dict:
    return {k: v for k, v in kw.items() if v not in (None, "", [], {}) and not (isinstance(v, float) and np.isnan(v))}


def link_rows(T: Tables, L: pd.DataFrame) -> list[dict]:
    if not len(L):
        return []
    ev = T.evidence[T.evidence.target.isin(set(L.target)) & T.evidence.trait.isin(set(L.trait))]
    pubs = ev.groupby(["target", "trait"]).pmid.apply(lambda s: sorted({f"PMID:{int(p)}" for p in s})).to_dict()
    studies = ev.groupby(["target", "trait"]).study.apply(lambda s: sorted(set(s))).to_dict()
    cond = graph.CONDITION
    out = []
    for r in L.itertuples(index=False):
        cat = T.nodes.category.get(r.trait)
        rule = T.link_rule["gwas_condition" if cat in cond else "gwas_other"]
        out.append(_edge(id=r.edge_id, subject=r.target, predicate=rule.predicate, object=r.trait, category=[rule.category],
                         primary_knowledge_source=T.link_source, aggregator_knowledge_source=[graph.AGGREGATOR],
                         knowledge_level=rule.knowledge_level, agent_type=rule.agent_type,
                         publications=pubs.get((r.target, r.trait), []), **{f"{PREFIX}:studies": studies.get((r.target, r.trait), []),
                                                                            f"{PREFIX}:first_year": int(r.first_year), f"{PREFIX}:n_studies": int(r.n_studies)}))
    return out


def assertion_rows(T: Tables, A: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    edges, gaps = [], []
    for r in A.itertuples(index=False):
        ks, pred, level, agent = T.sources[r.source]
        fell = pred in T.gene_only and not str(r.gene).startswith("NCBIGene:")
        pred = T.no_gene if fell else pred
        if not r.withheld:
            edges.append(_edge(id=r.edge, subject=r.gene, predicate=pred, object=r.condition, category=[T.assoc if fell else T.category.get(r.source, T.assoc)], primary_knowledge_source=ks,
                               aggregator_knowledge_source=[graph.AGGREGATOR], knowledge_level=level, agent_type=agent,
                               **{T.GENE_PRODUCT: r.target, f"{PREFIX}:source_edge_type": r.source, f"{PREFIX}:source_evidence": r.evidence or None}))
        else:
            rules = str(r.rule).split("|") if r.rule else ["uncorroborated"]
            counted = "literature_lift" in rules
            gaps.append(_edge(id=r.edge, subject=r.gene, predicate=GAP_PREDICATE, object=f"{GAP_NODE}precondition_unmet", category=[GAP_ASSOC],
                              primary_knowledge_source=ks, aggregator_knowledge_source=[graph.AGGREGATOR], knowledge_level="logical_entailment",
                              agent_type="data_analysis_pipeline", **{f"{PREFIX}:withheld_from": pred, f"{PREFIX}:withheld_object": r.condition, T.GENE_PRODUCT: r.target,
                                                                       f"{PREFIX}:source_edge_type": r.source, f"{PREFIX}:rule": "|".join(rules),
                                                                       f"{PREFIX}:gap_reason": "; ".join(T.rules[x][0] for x in rules),
                                                                       f"{PREFIX}:restored_by": "; or ".join(T.rules[x][1] for x in rules),
                                                                       f"{PREFIX}:papers_mentioning_both": int(r.papers) if counted and pd.notna(r.papers) else None,
                                                                       f"{PREFIX}:lift": float(r.lift) if counted and pd.notna(r.lift) else None,
                                                                       f"{PREFIX}:lift_cut": float(r.cut) if counted and pd.notna(r.cut) else None}))
    return edges, gaps


def expression_rows(T: Tables, X: pd.DataFrame) -> list[dict]:
    r0 = T.expr_rule
    return [_edge(id=r.edge, subject=r.target, predicate=r0.predicate, object=r.condition, category=[r0.category], primary_knowledge_source=T.expr_source,
                  aggregator_knowledge_source=[graph.AGGREGATOR], knowledge_level=r0.knowledge_level, agent_type=r0.agent_type,
                  subject_aspect_qualifier=r0.subject_aspect_qualifier, **{T.DIRECTION: "increased" if r.log2fc > 0 else "decreased",
                                                                           f"{PREFIX}:log2_fold_change": float(r.log2fc), f"{PREFIX}:p_value": float(r.p),
                                                                           f"{PREFIX}:condition_as_named": r.label})
            for r in X.itertuples(index=False)]


def drug_rows(T: Tables, D: pd.DataFrame, I: pd.DataFrame) -> list[dict]:
    t, i = T.drug_rules[T.DT], T.drug_rules[T.DI]
    out = [_edge(id=r.edge, subject=r.drug, predicate=t.predicate, object=r.target, category=[t.category], primary_knowledge_source=T.drug_source,
                 aggregator_knowledge_source=[graph.AGGREGATOR], knowledge_level=t.knowledge_level, agent_type=t.agent_type,
                 **{f"{PREFIX}:action_type": r.action or None, f"{PREFIX}:mechanism_of_action": bool(r.mechanism), f"{PREFIX}:activity": r.activity or None})
           for r in D.itertuples(index=False)]
    out += [_edge(id=r.edge, subject=r.drug, predicate=i.predicate, object=r.condition, category=[i.category], primary_knowledge_source=T.drug_source,
                  aggregator_knowledge_source=[graph.AGGREGATOR], knowledge_level=i.knowledge_level, agent_type=i.agent_type)
            for r in I.itertuples(index=False)]
    return out


def gene_rows(T: Tables, pairs: dict) -> list[dict]:
    """pairs: protein -> its gene. The edges that make each protein its gene's product, as the export writes them."""
    if not pairs:
        return []
    _, e, _, _ = T.gene_frames(pairs, T.gene_names)
    keep = ["id", "subject", "predicate", "object", "primary_knowledge_source", "knowledge_level", "agent_type", f"{PREFIX}:source_edge_type"]
    return [_edge(**{k: x[k] for k in keep}, category=[x["category"]]) for x in e.to_dict("records")]


def target_gap_rows(T: Tables, G: pd.DataFrame, cited: set[str]) -> list[dict]:
    out = []
    for r in G.itertuples(index=False):
        for col, kind in (("dark_gap_id", T.DARK), ("ligand_gap_id", T.LIGAND)):
            gid = getattr(r, col)
            if isinstance(gid, str) and gid in cited:
                out.append(_edge(id=gid, subject=r.id, predicate=GAP_PREDICATE, object=f"{GAP_NODE}{kind}", category=[GAP_ASSOC],
                                 primary_knowledge_source=T.tgap_source, knowledge_level="logical_entailment", agent_type="data_analysis_pipeline",
                                 **{f"{PREFIX}:gap_type": kind, f"{PREFIX}:level": r.level,
                                    f"{PREFIX}:proposal": r.proposal if kind == T.DARK and isinstance(r.proposal, str) else None,
                                    f"{PREFIX}:potency_cutoff": r.cutoff if kind == T.LIGAND else None}))
    return out


def of(answer: dict, bundle_dir: Path | None = None, limit: int | None = LIMIT, question: dict | None = None) -> dict:
    """The answer's subgraph as KGX JSON. `limit` caps the edges and gaps returned (cited ones first); None returns all."""
    T = tables(str(bundle_dir or ask.default_bundle()))
    found = ids_in(answer)
    named = {i for i in found if i in T.nodes.index}
    cited = {i for i in found if i.startswith(f"{PREFIX}:e") or i.startswith(f"{PREFIX}:g")}
    as_of = answer.get("as_of") if isinstance(answer.get("as_of"), int) else None
    L = T.links
    Lm = L[(L.edge_id.isin(cited)) | (L.target.isin(named) & L.trait.isin(named) & ((L.first_year <= as_of) if as_of else True))]
    A = T.assertions
    Am = A[A.edge.isin(cited) | (A.target.isin(named) & A.condition.isin(named))]
    X = T.expression
    Xm = X[X.edge.isin(cited) | (X.target.isin(named) & X.condition.isin(named))]
    D, I = T.drug_targets, T.drug_indications
    Dm = D[D.edge.isin(cited) | (D.drug.isin(named) & D.target.isin(named))]
    Im = I[I.edge.isin(cited) | (I.drug.isin(named) & I.condition.isin(named))]
    G = T.target_gaps
    Gm = G[G.dark_gap_id.isin(cited) | G.ligand_gap_id.isin(cited)]
    a_edges, a_gaps = assertion_rows(T, Am)
    # a gene joins the targets the answer names, and the proteins its assertions name
    genes = {T.gene_of[t] for t in named if t in T.gene_of} | {e["subject"] for e in a_edges + a_gaps if str(e["subject"]).startswith("NCBIGene:")}
    products = {t for t in named if t in T.gene_of} | set(Am.target)
    g_edges = gene_rows(T, {p: T.gene_of[p] for p in sorted(products) if p in T.gene_of and T.gene_of[p] in genes})
    edges = link_rows(T, Lm) + a_edges + g_edges + expression_rows(T, Xm) + drug_rows(T, Dm, Im)
    gaps = a_gaps + target_gap_rows(T, Gm, cited)
    order = lambda rows: sorted(rows, key=lambda e: (e["id"] not in cited, e["id"]))
    edges, gaps = order(edges), order(gaps)
    n_edges, n_gaps = len(edges), len(gaps)
    if limit is not None:
        gaps = gaps[:max(0, limit - min(len(edges), limit))] if len(edges) < limit else []
        edges = edges[:limit]
    ends = {e["subject"] for e in edges + gaps} | {e["object"] for e in edges + gaps} | named
    nodes = []
    for n in sorted(ends):
        if n in T.nodes.index:
            row = T.nodes.loc[n]
            nodes.append(_edge(id=n, category=[row.category], name=row["name"] or None))
        elif n.startswith("NCBIGene:") and n in T.gene_names:
            nodes.append(_edge(id=n, category=["biolink:Gene"], name=T.gene_names[n]))
        elif n.startswith(GAP_NODE):
            nodes.append(_edge(id=n, category=[GAP_CATEGORY], name=T.gap_names.get(n[len(GAP_NODE):], n[len(GAP_NODE):].replace("_", " "))))
    m = T.manifest
    meta = {"release": m.get("release"), "source_file": m.get("source_file"), "source_sha256": m.get("source_sha256"), "as_of": as_of,
            "question": question, "rule": ("every edge and gap the answer cites, every edge and gap between two nodes it names (a GWAS link only "
                                           "if first reported by as_of; the other sources as of the release), and the nodes at their ends"),
            "counts": {"nodes": len(nodes), "edges": len(edges), "gaps": len(gaps), "edges_in_full": n_edges, "gaps_in_full": n_gaps,
                       "cited_ids": len(cited), "cited_ids_not_found": sorted(cited - {e["id"] for e in order(link_rows(T, Lm) + a_edges + g_edges + expression_rows(T, Xm)
                                                                                                             + drug_rows(T, Dm, Im) + a_gaps + target_gap_rows(T, Gm, cited))})},
            "omitted": (n_edges + n_gaps) - (len(edges) + len(gaps))}
    return {"graph": meta, "nodes": nodes, "edges": edges, "gaps": gaps}


def answer_for(G: ask.Graph, question: str, args: list, as_of: int | None = None, **kw) -> dict:
    """Ask a question as the CLI and server do."""
    if question == "known":
        return G.known(args[0], as_of, kw.get("top") or 10)
    if question == "why":
        return G.why(args[0], args[1], as_of)
    if question == "connect":
        return G.connect(args[0], args[1], as_of)
    if question == "predict":
        return G.predict(args[0], as_of, kw.get("top") or 10, bool(kw.get("unexplained_only")), bool(kw.get("unlinked")))
    if question == "whats_new":
        return G.whats_new(args[0], int(args[1]), kw.get("top") or 10)
    if question == "gaps":
        return G.gaps(args[0] if args and not kw.get("disease") else None, kw.get("family"), kw.get("top"), kw.get("disease"), as_of, kw.get("kind"))
    raise ValueError(f"no subgraph for the question {question!r}")


# ---- the check against the export ------------------------------------------------------------------

CHECK_SET = [("why", ["APOE", "Alzheimer disease"], {}), ("why", ["APOE", "Alzheimer disease"], {"as_of": 2010}),
             ("known", ["PCSK9"], {}), ("connect", ["LDLR", "PCSK9"], {}), ("predict", ["Parkinson disease"], {}),
             ("whats_new", ["type 2 diabetes mellitus", 2020], {}), ("gaps", [], {"disease": "Parkinson disease"}),
             ("gaps", ["POLR3E"], {}), ("known", ["atorvastatin"], {}), ("why", ["HMGCR", "hypercholesterolemia"], {})]


def flat(v) -> str:
    return "|".join(sorted(map(str, v))) if isinstance(v, list) else ("" if v is None else str(v))


def check(release: str, bundle_dir: Path | None = None) -> dict:
    """Every row of the check set's subgraphs against the graph's own export, column by column."""
    g = paths.GRAPH / release
    E = pd.read_csv(g / "edges.tsv", sep="\t", dtype=str, keep_default_na=False, usecols=CORE + ["publications", "subject_aspect_qualifier"]).set_index("id")
    Gp = pd.read_csv(g / "gaps.tsv", sep="\t", dtype=str, keep_default_na=False, usecols=CORE + [f"{PREFIX}:withheld_object"]).set_index("id")
    N = pd.read_csv(g / "nodes.tsv", sep="\t", dtype=str, keep_default_na=False, usecols=["id", "category", "name"]).set_index("id")
    G = ask.Graph(bundle_dir or ask.default_bundle())
    out, bad = [], []
    for q, args, kw in CHECK_SET:
        sg = of(answer_for(G, q, args, **kw), bundle_dir, limit=None, question={"question": q, "args": args, **kw})
        rows = 0
        for kind, ref, cols in (("edges", E, CORE[1:] + ["publications", "subject_aspect_qualifier"]), ("gaps", Gp, CORE[1:] + [f"{PREFIX}:withheld_object"])):
            for e in sg[kind]:
                rows += 1
                if e["id"] not in ref.index:
                    bad.append(f"{q} {args}: {kind} id {e['id']} is not in the export")
                    continue
                x = ref.loc[e["id"]]
                for c in cols:
                    mine, theirs = flat(e.get(c)), x[c]
                    if c == "publications":
                        theirs = "|".join(sorted(theirs.split("|"))) if theirs else ""
                    if mine != theirs:
                        bad.append(f"{q} {args}: {e['id']} {c}: {mine!r} here, {theirs!r} in the export")
        for n in sg["nodes"]:
            if n["id"] not in N.index or N.loc[n["id"], "category"] != n["category"][0]:
                bad.append(f"{q} {args}: node {n['id']} differs from the export")
        out.append({"question": q, "args": args, **kw, "nodes": len(sg["nodes"]), "edges": len(sg["edges"]), "gaps": len(sg["gaps"]),
                    "cited_ids_not_found": sg["graph"]["counts"]["cited_ids_not_found"], "rows_compared": rows})
    full = every_row(release, bundle_dir, E, Gp)
    bad += full.pop("differences")
    return {"release": release, "answers": out, "rows_compared": sum(o["rows_compared"] for o in out), "every_row": full,
            "differences": bad[:50], "n_differences": len(bad)}


def every_row(release: str, bundle_dir: Path | None, E: pd.DataFrame, Gp: pd.DataFrame) -> dict:
    """Every edge and gap the bundle can give, against the export: the core columns, row for row."""
    T = tables(str(bundle_dir or ask.default_bundle()))
    a_edges, a_gaps = assertion_rows(T, T.assertions)
    cited = set(T.target_gaps.dark_gap_id.dropna()) | set(T.target_gaps.ligand_gap_id.dropna())
    rows = {"edges": link_rows(T, T.links) + a_edges + gene_rows(T, T.gene_of) + expression_rows(T, T.expression) + drug_rows(T, T.drug_targets, T.drug_indications),
            "gaps": a_gaps + target_gap_rows(T, T.target_gaps, cited)}
    bad, counts = [], {}
    for kind, ref in (("edges", E), ("gaps", Gp)):
        mine = pd.DataFrame([{c: flat(e.get(c)) for c in CORE} for e in rows[kind]]).set_index("id")
        counts[kind] = {"here": int(len(mine)), "in_export": int(len(ref))}
        missing, extra = ref.index.difference(mine.index), mine.index.difference(ref.index)
        if kind == "edges":       # the ontology's same_as edges are not in the bundle, and no answer cites them
            missing = missing[~ref.loc[missing, "predicate"].eq("biolink:same_as")]
            counts[kind]["same_as_not_in_bundle"] = int((ref.predicate == "biolink:same_as").sum())
        bad += [f"every row: {kind} {i} is in the export and not here" for i in missing[:10]]
        bad += [f"every row: {kind} {i} is here and not in the export" for i in extra[:10]]
        both = mine.index.intersection(ref.index)
        for c in CORE[1:]:
            d = (mine.loc[both, c] != ref.loc[both, c])
            counts[kind][f"differ_{c}"] = int(d.sum())
            bad += [f"every row: {kind} {i} {c}: {mine.at[i, c]!r} here, {ref.at[i, c]!r} in the export" for i in d[d].index[:5]]
    return {"counts": counts, "differences": bad}


def kgx_tsv(sg: dict, out_dir: Path) -> None:
    """The subgraph as KGX TSV files, for the reference validator."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for kind in ("nodes", "edges"):
        rows = [{k: flat(v) for k, v in r.items()} for r in sg[kind]]
        pd.DataFrame(rows).to_csv(out_dir / f"{kind}.tsv", sep="\t", index=False)


def _selftest() -> int:
    bad = []

    def ck(name: str, cond: bool) -> None:
        print(("ok   " if cond else "FAIL ") + name)
        if not cond:
            bad.append(name)
    ck("identifiers are found anywhere in an answer, and prose is not",
       ids_in({"a": ["UniProtKB:P1", {"b": "MONDO:0000001"}], "c": "a sentence: with a colon", "d": 3}) == {"UniProtKB:P1", "MONDO:0000001"})
    ck("an empty field is left out of a row", _edge(id="x", a="", b=None, c=[], d=float("nan"), e=0) == {"id": "x", "e": 0})
    ck("lists flatten sorted for comparison", flat(["PMID:2", "PMID:10"]) == "PMID:10|PMID:2" and flat(None) == "")
    print(f"\n{len(bad)} failed" if bad else "\nall passed")
    return 1 if bad else 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if argv and argv[0] == "--tsv":
        # The check set's subgraphs, joined, as KGX TSV for the reference validator:
        #   .venv-standards/bin/python -I src/pharos_gap_kg/standards.py --dir data/graph/pharos319/subgraph_check docs/data/subgraph_standards_pharos319.json
        out = Path(argv[1]) if len(argv) > 1 else paths.GRAPH / "pharos319" / "subgraph_check"
        G = ask.Graph(ask.default_bundle())
        rows = {"nodes": {}, "edges": {}}
        for q, args, kw in CHECK_SET:
            s = of(answer_for(G, q, args, **kw), None, limit=None)
            for kind in rows:
                rows[kind].update({r["id"]: r for r in s[kind]})
        kgx_tsv({k: [v[i] for i in sorted(v)] for k, v in rows.items()}, out)
        print(f"{out}: {len(rows['nodes'])} nodes, {len(rows['edges'])} edges")
        return 0
    if argv and argv[0] == "--check":
        res = check(argv[1] if len(argv) > 1 else "pharos319")
        dest = paths.ROOT / "docs" / "data" / f"subgraph_check_{res['release']}.json"
        dest.write_text(json.dumps(res, indent=1, default=str) + "\n")
        print(f"{res['rows_compared']:,} rows compared over {len(res['answers'])} answers; {res['n_differences']} differ")
        for d in res["differences"][:20]:
            print("  " + d)
        return 1 if res["n_differences"] else 0
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("question", choices=["known", "why", "connect", "predict", "whats_new", "gaps"])
    ap.add_argument("args", nargs="*")
    ap.add_argument("--as-of", type=int)
    ap.add_argument("--top", type=int)
    ap.add_argument("--disease")
    ap.add_argument("--family")
    ap.add_argument("--kind")
    ap.add_argument("--unlinked", action="store_true")
    ap.add_argument("--unexplained-only", action="store_true")
    ap.add_argument("--bundle", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--limit", type=int, help="at most this many edges and gaps (default all)")
    a = ap.parse_args(argv)
    G = ask.Graph(a.bundle or ask.default_bundle())
    kw = {k: getattr(a, k) for k in ("top", "disease", "family", "kind", "unlinked", "unexplained_only")}
    sg = of(answer_for(G, a.question, a.args, a.as_of, **kw), a.bundle, limit=a.limit,
            question={"question": a.question, "args": a.args, "as_of": a.as_of, **{k: v for k, v in kw.items() if v}})
    text = json.dumps(sg, indent=1, default=str)
    if a.out:
        a.out.write_text(text + "\n")
        c = sg["graph"]["counts"]
        print(f"{a.out}: {c['nodes']} nodes, {c['edges']} edges, {c['gaps']} gaps")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
