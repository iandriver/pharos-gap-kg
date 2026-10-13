"""One graph for a release: every dated link, what each source asserts, drugs, observations, gaps.

Five parts, each exported under its own knowledge source and joined on shared node ids:

    links         every dated GWAS link of the release, with its first year, studies and papers
    assertions    the nine disease sources that assert a gene with a condition. A row a pruning rule
                  withheld is a gap that names the rule, what was counted and what restores it.
    expression    Expression Atlas: a gene's expression is higher or lower in a condition, with the
                  fold change, the smallest change the release shows, and the p-value. An
                  observation of a difference, not an assertion that the gene is involved.
    drugs         a drug acts on a target; a drug treats an indication. A DrugCentral indication row
                  of the release is the path through both, and is not an edge of its own: the
                  release makes it by joining a drug's targets to the drug's indications.
    ontology      a GWAS trait whose id the ontology has retired is the same as its current term.
                  The GWAS table names traits by EFO ids and the disease sources mostly by MONDO
                  ids, so this is what joins the dated links to the rest on the condition side.
    target gaps   what each target lacks for the next development level

Every input row is an edge or a gap, and nothing is written unless all ten conformance checks pass.
The definitions are in docs/PROTOCOL.md, "one graph, and each source exported as what it asserts".

Run:  python -m pharos_gap_kg.kg build pharos319
      python -m pharos_gap_kg.kg --selftest
"""
from __future__ import annotations

import hashlib
import json
import sys

import pandas as pd
import pyarrow.parquet as pq
from kgx_gaps import Mapping, Rule, conformance
from kgx_gaps.export import Exporter

from pharos_gap_kg import ctxgraph, ctxtables, graph, ontology, paths, pubyear, registry, tgapgraph

PREFIX = graph.PREFIX
EXPRESSION, INDICATION = "Expression Atlas", "DrugCentral Indication"
REMODELLED = (EXPRESSION, INDICATION)
FLOOR = 0.999                     # the release holds no fold change below 1.0 in size, at three decimals
AGG = "aggregator_knowledge_source"
EVIDENCE, ROW, P = f"{PREFIX}:source_evidence", f"{PREFIX}:source_row", f"{PREFIX}:p_value"
DIRECTION = "subject_direction_qualifier"
ACTION, MOA, ACTIVITY, ACT_SOURCE = (f"{PREFIX}:{c}" for c in ("action_type", "mechanism_of_action", "activity", "activity_source"))
DRUG_TARGET, DRUG_INDICATION = "drug_target", "drug_indication"
CONDITION = "biolink:DiseaseOrPhenotypicFeature"
INT_COLUMNS = graph.INT_COLUMNS + [ROW]


def graph_dir(release: str):
    return paths.GRAPH / release


def fmt(v) -> str:
    return format(float(v), ".4g")


# ---- assertions ------------------------------------------------------------------------------------

def evidence_text(ev: pd.DataFrame) -> dict:
    """(target, disease, source) -> the source's own evidence for the pair, as "field=value; field=value"."""
    out: dict = {}
    for r in ev.sort_values(["target", "disease", "source", "field"]).itertuples(index=False):
        v = r.value
        out.setdefault((r.target, r.disease, r.source), []).append(f"{r.field}={fmt(v) if isinstance(v, (int, float)) else v}")
    return {k: "; ".join(v) for k, v in out.items()}


def assertions(rows: pd.DataFrame, names: dict, evidence: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    """The nine sources that assert a gene with a condition. rows: subject, object, type, withheld (+ rules, papers, lift, cut)."""
    if set(rows.type) & set(REMODELLED):
        raise SystemExit(f"REFUSED: {sorted(set(rows.type) & set(REMODELLED))} are not assertions of a gene with a condition here")

    def resolver(label):
        return label, ("biolink:Protein" if label.startswith("UniProtKB:") else "biolink:Disease"), True, names.get(label, "")

    def extras(row, kind):
        if kind == "edge":
            return {AGG: graph.AGGREGATOR, EVIDENCE: evidence.get((row["subject"], row["object"], row["type"]), "")}
        return ctxgraph.gap_extras(row)

    N, E, G, merged = [], [], [], 0
    for ks in sorted({ctxgraph.SOURCES[s][0] for s in set(rows.type)}):
        mine = rows[rows.type.map(lambda s: ctxgraph.SOURCES[s][0]) == ks]
        ex = Exporter(mapping=ctxgraph.mapping(ks), resolver=resolver, extra_edge_columns=[AGG, EVIDENCE], extra_gap_columns=ctxgraph.EXTRA_GAP, extras=extras)
        for row in mine.to_dict("records"):
            ex.add(row)
        n, e, g = ex.frames()
        N.append(n), E.append(e), G.append(g)
        merged += ex.merged
    return pd.concat(N, ignore_index=True).drop_duplicates("id"), pd.concat(E, ignore_index=True), pd.concat(G, ignore_index=True), merged


# ---- expression ------------------------------------------------------------------------------------

def expression_mapping() -> Mapping:
    m = Mapping(prefix=PREFIX, knowledge_source="infores:expression-atlas", edge_identity_columns=[DIRECTION, ROW])
    m.rules[EXPRESSION] = Rule(predicate="biolink:correlated_with", category="biolink:GeneToDiseaseAssociation", knowledge_level="statistical_association",
                               agent_type="data_analysis_pipeline", subject_aspect_qualifier="expression")
    problems = m.validate()
    if problems:
        raise SystemExit("REFUSED mapping: " + "; ".join(problems))
    return m


def expression(rows: pd.DataFrame, names: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    """rows: subject, object, label (the source's name for the condition), log2fc, p, row (the release's row id)."""
    def resolver(label):
        return label, ("biolink:Protein" if label.startswith("UniProtKB:") else "biolink:Disease"), True, names.get(label, "")

    def extras(row, kind):
        return {AGG: graph.AGGREGATOR, DIRECTION: "increased" if row["effect"] > 0 else "decreased", P: fmt(row["p"]), ROW: int(row["row"])} if kind == "edge" else {}

    ex = Exporter(mapping=expression_mapping(), resolver=resolver, extra_edge_columns=[AGG, DIRECTION, P, ROW], extras=extras)
    for r in rows.itertuples(index=False):
        if ex.add({"type": EXPRESSION, "subject": r.subject, "object": r.object, "context": r.label, "effect": float(r.log2fc), "floor": FLOOR, "p": float(r.p), "row": int(r.row)}) != "edge":
            raise SystemExit(f"REFUSED: the expression row {r.row} did not clear the release's own threshold")
    n, e, g = ex.frames()
    return n, e, g, ex.merged


# ---- drugs -----------------------------------------------------------------------------------------

def drug_mapping() -> Mapping:
    m = Mapping(prefix=PREFIX, knowledge_source="infores:drugcentral")
    m.rules[DRUG_TARGET] = Rule(predicate="biolink:directly_physically_interacts_with", category="biolink:ChemicalGeneInteractionAssociation",
                                knowledge_level="knowledge_assertion", agent_type="manual_agent")
    m.rules[DRUG_INDICATION] = Rule(predicate="biolink:treats", category="biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation",
                                    knowledge_level="knowledge_assertion", agent_type="manual_agent")
    problems = m.validate()
    if problems:
        raise SystemExit("REFUSED mapping: " + "; ".join(problems))
    return m


def indication_id(name: str) -> str:
    """An indication the release names without an identifier: a node under the local prefix, from the name."""
    return f"{PREFIX}:indication_{hashlib.sha256(' '.join(name.lower().split()).encode()).hexdigest()[:12]}"


def joined(values) -> str:
    return "|".join(sorted({str(v) for v in values if isinstance(v, str) and v}))


def drug_tables(activity: pd.DataFrame, indications: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """activity: drug, dcid, target (a CURIE), action_type, has_moa, act_type, act_value, source, chembl.
    indications: target (a CURIE), drug, name, disease (an id or None). Returns the drugs, one row per
    (drug, target), one row per (drug, indication), and the counts the protocol quotes."""
    a = activity.copy()
    a["drug_id"] = "DrugCentral:" + a.dcid.astype(int).astype(str)
    drugs = a.sort_values("drug").groupby("drug_id").agg(name=("drug", "first"), chembl=("chembl", lambda s: next((x for x in s if isinstance(x, str) and x), ""))).reset_index()
    a["act"] = [f"{t}={fmt(v)}" if isinstance(t, str) and t and pd.notna(v) else "" for t, v in zip(a.act_type, a.act_value)]
    dt = a.groupby(["drug_id", "target"]).agg(action=("action_type", joined), moa=("has_moa", lambda s: bool((s == 1).any())), activity=("act", joined),
                                               source=("source", joined), rows=("drug", "size")).reset_index()
    by_name = {n.lower(): i for n, i in zip(a.drug, a.drug_id)}
    ind = indications.copy()
    ind["drug_id"] = ind.drug.str.lower().map(by_name)
    if ind.drug_id.isna().any():
        raise SystemExit(f"REFUSED: {int(ind.drug_id.isna().sum()):,} indication rows name a drug the drug-activity table does not hold")
    ind["name"] = ind.name.map(graph.repair)
    # one indication, one node: an id where any row of that name gives one, else a local id from the name
    known = ind.dropna(subset=["disease"]).groupby("name").disease.agg(lambda s: sorted(set(s))[0]).to_dict()
    ind["object"] = [known.get(n) or indication_id(n) for n in ind.name]
    ind["grounded"] = ind.name.isin(known)
    acts = dt.set_index(["drug_id", "target"])
    key = list(zip(ind.drug_id, ind.target))
    missing = [k for k in set(key) if k not in acts.index]
    if missing:
        raise SystemExit(f"REFUSED: {len(missing):,} (drug, target) pairs of the indication rows are not in the drug-activity table, e.g. {missing[:2]}")
    via_moa = acts.moa.reindex(pd.MultiIndex.from_tuples(key)).to_numpy()
    di = ind.sort_values("name").groupby(["drug_id", "object"]).agg(name=("name", "first"), grounded=("grounded", "first"), rows=("drug", "size")).reset_index()
    counts = {"indication_rows": len(ind), "through_the_mechanism_target": int(via_moa.sum()), "through_another_activity": int((~via_moa).sum()),
              "indication_rows_without_an_id": int(indications.disease.isna().sum()), "indications": int(di.object.nunique()),
              "indications_not_grounded": int(di[~di.grounded].object.nunique()), "drugs": len(drugs), "drugs_with_an_indication": int(di.drug_id.nunique()),
              "activity_rows": len(a), "mechanism_pairs": int(dt.moa.sum())}
    return drugs, dt, di, counts


def drug_frames(drugs: pd.DataFrame, dt: pd.DataFrame, di: pd.DataFrame, names: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    dn = dict(zip(drugs.drug_id, drugs.name))
    grounded = dict(zip(di.object, di.grounded))
    iname = dict(zip(di.object, di.name))

    def resolver(label):
        if label.startswith("DrugCentral:"):
            return label, "biolink:Drug", True, dn.get(label, "")
        if label.startswith("UniProtKB:"):
            return label, "biolink:Protein", True, names.get(label, "")
        ok = bool(grounded.get(label, True))
        return label, ("biolink:Disease" if ok else CONDITION), ok, names.get(label) or iname.get(label, "")

    def extras(row, kind):
        if kind != "edge":
            return {}
        if row["type"] == DRUG_TARGET:
            return {AGG: graph.AGGREGATOR, ACTION: row["action"], MOA: "true" if row["moa"] else "false", ACTIVITY: row["activity"], ACT_SOURCE: row["source"]}
        return {AGG: graph.AGGREGATOR}

    ex = Exporter(mapping=drug_mapping(), resolver=resolver, extra_edge_columns=[AGG, ACTION, MOA, ACTIVITY, ACT_SOURCE], extras=extras)
    for r in drugs.itertuples(index=False):
        ex.node(r.drug_id, "biolink:Drug", r.name, True, f"CHEMBL.COMPOUND:{r.chembl}" if r.chembl else "")
    for r in dt.itertuples(index=False):
        ex.add({"type": DRUG_TARGET, "subject": r.drug_id, "object": r.target, "action": r.action, "moa": r.moa, "activity": r.activity, "source": r.source})
    for r in di.itertuples(index=False):
        ex.add({"type": DRUG_INDICATION, "subject": r.drug_id, "object": r.object})
    n, e, g = ex.frames()
    return n, e, g, ex.merged


# ---- ontology: a retired trait id is the same as its current term ---------------------------------

def current_terms(traits, terms: dict) -> tuple[list[tuple[str, str]], int]:
    """(trait as published, its current term) where the ontology retired the id and names a replacement
    whose kind is declared; and how many replacements were left out for having no declared kind."""
    pairs, undeclared = [], 0
    for t in sorted(traits):
        c = ontology.canonical(t, terms)
        if c and c != t:
            if c.split(":", 1)[0] in graph.TRAIT_CATEGORY:
                pairs.append((t, c))
            else:
                undeclared += 1
    return pairs, undeclared


def ontology_frames(pairs: list[tuple[str, str]], names: dict, terms: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ex = Exporter(mapping=Mapping(prefix=PREFIX, knowledge_source="infores:efo"), extra_edge_columns=[AGG])
    for a, b in pairs:
        ex.node(a, graph.trait_category(a), names.get(a, ""), True)
        ex.node(b, graph.trait_category(b), names.get(b) or graph.repair(terms.get(b, {}).get("name", "")), True)
        ex.structural(a, "biolink:same_as", b, knowledge_level="knowledge_assertion", agent_type="manual_agent")
    return ex.frames()


def cross_references(path) -> dict:
    """id -> the ids its `xref:` lines name. A cross-reference is not a statement that two terms are one."""
    out, cur = {}, None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("id: "):
                cur = line[4:].strip()
            elif line.startswith("xref: ") and cur:
                out.setdefault(cur, set()).add(line[6:].split()[0])
            elif not line.strip():
                cur = None
    return out


# ---- the release's tables --------------------------------------------------------------------------

def release_tables(release: str) -> dict:
    d = paths.tables(release)
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot"]).to_pandas()
    target_of = dict(zip(prot.id, "UniProtKB:" + prot.uniprot))
    dz = pq.read_table(d / "disease.parquet", columns=["id", "dtype", "protein_id", "name", "did", "mondoid", "drug_name", "log2foldchange", "pvalue"]).to_pandas()
    dz["target"] = dz.protein_id.map(target_of)
    dz["disease"] = [ctxtables.disease_id(m, x) for m, x in zip(dz.mondoid, dz.did)]
    ea = dz[dz.dtype == EXPRESSION]
    if ea.target.isna().any() or ea.disease.isna().any():
        raise SystemExit("REFUSED: an Expression Atlas row has no target or no condition id")
    exp = pd.DataFrame({"subject": ea.target, "object": ea.disease, "label": ea.name.map(graph.repair), "log2fc": ea.log2foldchange.astype(float),
                        "p": pd.to_numeric(ea.pvalue), "row": ea.id.astype(int)}).reset_index(drop=True)
    dc = dz[dz.dtype == INDICATION]
    ind = pd.DataFrame({"target": dc.target, "drug": dc.drug_name, "name": dc.name, "disease": dc.disease}).reset_index(drop=True)
    act = pq.read_table(d / "drug_activity.parquet").to_pandas()
    t2 = pq.read_table(d / "t2tc.parquet", columns=["target_id", "protein_id"]).to_pandas().dropna()
    act["target"] = act.target_id.map(dict(zip(t2.target_id, t2.protein_id.map(target_of))))
    if act.target.isna().any():
        raise SystemExit(f"REFUSED: {int(act.target.isna().sum()):,} drug-activity rows are on a target with no protein")
    act = act.rename(columns={"cmpd_chemblid": "chembl"})[["drug", "dcid", "target", "action_type", "has_moa", "act_type", "act_value", "source", "chembl"]]
    return {"expression": exp, "indications": ind, "activity": act, "source_sha256": json.loads((d / "_manifest.json").read_text())["sha256"]}


# ---- the merged graph ------------------------------------------------------------------------------

def merge(parts: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """parts: name -> (nodes, edges, gaps), in the order that decides a shared node's name and category."""
    def cat(i):
        frames = [p[i].assign(**{f"{PREFIX}:part": k}) if i else p[i] for k, p in parts.items() if len(p[i])]
        cols = list(dict.fromkeys(c for f in frames for c in f.columns))
        return pd.concat([f.reindex(columns=cols) for f in frames], ignore_index=True).fillna("") if frames else pd.DataFrame()
    nodes = cat(0).drop_duplicates("id").reset_index(drop=True)
    edges, gaps = cat(1), cat(2)
    for c in INT_COLUMNS:
        if c in edges:
            edges[c] = [("" if v == "" else str(int(float(v)))) for v in edges[c]]
    return nodes, edges, gaps


def build(release: str) -> dict:
    T = release_tables(release)
    # links: every dated link of the release
    years = pubyear.load().set_index("pmid").year
    ev = graph.evidence(release, years, cutoff=10**6)
    known, xref = graph.nodes_of(release)
    ex = graph.export(ev, known, xref)
    links = graph.frames(ex)
    # assertions: the nine sources, kept or withheld
    rows, names = ctxgraph.load_rows(release)
    if rows[rows.type.isin(REMODELLED)].withheld.any():
        raise SystemExit("REFUSED: a pruning rule withheld a row of a remodelled source; the protocol says none is")
    nine = rows[~rows.type.isin(REMODELLED)]
    names.update({k: v[3] for k, v in known.items()})          # a trait or target named by the dated links keeps that name
    evt = evidence_text(pd.read_parquet(paths.context(release) / "disease_evidence.parquet").query("source not in @REMODELLED"))
    an, ae, ag, amerged = assertions(nine, names, evt)
    # expression, drugs, target gaps
    xn, xe, xg, xmerged = expression(T["expression"], names)
    drugs, dt, di, counts = drug_tables(T["activity"], T["indications"])
    dn, de, dg, dmerged = drug_frames(drugs, dt, di, names)
    gap_rows, tnames, summary, summary_file = tgapgraph.load_rows(release)
    tn, te, tg, tmerged = tgapgraph.export(gap_rows, tnames)
    # ontology: tie a retired trait id to its current term
    efo = registry.require("efo_obo")
    terms = ontology.read(efo)
    traits = sorted(set(ev.object))
    pairs, undeclared = current_terms(traits, terms)
    on, oe, og = ontology_frames(pairs, names, terms)
    parts = {"links": links, "assertions": (an, ae, ag), "expression": (xn, xe, xg), "drugs": (dn, de, dg), "ontology": (on, oe, og), "target_gaps": (tn, te, tg)}
    nodes, edges, gaps = merge(parts)
    n_rows = len(ev) + len(nine) + len(T["expression"]) + len(dt) + len(di) + len(gap_rows)
    merged = ex.merged + amerged + xmerged + dmerged + tmerged
    rep = conformance.check(nodes, edges.drop(columns=[f"{PREFIX}:part"]), gaps.drop(columns=[f"{PREFIX}:part"]),
                            Mapping(prefix=PREFIX, knowledge_source=graph.AGGREGATOR), n_input_rows=n_rows, n_merged=merged)
    if not rep.passed or rep.n_skipped:
        raise SystemExit(f"REFUSED: the merged graph does not pass conformance\n{rep}\nNothing was written.")
    part = f"{PREFIX}:part"
    out = graph_dir(release)
    files = graph.write(out, nodes, edges, gaps)
    lo = lambda p: set(edges[edges[part] == p].object) | set(edges[edges[part] == p].subject)
    conditions = set(nodes[~nodes.category.isin(["biolink:Protein", "biolink:Drug"])].id)
    elsewhere = (lo("assertions") | lo("expression") | lo("drugs")) & conditions
    now = dict(pairs)
    by_id = {t for t in traits if t in elsewhere}
    by_term = {t for t in traits if t not in by_id and now.get(t) in elsewhere}
    xr = cross_references(efo)
    by_xref = {t for t in traits if t not in by_id and t not in by_term and (xr.get(t, set()) | xr.get(now.get(t, ""), set())) & elsewhere}
    ea = T["expression"]
    man = {"release": release, "source_sha256": T["source_sha256"], "target_gaps_summary_sha256": graph.file_sha(summary_file),
           "rows": n_rows, "merged": merged, "nodes": len(nodes), "edges": len(edges), "gaps": len(gaps),
           "node_categories": {k: int(v) for k, v in nodes.category.value_counts().items()},
           "nodes_not_grounded": int((nodes.id_grounded.astype(str).str.lower() == "false").sum()),
           "edges_by_part": {k: int(v) for k, v in edges[part].value_counts().items()}, "gaps_by_part": {k: int(v) for k, v in gaps[part].value_counts().items()},
           "predicates": {k: int(v) for k, v in edges.predicate.value_counts().items()},
           "edge_categories": {k: int(v) for k, v in edges.category.value_counts().items()},
           "knowledge_sources": {k: int(v) for k, v in edges.primary_knowledge_source.value_counts().items()},
           "links": {"rows": len(ev), "first_year": int(ev.first_year.min()), "latest_year": int(ev.first_year.max())},
           "assertions": {"rows": len(nine), "edges": len(ae), "gaps": len(ag), "edges_by_source": {k: int(v) for k, v in ae[f"{PREFIX}:source_edge_type"].value_counts().items()},
                          "gaps_by_source": {k: int(v) for k, v in ag[f"{PREFIX}:source_edge_type"].value_counts().items()},
                          "gaps_by_rule": {k: int(v) for k, v in ag[f"{PREFIX}:rule"].value_counts().items()},
                          "edges_with_the_source's_own_evidence": int((ae[EVIDENCE] != "").sum())},
           "expression": {"rows": len(ea), "edges": len(xe), "increased": int((ea.log2fc > 0).sum()), "decreased": int((ea.log2fc < 0).sum()), "conditions_as_named": int(ea.label.nunique()),
                          "conditions": int(ea.object.nunique()), "smallest_fold_change": float(ea.log2fc.abs().min()), "largest_p": float(ea.p.max()), "floor": FLOOR},
           "drugs": {**counts, "drug_target_edges": int((de[f"{PREFIX}:source_edge_type"] == DRUG_TARGET).sum()), "drug_indication_edges": int((de[f"{PREFIX}:source_edge_type"] == DRUG_INDICATION).sum())},
           "target_gaps": {"gaps": len(tg), "by_type": {k: int(v) for k, v in tg[f"{PREFIX}:gap_type"].value_counts().items()}},
           "ontology": {"sha256": graph.file_sha(efo), "traits": len(traits), "traits_retired_and_replaced": len(pairs) + undeclared, "same_as_edges": len(oe),
                        "replacements_with_no_declared_kind": undeclared},
           "conditions": len(conditions),
           "traits_that_reach_another_part": {"under_their_own_id": len(by_id), "through_their_current_term": len(by_term), "by_a_cross_reference_only_and_not_joined": len(by_xref),
                                              "joined": len(by_id) + len(by_term)},
           "targets_in_more_than_one_part": int(sum(len({p for p in parts if t in sets[p]}) > 1 for sets in [{p: lo(p) | set(gaps[gaps[part] == p].subject) for p in parts}]
                                                    for t in set(nodes[nodes.category == "biolink:Protein"].id))),
           "conformance": "10 of 10 checks passed", "files": {k: graph.file_sha(p) for k, p in files.items()}}
    (out / "graph.json").write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"graph_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    names = {"UniProtKB:A": "AAA", "UniProtKB:B": "BBB", "MONDO:1": "a disease", "MONDO:2": "another"}
    # assertions carry the source's own evidence; a withheld row is a gap
    rows = pd.DataFrame([("UniProtKB:A", "MONDO:1", "DisGeNET", False, None, None, None, None), ("UniProtKB:A", "MONDO:1", "CTD", False, None, None, None, None),
                         ("UniProtKB:B", "MONDO:2", "CTD", True, "literature_lift", 3, 1.2, 2.6)], columns=["subject", "object", "type", "withheld", "rules", "papers", "lift", "cut"])
    ev = pd.DataFrame([("UniProtKB:A", "MONDO:1", "DisGeNET", "score", 0.3), ("UniProtKB:A", "MONDO:1", "DisGeNET", "pubmed_ids", 4.0), ("UniProtKB:A", "MONDO:1", "CTD", "therapeutic", 0.0)],
                      columns=["target", "disease", "source", "field", "value"])
    an, ae, ag, _ = assertions(rows, names, evidence_text(ev))
    ck("an assertion carries the source's own evidence, and a withheld row is a gap with what was counted",
       sorted(ae[EVIDENCE]) == ["pubmed_ids=4; score=0.3", "therapeutic=0"] and len(ag) == 1 and ag[ctxgraph.LIFT].tolist() == ["1.2"])
    try:
        assertions(rows.assign(type=EXPRESSION), names, {})
        ck("a remodelled source is refused as an assertion", False)
    except SystemExit:
        ck("a remodelled source is refused as an assertion", True)
    # expression: direction, size, floor, p, the source's label, one edge per release row
    x = pd.DataFrame([("UniProtKB:A", "MONDO:1", "active disease", 1.0, 0.04, 11), ("UniProtKB:A", "MONDO:1", "active disease", -2.5, 0.001, 12), ("UniProtKB:B", "MONDO:1", "a disease", 3.0, 0.01, 13)],
                     columns=["subject", "object", "label", "log2fc", "p", "row"])
    xn, xe, xg, xm = expression(x, names)
    ck("an expression row is an observation with a direction, a size above its floor and a p-value",
       len(xe) == 3 and not len(xg) and xe[DIRECTION].tolist() == ["increased", "decreased", "increased"] and set(xe.subject_aspect_qualifier) == {"expression"}
       and xe[f"{PREFIX}:effect_size"].astype(float).tolist() == [1.0, -2.5, 3.0] and set(xe[f"{PREFIX}:detection_floor"].astype(float)) == {FLOOR} and set(xe[f"{PREFIX}:detected"]) == {"true"}
       and xe[P].astype(float).tolist() == [0.04, 0.001, 0.01])
    ck("two rows of one pair and condition are two edges, told apart by the release's row", xe.id.is_unique and xe[f"{PREFIX}:context"].tolist()[:2] == ["active disease"] * 2)
    # drugs: the indication row is a path
    act = pd.DataFrame([("druga", 7, "UniProtKB:A", "INHIBITOR", 1, "IC50", 8.4, "CHEMBL", "CHEMBL1"), ("druga", 7, "UniProtKB:A", "INHIBITOR", 1, "Ki", 7.0, "IUPHAR", "CHEMBL1"),
                        ("druga", 7, "UniProtKB:B", None, 0, None, None, "DRUG LABEL", "CHEMBL1"), ("drugb", 9, "UniProtKB:B", "AGONIST", 1, None, None, "CHEMBL", None)],
                       columns=["drug", "dcid", "target", "action_type", "has_moa", "act_type", "act_value", "source", "chembl"])
    ind = pd.DataFrame([("UniProtKB:A", "DrugA", "a disease", "MONDO:1"), ("UniProtKB:B", "druga", "a disease", "MONDO:1"), ("UniProtKB:B", "drugb", "Pain after surgery", None),
                        ("UniProtKB:B", "drugb", "pain  after surgery", None)], columns=["target", "drug", "name", "disease"])
    drugs, dt, di, c = drug_tables(act, ind)
    ck("one edge per drug and target, with every action, activity and source the release gives, and whether it is the mechanism",
       len(dt) == 3 and dt.set_index(["drug_id", "target"]).loc[("DrugCentral:7", "UniProtKB:A")].tolist() == ["INHIBITOR", True, "IC50=8.4|Ki=7", "CHEMBL|IUPHAR", 2])
    ck("an indication row is counted by whether it runs through the drug's mechanism target",
       (c["indication_rows"], c["through_the_mechanism_target"], c["through_another_activity"], c["indication_rows_without_an_id"]) == (4, 3, 1, 2))
    ck("an indication with no identifier is one node under the local prefix, however its name is spaced or cased",
       c["indications"] == 2 and c["indications_not_grounded"] == 1 and di[~di.grounded].object.tolist() == [indication_id("Pain after surgery")] and indication_id("pain  after SURGERY") == indication_id("Pain after surgery"))
    dn, de, dg, dm = drug_frames(drugs, dt, di, names)
    ck("a drug is a node with its ChEMBL id; it acts on targets and treats indications",
       dn.set_index("id").loc["DrugCentral:7"].tolist()[:2] == ["biolink:Drug", "druga"] and dn.set_index("id").xref["DrugCentral:7"] == "CHEMBL.COMPOUND:CHEMBL1"
       and de.predicate.value_counts().to_dict() == {"biolink:directly_physically_interacts_with": 3, "biolink:treats": 2}
       and de[de.predicate == "biolink:directly_physically_interacts_with"][MOA].tolist() == ["true", "false", "true"])
    ck("the indication with no identifier is flagged as not grounded", dn.set_index("id").id_grounded[indication_id("Pain after surgery")] == "false")
    try:
        drug_tables(act, pd.DataFrame([("UniProtKB:A", "drugb", "a disease", "MONDO:1")], columns=["target", "drug", "name", "disease"]))
        ck("an indication row whose drug has no activity on that target is refused", False)
    except SystemExit:
        ck("an indication row whose drug has no activity on that target is refused", True)
    # ontology: a retired id is the same as its current term
    terms = ontology.load("""[Term]
id: EFO:1
name: old
is_obsolete: true
replaced_by: http://purl.obolibrary.org/obo/MONDO_1

[Term]
id: MONDO:1
name: a disease

[Term]
id: EFO:2
name: live

[Term]
id: EFO:3
name: gone
is_obsolete: true
replaced_by: OBA:9

[Term]
id: OBA:9
name: an attribute
""".splitlines())
    pairs, undeclared = current_terms(["EFO:1", "EFO:2", "EFO:3"], terms)
    on, oe, og = ontology_frames(pairs, {"EFO:1": "old name"}, terms)
    ck("a retired trait id is tied to its current term by same_as; a live id is not; a replacement of no declared kind is counted",
       pairs == [("EFO:1", "MONDO:1")] and undeclared == 1 and oe.predicate.tolist() == ["biolink:same_as"] and (oe.subject[0], oe.object[0]) == ("EFO:1", "MONDO:1")
       and set(oe[f"{PREFIX}:source_edge_type"]) == {""} and on.set_index("id").name.to_dict() == {"EFO:1": "old name", "MONDO:1": "a disease"})
    # the merged graph passes conformance, every row accounted for, one node per id
    nodes, edges, gaps = merge({"assertions": (an, ae, ag), "expression": (xn, xe, xg), "drugs": (dn, de, dg), "ontology": (on, oe, og)})
    part = f"{PREFIX}:part"
    rep = conformance.check(nodes, edges.drop(columns=[part]), gaps.drop(columns=[part]), Mapping(prefix=PREFIX, knowledge_source=graph.AGGREGATOR),
                            n_input_rows=len(rows) + len(x) + len(dt) + len(di), n_merged=0)
    ck("the parts join on shared ids into one graph that passes all ten checks", rep.passed and not rep.n_skipped and nodes.id.is_unique and edges.id.is_unique
       and set(edges[part]) == {"assertions", "expression", "drugs", "ontology"} and nodes.set_index("id").name["MONDO:1"] == "a disease")
    ck("a whole-number column stays a whole number", edges[edges[part] == "expression"][ROW].tolist() == ["11", "12", "13"])
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 2 or argv[0] != "build":
        raise SystemExit(__doc__)
    m = build(argv[1])
    print(f"{m['release']}: {m['rows']:,} rows -> {m['edges']:,} edges + {m['gaps']:,} gaps, {m['nodes']:,} nodes; {m['conformance']}")
    print(f"  edges by part: {m['edges_by_part']}; gaps by part: {m['gaps_by_part']}")
    print(f"  drugs: {m['drugs']}")
    print(f"  expression: {m['expression']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
