"""Everything the question-answering front end reads, built once into one directory with a manifest.

The front end answers from the bundle and from nothing else, so what a user is told can be traced to
these files and their digests:

    nodes.parquet        every Pharos target and every GWAS trait, with what names and grades them
    links.parquet        one row per target-trait link: its edge id, first and last year, study counts
    evidence.parquet     one row per (link, study): the study, its paper and the paper's year
    lineage.parquet      pairs of traits that are one term, or broader and narrower
    neighbours.parquet   pairs of targets within 500,000 bases, and the bases between them
    backtest.json        the scorer the front end ranks by, how it held on the prospective test, how often a
                         candidate at each rank was later reported, and the baselines on the 2018 split
    target_gaps.*        what each target lacks for the next development level, and how often targets in
                         each band were promoted in the back-test (where `tgaps` has been run)
    assertions.parquet, expression.parquet, drug_targets.parquet, drug_indications.parquet,
    context_nodes.parquet
                         the merged graph's other parts (where `kg` has built it): which sources assert
                         a target with a condition and which rows a rule withheld, expression
                         observations, drugs, and the conditions and drugs those name
    manifest.json        the inputs' digests, each file's digest, and the counts

The links are every link in the release, whatever its year. "As of 2018" is a filter the front end
applies; the bundle does not choose a date. Edge ids come from the same exporter as the KGX graph,
and where a dated graph has been built its edges must be among these under the same ids.

Run:  python -m pharos_gap_kg.bundle pharos319
      python -m pharos_gap_kg.bundle pharos319 --backtest-only     # refresh backtest.json in a built bundle
      python -m pharos_gap_kg.bundle pharos319 --context-only      # add the merged graph's other parts to a built bundle
      python -m pharos_gap_kg.bundle pharos319 --sizes-only        # add the size of each study behind the links
      python -m pharos_gap_kg.bundle --selftest
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from kgx_gaps import conformance

from pharos_gap_kg import context, ctxtables, graph, ontology, paths, position, pubyear, registry

FILES = ["nodes.parquet", "links.parquet", "evidence.parquet", "lineage.parquet", "neighbours.parquet"]
GAP_FILES = ["target_gaps.parquet", "target_gaps.json"]
SIZE_FILES = ["study_sizes.parquet"]
CONTEXT_FILES = ["context_nodes.parquet", "assertions.parquet", "expression.parquet", "drug_targets.parquet", "drug_indications.parquet"]


def context_tables(nodes: pd.DataFrame, edges: pd.DataFrame, gaps: pd.DataFrame, traits: set) -> dict[str, pd.DataFrame]:
    """The merged graph's other parts as the tables the questions read, each row with the id of its edge
    or gap. `traits` are the GWAS traits the bundle's node table already holds."""
    P = graph.PREFIX
    part, src = f"{P}:part", f"{P}:source_edge_type"
    num = lambda s: pd.to_numeric(s, errors="coerce")
    a, w = edges[edges[part] == "assertions"], gaps[gaps[part] == "assertions"]
    kept = pd.DataFrame({"target": a.subject, "condition": a.object, "source": a[src], "edge": a.id, "withheld": False, "rule": "",
                         "papers": np.nan, "lift": np.nan, "cut": np.nan, "evidence": a[f"{P}:source_evidence"]})
    held = pd.DataFrame({"target": w.subject, "condition": w[f"{P}:withheld_object"], "source": w[src], "edge": w.id, "withheld": True, "rule": w[f"{P}:rule"],
                         "papers": num(w[f"{P}:papers_mentioning_both"]), "lift": num(w[f"{P}:lift"]), "cut": num(w[f"{P}:lift_cut"]), "evidence": ""})
    x = edges[edges[part] == "expression"]
    d = edges[edges[part] == "drugs"]
    dt, di = d[d[src] == "drug_target"], d[d[src] == "drug_indication"]
    other = nodes[~nodes.category.isin(["biolink:Protein", f"{P}:KnowledgeGap"]) & ~nodes.id.isin(traits)]
    return {
        "context_nodes": pd.DataFrame({"id": other.id, "kind": np.where(other.category == "biolink:Drug", "drug", "condition"), "category": other.category, "name": other.name,
                                       "grounded": other.id_grounded.astype(str).str.lower() == "true", "xref": other.xref}).sort_values("id"),
        "assertions": pd.concat([kept, held], ignore_index=True).sort_values(["target", "condition", "source"]),
        "expression": pd.DataFrame({"target": x.subject, "condition": x.object, "label": x[f"{P}:context"], "log2fc": num(x[f"{P}:effect_size"]), "p": num(x[f"{P}:p_value"]),
                                    "edge": x.id}).sort_values(["target", "condition", "label"]),
        "drug_targets": pd.DataFrame({"drug": dt.subject, "target": dt.object, "action": dt[f"{P}:action_type"], "mechanism": dt[f"{P}:mechanism_of_action"] == "true",
                                      "activity": dt[f"{P}:activity"], "source": dt[f"{P}:activity_source"], "edge": dt.id}).sort_values(["drug", "target"]),
        "drug_indications": pd.DataFrame({"drug": di.subject, "condition": di.object, "edge": di.id}).sort_values(["drug", "condition"]),
    }


def context_files(release: str, out: Path) -> tuple[list[str], dict | None]:
    """Write what the questions read of the merged graph's other parts, where `kg` has built it, and
    return the files written and what they hold. The graph is read only if it is the one its record describes."""
    g = paths.GRAPH / release
    rec = g / "graph.json"
    if not rec.is_file():
        return [], None
    man = json.loads(rec.read_text())
    for k, want in man["files"].items():
        if graph.file_sha(g / f"{k}.tsv") != want:
            raise SystemExit(f"REFUSED: {g / (k + '.tsv')} is not the file {rec} describes; build the merged graph again:  python -m pharos_gap_kg.kg build {release}")
    rd = lambda k: pd.read_csv(g / f"{k}.tsv", sep="\t", dtype=str, keep_default_na=False)
    traits = set(pd.read_parquet(out / "nodes.parquet", columns=["id", "kind"]).query("kind == 'trait'").id)
    T = context_tables(rd("nodes"), rd("edges"), rd("gaps"), traits)
    for name, df in T.items():
        df.reset_index(drop=True).to_parquet(out / f"{name}.parquet", index=False)
    a = T["assertions"]
    what = {"merged_graph": man["files"], "conditions": int((T["context_nodes"].kind == "condition").sum()), "drugs": int((T["context_nodes"].kind == "drug").sum()),
            "assertions": int((~a.withheld).sum()), "withheld": int(a.withheld.sum()), "sources": int(a.source.nunique()), "expression": len(T["expression"]),
            "drug_targets": len(T["drug_targets"]), "drug_indications": len(T["drug_indications"]), "graph_edges": man["edges"], "graph_gaps": man["gaps"], "graph_nodes": man["nodes"]}
    return list(CONTEXT_FILES), what


def size_files(release: str, out: Path) -> tuple[list[str], dict | None]:
    """Write the size of each study behind the bundle's links, where `studysize` has read them from
    the GWAS Catalog, and return the file written and what it holds."""
    rec = paths.DOCS_DATA / f"study_sizes_{release}.json"
    src = paths.TABLES / "gwas_catalog_20261005" / "study_sizes.parquet"
    if not (rec.is_file() and src.is_file()):
        return [], None
    r = json.loads(rec.read_text())
    if graph.file_sha(src) != r["file_sha256"]:
        raise SystemExit(f"REFUSED: {src} is not the table {rec.name} describes; run:  python -I -m pharos_gap_kg.studysize build {release}")
    mine = set(pd.read_parquet(out / "evidence.parquet", columns=["study"]).study)
    t = pd.read_parquet(src)
    t = t[t.study.isin(mine) & t.cases.notna()][["study", "cases", "controls", "individuals"]].sort_values("study").reset_index(drop=True)
    t.astype({"cases": "int64", "controls": "int64", "individuals": "int64"}).to_parquet(out / SIZE_FILES[0], index=False)
    return list(SIZE_FILES), {"studies": len(mine), "with_a_size": len(t), "catalog": r["catalog"], "catalog_sha256": r["catalog_sha256"]}


def add_sizes(release: str) -> dict:
    """Add the study sizes to a bundle that is already built, and record them in its manifest."""
    out = paths.bundle(release)
    mf = out / "manifest.json"
    if not mf.is_file():
        raise SystemExit(f"no bundle at {out}; build one with:  python -m pharos_gap_kg.bundle {release}")
    files, what = size_files(release, out)
    if not files:
        raise SystemExit(f"no study sizes to add; run:  python -I -m pharos_gap_kg.studysize build {release}")
    man = json.loads(mf.read_text())
    man["files"].update({f: graph.file_sha(out / f) for f in files})
    man["study_sizes"] = what
    mf.write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"bundle_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def add_context(release: str) -> dict:
    """Add the merged graph's other parts to a bundle that is already built, and record them in its manifest."""
    out = paths.bundle(release)
    mf = out / "manifest.json"
    if not mf.is_file():
        raise SystemExit(f"no bundle at {out}; build one with:  python -m pharos_gap_kg.bundle {release}")
    files, what = context_files(release, out)
    if not files:
        raise SystemExit(f"no merged graph to read; build it with:  python -m pharos_gap_kg.kg build {release}")
    man = json.loads(mf.read_text())
    man["files"].update({f: graph.file_sha(out / f) for f in files})
    man["context"] = what
    mf.write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"bundle_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


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
    exported = gap_edge_ids(release, t)
    for kind in ("dark", "ligand"):
        t[f"{kind}_gap_id"] = t.id.map(exported[kind]) if exported else None
    t.to_parquet(out / GAP_FILES[0], index=False)
    dg, lg = s["dark_gap"], s["ligand_gap"]
    m = re.fullmatch(r"pharos(\d)(\d+)", release)
    rates = {"earlier": f"Pharos {m.group(1)}.{m.group(2)}" if m else release, "later": "Pharos 4.0.0 (May 2026)", "summary_sha256": graph.file_sha(summary),
             "graph": exported["graph"] if exported else None,
             "dark": {"targets": dg["with_a_gap"], "in_both_releases": dg["in_both_releases"], "promoted": dg["promoted"],
                      "bands": {b: {"promoted": v, "fell_in_the_mirror": dg["mirror"]["bands"][b]} for b, v in dg["bands"].items()},
                      "route_cost_auc": dg["scores"]["route cost"]["auc"], "best_single_count": dg["best_single_count"],
                      "best_single_count_auc": dg["scores"][dg["best_single_count"]]["auc"],
                      "named_criterion_cleared": dg["proposal"]["named_criterion_cleared"], "by_naming_at_random": dg["proposal"]["by_naming_at_random"]},
             "ligand": {"candidates": lg["candidates"], "gained": lg["gained"]}}
    (out / GAP_FILES[1]).write_text(json.dumps(rates, indent=1) + "\n")
    return list(GAP_FILES)


def gap_edge_ids(release: str, t: pd.DataFrame) -> dict | None:
    """Where the gaps have been exported as graph edges (`python -m pharos_gap_kg.tgapgraph <release>`),
    the id of each target's gap edge of each kind. The export and this table must hold the same gaps,
    or an answer would cite an edge that says something else."""
    from pharos_gap_kg import tgapgraph
    d, record = paths.GRAPH / f"{release}_target_gaps", paths.DOCS_DATA / f"graph_{release}_target_gaps.json"
    if not (d / "gaps.tsv").is_file() or not record.is_file():
        return None
    man = json.loads(record.read_text())
    if graph.file_sha(d / "gaps.tsv") != man["files"]["gaps"]:
        raise SystemExit(f"REFUSED: {d / 'gaps.tsv'} is not the file {record.name} describes; export the target-level gaps again")
    g = pd.read_csv(d / "gaps.tsv", sep="\t", dtype=str, keep_default_na=False, usecols=["id", "subject", f"{graph.PREFIX}:gap_type"])
    out = {kind: dict(zip(x.subject, x.id)) for kind, name in (("dark", tgapgraph.DARK), ("ligand", tgapgraph.LIGAND))
           for x in [g[g[f"{graph.PREFIX}:gap_type"] == name]]}
    want = {"dark": set(t.id[t.flips.notna()]), "ligand": set(t.id[t.lacks_ligand == True])}        # noqa: E712
    for kind in want:
        if set(out[kind]) != want[kind]:
            raise SystemExit(f"REFUSED: the exported graph holds {len(out[kind]):,} {kind} gaps and the table {len(want[kind]):,}, or they are for other targets; "
                             f"run:  python -m pharos_gap_kg.tgapgraph {release}")
    out["graph"] = {"directory": d.name, "gaps_sha256": man["files"]["gaps"], "gaps": man["gaps"], "conformance": man["conformance"]}
    return out


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
    out = {"headline": head, "calibration": {k_: c[k_] for k_ in ("cutoff", "data_ends", "base_rate", "by_trait", "by_target")}}
    out.update(ranking() or {})
    named = sources_then_gwas(c["cutoff"], c["split"])
    if named:
        out["sources_then_gwas"] = named
    return out


def sources_then_gwas(cutoff: int, split: dict) -> dict | None:
    """For each disease source: of the pairs it names that had no GWAS link by the cutoff, the share that
    were GWAS-linked later. Read from the context audit of the same split (findings/CONTEXT.md). The
    source that is the GWAS data again is left out."""
    f = paths.DOCS_DATA / f"context_gene_trait_{cutoff}.json"
    if not f.is_file():
        return None
    a = json.loads(f.read_text())
    if a["split"] != split:
        raise SystemExit("REFUSED: the context audit and the calibration are not from one split")
    a = a["audit"]
    by = {k: {"named_without_a_link": v["candidates_in_source"], "linked_later": v["positives_in_source"], "share": v["precision"]}
          for k, v in a["sources"].items() if k not in ("any source", ctxtables.GWAS_AGAIN)}
    return {"cutoff": cutoff, "disease_traits": a["traits"], "base_rate": a["base_rate"], "by_source": by,
            "reading": "The sources carry no dates, so a source may have named a pair because of the very studies that linked it later."}


def ranking() -> dict | None:
    """The scorer the front end ranks by since 11 Oct 2026, how it did on the prospective test, and its
    calibration by rank there. Read from the committed summaries; None where they are absent, and the
    front end then ranks by neighbours alone."""
    cp, pp = paths.DOCS_DATA / "calibration_prospect_2022.json", paths.DOCS_DATA / "prospect_gene_trait_2022.json"
    if not (cp.is_file() and pp.is_file()):
        return None
    c, p = json.loads(cp.read_text()), json.loads(pp.read_text())
    if c["split"] != p["split"]:
        raise SystemExit("REFUSED: the calibration and the scoring of the prospective test are not from one split")
    name = c["scorer"]
    part = lambda term: re.fullmatch(r"log\(1 \+ (B\d)\)", term).group(1)
    d = p["paired"][f"{name}-B3"]
    three = lambda v: [d[v]["difference"], d[v]["low"], d[v]["high"]]
    scorer = {"name": name, "weights": {part(k): v for k, v in zip(c["terms"], c["weights"])}, "each_part_enters_as": "log(1 + part)",
              "parts": {"B3": "shared-trait neighbours", "B1": "traits the target is already linked to", "B5": "genes within 500,000 bases already linked to the trait"},
              "weights_fitted_on_links_to": c["cutoff"], "served_by_target": c["served_by_target"]}
    head = {"cutoff": c["cutoff"], "first_year": c["first_year"], "data_ends": c["data_ends"], "measure": "mean average precision over traits",
            "traits_scored": c["traits_scored"], "test_links": c["test_links"], "combiner": p["results"][name]["full"][0], "neighbours": p["results"]["B3"]["full"][0],
            "degree": p["results"]["B1"]["full"][0], "combiner_minus_neighbours_all": three("full"), "combiner_minus_neighbours_unexplained_lower_bound": three("low"),
            "combiner_minus_neighbours_unexplained_neutral": three("neutral"), "combiner_minus_neighbours_by_locus": three("locus"),
            "share_of_later_links_explained_by_wording_or_position": p["counts"]["positives_lineage_or_position_linked"] / p["counts"]["test_links"],
            "by_target": c["by_target_against_neighbours"], "split": c["split"]}
    cal = {**{k_: c[k_] for k_ in ("cutoff", "first_year", "data_ends", "base_rate", "by_trait", "by_target")}, "prospective": True, "scorer": name}
    out = {"scorer": scorer, "prospective": head, "calibration": cal}
    if "pooled_over" in c:
        cal["pooled_over"] = c["pooled_over"]          # every trait and target ranked, not only those that had a later link
    nl = no_link(c)
    if nl is not None:
        out["no_link"] = nl
    return out


def no_link(c: dict) -> dict | None:
    """How the position rankings for targets with no link have held: the cold-start scoring of the
    prospective split and of the two earlier periods, and the calibration by rank. None where the
    summaries are absent, and the front end then ranks no target that has no link."""
    pp, ep = paths.DOCS_DATA / "prospect_coldstart_2022.json", paths.DOCS_DATA / "coldstart_gene_trait_2018.json"
    if "no_link" not in c or not (pp.is_file() and ep.is_file()):
        return None
    p, e = json.loads(pp.read_text()), json.loads(ep.read_text())
    if p["split"] != c["split"]:
        raise SystemExit("REFUSED: the cold-start scoring and the calibration of the prospective test are not from one split")
    three = lambda d: [d["difference"], d["low"], d["high"]]
    g, b, k = p["new_targets"], p["by_new_target"], c["no_link"]
    head = {"cutoff": p["cutoff"], "first_year": k["first_year"], "data_ends": k["data_ends"], "later_links": g["positives"], "later_links_listed": k["positives_listed"],
            "targets_with_no_link": k["targets_with_no_link"], "of_them_listed_for_some_trait": k["of_them_listed_for_some_trait"],
            "by_trait": {"traits": g["units"], "position": g["results"]["P3"]["all"]["value"][0], "tie_breaker_alone": g["results"]["P1"]["all"]["value"][0],
                         "chance": g["results"]["P0"]["all"]["value"][0], "position_minus_tie_breaker": three(g["paired"]["P3-P1"]["all"])},
            "by_target": {"targets": b["units"], "position": b["results"]["Q3"]["all"]["value"][0], "tie_breaker_alone": b["results"]["Q1"]["all"]["value"][0],
                          "chance": b["results"]["Q0"]["all"]["value"][0], "position_minus_tie_breaker": three(b["paired"]["Q3-Q1"]["all"])},
            "earlier_periods": {y: {"by_trait": three(e["periods"][y]["new_targets"]["paired"]["P3-P1"]["all"]), "by_target": three(e["periods"][y]["by_new_target"]["paired"]["Q3-Q1"]["all"])}
                                for y in sorted(e["periods"])}}
    cal = {k_: k[k_] for k_ in ("bins", "by_trait", "by_target", "first_year", "data_ends", "base_rate_among_listed", "base_rate_among_all", "listed_pairs")}
    return {"headline": head, "calibration": cal}


def refresh_backtest(release: str) -> dict:
    """Write backtest.json again in a bundle that is already built, and record it in the manifest."""
    out = paths.bundle(release)
    mf = out / "manifest.json"
    if not mf.is_file():
        raise SystemExit(f"no bundle at {out}; build one with:  python -m pharos_gap_kg.bundle {release}")
    bt = backtest()
    if bt is None:
        raise SystemExit("no backtest summaries to write")
    (out / "backtest.json").write_text(json.dumps(bt, indent=1) + "\n")
    man = json.loads(mf.read_text())
    man["files"]["backtest.json"] = graph.file_sha(out / "backtest.json")
    man["backtest"] = backtest_label(bt)
    mf.write_text(json.dumps(man, indent=1) + "\n")
    (paths.DOCS_DATA / f"bundle_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def backtest_label(bt: dict | None) -> str:
    if not bt:
        return "none"
    s = f"split frozen at {bt['headline']['cutoff']}"
    if "scorer" in bt:
        s += f"; ranked by {bt['scorer']['name']}, tested on links first reported {bt['prospective']['first_year']} to {bt['prospective']['data_ends']}"
    return s


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
        "backtest": backtest_label(bt),
        "files": {f: graph.file_sha(out / f) for f in files},
    }
    if GAP_FILES[0] in files:
        man["target_gaps"] = int(len(pd.read_parquet(out / GAP_FILES[0])))
    cfiles, what = context_files(release, out)
    if cfiles:
        man["files"].update({f: graph.file_sha(out / f) for f in cfiles})
        man["context"] = what
    sfiles, what = size_files(release, out)
    if sfiles:
        man["files"].update({f: graph.file_sha(out / f) for f in sfiles})
        man["study_sizes"] = what
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
    ck("a bundle says what it ranks by only when a scorer is recorded",
       backtest_label(None) == "none" and backtest_label({"headline": {"cutoff": 2018}}) == "split frozen at 2018"
       and "ranked by x, tested on links first reported 2023 to 2026" in backtest_label({"headline": {"cutoff": 2018}, "scorer": {"name": "x"}, "prospective": {"first_year": 2023, "data_ends": 2026}}))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) == 2 and argv[1] == "--gaps-only":
        man = add_gaps(argv[0])
        print(f"{man['release']}: {man['target_gaps']:,} target-level gaps added to {paths.bundle(man['release'])}")
        return 0
    if len(argv) == 2 and argv[1] == "--sizes-only":
        man = add_sizes(argv[0])
        print(f"{man['release']}: the size of {man['study_sizes']['with_a_size']:,} of {man['study_sizes']['studies']:,} studies added to {paths.bundle(man['release'])}")
        return 0
    if len(argv) == 2 and argv[1] == "--context-only":
        man = add_context(argv[0])
        c = man["context"]
        print(f"{man['release']}: {c['assertions']:,} assertions and {c['withheld']:,} withheld rows from {c['sources']} sources, {c['expression']:,} expression observations, "
              f"{c['drugs']:,} drugs ({c['drug_targets']:,} target and {c['drug_indications']:,} indication edges), {c['conditions']:,} more conditions added to {paths.bundle(man['release'])}")
        return 0
    if len(argv) == 2 and argv[1] == "--backtest-only":
        man = refresh_backtest(argv[0])
        print(f"{man['release']}: backtest.json written again; {man['backtest']}")
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
