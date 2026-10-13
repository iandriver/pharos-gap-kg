"""The merged graph, the bundle the questions read and the earlier records say the same thing.

Four artifacts hold the same rows: the merged graph's files, its record, the bundle's tables, and the
records of the three earlier exports and of the pruning round. Where two of them can be counted
against each other they are, and an id the questions cite must be an id in the graph.

The record-against-record checks run anywhere. The file checks are skipped where the graph or the
bundle has not been built."""
import json

import pandas as pd
import pytest

from pharos_gap_kg import graph, kg, paths

R = "pharos319"
REC = paths.DOCS_DATA / f"graph_{R}.json"
G_DIR, B_DIR = paths.GRAPH / R, paths.bundle(R)
J = lambda f: json.loads((paths.DOCS_DATA / f).read_text())
built = pytest.mark.skipif(not ((G_DIR / "graph.json").is_file() and (B_DIR / "assertions.parquet").is_file()), reason="the merged graph or the bundle's tables are not built")
pytestmark = pytest.mark.skipif(not REC.is_file(), reason="no record of a merged graph")


def test_the_record_agrees_with_the_earlier_records():
    m, ctx, tg, applied, b = J(f"graph_{R}.json"), J(f"graph_{R}_context_disease.json"), J(f"graph_{R}_target_gaps.json"), J(f"pruning2_applied_{R}.json"), J(f"bundle_{R}.json")
    assert m["conformance"] == "10 of 10 checks passed" and m["merged"] == 0 and m["rows"] == m["edges"] - m["edges_by_part"]["ontology"] + m["gaps"]
    assert m["edges_by_part"]["links"] == b["links"] == m["links"]["rows"] and (m["links"]["first_year"], m["links"]["latest_year"]) == (b["first_year"], b["latest_year"])
    # the nine sources: the disease-source export less the two remodelled ones, row for row
    two = ctx["edges_by_source"]["Expression Atlas"] + ctx["edges_by_source"]["DrugCentral Indication"]
    assert m["assertions"]["rows"] == ctx["rows"] - two and m["assertions"]["edges"] == ctx["edges"] - two
    assert m["assertions"]["gaps"] == ctx["gaps"] == applied["withheld"] and m["assertions"]["gaps_by_source"] == ctx["gaps_by_source"] and m["assertions"]["gaps_by_rule"] == ctx["gaps_by_rule"]
    assert {k: v for k, v in ctx["edges_by_source"].items() if k not in kg.REMODELLED} == m["assertions"]["edges_by_source"]
    assert m["target_gaps"]["gaps"] == tg["gaps"] and m["target_gaps"]["by_type"] == tg["gaps_by_type"] and m["gaps"] == ctx["gaps"] + tg["gaps"]
    # what the protocol's note of 13 Oct quotes
    d, x = m["drugs"], m["expression"]
    assert (d["indication_rows"], d["through_the_mechanism_target"], d["through_another_activity"], d["indication_rows_without_an_id"]) == (13919, 10327, 3592, 7738)
    assert d["through_the_mechanism_target"] + d["through_another_activity"] == d["indication_rows"]
    assert (x["rows"], x["conditions_as_named"]) == (159846, 142) and x["smallest_fold_change"] >= 1.0 > x["floor"] and x["largest_p"] < 0.05 and x["edges"] == x["rows"] == x["increased"] + x["decreased"]
    t = m["traits_that_reach_another_part"]
    assert (m["ontology"]["traits"], t["under_their_own_id"]) == (b["traits"], 41) and t["joined"] == t["under_their_own_id"] + t["through_their_current_term"]
    assert m["ontology"]["same_as_edges"] + m["ontology"]["replacements_with_no_declared_kind"] == m["ontology"]["traits_retired_and_replaced"]


@built
def test_the_files_are_the_ones_the_record_describes_and_the_bundle_was_read_from_them():
    m, b = J(f"graph_{R}.json"), json.loads((B_DIR / "manifest.json").read_text())
    assert {k: graph.file_sha(G_DIR / f"{k}.tsv") for k in ("nodes", "edges", "gaps")} == m["files"] == b["context"]["merged_graph"]
    c = b["context"]
    assert (c["assertions"], c["withheld"], c["expression"]) == (m["assertions"]["edges"], m["assertions"]["gaps"], m["expression"]["edges"])
    assert (c["drug_targets"], c["drug_indications"], c["drugs"], c["sources"]) == (m["drugs"]["drug_target_edges"], m["drugs"]["drug_indication_edges"], m["drugs"]["drugs"], 9)


@built
def test_every_id_the_questions_can_cite_is_in_the_graph():
    part = f"{graph.PREFIX}:part"
    e = pd.read_csv(G_DIR / "edges.tsv", sep="\t", dtype=str, keep_default_na=False, usecols=["id", "subject", "object", "predicate", part])
    g = pd.read_csv(G_DIR / "gaps.tsv", sep="\t", dtype=str, keep_default_na=False, usecols=["id", part])
    n = pd.read_csv(G_DIR / "nodes.tsv", sep="\t", dtype=str, keep_default_na=False, usecols=["id", "category"])
    assert e.id.is_unique and g.id.is_unique and n.id.is_unique and not set(e.id) & set(g.id)
    rd = lambda f: pd.read_parquet(B_DIR / f)
    a, x, dt, di, links = rd("assertions.parquet"), rd("expression.parquet"), rd("drug_targets.parquet"), rd("drug_indications.parquet"), rd("links.parquet")
    by = lambda p: set(e[e[part] == p].id)
    assert set(links.edge_id) == by("links"), "the dated links of the merged graph are the bundle's links, under the same ids"
    assert set(a[~a.withheld].edge) == by("assertions") and set(a[a.withheld].edge) == set(g[g[part] == "assertions"].id)
    assert set(x.edge) == by("expression") and set(dt.edge) | set(di.edge) == by("drugs")
    nodes = set(n.id)
    assert set(a.target) <= nodes and set(x.condition) <= nodes and set(dt.drug) <= nodes and set(di.condition) <= nodes
    assert set(a[~a.condition.str.startswith("NCBIGene:")].condition) <= nodes, "a condition named only by a withheld row is still a node"
    # a DrugCentral indication is a path, and is no longer an edge from target to disease
    assert not (e.predicate == "biolink:target_for").any() and set(e[e.predicate == "biolink:treats"].subject) <= set(n[n.category == "biolink:Drug"].id)
    same = e[e.predicate == "biolink:same_as"]
    assert len(same) == J(f"graph_{R}.json")["ontology"]["same_as_edges"] and set(same.subject) <= set(links.trait)


@built
def test_names_still_lead_to_one_thing():
    """Naming conditions and drugs must not make a gene symbol or a trait name ambiguous. The one symbol
    that is also, letter for letter, the name of a condition is listed here; a second one fails the build."""
    from pharos_gap_kg import ask, others
    G = ask.Graph(B_DIR)
    O = G.others
    clash = sorted({s for s in G.t["name"] if others.norm(s) in O._cname or others.norm(s) in O._dname})
    lost = [s for s in clash if not G.resolve(s)["ok"] or G.resolve(s)["kind"] != "target"]
    assert lost == ["HBD"], lost
    assert G.resolve("ACHE")["kind"] == "target" and G.resolve("Ache")["kind"] == "condition"
    exact = [n for n in G.d["name"] if others.norm(n) in O._cname]
    assert not exact, "a condition with a GWAS trait's name is read as that trait, not offered beside it"
    c, b = G.status()["other_sources"], json.loads((B_DIR / "manifest.json").read_text())["context"]
    assert (c["assertions"], c["withheld"], c["drugs"], c["expression_observations"], c["conditions"]) == (b["assertions"], b["withheld"], b["drugs"], b["expression"], b["conditions"])
