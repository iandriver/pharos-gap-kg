"""The subgraph behind an answer: its rows are the export's, and the reference validator finds in it only
what it finds in the whole graph."""
import json

import pytest

from pharos_gap_kg import ask, paths, subgraph

R = json.loads((paths.DOCS_DATA / "subgraph_check_pharos319.json").read_text())
G = json.loads((paths.DOCS_DATA / "graph_pharos319.json").read_text())
S = json.loads((paths.DOCS_DATA / "subgraph_standards_pharos319.json").read_text())


def test_every_row_the_bundle_gives_is_the_exports():
    c = R["every_row"]["counts"]
    assert R["n_differences"] == 0 and R["rows_compared"] > 100
    assert c["edges"]["here"] == G["edges"] - G["predicates"]["biolink:same_as"] and c["edges"]["in_export"] == G["edges"]
    assert c["gaps"]["here"] == c["gaps"]["in_export"] == G["gaps"]
    assert all(v == 0 for k, v in {**c["edges"], **c["gaps"]}.items() if k.startswith("differ_"))


def test_the_reference_validator_finds_only_the_graphs_own_causes():
    assert S["checked"]["edges"] > 100 and S["biolink_model"] == "4.4.3"
    assert {e["prefix"] for e in S["errors"]} <= {"pharosgapkg", None}
    assert {e["type"] for e in S["errors"]} <= {"INVALID_NODE_PROPERTY_VALUE", "INVALID_EDGE_PROPERTY_VALUE", "INVALID_CATEGORY"}


def test_an_answer_gives_the_edges_it_cites():
    bundle = ask.default_bundle()
    if not (bundle / "assertions.parquet").is_file():
        pytest.skip("the bundle is not built")
    a = ask.Graph(bundle).why("APOE", "Alzheimer disease")
    s = subgraph.of(a, bundle, limit=None)
    ids = {e["id"] for e in s["edges"]}
    assert a["edge"] in ids and not s["graph"]["counts"]["cited_ids_not_found"]
    assert {e["subject"] for e in s["edges"]} | {e["object"] for e in s["edges"]} <= {n["id"] for n in s["nodes"]}
