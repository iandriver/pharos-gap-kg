"""The reference validator's record against the graph it read: every node and edge was read, and each
error is one of the four causes findings/STANDARDS.md names, in the numbers the graph's own files give."""
import json

import pandas as pd
import pytest

from pharos_gap_kg import paths

R = json.loads((paths.DOCS_DATA / "standards_pharos319.json").read_text())
G = json.loads((paths.DOCS_DATA / "graph_pharos319.json").read_text())


def test_every_node_and_edge_was_read_and_only_identifiers_fail():
    assert R["checked"] == {"nodes": G["nodes"], "edges": G["edges"]} and R["biolink_model"] == "4.4.3"
    assert {e["type"] for e in R["errors"]} == {"INVALID_NODE_PROPERTY_VALUE", "INVALID_EDGE_PROPERTY_VALUE", "INVALID_CATEGORY"}
    assert {e["prefix"] for e in R["errors"]} == {"pharosgapkg", "Orphanet", "AmyCo", None}
    assert R["nodes_with_an_error"] == 1042 and R["subject_object_pairs_with_an_error"] == 1062


def test_the_failing_identifiers_are_the_ones_the_graph_holds():
    f = paths.GRAPH / "pharos319" / "nodes.tsv"
    if not f.is_file():
        pytest.skip("the graph is not built")
    n = pd.read_csv(f, sep="\t", usecols=["id", "category"], dtype=str)
    pre = n.id.str.split(":").str[0]
    gap = n.category.str.contains("KnowledgeGap")
    node_errors = {e["prefix"]: e["records"] for e in R["errors"] if e["type"] == "INVALID_NODE_PROPERTY_VALUE" and e["prefix"]}
    assert node_errors == {"pharosgapkg": int(((pre == "pharosgapkg") & ~gap).sum()), "Orphanet": int((pre == "Orphanet").sum()), "AmyCo": int((pre == "AmyCo").sum())}
    assert int(gap.sum()) == 3 and node_errors["pharosgapkg"] + int(gap.sum()) == G["nodes_not_grounded"]
