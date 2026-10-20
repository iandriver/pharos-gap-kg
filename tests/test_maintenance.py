"""The maintenance demonstration and the reproduction check, read against the traced map of the build."""
import json

from pharos_gap_kg import maintain, paths

J = lambda f: json.loads((paths.DOCS_DATA / f).read_text())


def test_the_stages_that_ran_are_the_stale_ones_and_wrote_every_changed_file():
    g, d = J("build_graph_pharos319.json")["stages"], J("maintenance_demo_pharos319.json")
    assert len(g) == d["stages"] == 20 and d["ran"] == list(d["stale"])
    writes = {s["name"]: set(s["writes"]) for s in g}
    outs = set().union(*(writes[n] for n in d["ran"]))
    assert set(d["files_changed"]) - outs == {d["change"]["file"]} and d["every_changed_file_is_an_output_of_a_stage_that_ran_or_the_change"]
    assert d["stale"] == {k: v for k, v in maintain.stale(g, {s["name"]: {"inputs": {f: "h" for f in s["reads"]}} for s in g},
                                                          digest=lambda f: "x" if f == d["change"]["file"] else "h").items()}
    assert set(maintain.whatif(g, d["change"]["file"])) == set(d["ran"])


def test_the_incremental_rebuild_is_the_full_build_and_a_clean_clone_is_this_checkout():
    c, r = J("maintenance_check_pharos319.json"), J("reproduce_check_pharos319.json")
    n = len({f for s in J("build_graph_pharos319.json")["stages"] for f in s["reads"] + s["writes"]})
    assert c["files_compared"] == r["files_compared"] == n == 112 and c["files_that_differ"] == []
    assert r["files_that_differ"] == ["data:tables/pharos319/_manifest.json"]


def test_a_clean_clone_after_the_remodel_is_this_checkout():
    r = J("reproduce_check_pharos319_e4e581f.json")
    assert r["files_compared"] == 112 and r["files_that_differ"] == ["data:tables/pharos319/_manifest.json"]
