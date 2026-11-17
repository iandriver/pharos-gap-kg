"""The walkthrough page's data, against the records it is built from: current, and holding no cell of a study's file."""
import json

import pytest

from pharos_gap_kg import add, paths, web

F = web.HERE / "walkthrough_pharos319.json"
W = json.loads(F.read_text())
HAVE = all(add.record_path("pharos319", s, r).is_file() for s, r, _ in web.RUNS) and (add.run_dir("pharos319", "D") / "card1.json").is_file()


@pytest.mark.skipif(not HAVE, reason="the records of the runs, or the runs' cards, are not here")
def test_the_committed_walkthrough_is_the_one_the_records_make():
    assert web.walkthrough("pharos319") == W, "rebuild it:  python -m pharos_gap_kg.web walkthrough pharos319"


def test_the_walkthrough_carries_sheet_names_and_sizes_and_no_cell():
    assert [s["id"] for s in W["studies"]] == ["A", "B", "B_rerun", "C", "C_rerun", "D"]
    for s in W["studies"]:
        assert all(set(x) == {"name", "rows", "columns"} for x in s["sheets"]), s["id"]
    d = W["studies"][-1]
    replay = json.loads((paths.DOCS_DATA / "cheng_reference_replay_pharos319.json").read_text())
    assert replay["full"]["rows"] == replay["bundle"]["accounted"] == 539
    assert replay["full"]["graph_genes"] == d["join"]["targets"]
    assert f"{replay['bundle']['graph_genes']} genes" in d["reference_note"]
    assert all(row["symbol"] in d["reference_note"] for row in replay["resolver_differences"])
    assert d["join"]["rows"] == 500 and "disagree" in d["outcome"] and "at most 68%" in d["outcome"] and not d["joined"]
    assert W["studies"][3]["outcome"].startswith("Nothing joins") and "95 of 102 rows, 93.1%" in W["studies"][3]["outcome"]
    assert W["studies"][1]["check"]["against"].startswith("the genes the GWAS Catalog reports")
