"""Adding a table: the records of the two studies and the control, read against the protocol, the tables they
describe, and the MCP tools' verdicts on the same proposals.

The files are kept outside the repository (study A is CC BY-NC-ND), so the tests that read them skip where they
are absent."""
import json
import re

import pytest

from pharos_gap_kg import add, paths, registry

R = {k: paths.DOCS_DATA / f"add_{s['key']}_pharos319.json" for k, s in add.STUDIES.items()}
CTRL = paths.DOCS_DATA / "add_uffelmann2026_supp_control_pharos319.json"
CMP = paths.DOCS_DATA / "add_compare_pharos319.json"
if not all(p.is_file() for p in (*R.values(), CTRL, CMP)) or paths.ROOT is None:
    pytest.skip("no records of the add stage", allow_module_level=True)
A, B, C, K = (json.loads(p.read_text()) for p in (R["A"], R["B"], CTRL, CMP))
PROTO = re.sub(r"\s+", " ", (paths.ROOT / "docs" / "PROTOCOL.md").read_text())
HAVE = all(registry.SOURCES[s["key"]].path.is_file() for s in add.STUDIES.values()) and (paths.bundle("pharos319") / "manifest.json").is_file()
need = pytest.mark.skipif(not HAVE, reason="the studies' files or the bundle are not here")


def test_the_descriptions_are_the_protocols_word_for_word():
    for k, s in add.STUDIES.items():
        assert add.text(s["about"]) in PROTO, f"the description of study {k} is not the protocol's"
    assert A["description"] == add.text(add.STUDIES["A"]["about"]) and B["description"] == add.text(add.STUDIES["B"]["about"])


def test_the_control_passes_and_joins_later_pys_table():
    assert C["stage1"]["accepted"] and C["stage2"]["accepted"] and C["equals_later_py"] and C["join"]["rows"] == 181


def test_study_a_is_reproduced():
    a = K["A"]
    assert a["reproduced"] and a["hand_pairs"] == a["model_pairs"] == 181 and not a["only_the_model"] and not a["only_by_hand"]
    assert A["sheet"] == "SD2" and A["model"] == {"sheet": add.MODEL, "mapping": add.FALLBACK}
    assert [t["accepted"] for t in A["stage2"]["tries"]] == [False, False, True]
    assert A["spent_usd"] <= add.CAP_USD and B["spent_usd"] <= add.CAP_USD


def test_study_b_chose_a_sheet_with_no_gene_and_joined_nothing():
    assert B["stage1"]["accepted"] and B["sheet"] == "Supplementary Table 5 - Genomic" and B["join"] is None
    assert B["verdict"] == {"accepted": True, "map": False, "reason": B["mapping"]["reason"]} and "no gene column" in B["mapping"]["reason"]
    b = K["B"]
    assert b["catalog_associations"] == 121 and b["catalog_reported_genes"] == 0 and b["reading"].startswith("not computable")


@need
def test_the_tools_give_the_commands_verdicts():
    from pharos_gap_kg import ask
    S = add.Session(ask.Graph(paths.bundle("pharos319")))
    for k, r in (("A", A), ("B", B)):
        s = add.STUDIES[k]
        f = str(registry.SOURCES[s["key"]].path)
        for t in r["stage1"]["tries"]:
            v = add.tool_check(S, f, s["about"], "sheet", t["proposal"])
            assert (v["accepted"], v["problems"]) == (t["accepted"], t["problems"]), (k, "sheet", t["attempt"])
        for t in r["stage2"]["tries"]:
            v = add.tool_check(S, f, s["about"], "mapping", t["proposal"], r["sheet"], r["header_row"])
            assert (v["accepted"], v["problems"]) == (t["accepted"], t["problems"]), (k, "mapping", t["model"], t["attempt"])
    v = add.tool_check(S, str(registry.SOURCES["uffelmann2026_supp"].path), add.STUDIES["A"]["about"], "mapping", A["mapping"], A["sheet"], A["header_row"])
    assert v["would_join"] == {"links": 181, "graph_genes": 181, "traits": ["MONDO:0004975"]}
