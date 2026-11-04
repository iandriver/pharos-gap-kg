"""Adding a table: the records of the two studies and the control, read against the protocol, the tables they
describe, and the MCP tools' verdicts on the same proposals.

The files are kept outside the repository (study A is CC BY-NC-ND), so the tests that read them skip where they
are absent."""
import json
import re

import pytest

from pharos_gap_kg import add, paths, registry

R = {k: paths.DOCS_DATA / f"add_{s['key']}_pharos319.json" for k, s in add.STUDIES.items()}
R["B_rerun"] = paths.DOCS_DATA / "add_kim2023_pd_supp_rerun_pharos319.json"
CTRL = paths.DOCS_DATA / "add_uffelmann2026_supp_control_pharos319.json"
CMP = paths.DOCS_DATA / "add_compare_pharos319.json"
if not all(p.is_file() for p in (*R.values(), CTRL, CMP)) or paths.ROOT is None:
    pytest.skip("no records of the add stage", allow_module_level=True)
A, B, C, K = (json.loads(p.read_text()) for p in (R["A"], R["B"], CTRL, CMP))
BR, SC = json.loads(R["B_rerun"].read_text()), json.loads(R["C"].read_text())
F = re.sub(r"\s+", " ", (paths.ROOT / "findings" / "ADD.md").read_text().split("## Added 4 Nov 2026")[1])
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


def test_study_b_again_finds_the_genes_and_is_not_joined():
    b, k = BR["back_to_the_sheet"], K["B_rerun"]
    assert BR["rerun"].startswith("a rerun after the registered result") and BR["back_to_the_sheet_rule"]
    assert [t["accepted"] for t in BR["stage1"]["tries"] + BR["stage2"]["tries"]] == [True, True] and BR["stage2"]["tries"][0]["proposal"]["map"] is False
    assert BR["sheet"] == "Supplementary Table 3 - Lead SN" and [t["accepted"] for t in b["stage2"]["tries"]] == [False, False, False, True]
    assert (BR["join"]["rows"], BR["join"]["targets"], k["rows"], k["a_table_of_the_studys_loci"], k["reading"]) == (191, 131, 303, False, "partial")
    assert "joined" not in BR
    s1, s2 = k["share_of_table_genes_the_catalog_maps"], k["share_of_catalog_genes_the_table_holds"]
    for w in (f"joins {BR['join']['rows']} links from {BR['join']['targets']} graph genes", f"({k['of_them_graph_genes']} graph genes), {k['shared']} genes are shared",
              f"{s1[0]:.0%} of the table's genes [{s1[1]:.0%}, {s1[2]:.0%}]", f"{s2[0]:.0%} of the Catalog's [{s2[1]:.0%}, {s2[2]:.0%}]",
              f"The sheet's {k['rows']} rows lie outside {k['rows_band'][0]} to {k['rows_band'][1]}"):
        assert w in F, f"findings/ADD.md (4 Nov) no longer says: {w}"


def test_study_c_meets_check_8_with_ungrounded_names_and_joins_nothing():
    k = K["C"]
    assert SC["sheet"] == "Supplementary Table 2" and k["a_table_of_the_studys_loci"] and k["rows"] == 102 and "back_to_the_sheet" not in SC
    t = SC["stage2"]["tries"]
    assert [x["accepted"] for x in t] == [False, False, False, True]
    assert t[2]["problems"] == ["8_identifiers: 95 of 102 rows (93.1%) resolve as they must (both sides, or the subject of a row the condition makes a gap); at least 95% must"]
    assert t[2]["proposal"]["subject"]["resolver"] == "gene_symbol" and t[3]["proposal"]["subject"]["resolver"] == "name:gene"
    assert SC["join"]["rows"] == 0 and not k["joined_a_table"] and k["catalog_mapped_genes"] == 0 and k["reading"].startswith("not computable")
    assert SC["spent_usd"] <= add.CAP_USD
    for w in ("95 of 102 rows (93.1%) resolve", "Code accepted it, because check 8 counts such a node as resolved", "`add join` refuses it"):
        assert w in F, f"findings/ADD.md (4 Nov) no longer says: {w}"


@need
def test_the_tools_give_the_commands_verdicts():
    from pharos_gap_kg import ask
    S = add.Session(ask.Graph(paths.bundle("pharos319")))
    for k, r in (("A", A), ("B", B), ("B", BR), ("C", SC)):
        s = add.STUDIES[k]
        f = str(registry.SOURCES[s["key"]].path)
        first = next(t["proposal"] for t in r["stage1"]["tries"] if t["accepted"])
        back = r.get("back_to_the_sheet") or {}
        for t in r["stage1"]["tries"] + (back.get("stage1") or {}).get("tries", []):
            v = add.tool_check(S, f, s["about"], "sheet", t["proposal"])
            assert (v["accepted"], v["problems"]) == (t["accepted"], t["problems"]), (k, "sheet", t["attempt"])
        for sheet, tries in ((first, r["stage2"]["tries"]), ({"sheet": r["sheet"], "header_row": r["header_row"]}, (back.get("stage2") or {}).get("tries", []))):
            for t in tries:
                v = add.tool_check(S, f, s["about"], "mapping", t["proposal"], sheet["sheet"], sheet["header_row"])
                assert (v["accepted"], v["problems"]) == (t["accepted"], t["problems"]), (k, t["stage"], t["model"], t["attempt"])
    v = add.tool_check(S, str(registry.SOURCES["uffelmann2026_supp"].path), add.STUDIES["A"]["about"], "mapping", A["mapping"], A["sheet"], A["header_row"])
    assert v["would_join"] == {"links": 181, "graph_genes": 181, "traits": ["MONDO:0004975"]}
