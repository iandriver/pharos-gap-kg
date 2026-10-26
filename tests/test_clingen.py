"""The ClinGen audit's findings and protocol, read against its record."""
import json
import re

import pytest

from pharos_gap_kg import clingen, paths

REC = paths.DOCS_DATA / "clingen_audit_pharos319.json"
if not REC.is_file() or paths.ROOT is None:
    pytest.skip("no audit record", allow_module_level=True)
D = json.loads(REC.read_text())
F = re.sub(r"\s+", " ", (paths.ROOT / "findings" / "CLINGEN.md").read_text())
PROTO = re.sub(r"\s+", " ", (paths.ROOT / "docs" / "PROTOCOL.md").read_text())
pct = lambda x: f"{x[0]:.1%} [{x[1]:.1%}, {x[2]:.1%}]"
sig = lambda x: f"{x[0]:+.3f} [{x[1]:+.3f}, {x[2]:+.3f}]".replace("+-", "-")
n = lambda x: f"{x:,}"


def has(*wanted):
    for w in wanted:
        assert w in F, f"findings/CLINGEN.md no longer says: {w}"


def test_the_selftest_passes():
    assert clingen._selftest() == 0


def test_the_file_and_the_scope():
    c = D["clingen"]
    assert c["sha256"].startswith("b1fc3a11") and c["file_created"] == "2026-10-26"
    has(f"{n(c['curations'])} curations of {n(c['genes'])} genes", f"{n(c['pairs'])} gene-disease pairs",
        f"{n(c['pairs_by_class']['supported'])} supported", f"{c['pairs_by_class']['limited']} limited", f"{c['pairs_by_class']['against']} against",
        f"{n(D['rows'])} rows of the nine disease sources, {n(D['in_scope'])} are in scope", f"{n(D['in_scope_withheld'])} of those are withheld",
        f"{n(D['rows'] - D['in_scope'])} of the {n(D['rows'])} rows are out of scope")
    for w in (f"{n(c['curations'])} curations of {n(c['genes'])} genes", f"{n(D['in_scope'])} rows, {n(D['in_scope_withheld'])} of them withheld",
              f"{n(D['rows'] - D['in_scope'])} of the {n(D['rows'])} rows are out of scope"):
        assert w in PROTO, f"docs/PROTOCOL.md no longer says: {w}"
    withheld = {b["source"]: b["exact"]["in_scope"] for b in D["by_source"] if b["side"] == "withheld"}
    assert withheld == {"CTD": 946, "Monarch": 353, "JensenLab Knowledge UniProtKB-KW": 110, "JensenLab Knowledge MedlinePlus": 48, "UniProt Disease": 46}


def test_the_kept_rows():
    k = D["kept_rows_of_all_sources"]["exact"]
    has(f"ClinGen has judged {n(k['judged'])} of the {n(k['in_scope'])} kept rows in scope", f"supports {pct(k['agreement'])}",
        f"is against {pct(k['against_rate'])}", f"the rest, {k['limited']}, it rates limited", f"The {k['against']} kept rows ClinGen is against")
    for b in D["by_source"]:
        if b["side"] == "kept":
            e = b["exact"]
            has(f"| {b['source']} | {n(e['in_scope'])} | {n(e['judged'])} | {pct(e['agreement'])} | {pct(e['against_rate'])} |")
    has(f"| all nine | {n(k['in_scope'])} | {n(k['judged'])} | {pct(k['agreement'])} | {pct(k['against_rate'])} |")
    s = D["kept_rows_clingen_is_against"]["in_the_pruned_sources"]
    assert (s["rows"], s["not_searched_and_asserted_by_another_source"], s["searched_by_the_cut_and_kept"]) == (141, 126, 14)
    assert s["papers_naming_both_when_searched"][0] == 1 and s["papers_naming_both_when_searched"][2] == 64
    assert s["not_searched_and_asserted_by_no_other_source"] == ["SCN9A MONDO:0005027 (JensenLab Knowledge UniProtKB-KW)"]
    has("141 of the 299 come from the five pruned sources", "did not test 126", "searched 14 and kept them, with 1 to 64 papers",
        "SCN9A with epilepsy", "126 of the 141 are asserted by two or more sources")
    ex = {(e["source"], e["symbol"]) for e in D["kept_rows_clingen_is_against"]["examples"]}
    assert {("CTD", g) for g in ("BRAF", "KRAS", "MAP2K1", "MYH6", "KCNE2")} <= ex


def test_kept_against_withheld():
    w = D["withheld_rows_of_all_sources"]["exact"]
    has(f"ClinGen has judged {w['judged']} of the {n(w['in_scope'])} withheld rows in scope",
        f"Of the {w['judged']} it supports {w['supported']}, rates {w['limited']} limited and is against {w['against']}")
    E, Q = D["pruned_sources"]["exact"], D["pruned_sources"]["near"]
    W = {b["source"]: b["exact"] for b in D["by_source"] if b["side"] == "withheld"}
    for s in clingen.PRUNED:
        e, q = W[s], Q[s]
        if e["judged"]:
            assert f"| {s} | {n(e['in_scope'])} | {e['judged']}: {e['supported']} / {e['limited']} / {e['against']} | {sig(E[s]['agreement_kept_minus_withheld'])}" in F
            assert f"{sig(q['agreement_kept_minus_withheld'])} on {q['judged_withheld']} |" in F
        assert E[s]["verdict_of_its_own"] == (s == "UniProt Disease")
    P, PN = E["pooled"], Q["pooled"]
    has(f"| pooled, weighted by judged withheld rows | {n(D['in_scope_withheld'])} | {w['judged']} | {sig(P['agreement_kept_minus_withheld'])} | {sig(PN['agreement_kept_minus_withheld'])} |",
        f"| pooled without Monarch | | | {sig(P['agreement_without_monarch'])} | {sig(PN['agreement_without_monarch'])} |",
        f"judges {D['withheld_rows_of_all_sources']['near']['judged']} withheld rows",
        f"is {pct(D['kept_rows_of_all_sources']['exact']['confirmation'])} for kept rows and {pct(w['confirmation'])} for withheld ones")
    assert P["agreement_kept_minus_withheld"][1] > 0 and P["agreement_without_monarch"][1] > 0, "the reading: the pooled interval lies above 0"
    u = W["UniProt Disease"]
    assert (u["supported"], u["judged"]) == (10, 11) and E["UniProt Disease"]["agreement_kept_minus_withheld"][1] < 0 < E["UniProt Disease"]["agreement_kept_minus_withheld"][2]
    has("ClinGen supports 10 of the 11 rows the rule withheld from it")


def test_what_the_rule_withheld_that_clingen_supports():
    s = D["withheld_rows_clingen_supports"]
    b = s["by_source"]
    has(f"{sum(b.values())} withheld rows are pairs ClinGen supports: UniProt Disease {b['UniProt Disease']}, Monarch {b['Monarch']}, "
        f"JensenLab Knowledge UniProtKB-KW {b['JensenLab Knowledge UniProtKB-KW']}", f"For {s['no_paper_found']} of them the PubMed search found no paper")
    assert s["by_rule"] == {"literature_lift": 12, "uncorroborated": 2, "literature_lift|uncorroborated": 1}
    a = D["withheld_rows_with_no_paper_found"]
    has(f"found no paper for {n(a['all'])} of the {n(a['of'])} withheld rows, {a['in_scope']} of them in scope")
    ex = {e["symbol"] for e in s["examples"]}
    assert {"CFAP43", "KIF5A", "OGT", "PRKG1", "TRAPPC9"} <= ex
