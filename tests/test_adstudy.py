"""The Alzheimer study's findings and protocol, read against its record."""
import json
import re

import pytest

from pharos_gap_kg import adstudy, paths

REC = paths.DOCS_DATA / "adstudy_pharos319.json"
if not REC.is_file() or paths.ROOT is None:
    pytest.skip("no audit record", allow_module_level=True)
D = json.loads(REC.read_text())
F = re.sub(r"\s+", " ", (paths.ROOT / "findings" / "AD_STUDY.md").read_text())
PROTO = re.sub(r"\s+", " ", (paths.ROOT / "docs" / "PROTOCOL.md").read_text())
pct = lambda x, p=1: f"{x[0]:.{p}%} [{x[1]:.{p}%}, {x[2]:.{p}%}]"
sig = lambda x, p=3: f"{x[0]:+.{p}f} [{x[1]:+.{p}f}, {x[2]:+.{p}f}]".replace("+-", "-")
n = lambda x: f"{x:,}"


def has(*wanted):
    for w in wanted:
        assert w in F, f"findings/AD_STUDY.md no longer says: {w}"


def test_the_selftest_passes():
    assert adstudy._selftest() == 0


def test_the_study_the_trait_and_the_crosscheck():
    s, m = D["study"], D["matching"]
    assert s["sha256"].startswith("cb450c90") and (s["loci"], s["new_loci"]) == (127, 48)
    assert D["trait"]["id"] == "MONDO:0004975" and D["crosscheck"]["agree"]
    has(f"{s['loci']} loci, 48 of them new, at an effective size of {n(s['effective_size'])}",
        f"{m['effector_symbols']} symbols, of which {m['effector_symbols_matched']} name a graph gene",
        f"as of {D['trait']['as_of']}", f"{n(D['groups']['linked'])} genes with a GWAS link", f"{D['groups']['tied_without_a_link']} tied by sources",
        f"{D['crosscheck']['served_counts']['locus_not_resolved']} loci that carry the link")
    assert m["matched_through_ensembl_symbol"] == ["GBA1"] and {"FAM105B", "C6orf57", "HLA", "KIAA0125"} <= set(m["effector_symbols_not_in_the_graph"])
    assert len(m["effector_symbols_not_in_the_graph"]) == m["effector_symbols"] - m["effector_symbols_matched"] == 10
    for w in ("Defined on 28 Oct 2026, on Ian's word, before any gene of the study was matched to the graph", "registry key `uffelmann2026_supp`, 897,272 bytes"):
        assert w in PROTO, f"docs/PROTOCOL.md no longer says: {w}"


def test_gaps_that_close():
    a, w = D["gaps_that_close"]["primary:all_loci"], D["gaps_that_close"]["primary:new_loci"]
    for g, label in (("tied_without_a_link", None), ("one_source", "| one source |")):
        x, y = a["groups"][g], w["groups"][g]
        has(f"{x['genes']} | {x['closed']}, {pct(x['share'])} | {y['closed']}, {pct(y['share'])} |")
    z, zw = a["groups"]["nothing"], w["groups"]["nothing"]
    has(f"| nothing | {n(z['genes'])} | {z['closed']}, {pct(z['share'], 2)} | {zw['closed']}, {pct(zw['share'], 2)} |")
    has(f"Tied minus nothing is {sig(a['tied_minus_nothing'])} over all loci", f"the difference is {sig(w['tied_minus_nothing'])}",
        f"({w['detectable_at_80_percent_power']:.3f})", f"-0.030 [-0.130, +0.054]", f"The study names {a['linked_named']} of the {a['linked']} genes")
    assert a["tied_minus_nothing"][1] > 0, "the reading: the primary interval lies above 0"
    assert a["tied_genes_closed"] == ["ABI3", "APP", "PLCG2"] and w["tied_genes_closed"] == []
    has("ABI3, APP and PLCG2", "none of the 74 is named")


def test_the_order_of_a_locus():
    o, d = D["order_of_a_locus"]["primary"], D["order_of_a_locus"]["data_only"]
    has(f"{o['assessed']} of the graph's {o['loci']} loci", f"at {o['with_first_genes']} of them", f"at {o['graph']:.0f} of the {o['with_first_genes']}",
        f"chance gives {o['chance']:.1f}", f"gives {o['baseline']:.1f}", f"Graph minus chance is {sig(o['graph_minus_chance'], 1)}",
        f"that baseline {sig(o['graph_minus_baseline'], 1)}", f"at {o['singled_out_named']} of {o['singled_out_assessed']} loci",
        f"gives {d['graph']:.0f} against a chance of {d['chance']:.1f} at {d['with_first_genes']} loci ({sig(d['graph_minus_chance'], 1)})")
    assert o["graph_minus_chance"][1] > 0 and o["graph_minus_baseline"][1] > 0
    miss = {(tuple(r["first"]), tuple(r["named"])) for r in o["rows"] if r["graph"] == 0}
    assert miss == {(("MPO",), ("TSPOAP1",)), (("NDUFAF6",), ("TP53INP1",)), (("ACTR2",), ("SPRED2",))}


def test_the_largest_study():
    f = D["largest_study"]
    L, c = f["largest"], f["check_with_individuals_as_cases"]
    has(f"{L['study']} ({L['year']}), held as {n(L['cases'])} cases and {n(L['controls'])}", f"odds ratio of {L['smallest_odds_ratio_it_could_find']:.3f}",
        f"none of the {f['new_loci']} should have been seen", f"from {f['odds_ratio_above_1_range'][0]:.3f} to {f['odds_ratio_above_1_range'][1]:.3f}",
        f"effective size of {n(f['needed_at_1_05']['effective'])}", f"{n(c['individuals'])} proxy cases", f"the floor at 20% is {c['at_20_percent']:.3f}",
        f"{c['should_have_been_seen']} of the 48 new loci could have been found")
    assert f["should_have_been_seen"] == 0 and f["check_at_20_percent"]
    for s in c["seen"]:
        has(f"{s['genes']} (odds ratio {s['odds_ratio_above_1']:.3f} against {s['could_have_found']:.3f}" if s["genes"] == "DOK3" else
            f"{s['genes']} ({s['odds_ratio_above_1']:.3f} against {s['could_have_found']:.3f})")
    assert "28 Oct 2026, after the result of \"A study published after the graph: Alzheimer disease\"" in PROTO


def test_the_ranking():
    R, a = D["ranking"]["primary:new_loci"], D["ranking"]["primary:all_loci"]
    has(f"ranks {n(D['ranking']['candidates'])} candidate genes", f"names {R['positives']} of them",
        f"{R['roc_area']['served']:.3f} under the served score, {R['roc_area']['neighbours']:.3f} under the neighbour score and {R['roc_area']['degree']:.3f}",
        f"Served minus neighbours is {sig(R['served_minus_neighbours'])}", f"served minus degree {sig(R['served_minus_degree'])}",
        f"a gain of more than {R['served_minus_neighbours'][2]:.3f}", f"({a['positives']} positives) served minus neighbours is {sig(a['served_minus_neighbours'])}")
    assert R["served_minus_neighbours"][1] < 0 < R["served_minus_neighbours"][2], "the reading: no difference detected"
