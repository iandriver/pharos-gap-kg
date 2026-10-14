"""The figures in the progress report are the figures in the committed records.

The report is read by people who were not there when a number was measured. Each figure it gives is
looked up here in the record that holds it, so a record that moves after the report was written, or a
figure retyped wrongly, fails the build and names the line.

Skipped where the report is absent."""
import json
import re

import pytest

from pharos_gap_kg import paths

REPORT = paths.ROOT / "docs" / "PROGRESS_20261027.md" if paths.ROOT else None
pytestmark = pytest.mark.skipif(REPORT is None or not REPORT.is_file(), reason="no progress report in this checkout")

J = lambda f: json.loads((paths.DOCS_DATA / f).read_text())
pc = lambda x: f"{100 * x:.1f}%"
iv = lambda d: f"{d['difference']:+.4f} [{d['low']:+.4f}, {d['high']:+.4f}]"
iv3 = lambda d: f"{d['difference']:+.3f} [{d['low']:+.3f}, {d['high']:+.3f}]"
MODELS = ("inception/mercury-2.5", "deepseek/deepseek-v4-flash", "z-ai/glm-5.3-flash")


@pytest.fixture(scope="module")
def T():
    return re.sub(r"\s+", " ", REPORT.read_text())


def has(T, *wanted):
    for w in wanted:
        assert w in T, f"the report no longer says: {w}"


def test_no_em_dash(T):
    assert "—" not in T


def test_the_graphs(T):
    b, c, m = J("bundle_pharos319.json"), J("graph_pharos319_context_disease.json"), J("graph_pharos319.json")
    d = m["drugs"]
    has(T, f"{m['edges']:,} edges and {m['gaps']:,} gaps over {m['nodes']:,} nodes",
        f"{b['links']:,} GWAS links between {b['targets_with_a_link']:,} of Pharos's {b['targets']:,} protein targets and {b['traits']:,} traits, from {b['studies']:,} studies in {b['papers']:,} papers",
        f"{m['assertions']['edges']:,} assertions of a target with a condition from nine disease sources, and the {m['assertions']['gaps']:,} rows",
        f"{m['expression']['edges']:,} expression observations", f"{d['drugs']:,} drugs with {d['drug_target_edges']:,} target and {d['drug_indication_edges']:,} indication edges",
        f"{d['through_the_mechanism_target']:,} of {d['indication_rows']:,} rows", f"{m['target_gaps']['gaps']:,} target-level gaps",
        f"{m['ontology']['same_as_edges']} `biolink:same_as` edges", f"joins {m['traits_that_reach_another_part']['joined']} GWAS traits",
        f"{c['gaps']:,} of {c['rows']:,} rows are withheld", f"{d['indications_not_grounded']} indications the release names without an identifier")
    assert m["conformance"] == "10 of 10 checks passed" and len(m["assertions"]["edges_by_source"]) == 9 and m["edges_by_part"]["links"] == b["links"]
    has(T, "The merged graph passes 10 of 10 checks")
    cats = [k for k in m["node_categories"] if k.startswith("biolink:")]
    for k in cats + list(m["predicates"]) + [J("graph_pharos319_target_gaps.json")["gap_predicate"]]:
        has(T, f"`{k}`")
    assert len(cats) == 7 and len(m["predicates"]) == 6 and "seven BioLink categories" in T and "six BioLink predicates" in T
    assert m["nodes_not_grounded"] == d["indications_not_grounded"] + m["node_categories"]["pharosgapkg:KnowledgeGap"], "the only nodes not grounded are those indications and the kinds of gap"


def test_usability(T):
    rows = {"usability_results_eval.json": "the first seven questions", "usability_results_gaps.json": "target-level gaps", "usability_results_rank.json": "rankings and their rates",
            "usability_results_new.json": "other sources, drugs and the gaps of a disease"}
    typed = []
    for f, name in rows.items():
        u = J(f)
        assert u["questions"] == 35 and u["verdicts"]["usable_by_a_small_model"]
        cell = lambda cond: ", ".join(str(u["cells"][f"{m}|{cond}"]["correct"]) for m in MODELS)
        has(T, f"| {name} | {cell('typed')} | {cell('sql')} | {cell('none')} |")
        typed += [u["cells"][f"{m}|typed"]["correct"] for m in MODELS]
        assert not any(v["fabricated"] for v in u["cells"].values())
    assert set(typed) == {34, 35}
    e, g, r, n = (J(f)["pooled_typed_minus_sql"] for f in rows)
    assert [J(f)["verdicts"]["typed_better_than_sql"] for f in rows] == [True, False, False, True]
    has(T, f"beat SQL on the first set ({iv3(e)}) and on the fourth ({iv3(n)})", iv3(g), iv3(r), "on four separate frozen task sets")
    by_model = J("usability_results_new.json")["typed_minus_sql"]
    assert sorted(v["a_only"] for v in by_model.values()) == [0, 1, 6] and not any(v["b_only"] for v in by_model.values()) and "six of the seven answers that differ are one model's" in T


def test_pruning(T):
    p1, p2, ap = J("planted_pharos319.json"), J("pruning2_pharos319.json"), J("pruning2_applied_pharos319.json")
    assert (p1["bar"]["planted_withheld_at_least"], p1["bar"]["real_withheld_at_most"]) == (0.8, 0.5) == (p2["bar"]["planted_withheld_at_least"], p2["bar"]["real_withheld_at_most"])
    passing = lambda d: sum(bool(v.get("passes")) for k, v in d.items() if k != "all sources")
    assert len(p1["sources"]) - 1 == 11 and len(p1["rule_applied_to"]) == 1
    assert (passing(p2["code_only_L0"]), passing(p2["lift_cut_L0b"]), passing(p2["reader"]["test_L1"])) == (2, 5, 3)
    has(T, "| in-graph corroboration: another source, a related disease or an interaction partner has the pair | 1 |",
        "| no PubMed abstract mentions the gene with the disease | 2 |", "scored on the other | 5 |", "must quote the passage it relies on | 3 |",
        f"withheld {pc(p1['sources']['all sources']['planted_withheld'])} of {p1['planted_pairs']:,} planted pairs")
    a = p2["reader"]["against_lift"]
    assert abs(a["l1_real_withheld"] - a["lift_cut_real_withheld"]) < 1e-12 and not a["reader_adds"]
    has(T, f"withholds {pc(a['l1_planted_withheld'])} of planted pairs and the lift cut {pc(a['lift_cut_planted_withheld'])} ({iv3(a)})")
    ctd = ap["sources"]["CTD"]
    assert len(ap["sources"]) == 5 and ap["withheld"] == sum(v["withheld_in_all"] for v in ap["sources"].values())
    has(T, f"Of the {ap['withheld']:,} rows withheld, {ctd['withheld_in_all']:,} are CTD's", f"is {pc(ctd['share_of_one_source_pairs_withheld_in_all'])} against the bar of 50%",
        f"rests on {p2['lift_cut_L0b']['CTD']['planted_pairs_the_cut_was_set_on']} planted pairs")
    out = [k for k, v in p2["lift_cut_L0b"].items() if not v.get("passes") and not p2["reader"]["test_L1"][k].get("passes") and not p2["code_only_L0"][k].get("passes")]
    assert {"Expression Atlas", "DrugCentral Indication", "JensenLab Knowledge AmyCo"} <= set(out)


def test_the_key_question(T):
    s18, s22 = J("split_gene_trait_2018.json"), J("split_gene_trait_2022.json")
    b, e, c, pr = J("baselines_gene_trait_2018.json"), J("enrich_gene_trait_2018.json"), J("combine_gene_trait_2018_scored.json"), J("prospect_gene_trait_2022.json")
    has(T, f"the {s18['train']['links']:,} links reported by 2018 and scores {s18['test']['links']:,} later links over {s18['test']['traits']} traits",
        f"scores {s22['test']['links']:,} links first reported from 2023 in the GWAS Catalog, over {s22['test']['traits']:,} traits")
    f4 = lambda x: f"{x:.4f}"
    has(T, f"| chance | {f4(b['results']['B0']['primary'][0])} |", f"| target degree | {f4(b['results']['B1']['primary'][0])} |",
        f"| neighbours: parameter-free, Pharos links only | {f4(b['results']['B3']['primary'][0])} |")
    blends = {k: v["primary"][0] for k, v in e["results"].items() if k.startswith("B3+X")}
    best = max(blends, key=blends.get)
    assert len(blends) == 5 and all(e["paired"][f"{k}-B3"]["primary"]["high"] < 0 for k in blends) and not e["context_adds_to_neighbours"]
    has(T, f"best of five kinds | {f4(blends[best])} | {iv(e['paired'][best + '-B3']['primary'])} |")
    alone = sorted(e["results"][f"X{i}"]["primary"][0] for i in range(1, 6))
    has(T, f"scores {alone[0]:.3f} to {alone[-1]:.3f} alone")
    for name, m in (("matrix factorisation", "mf"), ("LightGCN", "lightgcn"), ("LightGCN with lineage and position edges", "lightgcn_lin_pos")):
        r = J(f"model_{m}_2018_scored.json")
        has(T, f"| {name} | {f4(r['results']['model']['full'][0])} | {iv(r['paired']['model-B3']['full'])} |")
    has(T, f"six fitted weights | {f4(c['results']['combine_pos']['full'][0])} | {iv(c['paired']['combine_pos-B3']['full'])} |")
    assert len(pr["scorers"]["combine_pos"]["weights"]) == 6
    served, graph = "combine[B1+B3+B5]:F1", "lightgcn_lin_pos:L05:steps"
    has(T, f"| links from 2023 | neighbours | {f4(pr['results']['B3']['full'][0])} |",
        f"retuned | {f4(pr['results'][graph]['full'][0])} | {iv(pr['paired'][graph + '-B3']['full'])} |",
        f"two fitted weights | {f4(pr['results'][served]['full'][0])} | {iv(pr['paired'][served + '-B3']['full'])} |")
    w = J("prospect_weights_2022.json")["combiners"][served]["weights"]
    assert sum(1 for v in (w.values() if isinstance(w, dict) else w) if v != 1) == 2, "the served combiner has two fitted weights beside the fixed one"
    d = pr["paired"][served + "-B3"]
    assert d["full"]["low"] > 0 and d["low"]["high"] < 0 and d["neutral"]["low"] < 0 < d["neutral"]["high"], "ahead on all links, not ahead on links nothing earlier explains"
    n = pr["counts"]
    has(T, f"{pc(n['positives_lineage_or_position_linked'] / n['test_links'])} of the {n['test_links']:,} test links follow from an earlier link")
    arms = J("rehearsal_gene_trait_2018.json")["arms"]
    ahead = sorted(k for k, v in arms.items() if v["against_neighbours"]["low"] > 0)
    assert ahead == [served, "combine[B1+B3+B5]:F3"] and len(arms) == 28 and len({k.split(":")[0] for k in arms}) == 8
    has(T, "the only one of eight candidates, tried in 28 variants", f"({iv(arms[served]['against_neighbours'])})")
    worse = [k for k in ("mf", "lightgcn", "lightgcn_pos", "lightgcn_lin_pos") if J(f"model_{k}_2018_scored.json")["paired"]["model-B3"]["full"]["low"] > 0]
    assert not worse and pr["paired"][graph + "-B3"]["full"]["low"] < 0 and pr["paired"]["mf-B3"]["full"]["high"] < 0 and pr["paired"]["lightgcn_lin_pos-B3"]["full"]["high"] < 0


def test_explainability(T):
    c = J("calibration_prospect_2022.json")
    top = c["by_trait"]["1-10"]
    assert c["pooled_over"]["traits"] == 2225 and (c["first_year"], c["data_ends"]) == (2023, 2026)
    has(T, f"{pc(top['all']['share'])} were reported in 2023 to 2026, against {100 * c['base_rate']:.2f}% of all candidates", f"the rate is {pc(top['not_linked']['share'])}")
    ask_doc = (paths.ROOT / "docs" / "ASK.md").read_text()
    assert "4 GWAS studies in 4 paper(s), the first in 2007" in ask_doc and "25 more studies came after 2010" in ask_doc and "share 45 target(s)" in ask_doc
    has(T, "the four GWAS studies behind that edge", "the first in 2007", "25 more came later", "the 45 targets that Parkinson disease and Crohn's disease share")


def test_gaps(T):
    t = J("target_gaps_pharos319.json")["dark_gap"]
    r, d = t["scores"]["route cost"], t["route_cost_against_best_single_count"]
    assert t["best_single_count"] == "pubmed score"
    has(T, f"For {t['with_a_gap']:,} Tdark targets, Pharos's least-studied level", f"the {t['promoted']} of {t['in_both_releases']:,} targets that were promoted",
        f"{r['auc']:.3f} [{r['auc_low']:.3f}, {r['auc_high']:.3f}]", f"{d['difference']:+.3f} [{d['low']:+.3f}, {d['high']:+.3f}] over the PubMed score alone",
        f"{pc(t['bands']['1 flip, nearest third']['share'])} were promoted")
    if (paths.bundle("pharos319") / "study_sizes.parquet").is_file():       # the worked disease, as the question layer answers it
        from pharos_gap_kg import ask
        o = ask.Graph(paths.bundle("pharos319")).gaps("Parkinson disease", kind="no_genetic_support", top=2)
        h, k = o["holds"], {g["kind"]: g for g in o["gaps"]}
        f = k["no_genetic_support"]["what_an_absent_link_could_have_seen"]["largest"]
        has(T, f"For Parkinson disease the graph ties {h['tied_genes']} genes: {h['gwas_linked_genes']} by a GWAS link, {h['genes_asserted_by_two_or_more_sources']} by two or more sources and "
               f"{h['genes_that_are_the_mechanism_target_of_one']} as the mechanism target", f"{o['counts']['no_compound']} have no drug and no compound", f"and {o['counts']['dark']} are still Tdark",
            f"{k['locus_not_resolved']['genes']} of the {h['gwas_linked_genes']} GWAS-linked genes sit in {o['counts']['locus_not_resolved']} loci",
            f"{o['counts']['no_genetic_support']} genes are tied with no GWAS link", f"an odds ratio of about {f['smallest_odds_ratio_it_could_find']:.2f} to find a variant")
        top = k["no_genetic_support"]["top"]
        assert [x["symbol"] for x in top] == ["ATP13A2", "PARK7"] and all(len(x["tied_by"]["sources"]) == 6 for x in top) and "ATP13A2 and PARK7, are each asserted by six sources" in T
    c = J("disease_gaps_pharos319.json")
    o = c["smallest_odds_ratio_the_largest_study_could_find"]
    has(T, f"Over the {c['traits']} disease traits", f"{pc(c['share_of_gwas_linked_genes_in_a_locus_with_another'])} of GWAS-linked genes share their locus",
        f"singles out one gene in {pc(c['share_of_such_loci_where_other_evidence_singles_out_one_gene'])} of those loci", f"{pc(c['share_of_tied_genes_with_no_compound'])} of tied genes have no drug",
        f"for the {c['of_them_cases_and_controls']} traits whose largest study has cases and controls, the median is an odds ratio of {o['median']:.2f}")
    assert c["min_sources"] == 2 and "two or more disease sources" in T


def test_findings_about_the_data(T):
    s, c, k = J("snapshots_pharos316_pharos319.json"), J("context_gene_trait_2018.json")["audit"], J("controls_gene_trait_2018.json")
    assert {"gwas", "tiga", "disease", "ncats_ppi", "pathway", "ncats_ligand_activity"} <= set(s["identical"])
    has(T, f"{len(s['identical'])} of the {len(s['identical']) + len(s['differ'])} tables they share are identical")
    has(T, f"it scores {c['sources']['JensenLab Experiment TIGA']['primary'][0]:.2f} where the neighbour score reaches {c['baselines']['B3'][0]:.2f}")
    sh, res = k["shares"], k["results"]
    has(T, f"({sh['either']:,} of {sh['positives']:,})", f"links are {res['B1']['loci']:,} loci")
    assert round(sh["either"] / sh["positives"], 1) == 0.2
    assert round(res["B1"]["primary"][0] / res["B0"]["primary"][0]) == 4 and round(res["B1"]["locus_collapsed"][0] / res["B0"]["locus_collapsed"][0]) == 2
    controls = re.sub(r"\s+", " ", (paths.ROOT / "findings" / "CONTROLS.md").read_text())
    assert "505 retired EFO disease ids" in controls and "505 of the 1,404 trait ids" in controls and k["coverage"]["traits"] == 1404
    has(T, "Of the 1,404 traits with a link by 2018, 505 carry an ontology id that has since been retired and replaced")
