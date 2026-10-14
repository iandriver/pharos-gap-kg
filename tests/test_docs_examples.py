"""The answers quoted in docs/ASK.md are the answers the code gives on the real bundle.

Skipped where the bundle has not been built. A quoted example that has drifted from what a reader
would get by running the command is a claim nobody checked."""
import pytest

from pharos_gap_kg import ask, paths

BUNDLE = paths.bundle("pharos319")
pytestmark = pytest.mark.skipif(not (BUNDLE / "manifest.json").is_file(), reason="the pharos319 bundle is not built")


@pytest.fixture(scope="module")
def G():
    return ask.Graph(BUNDLE)


def test_status(G):
    assert G.status()["answer"].startswith("A graph of 295,327 GWAS links between 18,005 of Pharos's 20,412 protein targets and 2,247 traits")


def test_why_with_a_link(G):
    w = G.why("APOE", "Alzheimer disease", as_of=2010)
    assert w["edge"] == "pharosgapkg:e2e3b1a3f2814" and w["studies_after_as_of"] == 25
    assert [s["study"] for s in w["studies"]] == ["GCST000018", "GCST000479", "GCST000682", "GCST000808"]
    assert w["source"]["sha256"] == "86ebfffdd5fd03bcc23ef43a31251d40572ccd7feec62dadd3215f276f8b0953"


def test_why_without_a_link(G):
    assert G.why("TREM2", "Alzheimer disease", as_of=2012)["answer"].startswith(
        "As of 2012 there was no link between TREM2 and Alzheimer disease. One was first reported in 2013.")
    h = G.why("HLA-B", "Parkinson disease")
    assert "ranks 16 of 20,152 candidates" in h["answer"] and "MHC region" in h["answer"] and "not evidence" in h["answer"]
    assert "by this graph's ranking (shared-trait neighbours, traits already linked, and linked genes nearby) HLA-B ranks" in h["answer"]


def test_connect(G):
    assert G.connect("Parkinson disease", "Crohn's disease")["answer"] == \
        "As of 2023, Parkinson disease and Crohn's disease share 45 target(s) (of 260 and 650)."


def test_predict(G):
    q = G.predict("type 2 diabetes mellitus", top=5)
    assert q["answer"] == ("The 5 highest-ranked candidate targets for type 2 diabetes mellitus as of 2023, out of 17,891. Genes within 500,000 bases of "
                           "each other are shown as one lead. In the MHC region of chromosome 6, where neighbouring candidates are best read as one lead: 2 of them. "
                           "A target with no link is not in this list; 431 of them lie beside a gene with this link, and unlinked=true ranks those.")
    got = [(c["rank"], c["symbol"], c["score_parts"]["shared_trait_neighbours"], c["score_parts"]["traits_already_linked"], c["score_parts"]["genes_nearby_with_this_link"])
           for c in q["candidates"]]
    assert got == [(1, "GPANK1", 19.214, 104, 90), (16, "HIST1H1T", 19.598, 103, 11), (17, "MYBPC3", 19.181, 97, 12), (20, "MARK4", 15.893, 73, 20), (23, "CTDNEP1", 14.547, 74, 23)]
    assert q["candidates"][0]["same_locus"][:3] == ["C6orf47", "CSNK2B", "LY6G5B"] and q["candidates"][1]["same_locus"] == ["HFE"] and q["candidates"][4]["same_locus"] == ["ELP5", "GABARAP", "PHF23", "CLDN7"]
    assert all(c["explained_by"] for c in q["candidates"]) and [bool(c["caution"]) for c in q["candidates"]] == [True, True, False, False, False]
    assert q["candidates"][0]["how_often_this_held"] == ("In the test on links first reported from 2023 to 2026, 11% of candidates at ranks 1-10 that an earlier link "
                                                         "explains by trait wording or by position were reported by 2026, counted over every trait ranked (base rate "
                                                         "0.5%). A candidate not reported is not known to be absent.")
    scored_only = G.backtest["calibration"]["pooled_over"]["traits"] == 2225
    assert scored_only, "the bundle's calibration must be counted over every trait ranked"
    m = q["method"]
    assert "weighted 1, 1/32 and 1/4" in m and "fitted on links to 2022" in m and "illustration" not in m
    assert "On 194,017 links first reported from 2023 to 2026" in m and "scored 0.101 (mean average precision over 1,110 traits) against 0.095" in m
    assert "(+0.0066, 95% interval +0.0031 to +0.0102)" in m and "level with neighbours alone" in m
    t = G.predict("TREM2", top=4)
    assert "it scored 0.063 against 0.062 on the same links" in t["method"] and "not declared before it" in t["method"]
    p = G.predict("type 2 diabetes mellitus", as_of=2018, top=8, unexplained_only=True)
    assert "out of 19,260" in p["answer"] and "Reported between 2019 and 2023: 3 of these 8." in p["answer"] and "588 of them lie beside a gene with this link" in p["answer"]
    assert [(c["rank"], c["symbol"]) for c in p["candidates"]] == [(41, "ZPR1"), (49, "AIM2"), (54, "HIST1H4A"), (60, "TMEM258"),
                                                                 (65, "HLA-A"), (101, "SLC22A5"), (116, "MYB"), (132, "CTCFL")]
    later = {c["symbol"]: c["any_gene_at_this_locus_reported_after_as_of"] for c in p["candidates"]}
    assert {k for k, v in later.items() if v} == {"HIST1H4A", "TMEM258", "HLA-A"} and later["TMEM258"] == 2020
    assert "illustration and not a test" in p["method"] and all(c["score_parts"]["genes_nearby_with_this_link"] == 0 for c in p["candidates"])


def test_targets_with_no_link(G):
    assert int((G.t.shape[0])) == 20412
    u = G.predict("Parkinson disease", top=5, unlinked=True)
    assert u["answer"] == ("The 5 highest-ranked targets with no link as of 2023 that lie beside a gene linked to Parkinson disease, out of 38. "
                           "Targets within 500,000 bases of each other are shown as one lead.")
    got = [(c["rank"], c["symbol"], c["pharos_level"], c["genes_nearby_with_this_link"], [(n["symbol"], n["bases_apart"]) for n in c["nearest"]], c["same_locus"], bool(c["caution"]))
           for c in u["candidates"]]
    assert got == [(1, "RUSC1-AS1", "Tdark", 27, [("FDPS", 0), ("RUSC1", 0), ("ASH1L", 10883)], [], False),
                   (2, "PPP1R2P1", "Tdark", 20, [("HLA-DQB2", 115684), ("HLA-DQA2", 132018), ("HLA-DQB1", 210833)], [], True),
                   (3, "FBXL19-AS1", "Tdark", 16, [("FBXL19", 0), ("ORAI3", 25789), ("SETD1A", 34015)], ["PYCARD-AS1"], False),
                   (5, "FAM223B", "Tdark", 13, [("SLC10A3", 141780), ("UBL4A", 145855), ("LAGE3", 153189)], ["FAM223A"], False),
                   (7, "TMEM271", "Tdark", 9, [("PCGF3", 129453), ("CPLX1", 208662), ("GAK", 272981)], [], False)]
    assert u["candidates"][0]["how_often_this_held"].startswith("In the test on links first reported from 2023 to 2026, 5.0% of targets with no link listed at ranks 1-5 for a trait")
    assert "Nothing is fitted" in u["method"] and "0.143 (mean average precision over 714 traits) against 0.118" in u["method"] and "1,199 of those 6,729 links" in u["method"]
    d = G.predict("TMEM271", top=4)
    assert d["answer"] == ("TMEM271 has no link as of 2023. The 4 highest-ranked traits that a gene within 500,000 bases of it is linked to, out of 65. "
                           "Traits the ontology files under one another are shown as one lead.")
    assert [(c["rank"], c["name"], c["genes_nearby_with_this_link"], [n["symbol"] for n in c["nearest"]]) for c in d["candidates"]] == [
        (1, "appendicular lean mass", 13, ["PDE6B", "ATP5ME", "MYL5"]), (2, "heel bone mineral density", 13, ["PDE6B", "ATP5ME", "MYL5"]),
        (3, "BMI-adjusted waist-hip ratio", 13, ["PDE6B", "ATP5ME", "MYL5"]), (4, "high density lipoprotein cholesterol measurement", 9, ["PCGF3", "CPLX1", "GAK"])]
    assert d["candidates"][0]["how_often_this_held"] == ("In the test on links first reported from 2023 to 2026, 3.3% of traits listed at ranks 1-5 for a target with no link "
                                                         "were linked to it by 2026. Of all listed pairs 2.8% were, and of all pairs of a target with no link and a trait "
                                                         "0.13%. A pair not reported is not known to be absent.")
    w = G.why("TMEM271", "Parkinson disease")
    assert w["answer"] == ("There is no link between TMEM271 and Parkinson disease in this graph (data to 2023). No link here means no study in this source had reported one by "
                           "that year. It is not evidence that there is no association. What nearby evidence exists: PCGF3, 129,453 bases away, was linked to this trait; TMEM271 "
                           "has no link of its own, and by position it ranks 7 of the 38 targets with no link that lie beside a gene with this one. The largest study behind this "
                           "trait's links (GCST010991, 2021): it had 20,184 cases and 397,324 controls, so a variant carried by 20% of chromosomes needed an odds ratio of about "
                           "1.08 to be found. An absent link is no effect of about that size found at a common variant; smaller effects and rarer variants were not visible. This "
                           "is a rough floor, and a study that reported nothing for the trait is not in the source. The study an absent link calls for: for a common variant, "
                           "about 55,000 cases for an odds ratio of 1.05, 14,000 for 1.1 and 3,900 for 1.2, with the largest study's 19.7 controls per case; the largest study had"
                           " enough for 1.1 and 1.2. For rare variants (1 in 1,000 chromosomes combined, odds ratio 3), a gene-based test on sequence data would need about 13,000"
                           " cases. Apart from GWAS, nothing in this graph ties TMEM271 to this condition: no disease source asserts the pair, no drug path runs between them, and"
                           " no expression observation is on record. That is not evidence that there is no association.")
    assert "same association signal" in u["caveat"] and "same association signal" in d["caveat"]


def test_status_names_the_ranking(G):
    r = G.status()["ranking"]
    assert r["scorer"]["name"] == "combine[B1+B3+B5]:F1" and r["scorer"]["weights"] == {"B1": 1 / 32, "B3": 1.0, "B5": 0.25}
    assert r["tested_on_later_links"]["test_links"] == 194017 and r["tested_on_later_links"]["first_year"] == 2023


def test_gaps(G):
    if G.tg is None:
        pytest.skip("the bundle was built without target-level gaps")
    d = G.gaps("LSM14B")
    assert d["answer"] == ("LSM14B is Tdark. To leave it needs 2 criteria to clear; the nearest route is: PubMed score up by 0.085 (now 4.915, clears at 5); "
                           "and 1 more antibody (now 50, clears at 51). Of 1,021 targets in this band in the Pharos 3.19 database, 11% (10% to 13%) were no "
                           "longer Tdark in Pharos 4.0.0 (May 2026). In the mirror band, 0% of 1,589 targets fell the other way, which is what the pipeline "
                           "alone moves at this distance.")
    assert G.gaps("TREM2")["answer"] == ("TREM2 is Tbio: it has biological evidence and no qualifying compound. It lacks a compound at 1 uM or better, the "
                                         "potency cutoff for its family. Of 17,517 targets with no qualifying compound in the Pharos 3.19 database, 1,359 (8%) "
                                         "held one in Pharos 4.0.0 (May 2026).")
    k = G.gaps(family="Kinase", top=3)
    assert k["answer"] == "The 3 Tdark targets in the family Kinase nearest to leaving Tdark, of 5,499 with a gap. Nearest first."
    assert [(t["symbol"], t["nearest_route"].split(" (")[0]) for t in k["targets"]] == [("SBK2", "6 more antibodies"), ("NRK", "12 more antibodies"), ("RPS6KC1", "1 more GeneRIF")]
    assert "not the same as learning its function" in d["reading"]
    if G.tg_rates.get("graph"):
        assert d["gaps"][0]["graph_edge"] == "pharosgapkg:g46b1570d5751" and G.tg_rates["graph"]["gaps"] == 23236


def test_the_other_sources_as_the_guide_quotes_them(G):
    """Every answer the guide quotes for the merged graph's other parts is the answer the code gives."""
    import re
    if not G.has_others:
        pytest.skip("the bundle does not hold the merged graph's other parts")
    doc = re.sub(r"\s+", " ", (paths.ROOT / "docs" / "ASK.md").read_text().replace("\n> ", " "))
    said = [G.status()["answer"], G.why("TREM2", "Alzheimer disease", as_of=2012)["answer"], G.why("HLA-B", "Parkinson disease")["answer"], G.why("TMEM271", "Parkinson disease")["answer"],
            G.why("ACHE", "Alzheimer disease")["answer"], G.why("MT-RNR1", "Tinnitus")["answer"], G.why("SH3D21", "malignant mesothelioma")["answer"],
            G.known("donepezil")["answer"], G.connect("donepezil", "ACHE")["answer"], G.known("Anesthesia for cesarean section")["answer"]]
    if G.sizes is not None:
        pd_ = G.gaps("Parkinson disease", top=2)
        said += [pd_["answer"], pd_["reading"], G.gaps(disease="malignant mesothelioma")["answer"]]
        import json
        k = {g["kind"]: g for g in pd_["gaps"]}
        for shown in (k["no_compound"]["top"][1], k["locus_not_resolved"]["top"][0]):
            assert re.sub(r"\s+", " ", json.dumps(shown, indent=1)) in doc, "the guide no longer shows this row as the code gives it"
    for a in said:
        assert re.sub(r"\s+", " ", a) in doc, f"the guide no longer quotes: {a[:90]}"
    w = G.why("ACHE", "Alzheimer disease")["other_sources"]
    assert [p["drug"] for p in w["through_a_drug"]] == ["donepezil", "galantamine", "rivastigmine", "tacrine"] and all(p["target_is_the_drugs_mechanism"] and len(p["edges"]) == 2 for p in w["through_a_drug"])
    assert [t["tie"] for t in w["read_under"]] == ["as published", "the sources give it the same name"]
    h = G.why("MT-RNR1", "Tinnitus")["other_sources"]["withheld"]
    assert len(h) == 1 and h[0]["source"] == "CTD" and h[0]["rule"] == "literature_lift" and h[0]["papers_mentioning_both"] == 0 and h[0]["gap"].startswith("pharosgapkg:g")
    x = G.why("SH3D21", "malignant mesothelioma")
    assert x["found"] is False and x["other_sources"]["expression_observations"][0]["direction"] == "higher" and not x["other_sources"]["asserted_by"]
    assert G.predict("donepezil")["problem"] == "no dated links" and G.whats_new("donepezil", 2015)["problem"] == "no dated links"
    assert G.known("ACHE")["node"]["symbol"] == "ACHE" and "Ache" in G.known("ACHE")["note"] and G.known("Ache")["node"]["kind"] == "condition"
