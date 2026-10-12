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
                           "each other are shown as one lead. In the MHC region of chromosome 6, where neighbouring candidates are best read as one lead: 2 of them.")
    got = [(c["rank"], c["symbol"], c["score_parts"]["shared_trait_neighbours"], c["score_parts"]["traits_already_linked"], c["score_parts"]["genes_nearby_with_this_link"])
           for c in q["candidates"]]
    assert got == [(1, "GPANK1", 19.214, 104, 90), (16, "HIST1H1T", 19.598, 103, 11), (17, "MYBPC3", 19.181, 97, 12), (20, "MARK4", 15.893, 73, 20), (23, "CTDNEP1", 14.547, 74, 23)]
    assert q["candidates"][0]["same_locus"][:3] == ["C6orf47", "CSNK2B", "LY6G5B"] and q["candidates"][1]["same_locus"] == ["HFE"] and q["candidates"][4]["same_locus"] == ["ELP5", "GABARAP", "PHF23", "CLDN7"]
    assert all(c["explained_by"] for c in q["candidates"]) and [bool(c["caution"]) for c in q["candidates"]] == [True, True, False, False, False]
    assert q["candidates"][0]["how_often_this_held"] == ("In the test on links first reported from 2023 to 2026, 21% of candidates at ranks 1-10 that an earlier link "
                                                         "explains by trait wording or by position were reported by 2026 (base rate 0.5%). A candidate not reported "
                                                         "is not known to be absent.")
    m = q["method"]
    assert "weighted 1, 1/32 and 1/4" in m and "fitted on links to 2022" in m and "illustration" not in m
    assert "On 194,017 links first reported from 2023 to 2026" in m and "scored 0.101 (mean average precision over 1,110 traits) against 0.095" in m
    assert "(+0.0066, 95% interval +0.0031 to +0.0102)" in m and "level with neighbours alone" in m
    t = G.predict("TREM2", top=4)
    assert "it scored 0.063 against 0.062 on the same links" in t["method"] and "not declared before it" in t["method"]
    p = G.predict("type 2 diabetes mellitus", as_of=2018, top=8, unexplained_only=True)
    assert "out of 19,260" in p["answer"] and "Reported between 2019 and 2023: 3 of these 8." in p["answer"]
    assert [(c["rank"], c["symbol"]) for c in p["candidates"]] == [(41, "ZPR1"), (49, "AIM2"), (54, "HIST1H4A"), (60, "TMEM258"),
                                                                 (65, "HLA-A"), (101, "SLC22A5"), (116, "MYB"), (132, "CTCFL")]
    later = {c["symbol"]: c["any_gene_at_this_locus_reported_after_as_of"] for c in p["candidates"]}
    assert {k for k, v in later.items() if v} == {"HIST1H4A", "TMEM258", "HLA-A"} and later["TMEM258"] == 2020
    assert "illustration and not a test" in p["method"] and all(c["score_parts"]["genes_nearby_with_this_link"] == 0 for c in p["candidates"])


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
