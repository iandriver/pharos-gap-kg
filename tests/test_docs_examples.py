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
    assert "ranks 1 of 20,152 candidates" in h["answer"] and "MHC region" in h["answer"] and "not evidence" in h["answer"]


def test_connect(G):
    assert G.connect("Parkinson disease", "Crohn's disease")["answer"] == \
        "As of 2023, Parkinson disease and Crohn's disease share 45 target(s) (of 260 and 650)."


def test_predict(G):
    p = G.predict("type 2 diabetes mellitus", as_of=2018, top=8, unexplained_only=True)
    assert "out of 19,260" in p["answer"] and "Reported between 2019 and 2023: 3 of these 8." in p["answer"]
    assert [(c["rank"], c["symbol"]) for c in p["candidates"]] == [(2, "ZPR1"), (5, "AIM2"), (6, "HIST1H4A"), (11, "TMEM258"),
                                                                 (13, "HLA-A"), (50, "SLC22A5"), (60, "MYB"), (69, "CTCFL")]
    later = {c["symbol"]: c["any_gene_at_this_locus_reported_after_as_of"] for c in p["candidates"]}
    assert {k for k, v in later.items() if v} == {"HIST1H4A", "TMEM258", "HLA-A"} and later["TMEM258"] == 2020
    assert p["candidates"][0]["how_often_this_held"].startswith("In the 2018 backtest, 8% of candidates at ranks 1-10 for a disease")
