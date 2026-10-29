"""A study added after the release: its record, and the answers that read it against the test of 28 Oct.

The table is built from a licensed file kept outside the repository, so the tests that read it skip where
it has not been built (python -I -m pharos_gap_kg.later build pharos319)."""
import json
import re

import pytest

from pharos_gap_kg import ask, paths

B = paths.bundle("pharos319")
REC = paths.DOCS_DATA / "later_pharos319.json"
AD = paths.DOCS_DATA / "adstudy_pharos319.json"
if not (REC.is_file() and AD.is_file() and (B / "manifest.json").is_file()):
    pytest.skip("no record of the added study, or no bundle", allow_module_level=True)
R, T = json.loads(REC.read_text()), json.loads(AD.read_text())
HAVE = (ask.later_dir(B) / "manifest.json").is_file()
need = pytest.mark.skipif(not HAVE, reason="the added study's table is not built here")


@pytest.fixture(scope="module")
def both():
    return ask.Graph(B), ask.Graph(B, later=ask.later_dir(B))


def test_the_record_matches_the_test_of_28_oct():
    m = T["matching"]
    assert (R["loci_read"], R["new_loci_read"]) == (T["study"]["loci"], T["study"]["new_loci"]) == (127, 48)
    assert (R["effector_symbols"], R["effector_symbols_naming_a_graph_gene"]) == (m["effector_symbols"], m["effector_symbols_matched"]) == (191, 181)
    assert R["source_sha256"] == T["study"]["sha256"] and R["trait"] == T["trait"]["id"] and "not committed" in R["kept_in"]


def test_the_table_is_beside_the_bundle_and_not_in_the_archive():
    d = ask.later_dir(B)
    assert d.parent.name == "later" and d.name == B.name and B not in d.parents
    arc = json.loads((paths.DOCS_DATA / "bundle_archive_pharos319.json").read_text())
    assert not any("uffelmann" in f or "later" in f for f in arc["members"])


@need
def test_gaps_mark_what_the_test_of_28_oct_found(both):
    _, G = both
    o = G.gaps(disease="Alzheimer disease", kind="no_genetic_support", top=25)
    k = {g["kind"]: g for g in o["gaps"]}
    closed = T["gaps_that_close"]["primary:all_loci"]["tied_genes_closed"]
    assert k["no_genetic_support"]["named_by_a_later_study"] == len(closed) == 3
    assert f"({', '.join(closed)})" in o["answer"]
    lr = {g["kind"]: g for g in G.gaps(disease="Alzheimer disease")["gaps"]}["locus_not_resolved"]
    assert lr["loci_with_a_gene_named_by_a_later_study"] == T["order_of_a_locus"]["primary"]["assessed"] == 28 and lr["loci"] == T["order_of_a_locus"]["primary"]["loci"]


@need
def test_no_list_or_count_moves(both):
    G0, G = both
    for d in ("Alzheimer disease", "Parkinson disease", "malignant mesothelioma"):
        a, b = G0.gaps(disease=d), G.gaps(disease=d)
        assert a["counts"] == b["counts"] and b["answer"].startswith(a["answer"].split(" This was also read")[0])
    assert G.predict("Alzheimer disease")["candidates"] == G0.predict("Alzheimer disease")["candidates"]


@need
def test_answers_that_do_not_involve_the_trait_are_unchanged(both):
    G0, G = both
    doc = re.sub(r"\s+", " ", (paths.ROOT / "docs" / "ASK.md").read_text().replace("\n> ", " "))
    for q in (lambda g: g.why("HLA-B", "Parkinson disease"), lambda g: g.why("TMEM271", "Parkinson disease"), lambda g: g.known("donepezil"),
              lambda g: g.why("MT-RNR1", "Tinnitus"), lambda g: g.whats_new("Parkinson disease", 2019), lambda g: g.known("SNCA")):
        assert q(G) == q(G0)
    assert re.sub(r"\s+", " ", G.why("HLA-B", "Parkinson disease")["answer"]) in doc


@need
def test_why_says_what_the_study_added(both):
    _, G = both
    a = G.why("ABI3", "Alzheimer disease")
    assert not a["found"] and a["after_the_release"]["named"] and "reports a link the release does not hold" in a["answer"]
    assert a["answer"].endswith("The absence above stands as of the release.") or "The absence above stands as of the release." in a["answer"]
    doc = re.sub(r"\s+", " ", (paths.ROOT / "docs" / "ASK.md").read_text().replace("\n> ", " "))
    assert a["after_the_release"]["answer"] in doc and a["answer"].endswith(a["after_the_release"]["answer"]), "the guide no longer quotes the added study's sentence for ABI3"
    t = G.why("TREM2", "Alzheimer disease")
    assert t["found"] and "also names TREM2" in t["answer"]
    n = G.why("ACHE", "Alzheimer disease")
    assert not n["after_the_release"]["named"] and n["answer"] == ask.Graph(B).why("ACHE", "Alzheimer disease")["answer"]
    k = G.known("Alzheimer disease")["after_the_release"]
    assert k["links"] == R["rows"] and f"names {R['pharos_symbols']} genes at {R['loci_with_a_row']} loci" in k["answer"]
    linked = T["gaps_that_close"]["primary:all_loci"]["linked_named"]
    assert f"{R['pharos_symbols'] - linked} of the {R['pharos_symbols']} have no link to it in the release" in k["answer"]


def test_the_findings_quote_the_record():
    F = re.sub(r"\s+", " ", (paths.ROOT / "findings" / "AD_STUDY.md").read_text().split("## Added 29 Oct 2026: the study as a dated source")[1])
    by, linked = R["rows_also_named_by"], T["gaps_that_close"]["primary:all_loci"]["linked_named"]
    for w in (f"holds {R['rows']} links", f"({R['effector_symbols_naming_a_graph_gene']} of the {R['effector_symbols']} symbols)", f"at {R['loci_with_a_row']} of the {R['loci']} loci",
              f"{R['pharos_symbols_at_new_loci']} of the genes sit at loci the study calls new", f"{linked} of the {R['pharos_symbols']} already carry", f"and {R['pharos_symbols'] - linked} do not",
              f"FLAMES column also names the gene in {by['FLAMES']} rows, Coloc_MR in {by['Coloc_MR']}, its literature column in {by['Literature']} and the exome analysis in {by['Exome']}",
              f"marks 3 of the {T['crosscheck']['served_counts']['no_genetic_support']} genes", f"and 28 of the {T['crosscheck']['served_counts']['locus_not_resolved']} loci"):
        assert w in F, f"findings/AD_STUDY.md (29 Oct) no longer says: {w}"
