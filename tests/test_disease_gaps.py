"""The gaps of a disease, counted twice: by the question, and from the bundle's tables without it.

The question layer builds its lists from indexed lookups. Here every count is recomputed with plain
table operations and a different way of finding loci (connected components of a sparse matrix), for
several diseases, an earlier year and a condition that is not a GWAS trait. Where the two disagree,
one of them is wrong. Skipped where the bundle does not hold the merged graph's other parts."""
import json
import re

import numpy as np
import pandas as pd
import pytest

from pharos_gap_kg import ask, paths, studysize

B = paths.bundle("pharos319")
pytestmark = pytest.mark.skipif(not ((B / "assertions.parquet").is_file() and (B / "study_sizes.parquet").is_file() and (B / "target_gaps.parquet").is_file()),
                                reason="the bundle does not hold the other parts, the study sizes and the target gaps")
norm = lambda s: " ".join(re.sub(r"[^a-z0-9]+", " ", str(s).casefold()).split())


@pytest.fixture(scope="module")
def G():
    return ask.Graph(B)


@pytest.fixture(scope="module")
def T():
    rd = lambda f: pd.read_parquet(B / f)
    return {k: rd(f"{k}.parquet") for k in ("nodes", "links", "neighbours", "evidence", "context_nodes", "assertions", "expression", "drug_targets", "drug_indications", "target_gaps", "study_sizes")}


def by_tables(T, name: str, y: int) -> dict:
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components
    n = T["nodes"]
    t, d = n[n.kind == "target"].reset_index(drop=True), n[n.kind == "trait"]
    sym = dict(zip(t.id, t.name))
    genes = lambda ts: len({sym[x] for x in ts if x in sym})
    cond = T["context_nodes"][T["context_nodes"].kind == "condition"]
    same_name = set(cond[cond.name.map(norm) == norm(name)].id)
    row = d[d.name == name]
    if len(row):
        assert len(row) == 1
        tid = row.id.iat[0]
        ids = {tid} | ({row.live_term.iat[0]} if row.live_term.iat[0] else set()) | same_name
        L = T["links"]
        gw = set(L[(L.trait == tid) & (L.first_year <= y)].target)
    else:
        tid, ids, gw = None, same_name, set()
    a = T["assertions"][T["assertions"].condition.isin(ids)]
    per = a[~a.withheld].groupby("target").source.nunique()
    two, one = set(per[per >= 2].index), set(per[per < 2].index)
    drugs = set(T["drug_indications"][T["drug_indications"].condition.isin(ids)].drug)
    DT = T["drug_targets"]
    moa = set(DT[DT.drug.isin(drugs) & DT.mechanism].target)
    tied = (gw | two | moa) & set(sym)
    tg = T["target_gaps"]
    lig, dark = set(tg[tg.lacks_ligand == True].id), set(tg[tg.flips.notna()].id)   # noqa: E712
    x = T["expression"][T["expression"].condition.isin(ids)]
    out = {"tied": genes(tied), "gwas_entries": len(gw), "gwas": genes(gw), "two": genes(two), "one": genes(one), "moa": genes(moa), "drugs": len(drugs),
           "no_compound": genes(tied & lig - set(DT.target)), "dark": genes(tied & dark), "withheld": int(a.withheld.sum()), "observed_only": genes(set(x.target) - tied - one)}
    if tid is not None:
        out["no_genetic_support"] = genes((two | moa) & set(sym) - gw)
        idx = {v: i for i, v in enumerate(sorted(gw))}
        nb = T["neighbours"]
        nb = nb[nb.a.isin(idx) & nb.b.isin(idx)]
        M = sparse.coo_matrix((np.ones(len(nb)), (nb.a.map(idx), nb.b.map(idx))), shape=(len(idx), len(idx)))
        _, comp = connected_components(M, directed=False) if len(idx) else (0, np.array([], int))
        groups = pd.Series(sorted(gw)).map(sym).groupby(comp).agg(lambda s: len(set(s))) if len(idx) else pd.Series(dtype=int)
        out["loci"], out["locus_genes"] = int((groups >= 2).sum()), int(groups[groups >= 2].sum())
        e = T["evidence"]
        e = e[(e.trait == tid) & (e.year <= y)].drop_duplicates("study").merge(T["study_sizes"], on="study")
        if len(e):
            e["eff"] = [studysize.effective(c, k, i) for c, k, i in zip(e.cases, e.controls, e.individuals)]
            out["largest_study"] = e.sort_values(["eff", "study"]).study.iat[-1]
    return out


CASES = [("Parkinson disease", None), ("Alzheimer disease", None), ("Alzheimer disease", 2015), ("type 2 diabetes mellitus", None), ("schizophrenia", 2018), ("asthma", None),
         ("body height", None), ("malignant mesothelioma", None)]


@pytest.mark.parametrize("name,year", CASES)
def test_the_question_and_the_tables_count_the_same(G, T, name, year):
    y = year or G.latest
    want = by_tables(T, name, y)
    o = G.gaps(disease=name, as_of=year, top=3)
    assert o["found"], o
    h, k = o["holds"], {g["kind"]: g for g in o["gaps"]}
    got = {"tied": h["tied_genes"], "two": h["genes_asserted_by_two_or_more_sources"], "one": h["genes_asserted_by_one_source_only"], "moa": h["genes_that_are_the_mechanism_target_of_one"],
           "drugs": h["drugs_indicated"], "no_compound": k["no_compound"]["targets"], "dark": k["dark"]["targets"], "withheld": k["withheld"]["rows"], "observed_only": k["observed_only"]["targets"]}
    if "no_genetic_support" in want:
        got.update({"gwas_entries": h["gwas_linked_targets"], "gwas": h["gwas_linked_genes"], "no_genetic_support": k["no_genetic_support"]["targets"],
                    "loci": k["locus_not_resolved"]["loci"], "locus_genes": k["locus_not_resolved"]["genes"]})
        floor = k["no_genetic_support"].get("what_an_absent_link_could_have_seen")
        if "largest_study" in want:
            got["largest_study"] = floor["largest"]["study"]
        assert o["not_listed"] == []
    else:
        got.update({"gwas_entries": 0, "gwas": 0})
        assert o["not_listed"] == ["locus_not_resolved", "no_genetic_support"] and h["gwas_linked_targets"] is None and "not a GWAS trait" in o["answer"]
    assert got == want
    for g in o["gaps"]:
        assert g["closes_with"] and g["how_often_this_closed"] and len(g["top"]) <= 3 and g["more"] >= 0


def test_what_the_rows_say(G):
    o = G.gaps("Parkinson disease", top=5)
    k = {g["kind"]: g for g in o["gaps"]}
    r = k["no_compound"]["top"][0]
    assert r["graph_edge"].startswith("pharosgapkg:g") and r["tied_by"]["gwas_link"]["edge"].startswith("pharosgapkg:e") and "potency cutoff" in r["lacks"]
    assert G.gaps(r["symbol"])["gaps"][-1]["graph_edge"] == r["graph_edge"], "the same gap edge the target-level question cites"
    studies = [x["tied_by"].get("gwas_link", {}).get("studies", 0) for x in k["no_compound"]["top"]]
    assert studies == sorted(studies, reverse=True), "strongest first"
    assert all(x["pharos_level"] == "Tdark" and x["nearest_route"] for x in k["dark"]["top"])
    assert all("gwas_link" not in x["tied_by"] and (len(x["tied_by"].get("sources", [])) >= 2 or x["tied_by"].get("mechanism_target_of")) for x in k["no_genetic_support"]["top"])
    lr = k["locus_not_resolved"]
    assert all(x["genes_with_the_link"] >= 2 for x in lr["top"]) and lr["loci_where_other_evidence_singles_out_one_gene"] <= lr["loci"]
    f = k["no_genetic_support"]["what_an_absent_link_could_have_seen"]["largest"]
    assert f == {"study": f["study"], "year": f["year"], **studysize.smallest_effect(f["cases"], f["controls"], 0)} and f["smallest_odds_ratio_it_could_find"] > 1
    assert "rough floor" in o["answer"] and "not a finding that the target matters" in o["reading"]
    rate = json.loads((B / "backtest.json").read_text())["sources_then_gwas"]
    assert f"{rate['base_rate']:.1%} of all pairs" in k["no_genetic_support"]["how_often_this_closed"] and "carry no dates" in k["no_genetic_support"]["how_often_this_closed"]


def test_what_it_refuses_and_what_it_leaves_alone(G):
    assert G.gaps("donepezil")["problem"] == "not a target or a condition"
    assert G.gaps(disease="QXJV5")["problem"] == "not in the graph" and G.gaps(disease="amyloid")["problem"] == "ambiguous"
    assert G.gaps("LSM14B")["answer"].startswith("LSM14B is Tdark. To leave it needs 2 criteria to clear")
    assert G.gaps(family="Kinase", top=3)["answer"].startswith("The 3 Tdark targets in the family Kinase")
    assert "as of 2015" in G.gaps("Alzheimer disease", as_of=2015)["answer"]
    # the answer leads with its summary and counts, and stays short enough to be read whole; one list can be asked for in full
    o = G.gaps("Parkinson disease")
    assert list(o)[:2] == ["found", "answer"] and o["rows_shown_per_list"] == 2 and all(len(g["top"]) <= 2 for g in o["gaps"]) and len(json.dumps(o)) < 10_000
    assert o["counts"] == {g["kind"]: g.get("targets", g.get("loci", g.get("rows"))) for g in o["gaps"]}
    one = G.gaps("Parkinson disease", kind="no_genetic_support", top=4)
    assert [("top" in g) for g in one["gaps"]] == [False, False, False, True, False, False] and len(one["gaps"][3]["top"]) == 4 and one["counts"] == o["counts"]
    # a disease is found whether it is a GWAS trait or a condition of the other sources, whichever of the two kinds is asked
    for name in ("asthma", "mucocutaneous lymph node syndrome", "Alopecia", "Dengue Hemorrhagic Fever"):
        for kind in ("any", "trait", "condition"):
            f = G.find(name, kind)
            assert f["found"] and f["exact_match"] and f["matches"][0]["name"] == name and f["matches"][0]["kind"] == "disease or phenotype", (name, kind)
    assert G.find("malignant mesothelioma", "trait")["matches"][0]["kind"] == "condition" and not G.find("donepezil", "condition")["found"]
    # a count that cannot be made is null, not missing and not zero
    m = G.gaps(disease="malignant mesothelioma")
    assert m["counts"]["locus_not_resolved"] is None and m["counts"]["no_genetic_support"] is None and m["counts"]["no_compound"] > 0 and "null, not zero" in m["answer"]
    assert G.gaps("Parkinson disease", kind="nope")["problem"] == "no such kind" and "does not apply" in G.gaps(disease="malignant mesothelioma", kind="locus_not_resolved")["answer"]
    later = [x for x in G.gaps("Alzheimer disease", as_of=2010, top=25)["gaps"][3]["top"] if "gwas_link_first_reported_later" in x]
    assert all(x["gwas_link_first_reported_later"] > 2010 for x in later)


def test_the_tables_built_without_the_question_layer_agree_with_it(G):
    """`gaptables` computes every disease trait and condition at once and imports none of the question
    layer. Forty disease traits and forty conditions, drawn with a fixed seed, must read the same."""
    import random

    from pharos_gap_kg import diseasegaps, gaptables
    rng = random.Random(20261013)
    traits = rng.sample([int(i) for i in np.flatnonzero(G.d.disease_like.to_numpy())], 40)
    conds = rng.sample(sorted(G.others.groups), 40)
    s = gaptables.build(gaptables.read(B), only={G.d.id[i] for i in traits} | set(conds))["disease_gap_summary"].set_index("condition_id")
    same = lambda a, b: (a is None and (b is None or pd.isna(b))) or (a is not None and b is not None and not pd.isna(b) and (a == b or float(a) == float(b)))
    for r in [{"kind": "trait", "i": i} for i in traits] + [{"kind": "condition", "key": c} for c in conds]:
        o = diseasegaps.for_condition(G, r, G.latest, top=1)
        h, n, row = o["holds"], o["counts"], s.loc[G.d.id[r["i"]] if r["kind"] == "trait" else r["key"]]
        f = (next((g for g in o["gaps"] if g["kind"] == "no_genetic_support"), {}).get("what_an_absent_link_could_have_seen") or {}).get("largest") or {}
        pairs = [(h["tied_genes"], row.genes_tied), (h["gwas_linked_genes"], row.genes_with_a_gwas_link), (h["genes_asserted_by_two_or_more_sources"], row.genes_asserted_by_two_or_more_sources),
                 (h["genes_asserted_by_one_source_only"], row.genes_asserted_by_one_source_only), (h["drugs_indicated"], row.drugs_indicated),
                 (h["genes_that_are_the_mechanism_target_of_one"], row.genes_that_are_a_mechanism_target), (n["no_compound"], row.no_compound), (n["dark"], row.dark),
                 (n.get("locus_not_resolved"), row.loci_not_resolved), (n.get("no_genetic_support"), row.no_genetic_support), (n["withheld"], row.withheld_rows), (n["observed_only"], row.observed_only),
                 (f.get("study"), row.largest_study), (f.get("smallest_odds_ratio_it_could_find"), row.smallest_odds_ratio_it_could_find)]
        assert all(same(a, b) for a, b in pairs), (o["condition"], pairs)


def test_every_disease_trait_summed_two_ways():
    """The census asks the question layer about each of the disease traits in turn and records the
    sums. `gaptables` computes the same traits at once without it. The sums must be the same."""
    from pharos_gap_kg import gaptables
    rec = paths.DOCS_DATA / "disease_gaps_pharos319.json"
    if not rec.is_file():
        pytest.skip("no census on record")
    c = json.loads(rec.read_text())
    T = gaptables.read(B)
    d = T["nodes"][(T["nodes"].kind == "trait") & T["nodes"].disease_like]
    linked = set(T["links"].trait)
    s = gaptables.build(T, only=set(d[d.id.isin(linked)].id))["disease_gap_summary"]
    S = c["summed_over_traits"]
    got = {"traits": len(s), "tied": int(s.genes_tied.sum()), "gwas": int(s.genes_with_a_gwas_link.sum()), "no_compound": int(s.no_compound.sum()), "dark": int(s.dark.sum()),
           "loci": int(s.loci_not_resolved.sum()), "locus_genes": int(s.genes_in_those_loci.sum()), "loci_singled_out": int(s.loci_where_other_evidence_singles_out_one_gene.sum()),
           "no_genetic_support": int(s.no_genetic_support.sum()), "withheld": int(s.withheld_rows.sum()), "observed_only": int(s.observed_only.sum()),
           "two_sources": int(s.genes_asserted_by_two_or_more_sources.sum()), "mechanism": int(s.genes_that_are_a_mechanism_target.sum()),
           "sized": int(s.largest_study.notna().sum()), "cases_and_controls": int(s.smallest_odds_ratio_it_could_find.notna().sum())}
    want = {"traits": c["traits"], "tied": c["tied_genes_summed_over_traits"], "gwas": c["gwas_linked_genes_summed"], **{k: S[k] for k in ("no_compound", "dark", "loci", "locus_genes", "loci_singled_out",
            "no_genetic_support", "withheld", "observed_only", "two_sources", "mechanism")}, "sized": c["traits_whose_largest_study_has_a_size"], "cases_and_controls": c["of_them_cases_and_controls"]}
    assert got == want
    assert round(float(s.smallest_odds_ratio_it_could_find.median()), 3) == c["smallest_odds_ratio_the_largest_study_could_find"]["median"]
