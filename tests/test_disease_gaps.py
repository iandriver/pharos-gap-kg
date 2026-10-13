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
    later = [x for x in G.gaps("Alzheimer disease", as_of=2010, top=25)["gaps"][3]["top"] if "gwas_link_first_reported_later" in x]
    assert all(x["gwas_link_first_reported_later"] > 2010 for x in later)
