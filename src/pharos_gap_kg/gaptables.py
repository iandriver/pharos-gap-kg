"""What the questions say of the other sources and of the gaps of a disease, as tables, computed without them.

`ask`, `others` and `diseasegaps` answer one target, one condition or one pair at a time from indexed
lookups. This computes the same things for every disease trait and every condition at once from the
bundle's tables, a different way, so the usability test can read its expected answers from here, give
the same tables to the SQL condition, and check that the question layer agrees. Where the two
disagree, one of them is wrong. It imports none of the three.

    same_condition       the ids one trait or condition is read under: its own, the current term of a
                         retired id, and any id the sources give the same name
    disease_gap_summary  for each, the count of every kind of gap, and the largest study behind its
                         GWAS links with the smallest effect that study could have found
    disease_gap_genes    the genes in each list
    disease_loci         the loci that carry a trait's GWAS link on more than one gene

The definitions are in docs/PROTOCOL.md, "Gaps for a disease".

Run:  python -m pharos_gap_kg.gaptables --selftest
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MIN_SOURCES = 2
K, F = (5.4513 + 0.8416) ** 2, 0.32
FILES = ("nodes", "links", "neighbours", "evidence", "context_nodes", "assertions", "expression", "drug_targets", "drug_indications", "target_gaps", "study_sizes")
LISTS = ("tied", "no_compound", "dark", "no_genetic_support", "observed_only")


def norm(s) -> str:
    return " ".join(re.sub(r"[^a-z0-9]+", " ", str(s).casefold()).split())


def read(bundle: Path) -> dict[str, pd.DataFrame]:
    return {k: pd.read_parquet(Path(bundle) / f"{k}.parquet") for k in FILES}


def same_condition(T: dict) -> pd.DataFrame:
    """asked_id, read_id, tie. A trait is asked by its own id. A condition that no trait takes is asked
    by any of the ids the sources give its name, and read under all of them."""
    d = T["nodes"][T["nodes"].kind == "trait"]
    cond = T["context_nodes"][T["context_nodes"].kind == "condition"]
    named = cond.assign(key=cond.name.map(norm)).groupby("key").id.apply(list).to_dict()
    rows, taken = [], set()
    for t, name, live in zip(d.id, d.name, d.live_term):
        ids = [(t, "as published")] + ([(live, "its current ontology term")] if isinstance(live, str) and live and live != t else [])
        ids += [(c, "the sources give it the same name") for c in named.get(norm(name), []) if c not in {i for i, _ in ids}]
        rows += [(t, i, how) for i, how in ids]
        taken.update(i for i, _ in ids)
    for ids in named.values():
        left = sorted(c for c in ids if c not in taken)
        rows += [(a, b, "as the sources give it" if a == b else "the sources give it the same name") for a in left for b in left]
    return pd.DataFrame(rows, columns=["asked_id", "read_id", "tie"])


def floor(cases: int, controls: int, people: int) -> tuple[float | None, float | None]:
    """(smallest odds ratio, smallest effect in standard deviations per allele): one is given, by the study's design."""
    if cases > 0 and controls > 0:
        return round(math.exp(math.sqrt(K * (cases + controls) / (F * cases * controls))), 3), None
    n = cases + controls + people
    return (None, round(math.sqrt(K / (F * n)), 4)) if n > 0 else (None, None)


def size_for_comparison(cases: int, controls: int, people: int) -> float:
    return 4.0 * cases * controls / (cases + controls) if cases > 0 and controls > 0 else float(cases + controls + people)


def build(T: dict, only: set | None = None) -> dict[str, pd.DataFrame]:
    """T: the bundle's tables (`read`). only: asked ids to compute for, all disease traits and conditions if None."""
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components
    n = T["nodes"]
    t, d = n[n.kind == "target"], n[n.kind == "trait"]
    sym = dict(zip(t.id, t.name))
    is_target = set(t.id)
    genes = lambda ts: {sym[x] for x in ts if x in sym}
    S = same_condition(T)
    read_under = S.groupby("asked_id").read_id.apply(list).to_dict()
    A = T["assertions"]
    a_by = {k: g for k, g in A.groupby("condition")}
    x_by = {k: set(g.target) for k, g in T["expression"].groupby("condition")}
    i_by = {k: set(g.drug) for k, g in T["drug_indications"].groupby("condition")}
    DT = T["drug_targets"]
    moa_by = DT[DT.mechanism].groupby("drug").target.apply(set).to_dict()
    acted_on = set(DT.target)
    tg = T["target_gaps"]
    lig, dark = set(tg[tg.lacks_ligand == True].id), set(tg[tg.flips.notna()].id)   # noqa: E712
    L = T["links"]
    gw_by = L.groupby("trait").target.apply(set).to_dict()
    nb = T["neighbours"]
    ev = T["evidence"].drop_duplicates(["trait", "study"]).merge(T["study_sizes"], on="study")
    ev_by = {k: g for k, g in ev.groupby("trait")}
    names = dict(zip(d.id, d.name)) | dict(zip(T["context_nodes"].id, T["context_nodes"].name))
    traits = set(d.id)
    asked = sorted(set(d[d.disease_like].id) | (set(S.asked_id) - traits)) if only is None else sorted(only)
    summary, lists, loci = [], [], []
    for c in asked:
        ids = read_under.get(c, [c])
        a = pd.concat([a_by[i] for i in ids if i in a_by]) if any(i in a_by for i in ids) else A.iloc[0:0]
        per = a[~a.withheld].groupby("target").source.nunique()
        two, one = set(per[per >= MIN_SOURCES].index), set(per[per < MIN_SOURCES].index)
        drugs = set().union(*(i_by.get(i, set()) for i in ids))
        moa = set().union(*(moa_by.get(g, set()) for g in drugs)) if drugs else set()
        is_trait = c in traits
        gw = gw_by.get(c, set()) if is_trait else set()
        tied = (gw | two | moa) & is_target
        seen = set().union(*(x_by.get(i, set()) for i in ids))
        G = {"tied": genes(tied), "no_compound": genes(tied & lig - acted_on), "dark": genes(tied & dark), "observed_only": genes(seen - tied - one)}
        row = {"condition_id": c, "name": names.get(c, c), "is_a_gwas_trait": int(is_trait), "genes_tied": len(G["tied"]), "genes_with_a_gwas_link": len(genes(gw)) if is_trait else None,
               "genes_asserted_by_two_or_more_sources": len(genes(two)), "genes_asserted_by_one_source_only": len(genes(one)), "drugs_indicated": len(drugs),
               "genes_that_are_a_mechanism_target": len(genes(moa)), "no_compound": len(G["no_compound"]), "dark": len(G["dark"]), "loci_not_resolved": None, "genes_in_those_loci": None,
               "loci_where_other_evidence_singles_out_one_gene": None, "no_genetic_support": None, "withheld_rows": int(a.withheld.sum()), "observed_only": len(G["observed_only"]),
               "largest_study": None, "largest_study_year": None, "cases": None, "controls": None, "individuals": None, "smallest_odds_ratio_it_could_find": None,
               "smallest_effect_it_could_find_in_sd": None}
        if is_trait:
            G["no_genetic_support"] = genes((two | moa) & is_target - gw)
            row["no_genetic_support"] = len(G["no_genetic_support"])
            order = sorted(gw)
            idx = {v: i for i, v in enumerate(order)}
            e = nb[nb.a.isin(idx) & nb.b.isin(idx)]
            if len(order):
                M = sparse.coo_matrix((np.ones(len(e)), (e.a.map(idx).to_numpy(), e.b.map(idx).to_numpy())), shape=(len(order), len(order)))
                _, comp = connected_components(M, directed=False)
            else:
                comp = np.array([], int)
            found = []
            for k in np.unique(comp):
                members = [order[i] for i in np.flatnonzero(comp == k)]
                gs = sorted(genes(members))
                if len(gs) >= 2:
                    other = sorted(genes(set(members) & (two | moa)))
                    found.append((gs, other[0] if len(other) == 1 else None))
            for j, (gs, single) in enumerate(sorted(found, key=lambda r: (-len(r[0]), r[0])), 1):
                loci.append((c, j, len(gs), ", ".join(gs), single))
            row.update({"loci_not_resolved": len(found), "genes_in_those_loci": sum(len(g) for g, _ in found), "loci_where_other_evidence_singles_out_one_gene": sum(s is not None for _, s in found)})
            if c in ev_by:
                z = ev_by[c].assign(size=[size_for_comparison(*r) for r in zip(ev_by[c].cases, ev_by[c].controls, ev_by[c].individuals)]).sort_values(["size", "study"]).iloc[-1]
                o, sd = floor(int(z.cases), int(z.controls), int(z.individuals))
                cc = o is not None
                row.update({"largest_study": z.study, "largest_study_year": int(z.year), "cases": int(z.cases) if cc else None, "controls": int(z.controls) if cc else None,
                            "individuals": None if cc else int(z.cases + z.controls + z.individuals), "smallest_odds_ratio_it_could_find": o, "smallest_effect_it_could_find_in_sd": sd})
        summary.append(row)
        lists += [(c, k, g) for k in LISTS if k in G for g in sorted(G[k])]
    return {"same_condition": S, "disease_gap_summary": pd.DataFrame(summary), "disease_gap_genes": pd.DataFrame(lists, columns=["condition_id", "list", "symbol"]),
            "disease_loci": pd.DataFrame(loci, columns=["condition_id", "locus", "genes_with_the_link", "genes", "gene_singled_out_by_other_evidence"])}


# ---- one pair, one drug: what the usability questions ask ------------------------------------------

def pair(T: dict, S: pd.DataFrame, target: str, asked: str) -> dict:
    """What the other sources hold on one target and one condition."""
    ids = set(S[S.asked_id == asked].read_id) or {asked}
    A, X, DT, DI = T["assertions"], T["expression"], T["drug_targets"], T["drug_indications"]
    a = A[(A.target == target) & A.condition.isin(ids)]
    x = X[(X.target == target) & X.condition.isin(ids)]
    acts = DT[(DT.target == target) & DT.drug.isin(set(DI[DI.condition.isin(ids)].drug))]
    held = a[a.withheld]
    return {"sources": int(a[~a.withheld].source.nunique()), "withheld_rows": len(held), "papers": [int(p) for p in held.papers if pd.notna(p)], "rules": sorted(set(held.rule)),
            "drugs": int(acts.drug.nunique()), "mechanism": bool(acts.mechanism.any()), "observations": len(x), "directions": ["higher" if v > 0 else "lower" for v in x.log2fc]}


def drug(T: dict, drug_id: str) -> dict:
    t = T["drug_targets"][T["drug_targets"].drug == drug_id]
    return {"targets": int(t.target.nunique()), "mechanism_targets": int(t[t.mechanism].target.nunique()), "indications": int((T["drug_indications"].drug == drug_id).sum())}


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    N = lambda i, kind, name, live="", dis=True: {"id": i, "kind": kind, "name": name, "live_term": live, "disease_like": dis}
    nodes = pd.DataFrame([N("t:A", "target", "AAA"), N("t:B", "target", "BBB"), N("t:C", "target", "CCC"), N("t:D", "target", "DDD"), N("t:E", "target", "EEE"), N("t:E2", "target", "EEE"),
                          N("EFO:1", "trait", "alpha as published", "MONDO:1"), N("MONDO:2", "trait", "Beta"), N("EFO:9", "trait", "a height", "", False)])
    T = {"nodes": nodes,
         "links": pd.DataFrame([("t:A", "EFO:1", 2010), ("t:B", "EFO:1", 2012), ("t:E", "EFO:1", 2015), ("t:E2", "EFO:1", 2015), ("t:C", "MONDO:2", 2011)], columns=["target", "trait", "first_year"]),
         "neighbours": pd.DataFrame([("t:A", "t:B", 100), ("t:B", "t:A", 100), ("t:E", "t:E2", 1), ("t:E2", "t:E", 1)], columns=["a", "b", "gap"]),
         "evidence": pd.DataFrame([("t:A", "EFO:1", "S1", 1, 2010), ("t:B", "EFO:1", "S2", 2, 2012), ("t:C", "MONDO:2", "S3", 3, 2011)], columns=["target", "trait", "study", "pmid", "year"]),
         "study_sizes": pd.DataFrame([("S1", 1000, 1000, 0), ("S2", 10000, 10000, 0), ("S3", 0, 0, 100000)], columns=["study", "cases", "controls", "individuals"]),
         "context_nodes": pd.DataFrame([("MONDO:1", "condition", "Alpha disease"), ("UMLS:9", "condition", "BETA"), ("UMLS:7", "condition", "Gamma"), ("MESH:7", "condition", "gamma"),
                                        ("DrugCentral:1", "drug", "druga")], columns=["id", "kind", "name"]),
         "assertions": pd.DataFrame([("t:A", "MONDO:1", "S-one", False, "", np.nan), ("t:C", "MONDO:1", "S-one", False, "", np.nan), ("t:C", "MONDO:1", "S-two", False, "", np.nan),
                                     ("t:D", "MONDO:1", "S-one", False, "", np.nan), ("t:D", "MONDO:1", "CTD", True, "literature_lift", 2.0), ("t:B", "UMLS:7", "S-one", False, "", np.nan),
                                     ("t:B", "MESH:7", "S-two", False, "", np.nan)], columns=["target", "condition", "source", "withheld", "rule", "papers"]),
         "expression": pd.DataFrame([("t:D", "MONDO:1", 2.0), ("t:E", "MONDO:1", -1.5), ("t:A", "MONDO:1", 1.2)], columns=["target", "condition", "log2fc"]),
         "drug_targets": pd.DataFrame([("DrugCentral:1", "t:B", True), ("DrugCentral:1", "t:D", False)], columns=["drug", "target", "mechanism"]),
         "drug_indications": pd.DataFrame([("DrugCentral:1", "MONDO:1")], columns=["drug", "condition"]),
         "target_gaps": pd.DataFrame([("t:A", True, np.nan), ("t:B", True, np.nan), ("t:C", True, 1.0), ("t:D", True, np.nan), ("t:E", False, np.nan)], columns=["id", "lacks_ligand", "flips"])}
    S = same_condition(T)
    ck("a trait is read under its id, its current term and a condition of the same name",
       S[S.asked_id == "EFO:1"].read_id.tolist() == ["EFO:1", "MONDO:1"] and S[S.asked_id == "MONDO:2"].values.tolist() == [["MONDO:2", "MONDO:2", "as published"], ["MONDO:2", "UMLS:9", "the sources give it the same name"]])
    ck("two ids with one name are one condition, asked by either", sorted(S[S.asked_id == "UMLS:7"].read_id) == ["MESH:7", "UMLS:7"] == sorted(S[S.asked_id == "MESH:7"].read_id) and "MONDO:1" not in set(S.asked_id))
    R = build(T)
    s = R["disease_gap_summary"].set_index("condition_id")
    a = s.loc["EFO:1"]
    ck("ties: three by a GWAS link (two entries of one gene are one), one by two sources, one as a mechanism target", (a.genes_with_a_gwas_link, a.genes_asserted_by_two_or_more_sources, a.genes_that_are_a_mechanism_target, a.genes_tied) == (3, 1, 1, 4))
    ck("one source alone, and a target a drug acts on outside its mechanism, are not ties", a.genes_asserted_by_one_source_only == 2 and "DDD" not in set(R["disease_gap_genes"].query("condition_id == 'EFO:1' and list == 'tied'").symbol))
    ck("no compound: tied, the gap on record and no drug acting on it", sorted(R["disease_gap_genes"].query("condition_id == 'EFO:1' and list == 'no_compound'").symbol) == ["AAA", "CCC"] and a.no_compound == 2 and a.dark == 1)
    ck("a locus is linked genes that are neighbours, and other evidence singles out one", R["disease_loci"].values.tolist() == [["EFO:1", 1, 2, "AAA, BBB", "BBB"]] and (a.loci_not_resolved, a.genes_in_those_loci, a.loci_where_other_evidence_singles_out_one_gene) == (1, 2, 1))
    ck("tied with no GWAS link; a withheld row; an observation with nothing else", (a.no_genetic_support, a.withheld_rows, a.observed_only) == (1, 1, 0) and sorted(R["disease_gap_genes"].query("condition_id == 'EFO:1' and list == 'no_genetic_support'").symbol) == ["CCC"])
    ck("the largest study is the one with most cases and controls, and its floor is an odds ratio", (a.largest_study, a.cases, a.controls, a.smallest_odds_ratio_it_could_find) == ("S2", 10000, 10000, 1.17))
    b = s.loc["MONDO:2"]
    ck("a study of individuals gives a floor in standard deviations", (b.largest_study, b.individuals, b.smallest_effect_it_could_find_in_sd) == ("S3", 100000, 0.0352) and pd.isna(b.smallest_odds_ratio_it_could_find))
    g = s.loc["UMLS:7"]
    ck("a condition that is not a GWAS trait has no locus and no absent link to count", g.is_a_gwas_trait == 0 and pd.isna(g.loci_not_resolved) and pd.isna(g.no_genetic_support) and g.genes_asserted_by_two_or_more_sources == 1
       and s.loc["MESH:7"].genes_tied == 1)
    ck("a measurement is not among the diseases", "EFO:9" not in s.index)
    p = pair(T, S, "t:D", "EFO:1")
    ck("a pair: its sources, its withheld row with the papers counted, its drug path and whether through the mechanism",
       (p["sources"], p["withheld_rows"], p["papers"], p["drugs"], p["mechanism"], p["observations"], p["directions"]) == (1, 1, [2], 1, False, 1, ["higher"]))
    ck("a drug: its targets and how many are its mechanism", drug(T, "DrugCentral:1") == {"targets": 2, "mechanism_targets": 1, "indications": 1})
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else 0)
