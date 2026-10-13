"""Gaps for a disease: the targets the graph ties to a condition, what each lacks, and what would supply it.

Every part is read from the merged graph and every rate is one an earlier section measured. Nothing
here is fitted. The definitions are in docs/PROTOCOL.md, "Gaps for a disease".

Run:  python -m pharos_gap_kg.ask gaps "Parkinson disease"
      python -m pharos_gap_kg.diseasegaps census pharos319     # every disease trait, counted by kind of gap
      python -m pharos_gap_kg.diseasegaps --selftest

A target is TIED to the condition by a dated GWAS link, by assertions from two or more disease
sources, or by being the mechanism target of a drug indicated for it. One source alone, an expression
observation and a withheld row are counted and are not ties.

    no_compound          tied; no drug acts on it and it has no compound at its family's potency cutoff
    dark                 tied and Tdark
    locus_not_resolved   GWAS-linked genes within 500,000 bases of one another: one signal, several genes
    no_genetic_support   tied by sources or by a drug, with no GWAS link
    withheld             rows a pruning rule withheld for this condition
    observed_only        an expression observation and nothing else

A gap is evidence the graph does not hold. It is not a finding that the target matters in the
disease, and closing it may show that it does not.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from pharos_gap_kg import others, studysize

MIN_SOURCES = 2
READING = ("A gap is evidence this graph does not hold. It is not a finding that the target matters in the disease, and closing it may show that it does not. "
           "Two or more sources, a GWAS link or a drug's mechanism tie a target to the condition; one source alone does not.")
KINDS = ("no_compound", "dark", "locus_not_resolved", "no_genetic_support", "withheld", "observed_only")
NEEDS_A_DATED_LINK = ("locus_not_resolved", "no_genetic_support")
NOT_ASSERTIONS = ("Expression Atlas", "DrugCentral Indication")       # an observation and a path through a drug: neither asserts a gene with a condition


def loci(linked: list[int], neigh: dict, symbol) -> list[list[int]]:
    """Linked targets gathered into loci by the neighbour relation, one entry per gene symbol, loci of
    two or more genes only. `linked` are target indices; `neigh[i]` lists (neighbour, bases apart)."""
    inside = set(linked)
    parent = {i: i for i in linked}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in linked:
        for o, _ in neigh.get(i, ()):
            if o in inside:
                parent[find(i)] = find(o)
    groups: dict[int, list[int]] = {}
    for i in linked:
        groups.setdefault(find(i), []).append(i)
    out = []
    for members in groups.values():
        seen, one = set(), []
        for i in sorted(members):
            if symbol(i) not in seen:
                seen.add(symbol(i))
                one.append(i)
        if len(one) >= 2:
            out.append(one)
    return out


def study_floor(G, trait_id: str, y: int) -> dict | None:
    """The studies behind a trait's links as of the year, the largest of them, and the smallest effect it could have found."""
    if G.sizes is None:
        return None
    G._evidence(0, 0)                                              # reads the evidence table once
    e = G._ev[(G._ev.trait == trait_id) & (G._ev.year <= y)].drop_duplicates("study")
    if e.empty:
        return None
    sized = [(s, int(yr), G.sizes[s]) for s, yr in zip(e.study, e.year) if s in G.sizes]
    out = {"studies_behind_its_links": int(len(e)), "with_a_size": len(sized)}
    if not sized:
        return out
    s, yr, n = max(sized, key=lambda r: (studysize.effective(*r[2]), r[0]))
    floor = studysize.smallest_effect(*n)
    out["largest"] = {"study": s, "year": yr, **floor}
    out["reading"] = (f"The largest study behind its links ({s}, {yr}): {studysize.say(floor)}. An absent link is no effect of about that size found at a common variant; "
                      "smaller effects and rarer variants were not visible. This is a rough floor, and a study that reported nothing for the trait is not in the source.")
    return out


def for_condition(G, r: dict, y: int, top: int = 5) -> dict:
    """G: ask.Graph. r: a resolved trait or condition. y: the year the dated links are read as of."""
    O = G.others
    top = max(1, min(int(top), 25))
    sym = lambda i: G.t["name"][i]
    is_trait = r["kind"] == "trait"
    if is_trait:
        di = r["i"]
        lab, ids, mask = {"kind": "trait", **G._dlabel(di)}, O.ids("trait", G.d.id[di]), G._mask(y)
        gw = {G.t.id[G.lt[k]]: int(k) for k in np.flatnonzero(mask & (G.ld == di))}
    else:
        lab, ids, gw = {"kind": "condition", **O.label("condition", r["key"])}, O.ids("condition", r["key"]), {}
    name = lab["name"]
    pos = lambda idx: [p for c in ids for p in idx.get(c, ())]
    a = O.A.iloc[pos(O.a_c)]
    kept, held = a[~a.withheld], a[a.withheld]
    src = {t: sorted(set(s)) for t, s in kept.groupby("target").source}
    two = {t for t, s in src.items() if len(s) >= MIN_SOURCES}
    one = {t for t, s in src.items() if len(s) < MIN_SOURCES}
    drugs = list(dict.fromkeys(O.DI.iloc[pos(O.di_c)].drug))
    moa: dict[str, list[str]] = {}
    for d in drugs:
        for row in O.DT.iloc[O.dt_d.get(d, [])].itertuples():
            if row.mechanism:
                moa.setdefault(row.target, []).append(O.name.get(d, d))
    x = O.X.iloc[pos(O.x_c)]
    tied = (set(gw) | two | set(moa)) & set(G.ti)
    studies = lambda t: int(G.n_studies[gw[t]]) if t in gw else 0
    strength = lambda t: (t not in gw, -studies(t), -len(src.get(t, ())), t not in moa, sym(G.ti[t]), t)

    def label(t):
        return {k: v for k, v in G._tlabel(G.ti[t]).items() if k in ("id", "symbol", "pharos_level", "pharos_entries_with_this_symbol")}

    def tie(t):
        out = {}
        if t in gw:
            k = gw[t]
            out["gwas_link"] = {"studies": int(G.n_studies[k]), "first_reported": int(G.first[k]), "edge": G.edge_id[k]}
        if t in src:
            out["sources"] = src[t]
        if t in moa:
            out["mechanism_target_of"] = sorted(set(moa[t]))
        return out

    def listed(ts):
        """Targets strongest first, one per gene symbol."""
        seen, out = set(), []
        for t in sorted(ts, key=strength):
            if sym(G.ti[t]) not in seen:
                seen.add(sym(G.ti[t]))
                out.append(t)
        return out

    gaps = []
    # ---- no compound, and dark: read from the target-level gaps
    if G.tg is not None:
        g = {t: G._gap_row(t) for t in tied}
        nc = listed(t for t in tied if g[t] is not None and g[t]["lacks_ligand"] is True and t not in O.dt_t)
        rate = G._ligand_gap(g[nc[0]])["how_often_this_closed"] if nc else ""
        gaps.append({"kind": "no_compound", "what": "tied to the condition, with no drug acting on it and no compound at its family's potency cutoff", "targets": len(nc),
                     "top": [{**label(t), "tied_by": tie(t), "lacks": G._ligand_gap(g[t])["lacks"], **G._edge_of(g[t], "ligand")} for t in nc[:top]], "more": max(0, len(nc) - top),
                     "closes_with": "a compound active against the target at the potency cutoff of its family", "how_often_this_closed": rate})
        dk = listed(t for t in tied if g[t] is not None and pd.notna(g[t]["flips"]))
        bands = G.tg_rates["dark"]["bands"]
        gaps.append({"kind": "dark", "what": "tied to the condition and still Tdark, Pharos's least-studied level", "targets": len(dk),
                     "top": [{**label(t), "tied_by": tie(t), "nearest_route": g[t]["proposal"], "criteria_to_clear": int(g[t]["flips"]),
                              "promoted_in_its_band": f"{bands[g[t]['band']]['promoted']['share']:.0%}", **G._edge_of(g[t], "dark")} for t in dk[:top]], "more": max(0, len(dk) - top),
                     "closes_with": "the nearest route named for each target, or an experimental Gene Ontology annotation",
                     "how_often_this_closed": (f"By band, the share of such targets that were no longer Tdark in {G.tg_rates['later']} is given with each. Closing this gap raises the "
                                               "target's grade. It is not the same as learning its part in the disease.")})
    # ---- loci the link does not resolve to one gene
    floor = None
    if is_trait:
        L = loci([G.ti[t] for t in gw], G.neigh, sym)
        rows = []
        for members in L:
            ids_ = [G.t.id[i] for i in members]
            other = {sym(G.ti[t]): {k: v for k, v in tie(t).items() if k != "gwas_link"} for t in ids_ if t in two or t in moa}
            rows.append({"genes": [sym(i) for i in sorted(members, key=lambda i: (-studies(G.t.id[i]), sym(i)))], "genes_with_the_link": len(members),
                         "other_evidence_for": other, "singled_out_by_other_evidence": next(iter(other)) if len(other) == 1 else None,
                         "in_the_MHC": bool(G.mhc[members].any())})
        rows.sort(key=lambda r_: (r_["singled_out_by_other_evidence"] is not None, -r_["genes_with_the_link"], r_["genes"][0]))
        single = sum(r_["singled_out_by_other_evidence"] is not None for r_ in rows)
        gaps.append({"kind": "locus_not_resolved", "what": "GWAS-linked genes within 500,000 bases of one another: one association signal assigned to several genes",
                     "loci": len(rows), "genes": sum(r_["genes_with_the_link"] for r_ in rows), "loci_where_other_evidence_singles_out_one_gene": single,
                     "loci_in_the_MHC": sum(r_["in_the_MHC"] for r_ in rows), "top": [{**r_, "genes": r_["genes"][:8], "genes_more": max(0, len(r_["genes"]) - 8)} for r_ in rows[:top]],
                     "more": max(0, len(rows) - top),
                     "closes_with": ("evidence for one gene of the locus: fine-mapping, colocalisation with expression, or a functional result. In this graph that is a second "
                                     "source or a drug's mechanism on exactly one of them"),
                     "how_often_this_closed": "Not measured."})
        # ---- tied without a GWAS link
        ng = listed((two | set(moa)) & set(G.ti) - set(gw))
        sc, _, _ = G._score_trait(di, mask)
        cand = np.ones(len(G.t), bool)
        cand[G.lt[mask & (G.ld == di)]] = False
        floor = study_floor(G, G.d.id[di], y)

        def later(t):
            k = G._row(G.ti[t], di)
            return {"gwas_link_first_reported_later": int(G.first[k])} if k is not None and G.first[k] > y else {}

        rate = (G.backtest or {}).get("sources_then_gwas")
        how = "Not measured."
        if rate:
            asserting = {k: v["share"] for k, v in rate["by_source"].items() if k not in NOT_ASSERTIONS}
            by = {k: v for k, v in asserting.items() if k in set(kept.source)} or asserting
            lo, hi = min(by, key=by.get), max(by, key=by.get)
            how = (f"In the test on links after {rate['cutoff']}, over {rate['disease_traits']} disease traits, a pair one of these sources named with no GWAS link was "
                   f"GWAS-linked later between {by[lo]:.1%} of the time ({lo}) and {by[hi]:.1%} ({hi}), by source, against {rate['base_rate']:.1%} of all pairs. {rate['reading']}")
        gaps.append({"kind": "no_genetic_support", "what": "tied by two or more sources or by a drug's mechanism, with no GWAS link" + (f" as of {y}" if y < G.latest else ""),
                     "targets": len(ng),
                     "top": [{**label(t), "tied_by": tie(t), "rank_among_gwas_candidates": (int((sc[cand] > sc[G.ti[t]]).sum()) + 1) if sc[G.ti[t]] > 0 else None, **later(t)} for t in ng[:top]],
                     "more": max(0, len(ng) - top), "closes_with": "a genome-wide significant association at the gene, in a study of this trait",
                     "how_often_this_closed": how, **({"what_an_absent_link_could_have_seen": floor} if floor else {})})
    # ---- withheld rows, and observations with nothing else
    hrows = sorted(held.itertuples(), key=lambda r_: (sym(G.ti[r_.target]) if r_.target in G.ti else r_.target, r_.source))
    gaps.append({"kind": "withheld", "what": "rows a pruning rule withheld for this condition", "rows": len(held),
                 "top": [{**(label(r_.target) if r_.target in G.ti else {"id": r_.target}), **O._withheld(r_)} for r_ in hrows[:top]], "more": max(0, len(held) - top),
                 "closes_with": "what each row's rule says restores it", "how_often_this_closed": "Not measured."})
    only = x[~x.target.isin(tied | one) & x.target.isin(set(G.ti))]
    only = only.assign(symbol=[sym(G.ti[t]) for t in only.target])
    best = only.reindex(only.log2fc.abs().sort_values(ascending=False).index).drop_duplicates("symbol")      # one row per gene: its largest difference
    gaps.append({"kind": "observed_only", "what": "an expression observation in the condition and nothing else: no tie and no source", "targets": int(only.symbol.nunique()),
                 "top": [{**(label(r_.target) if r_.target in G.ti else {"id": r_.target}), "direction": "higher" if r_.log2fc > 0 else "lower", "log2_fold_change": round(float(r_.log2fc), 3),
                          "p": float(r_.p), "edge": r_.edge} for r_ in best[:top].itertuples()], "more": max(0, int(only.symbol.nunique()) - top),
                 "closes_with": "a source that asserts the gene with the condition, or a genetic association", "how_often_this_closed": "Not measured.", "reading": others.OBSERVATION})
    ties, tied_say = O.read_under(ids)
    genes = lambda ts: len({sym(G.ti[t]) for t in ts if t in G.ti})          # Pharos holds some genes as several entries
    held_ = {"gwas_linked_targets": len(gw) if is_trait else None, "gwas_linked_genes": genes(gw) if is_trait else None,
             "genes_asserted_by_two_or_more_sources": genes(two), "genes_asserted_by_one_source_only": genes(one),
             "drugs_indicated": len(drugs), "genes_that_are_the_mechanism_target_of_one": genes(moa), "genes_with_an_expression_observation": genes(set(x.target)),
             "withheld_rows": len(held), "tied_genes": genes(tied)}
    by = {g_["kind"]: g_ for g_ in gaps}
    n = lambda k, f="targets": by[k][f] if k in by else 0
    parts = []
    if is_trait:
        parts.append(f"{held_['gwas_linked_genes']:,} by a GWAS link" + (f" as of {y}" if y < G.latest else ""))
    parts += [f"{held_['genes_asserted_by_two_or_more_sources']:,} by two or more disease sources", f"{held_['genes_that_are_the_mechanism_target_of_one']:,} as the mechanism target of an indicated drug"]
    say = f"The graph ties {held_['tied_genes']:,} gene(s) to {name}, some in more than one way: " + ", ".join(parts) + ". Of them, "
    bits = []
    if "no_compound" in by:
        bits += [f"{n('no_compound'):,} have no drug and no compound at their family's cutoff", f"{n('dark'):,} {'is' if n('dark') == 1 else 'are'} still Tdark"]
    if is_trait:
        lr = by["locus_not_resolved"]
        if lr["loci"]:
            bits.append(f"{lr['genes']:,} GWAS-linked genes sit in {lr['loci']:,} {'locus' if lr['loci'] == 1 else 'loci'} that carry the link on more than one gene, and other evidence "
                        f"singles out one gene in {lr['loci_where_other_evidence_singles_out_one_gene']:,} of them")
        else:
            bits.append("no locus carries the link on more than one gene")
        bits.append(f"{n('no_genetic_support'):,} {'is' if n('no_genetic_support') == 1 else 'are'} tied with no GWAS link")
    say += "; ".join(bits) + f". {len(held):,} row(s) are withheld by a pruning rule, and {n('observed_only'):,} gene(s) have only an expression observation."
    if not is_trait:
        say += f" {name} is not a GWAS trait in this graph, so no dated link is present or absent for it and the kinds that need one are not listed."
    elif not G.d.disease_like[r["i"]]:
        say += f" {name} is a measurement and not a disease, so the disease sources and the drugs seldom name it."
    if floor and "reading" in floor:
        say += " " + floor["reading"]
    out = {"found": True, "condition": lab, "as_of": y, "read_under": ties, "holds": held_, "gaps": gaps,
           "not_listed": [] if is_trait else list(NEEDS_A_DATED_LINK), "answer": say + tied_say + " " + others.UNDATED.replace("this is", "all but the GWAS links are"), "reading": READING}
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    names = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "HLA-X", 6: "HLA-X", 7: "F"}
    neigh = {0: [(1, 10)], 1: [(0, 10), (2, 20)], 2: [(1, 20)], 3: [(4, 5)], 4: [(3, 5)], 5: [(6, 1)], 6: [(5, 1)], 7: [(0, 400)]}
    L = loci([0, 1, 2, 3, 5, 6], neigh, names.get)
    ck("genes chained by the neighbour relation are one locus", sorted(map(sorted, L)) == [[0, 1, 2]])
    ck("a linked gene whose neighbour is not linked is not a locus", not any(3 in m for m in L))
    ck("two entries of one gene symbol are one gene, and a locus needs two", not any(5 in m or 6 in m for m in L))
    ck("a neighbour that is linked joins the locus, from either side", sorted(map(sorted, loci([0, 1, 2, 7], {**neigh, 0: [(1, 10), (7, 400)]}, names.get))) == [[0, 1, 2, 7]])
    ck("no linked target, no locus", loci([], neigh, names.get) == [])
    ck("the kinds that need a dated link are among the kinds", set(NEEDS_A_DATED_LINK) < set(KINDS) and MIN_SOURCES == 2)
    return 0 if ok else 1


def census(release: str) -> dict:
    """Every disease-like GWAS trait with a link, asked once: how many gaps of each kind the graph holds.
    A gene tied to two traits is counted under each."""
    import json

    from pharos_gap_kg import ask, paths
    G = ask.Graph(paths.bundle(release))
    if not G.has_others or G.tg is None:
        raise SystemExit("the bundle does not hold the merged graph's other parts and the target-level gaps")
    deg = np.bincount(G.ld, minlength=len(G.d))
    rows = []
    for i in np.flatnonzero(G.d.disease_like.to_numpy() & (deg > 0)):
        o = for_condition(G, {"kind": "trait", "i": int(i)}, G.latest, top=1)
        k, h = {g["kind"]: g for g in o["gaps"]}, o["holds"]
        f = k["no_genetic_support"].get("what_an_absent_link_could_have_seen") or {}
        big = f.get("largest") or {}
        rows.append({"trait": G.d.id[i], "name": G.d["name"][i], "tied": h["tied_genes"], "gwas": h["gwas_linked_genes"], "two_sources": h["genes_asserted_by_two_or_more_sources"],
                     "mechanism": h["genes_that_are_the_mechanism_target_of_one"], "no_compound": k["no_compound"]["targets"], "dark": k["dark"]["targets"],
                     "loci": k["locus_not_resolved"]["loci"], "locus_genes": k["locus_not_resolved"]["genes"], "loci_singled_out": k["locus_not_resolved"]["loci_where_other_evidence_singles_out_one_gene"],
                     "no_genetic_support": k["no_genetic_support"]["targets"], "withheld": k["withheld"]["rows"], "observed_only": k["observed_only"]["targets"],
                     "odds_ratio": big.get("smallest_odds_ratio_it_could_find"), "sd": big.get("smallest_effect_it_could_find_in_standard_deviations_per_allele"), "sized": bool(big)})
    t = pd.DataFrame(rows)
    kinds = ["no_compound", "dark", "locus_genes", "no_genetic_support", "withheld", "observed_only"]
    cc = t.odds_ratio.dropna()
    man = {"release": release, "as_of": G.latest, "traits": len(t), "min_sources": MIN_SOURCES,
           "tied_genes_summed_over_traits": int(t.tied.sum()), "gwas_linked_genes_summed": int(t.gwas.sum()),
           "summed_over_traits": {c: int(t[c].sum()) for c in kinds + ["loci", "loci_singled_out", "two_sources", "mechanism"]},
           "traits_with_at_least_one": {c: int((t[c] > 0).sum()) for c in kinds + ["loci", "two_sources", "mechanism"]},
           "median_per_trait": {c: float(t[c].median()) for c in ["tied", "gwas", "no_compound", "locus_genes", "no_genetic_support"]},
           "share_of_gwas_linked_genes_in_a_locus_with_another": round(float(t.locus_genes.sum() / t.gwas.sum()), 4),
           "share_of_such_loci_where_other_evidence_singles_out_one_gene": round(float(t.loci_singled_out.sum() / t.loci.sum()), 4),
           "share_of_tied_genes_with_no_compound": round(float(t.no_compound.sum() / t.tied.sum()), 4),
           "traits_whose_largest_study_has_a_size": int(t.sized.sum()), "of_them_cases_and_controls": int(len(cc)),
           "smallest_odds_ratio_the_largest_study_could_find": {"median": round(float(cc.median()), 3), "lower_quartile": round(float(cc.quantile(0.25)), 3), "upper_quartile": round(float(cc.quantile(0.75)), 3),
                                                                "traits_above_1.5": int((cc > 1.5).sum())},
           "most_tied": [{"trait": r_.trait, "name": r_.name, "tied": int(r_.tied), "no_compound": int(r_.no_compound), "no_genetic_support": int(r_.no_genetic_support)}
                         for r_ in t.sort_values(["tied", "trait"], ascending=[False, True]).head(5).itertuples()]}
    (paths.DOCS_DATA / f"disease_gaps_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def main(argv: list[str] | None = None) -> int:
    import sys
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 2 or argv[0] != "census":
        raise SystemExit("Run:  python -m pharos_gap_kg.diseasegaps census pharos319\n      python -m pharos_gap_kg.diseasegaps --selftest")
    m = census(argv[1])
    print(f"{m['release']}: {m['traits']} disease traits; tied genes summed {m['tied_genes_summed_over_traits']:,}; by kind {m['summed_over_traits']}")
    print(f"  in a locus with another linked gene: {m['share_of_gwas_linked_genes_in_a_locus_with_another']:.1%} of GWAS-linked genes; other evidence singles one out in "
          f"{m['share_of_such_loci_where_other_evidence_singles_out_one_gene']:.1%} of those loci; smallest odds ratio the largest study could find: {m['smallest_odds_ratio_the_largest_study_could_find']}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
