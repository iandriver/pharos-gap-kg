"""Apply the lift cut to the sources it passed on, and write what is withheld.

The second pruning round (findings/PRUNING2.md) measured, on planted pairs, a rule that withholds a
(target, disease) pair PubMed mentions together no more often than a cut for its source. It passed on
five sources. This applies it to every one-source pair of those five, with the cut that was scored,
and writes each withheld pair with the count of papers, the lift and the cut. How it is applied was
fixed in docs/PROTOCOL.md before any pair outside the sample was searched.

Monarch also passed the in-graph rule of round one. Its list is recomputed here on the real table
alone, and a pair either rule withholds is withheld, naming which.

Run:  python -m pharos_gap_kg.applycut search pharos319    # hours at NCBI's limit; every answer is kept, so it resumes
      python -m pharos_gap_kg.applycut write pharos319     # the withheld list, the pruned table, the summary
      python -m pharos_gap_kg.applycut --selftest
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from pharos_gap_kg import graph, ontology, paths, planted as round1, prune2, registry

KEY = prune2.KEY
LIFT, INGRAPH = "literature_lift", "uncorroborated"
REASON = {LIFT: "PubMed mentions this gene with this disease no more often than the cut for its source",
          INGRAPH: "one source, and nothing else in the graph corroborates it"}
RESTORED_BY = {LIFT: "a lift above the cut for its source, or a second source",
               INGRAPH: "a second source, the same target with a broader or narrower disease, a GWAS link, or an interaction partner with this disease"}
UNION_AT_MOST = round1.LOSE_AT_MOST


def cuts(release: str) -> dict[str, float]:
    """source -> the cut that was scored for it, for the sources on which the lift cut passed."""
    f = paths.DOCS_DATA / f"pruning2_{release}.json"
    if not f.is_file():
        raise SystemExit(f"the second pruning round has no report; run:  python -m pharos_gap_kg.prune2 report {release}")
    rows = json.loads(f.read_text())["lift_cut_L0b"]
    return {s: float(r["cut"]) for s, r in rows.items() if s != "all sources" and r["passes"]}


def candidates(release: str, sources: list[str]) -> pd.DataFrame:
    links = pd.read_parquet(paths.context(release) / "disease_links.parquet")
    x = prune2.one_source(links)
    x = x[x.source.isin(sources)].reset_index(drop=True)
    name_of, symbol_of = prune2.names(release)
    x["symbol"] = x.target.map(symbol_of)
    x["name"] = [name_of.get(k) for k in zip(x.source, x.disease)]
    return x


def search(release: str) -> pd.DataFrame:
    c = cuts(release)
    x = candidates(release, sorted(c))
    found, n_searches, total = prune2.counts(x)
    out = pd.concat([x, found], axis=1)
    out["cut"] = out.source.map(c)
    out["lit_withheld"] = decide(out.lift, out.unsearchable, out.cut)
    W = prune2.work(release)
    out.to_parquet(W / "applied.parquet", index=False)
    (W / "applied.json").write_text(json.dumps({"pubmed_records": total, "searches": n_searches}, indent=1) + "\n")
    return out


def decide(lift: pd.Series, unsearchable: pd.Series, cut: pd.Series) -> np.ndarray:
    """Withheld: searchable, and a lift at or below the cut. An unsearchable pair is never withheld."""
    return (unsearchable.isna() & (lift.astype(float) <= cut.astype(float))).to_numpy()


def same_as_sample(applied: pd.DataFrame, sample_search: pd.DataFrame) -> int:
    """Each sampled real pair of an applied source must have the papers and the lift it had in the sample."""
    s = sample_search[~sample_search.planted & sample_search.source.isin(set(applied.source))]
    m = s.merge(applied, on=KEY, suffixes=("_sample", ""), how="left")
    if m.cut.isna().any():
        raise SystemExit(f"REFUSED: {int(m.cut.isna().sum())} sampled real pairs were not examined when the cut was applied")
    differ = (m.both_sample.fillna(-1) != m.both.fillna(-1)) | ~np.isclose(m.lift_sample.fillna(-1).astype(float), m.lift.fillna(-1).astype(float))
    if differ.any():
        raise SystemExit(f"REFUSED: {int(differ.sum())} sampled real pairs have other counts now than in the sample, e.g. {m[differ][KEY].iloc[0].tolist()}")
    return len(m)


def in_graph_list(release: str, sources: list[str]) -> pd.DataFrame:
    """The one-source pairs of `sources` that nothing in the REAL table corroborates (round one's four criteria)."""
    C = paths.context(release)
    links = pd.read_parquet(C / "disease_links.parquet")
    partners = pd.read_parquet(C / "partners.parquet")
    gwas = pd.read_parquet(paths.bundle(release) / "links.parquet", columns=["target", "trait"])
    terms = ontology.read(registry.require("efo_obo"))
    rel, onto = round1.lineage_of(sorted(set(links.disease)), terms)
    single = prune2.one_source(links)
    c = round1.corroboration(single[single.source.isin(sources)], links, rel, onto, gwas, partners, terms)
    return c[~c.corroborated][KEY].reset_index(drop=True)


def union_on_sample(release: str, source: str, real_ingraph: pd.DataFrame) -> dict:
    """On the sampled pairs of `source` that carry both verdicts: what each rule, either and both withhold."""
    W = prune2.work(release)
    se = pd.read_parquet(W / "search.parquet")
    se = se[(se.source == source) & se.unsearchable.isna()].assign(lit=lambda d: (d.both == 0).to_numpy())
    det = pd.read_parquet(W / "ingraph_detail.parquet")
    pl = se[se.planted].merge(det[det.planted & (det.source == source)][KEY + ["corroborated"]], on=KEY)
    re = se[~se.planted].merge(real_ingraph.assign(unc=True), on=KEY, how="left")
    re["corroborated"] = re.unc.isna()
    out = {}
    for name, g in (("planted", pl), ("real", re)):
        ing = (~g.corroborated).to_numpy()
        lit = g.lit.to_numpy()
        out[name] = {"n": int(len(g)), "in_graph": prune2.share(ing), "literature": prune2.share(lit), "either": prune2.share(ing | lit), "both": prune2.share(ing & lit)}
    out["union_applies"] = bool(out["real"]["either"]["share"] <= UNION_AT_MOST)
    return out


def combine(lit: pd.DataFrame, ingraph: pd.DataFrame) -> pd.DataFrame:
    """One row per withheld pair: its rules, why, what restores it, and the literature figures where they were read."""
    a = lit[lit.lit_withheld][KEY].assign(_lit=True)
    b = ingraph[KEY].assign(_ing=True)
    w = a.merge(b, on=KEY, how="outer")
    w["rules"] = [("|".join(r for r, on in ((LIFT, x), (INGRAPH, y)) if on is True)) for x, y in zip(w._lit, w._ing)]
    w["reason"] = w.rules.map(lambda r: "; ".join(REASON[k] for k in r.split("|")))
    w["restored_by"] = w.rules.map(lambda r: "; or ".join(RESTORED_BY[k] for k in r.split("|")))
    figs = lit[KEY + ["both", "lift", "cut"]].rename(columns={"both": "papers"})
    w = w.drop(columns=["_lit", "_ing"]).merge(figs, on=KEY, how="left")
    return w.sort_values(KEY).reset_index(drop=True)


def write(release: str) -> dict:
    C, W = paths.context(release), prune2.work(release)
    c = cuts(release)
    f = W / "applied.parquet"
    if not f.is_file():
        raise SystemExit(f"nothing has been searched; run:  python -m pharos_gap_kg.applycut search {release}")
    lit = pd.read_parquet(f)
    x = candidates(release, sorted(c))
    if len(lit) != len(x) or not lit[KEY].equals(x[KEY]):
        raise SystemExit("REFUSED: the searched pairs are not the one-source pairs of the sources the cut passed on; search again")
    checked = same_as_sample(lit, pd.read_parquet(W / "search.parquet"))
    rep = json.loads((paths.DOCS_DATA / f"pruning2_{release}.json").read_text())
    ing_sources = rep["in_graph_rule_on_matched_planting"]["rule_applied_to"]
    before = pd.read_parquet(C / "disease_links_withheld.parquet")
    ingraph = in_graph_list(release, ing_sources)
    unions = {s: union_on_sample(release, s, ingraph[ingraph.source == s]) for s in ing_sources if s in c}
    dropped = []
    for s, u in unions.items():
        if not u["union_applies"]:       # over the bar's cost: keep only the rule that withholds fewer real pairs
            if u["real"]["in_graph"]["share"] <= u["real"]["literature"]["share"]:
                lit.loc[lit.source == s, "lit_withheld"] = False
                dropped.append((s, LIFT))
            else:
                ingraph = ingraph[ingraph.source != s]
                dropped.append((s, INGRAPH))
    w = combine(lit, ingraph)
    links = pd.read_parquet(C / "disease_links.parquet")
    kept = links.merge(w[KEY].assign(_w=1), on=KEY, how="left")
    kept = kept[kept._w.isna()].drop(columns="_w")
    w.to_parquet(C / "disease_links_withheld.parquet", index=False)
    kept.to_parquet(C / "disease_links_pruned.parquet", index=False)

    per = {}
    for s in sorted(c):
        g = lit[lit.source == s]
        ok = g[g.unsearchable.isna()]
        exp = rep["lift_cut_L0b"][s]["real_withheld"]
        sh = prune2.share(ok.lit_withheld.to_numpy())
        ws = w[w.source == s]
        per[s] = {"cut": c[s], "one_source_pairs": int(len(g)), "unsearchable": int(g.unsearchable.notna().sum()), "withheld_by_the_cut": int(g.lit_withheld.sum()),
                  "share_of_searchable_withheld": sh["share"], "share_on_the_pairs_the_cut_was_scored_on": exp,
                  "inside_that_interval": bool(exp["low"] <= sh["share"] <= exp["high"]),
                  "withheld_in_all": int(len(ws)), "by_rule": {k: int(v) for k, v in ws.rules.value_counts().items()},
                  "share_of_one_source_pairs_withheld_in_all": float(len(ws) / len(g))}
    out = {"release": release, "pubmed": json.loads((W / "applied.json").read_text()), "sources": per,
           "sampled_real_pairs_checked_against_the_sample": checked,
           "in_graph_list": {"sources": ing_sources, "as_round_one_wrote_it": int(len(before)) if "rules" not in before else None, "on_the_real_table_alone": int(len(ingraph))},
           "two_rules_on_one_source": unions, "rule_dropped_for_cost": dropped,
           "withheld": int(len(w)), "rows_before": int(len(links)), "rows_kept": int(len(kept)),
           "files": {n: graph.file_sha(C / n) for n in ("disease_links_pruned.parquet", "disease_links_withheld.parquet")},
           "disease_links_sha256": json.loads((C / "manifest.json").read_text())["tables"]["disease_links"]["sha256"]}
    prev = paths.DOCS_DATA / f"pruning2_applied_{release}.json"
    if out["in_graph_list"]["as_round_one_wrote_it"] is None and prev.is_file():      # a re-run: keep the count round one wrote
        out["in_graph_list"]["as_round_one_wrote_it"] = json.loads(prev.read_text())["in_graph_list"]["as_round_one_wrote_it"]
    prev.write_text(prune2.dump(out))
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    lit = pd.DataFrame({"target": list("abcde"), "disease": ["x"] * 5, "source": ["s"] * 4 + ["m"], "both": pd.array([0, 3, 9, None, 0], "Int64"),
                        "lift": [0.0, 1.5, 4.0, None, 0.0], "unsearchable": [None, None, None, "PubMed did not find the disease name", None], "cut": [2.0, 2.0, 2.0, 2.0, 0.0]})
    lit["lit_withheld"] = decide(lit.lift, lit.unsearchable, lit.cut)
    ck("a pair at or below its source's cut is withheld, one above it is kept, and an unsearchable pair is never withheld", lit.lit_withheld.tolist() == [True, True, False, False, True])
    ing = pd.DataFrame({"target": ["e", "f"], "disease": ["x", "x"], "source": ["m", "m"]})
    w = combine(lit, ing).set_index("target")
    ck("a pair either rule withholds is withheld once, naming which", w.rules.to_dict() == {"a": LIFT, "b": LIFT, "e": f"{LIFT}|{INGRAPH}", "f": INGRAPH})
    ck("a gap from the cut carries the papers, the lift and the cut; one from the graph alone carries none",
       (int(w.papers["b"]), w.lift["b"], w.cut["b"]) == (3, 1.5, 2.0) and pd.isna(w.papers["f"]) and pd.isna(w.cut["f"]))
    ck("each rule's reason and what restores the pair are written out", w.reason["e"] == REASON[LIFT] + "; " + REASON[INGRAPH] and w.restored_by["f"] == RESTORED_BY[INGRAPH])
    se = pd.DataFrame({"target": ["a", "c", "z"], "disease": ["x"] * 3, "source": ["s", "s", "other"], "planted": [False] * 3, "both": pd.array([0, 9, 1], "Int64"), "lift": [0.0, 4.0, 1.0]})
    ck("sampled real pairs of an applied source are checked, others are not", same_as_sample(lit, se) == 2)
    for change, what in ((se.assign(lift=[0.0, 4.5, 1.0]), "another lift"), (pd.concat([se, se.iloc[[0]].assign(target="q")]), "a pair that was not examined")):
        try:
            same_as_sample(lit, change)
            ck(f"a sampled pair with {what} stops the write", False)
        except SystemExit:
            ck(f"a sampled pair with {what} stops the write", True)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 2 or argv[0] not in ("search", "write"):
        raise SystemExit(__doc__)
    stage, release = argv
    if stage == "search":
        out = search(release)
        print(f"  {len(out):,} one-source pairs examined; unsearchable {int(out.unsearchable.notna().sum()):,}; at or below the cut {int(out.lit_withheld.sum()):,}")
        return 0
    out = write(release)
    for s, r in out["sources"].items():
        e = r["share_on_the_pairs_the_cut_was_scored_on"]
        print(f"  {s:34s} cut {r['cut']:6.3f}  one-source {r['one_source_pairs']:>6,}  unsearchable {r['unsearchable']:>4,}  withheld by the cut {r['withheld_by_the_cut']:>5,} "
              f"({r['share_of_searchable_withheld']:.1%}; scored {e['share']:.1%} [{e['low']:.1%}, {e['high']:.1%}] {'inside' if r['inside_that_interval'] else 'OUTSIDE'})  in all {r['withheld_in_all']:>5,} {r['by_rule']}")
    g = out["in_graph_list"]
    print(f"  in-graph list {g['sources']}: {g['as_round_one_wrote_it']} as round one wrote it, {g['on_the_real_table_alone']} on the real table alone")
    for s, u in out["two_rules_on_one_source"].items():
        print(f"  {s}, both rules on the sampled pairs: planted n={u['planted']['n']} in-graph {u['planted']['in_graph']['share']:.3f} literature {u['planted']['literature']['share']:.3f} either {u['planted']['either']['share']:.3f} both {u['planted']['both']['share']:.3f}; "
              f"real n={u['real']['n']} in-graph {u['real']['in_graph']['share']:.3f} literature {u['real']['literature']['share']:.3f} either {u['real']['either']['share']:.3f} both {u['real']['both']['share']:.3f}; union applies {u['union_applies']}")
    print(f"  withheld {out['withheld']:,}; rows kept {out['rows_kept']:,} of {out['rows_before']:,}; sampled pairs checked {out['sampled_real_pairs_checked_against_the_sample']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
