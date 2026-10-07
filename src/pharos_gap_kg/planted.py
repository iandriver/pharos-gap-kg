"""Does a pruning rule remove false pairs? False pairs are planted, the rule is run, and both are counted.

The table is the (target, disease) pairs of the undated disease sources. Nobody knows which of its
real pairs are wrong, so a rule cannot be scored on them. It can be scored on pairs made to be wrong:
for each source, two of its pairs are drawn and crossed, giving a target and a disease that are each
as common in that source as a real pair's, and that no source puts together (docs/PROTOCOL.md).

THE RULE UNDER TEST keeps a pair that is corroborated and withholds one that is not. Corroborated is
any of: (a) a second source asserts the pair; (b) some source asserts the target with a disease in
this one's lineage; (c) a GWAS link joins the target to the disease or its lineage; (d) an
interaction partner of the target is asserted with the same disease.

WHAT IS COUNTED: the share of planted pairs withheld, and the share of a source's real single-source
pairs withheld. The second is the rule's cost. Those pairs are not known to be false.

Run:  python -m pharos_gap_kg.planted pharos319
      python -m pharos_gap_kg.planted --selftest
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np
import pandas as pd

from pharos_gap_kg import graph, ontology, paths, registry

SEED = 20261007
FRACTION = 0.05
KEEP_AT_LEAST, LOSE_AT_MOST = 0.80, 0.50      # the bar a source must clear for the rule to be applied to it


def plant(links: pd.DataFrame, fraction: float = FRACTION, seed: int = SEED) -> pd.DataFrame:
    """target, disease, source: false pairs, 5% of each source's pairs, by crossing two of its pairs."""
    rng = np.random.default_rng(seed)
    taken = set(zip(links.target, links.disease))
    out = []
    for s, g in sorted(links.groupby("source"), key=lambda kv: kv[0]):
        g = g.drop_duplicates(["target", "disease"])
        t, d = g.target.to_numpy(), g.disease.to_numpy()
        want, tries = math.ceil(fraction * len(g)), 0
        while want and tries < 200 * math.ceil(fraction * len(g)) + 1000:
            tries += 1
            i, j = rng.integers(0, len(g), 2)
            pair = (t[i], d[j])
            if pair not in taken:
                taken.add(pair)
                out.append((t[i], d[j], s))
                want -= 1
        if want:
            raise SystemExit(f"REFUSED: could not plant {want} more false pairs in {s}; its pairs cross into real ones too often")
    return pd.DataFrame(out, columns=["target", "disease", "source"])


def lineage_of(diseases: list[str], terms: dict) -> tuple[dict, dict]:
    """For each disease id: the other ids in the table in its lineage, and its live term's ancestors."""
    canon = {d: ontology.canonical(d, terms) for d in diseases}
    anc = {d: (ontology.ancestors(c, terms) if c else set()) for d, c in canon.items()}
    by_term: dict[str, list[str]] = {}
    for d, c in canon.items():
        if c:
            by_term.setdefault(c, []).append(d)
    rel: dict[str, set] = {d: set() for d in diseases}
    for d, c in canon.items():
        if not c:
            continue
        rel[d].update(x for x in by_term[c] if x != d)
        for a in anc[d]:
            for x in by_term.get(a, ()):
                rel[d].add(x)
                rel[x].add(d)
    return rel, {"canon": canon, "anc": anc}


def corroboration(pairs: pd.DataFrame, table: pd.DataFrame, rel: dict, onto: dict, gwas: pd.DataFrame, partners: pd.DataFrame, terms: dict) -> pd.DataFrame:
    """For each of `pairs` (target, disease): the four criteria, read against `table` (every pair the pruner holds)."""
    diseases_of = table.groupby("target").disease.agg(set).to_dict()
    holders_of = table.groupby("disease").target.agg(set).to_dict()
    n_sources = table.groupby(["target", "disease"]).source.nunique().to_dict()
    nb: dict[str, set] = {}
    for a, b in zip(partners.a, partners.b):
        nb.setdefault(a, set()).add(b)
        nb.setdefault(b, set()).add(a)
    g_terms: dict[str, set] = {}
    g_up: dict[str, set] = {}
    live = {t: ontology.canonical(t, terms) for t in set(gwas.trait)}
    up = {t: (ontology.ancestors(c, terms) if c else set()) for t, c in live.items()}
    for tg, tr in zip(gwas.target, gwas.trait):
        if live[tr]:
            g_terms.setdefault(tg, set()).add(live[tr])
            g_up.setdefault(tg, set()).update(up[tr])
    a_, b_, c_, d_ = [], [], [], []
    for t, d in zip(pairs.target, pairs.disease):
        a_.append(n_sources.get((t, d), 0) >= 2)
        b_.append(bool(rel.get(d, set()) & (diseases_of.get(t, set()) - {d})))
        c = onto["canon"].get(d)
        mine = g_terms.get(t, set())
        c_.append(bool(c) and (c in mine or c in g_up.get(t, set()) or bool(onto["anc"].get(d, set()) & mine)))
        d_.append(bool(nb.get(t, set()) & (holders_of.get(d, set()) - {t})))
    out = pairs[["target", "disease", "source"]].copy()
    out["second_source"], out["lineage"], out["gwas"], out["partner"] = a_, b_, c_, d_
    out["corroborated"] = out[["second_source", "lineage", "gwas", "partner"]].any(axis=1)
    return out


def neighbour_scores(pairs: pd.DataFrame, table: pd.DataFrame) -> np.ndarray:
    """For each pair, the neighbours score on the table with the pair itself left out: the sum, over the
    target's other diseases, of their cosine similarity to this disease by the targets they share."""
    from scipy import sparse
    tp = table.drop_duplicates(["target", "disease"])
    ti = {t: i for i, t in enumerate(sorted(set(tp.target)))}
    di = {d: i for i, d in enumerate(sorted(set(tp.disease)))}
    T = sparse.coo_matrix((np.ones(len(tp), np.float32), (tp.target.map(ti).to_numpy(), tp.disease.map(di).to_numpy())), shape=(len(ti), len(di))).tocsr()
    Tc = T.tocsc()
    deg = np.asarray(T.sum(axis=0)).ravel()
    out = np.zeros(len(pairs))
    pt, pd_ = pairs.target.map(ti).to_numpy(), pairs.disease.map(di).to_numpy()
    for d in np.unique(pd_):
        if deg[d] <= 1:
            continue
        shared = np.asarray(T.T @ Tc[:, d].toarray().ravel())       # targets each disease shares with d
        for k in np.flatnonzero(pd_ == d):
            mine = T.indices[T.indptr[pt[k]]:T.indptr[pt[k] + 1]]
            mine = mine[mine != d]
            if len(mine):
                out[k] = float(((shared[mine] - 1) / np.sqrt(deg[mine] * (deg[d] - 1))).sum())
    return out


def auc(real: np.ndarray, planted: np.ndarray) -> float:
    """The chance that a real pair scores above a planted one, ties counting a half."""
    from scipy.stats import rankdata
    if not len(real) or not len(planted):
        return float("nan")
    r = rankdata(np.r_[real, planted])
    return float((r[:len(real)].sum() - len(real) * (len(real) + 1) / 2) / (len(real) * len(planted)))


def evaluate(real: pd.DataFrame, planted: pd.DataFrame, gwas: pd.DataFrame, partners: pd.DataFrame, terms: dict) -> tuple[dict, pd.DataFrame]:
    table = pd.concat([real, planted], ignore_index=True)
    rel, onto = lineage_of(sorted(set(table.disease)), terms)
    n_src = real.groupby(["target", "disease"]).source.transform("nunique")
    single = real[n_src == 1].drop_duplicates(["target", "disease"]).assign(planted=False)
    both = pd.concat([single, planted.assign(planted=True)], ignore_index=True)
    c = corroboration(both, table, rel, onto, gwas, partners, terms)
    c["planted"] = both.planted.to_numpy()
    c["neighbours"] = neighbour_scores(both, table)
    dt = table.groupby(["source", "target"]).disease.transform("nunique")
    dd = table.groupby(["source", "disease"]).target.transform("nunique")
    key = dict(zip(zip(table.source, table.target, table.disease), ((dt - 1) * (dd - 1)).to_numpy()))
    c["degree_product"] = [key[(s, t, d)] for s, t, d in zip(c.source, c.target, c.disease)]
    res = {}
    for s, g in list(c.groupby("source")) + [("all sources", c)]:
        p, r = g[g.planted], g[~g.planted]
        row = {"pairs": int((real.source == s).sum()) if s != "all sources" else int(len(real.drop_duplicates(["target", "disease"]))),
               "real_single_source": len(r), "planted": len(p),
               "planted_withheld": float((~p.corroborated).mean()), "real_single_source_withheld": float((~r.corroborated).mean()),
               "corroborated_by": {k: {"planted": float(p[k].mean()), "real": float(r[k].mean())} for k in ("lineage", "gwas", "partner")},
               "auc_neighbours": auc(r.neighbours.to_numpy(), p.neighbours.to_numpy()),
               "auc_degree_product": auc(r.degree_product.to_numpy(), p.degree_product.to_numpy())}
        row["rule_applied"] = bool(row["planted_withheld"] >= KEEP_AT_LEAST and row["real_single_source_withheld"] <= LOSE_AT_MOST)
        res[s] = row
    return res, c


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    terms = ontology.load("""[Term]
id: MONDO:1
name: disease

[Term]
id: MONDO:2
name: heart disease
is_a: MONDO:1

[Term]
id: MONDO:3
name: arrhythmia
is_a: MONDO:2

[Term]
id: MONDO:4
name: lung disease
is_a: MONDO:1
""".splitlines())
    real = pd.DataFrame([("A", "MONDO:3", "s1"), ("A", "MONDO:3", "s2"), ("B", "MONDO:3", "s1"), ("B", "MONDO:2", "s1"), ("C", "MONDO:4", "s1"),
                         ("D", "MONDO:4", "s1"), ("E", "UMLS:9", "s1"), ("F", "MONDO:3", "s1"), ("G", "MONDO:4", "s2")], columns=["target", "disease", "source"])
    rel, onto = lineage_of(sorted(set(real.disease)), terms)
    ck("lineage runs both ways between a disease and its ancestor in the table, and siblings are apart",
       rel["MONDO:3"] == {"MONDO:2"} and rel["MONDO:2"] == {"MONDO:3"} and rel["MONDO:4"] == set() and rel["UMLS:9"] == set())
    gwas = pd.DataFrame({"target": ["C", "F"], "trait": ["MONDO:4", "MONDO:2"]})
    partners = pd.DataFrame({"a": ["D"], "b": ["C"]})
    q = pd.DataFrame([("A", "MONDO:3", "s1"), ("B", "MONDO:3", "s1"), ("C", "MONDO:4", "s1"), ("D", "MONDO:4", "s1"), ("E", "UMLS:9", "s1"),
                      ("F", "MONDO:3", "s1"), ("G", "MONDO:4", "s2")], columns=["target", "disease", "source"])
    c = corroboration(q, real, rel, onto, gwas, partners, terms).set_index("target")
    ck("(a) a second source corroborates", bool(c.loc["A", "second_source"]) and not c.loc["B", "second_source"])
    ck("(b) the same target with a broader disease corroborates, and a pair does not corroborate itself", bool(c.loc["B", "lineage"]) and not c.loc["A", "lineage"])
    ck("(c) a GWAS link to the disease, or to a broader form of it, corroborates", bool(c.loc["C", "gwas"]) and bool(c.loc["F", "gwas"]) and not c.loc["D", "gwas"])
    ck("(d) a partner asserted with the same disease corroborates, in either direction", bool(c.loc["D", "partner"]) and bool(c.loc["C", "partner"]) and not c.loc["G", "partner"])
    ck("a pair with none of the four is not corroborated, and a disease outside the ontology has no lineage", not c.loc["E", "corroborated"] and not c.loc["G", "corroborated"])
    big = pd.DataFrame([(f"T{i}", f"D{i % 7}", "s") for i in range(140)] + [(f"T{i}", f"D{(i + 1) % 7}", "s") for i in range(0, 140, 3)], columns=["target", "disease", "source"])
    pl = plant(big, fraction=0.05, seed=1)
    have = set(zip(big.target, big.disease))
    ck("planted pairs number 5% of the source's pairs, rounded up, and none is a real pair or a repeat",
       len(pl) == math.ceil(0.05 * len(big.drop_duplicates())) and not (set(zip(pl.target, pl.disease)) & have) and not pl.duplicated(["target", "disease"]).any())
    ck("a planted pair's target and disease both come from the source", set(pl.target) <= set(big.target) and set(pl.disease) <= set(big.disease))
    ck("planting is the same under the same seed", plant(big, 0.05, 1).equals(pl) and not plant(big, 0.05, 2).equals(pl))
    ck("the area is 1 when every real pair scores above every planted one, a half when they tie", auc(np.array([2.0, 3.0]), np.array([1.0])) == 1.0 and auc(np.array([1.0, 1.0]), np.array([1.0])) == 0.5)
    tb = pd.DataFrame([("A", "x", "s"), ("A", "y", "s"), ("B", "x", "s"), ("B", "y", "s"), ("C", "x", "s"), ("C", "z", "s")], columns=["target", "disease", "source"])
    ns = neighbour_scores(pd.DataFrame([("A", "x", "s"), ("C", "z", "s")], columns=["target", "disease", "source"]), tb)
    # (A,x) left out: A's other disease y shares B with x (shared 2 incl. A, minus A = 1); deg y = 2, deg x - 1 = 2 -> 1/sqrt(4) = 0.5
    ck("the neighbours score leaves the pair itself out", abs(ns[0] - 0.5) < 1e-6 and ns[1] == 0.0)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 1:
        raise SystemExit(__doc__)
    release = argv[0]
    C = paths.context(release)
    real = pd.read_parquet(C / "disease_links.parquet")
    partners = pd.read_parquet(C / "partners.parquet")
    gwas = pd.read_parquet(paths.bundle(release) / "links.parquet", columns=["target", "trait"])
    terms = ontology.read(registry.require("efo_obo"))
    planted = plant(real)
    res, detail = evaluate(real, planted, gwas, partners, terms)
    print(f"  {'source':34s} {'pairs':>8s} {'single':>8s} {'planted':>8s} {'planted withheld':>17s} {'real withheld':>14s} {'AUC neigh.':>11s} {'AUC degree':>11s}  applied")
    for s, r in sorted(res.items(), key=lambda kv: (kv[0] == "all sources", -kv[1]["pairs"])):
        print(f"  {s:34s} {r['pairs']:>8,} {r['real_single_source']:>8,} {r['planted']:>8,} {r['planted_withheld']:>17.3f} {r['real_single_source_withheld']:>14.3f} "
              f"{r['auc_neighbours']:>11.3f} {r['auc_degree_product']:>11.3f}  {r['rule_applied']}")
    a = res["all sources"]["corroborated_by"]
    print("  corroborated by each criterion alone, all sources (planted / real single-source): "
          + "; ".join(f"{k} {v['planted']:.3f} / {v['real']:.3f}" for k, v in a.items()))
    apply_to = sorted(s for s, r in res.items() if s != "all sources" and r["rule_applied"])
    single = detail[~detail.planted]
    withheld = single[single.source.isin(apply_to) & ~single.corroborated][["target", "disease", "source"]]
    kept = real.merge(withheld.assign(_w=1), on=["target", "disease", "source"], how="left")
    kept = kept[kept._w.isna()].drop(columns="_w")
    kept.to_parquet(C / "disease_links_pruned.parquet", index=False)
    withheld.assign(reason="one source, and nothing else in the graph corroborates it",
                    restored_by="a second source, the same target with a broader or narrower disease, a GWAS link, or an interaction partner with this disease"
                    ).to_parquet(C / "disease_links_withheld.parquet", index=False)
    man = json.loads((C / "manifest.json").read_text())
    out = paths.DOCS_DATA / f"planted_{release}.json"
    out.write_text(json.dumps({"release": release, "seed": SEED, "fraction": FRACTION, "bar": {"planted_withheld_at_least": KEEP_AT_LEAST, "real_withheld_at_most": LOSE_AT_MOST},
                               "disease_links_sha256": man["tables"]["disease_links"]["sha256"], "partners_sha256": man["tables"]["partners"]["sha256"],
                               "planted_pairs": len(planted), "sources": res, "rule_applied_to": apply_to,
                               "real_pairs_withheld": len(withheld), "rows_kept": len(kept), "rows_before": len(real),
                               "files": {f: graph.file_sha(C / f) for f in ("disease_links_pruned.parquet", "disease_links_withheld.parquet")}}, indent=1) + "\n")
    print(f"  rule applied to {len(apply_to)} source(s): {apply_to}; real pairs withheld {len(withheld):,}; rows kept {len(kept):,} of {len(real):,}")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
