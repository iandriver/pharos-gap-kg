"""Does the rest of Pharos add to the neighbours score? Five context scorers on the frozen split.

Each scorer reads one of the context tables and the training links, and has no fitted parameter
(docs/PROTOCOL.md). The context tables are not dated, so what they show is an upper bound on what
context known at the cutoff could offer.

    X1 partners          how many of the target's interaction partners are linked to the trait
    X2 pathways          over the target's pathways, the share of each pathway's other members linked to it
    X3 GO terms          the same over GO terms
    X4 shared compounds  the same over compounds active on two or more targets
    X5 disease sources   how many of the eleven sources assert the target with the trait's live term

Each is also combined with the neighbours score B3 by averaging percentile ranks within a trait.

Run:  python -m pharos_gap_kg.enrich pharos319 --cutoff 2018
      python -m pharos_gap_kg.enrich --selftest
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd

from pharos_gap_kg import baselines, context, controls, ontology, paths, position, registry, score, split

WEIGHTS = (0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 1.0)
NAMES = {"X1": "partners", "X2": "pathways", "X3": "GO terms", "X4": "shared compounds", "X5": "disease sources"}


def partner_count(pairs: pd.DataFrame, index: dict, A: np.ndarray) -> np.ndarray:
    """X1: for each (target, trait), the number of the target's partners with a training link to the trait."""
    from scipy import sparse
    p = pairs[pairs.a.isin(index) & pairs.b.isin(index)]
    i, j = p.a.map(index).to_numpy(), p.b.map(index).to_numpy()
    n = A.shape[0]
    P = sparse.coo_matrix((np.ones(2 * len(p), np.float32), (np.r_[i, j], np.r_[j, i])), shape=(n, n)).tocsr()
    P.data[:] = 1.0
    return np.asarray(P @ A.astype(np.float32))


def group_share(members: pd.DataFrame, index: dict, A: np.ndarray, chunk: int = 4000) -> np.ndarray:
    """X2 to X4: over the groups containing the target, the sum of the share of each group's OTHER
    members linked to the trait. A group's size counts every member, in the training set or not."""
    from scipy import sparse
    size = members.groupby("grp").target.nunique()
    m = members[members.target.isin(index)]
    groups = {g: k for k, g in enumerate(sorted(set(m.grp)))}
    G = sparse.coo_matrix((np.ones(len(m), np.float32), (m.target.map(index).to_numpy(), m.grp.map(groups).to_numpy())),
                          shape=(A.shape[0], len(groups))).tocsc()
    G.data[:] = 1.0
    other = np.array([size[g] - 1 for g in groups], np.float32)
    a = A.astype(np.float32)
    out = np.zeros(A.shape, np.float32)
    for k in range(0, len(groups), chunk):
        Gk = G[:, k:k + chunk]
        linked = np.asarray(Gk.T @ a)                          # groups x traits: members linked to the trait
        mine = np.asarray(Gk @ (linked / other[k:k + chunk, None]))
        own = np.asarray(Gk @ (1.0 / other[k:k + chunk])).ravel()   # what the target's own link would add
        out += mine - a * own[:, None]                          # a target does not vote for itself
    return out


def source_count(links: pd.DataFrame, targets: list[str], canon: list) -> np.ndarray:
    """X5: for each (target, trait), the number of disease sources asserting the target with the trait's live term."""
    C = context.source_matrices(links.rename(columns={"disease": "mondoid"}), targets, canon)
    out = np.zeros((len(targets), len(canon)), np.float32)
    for M in C.values():
        out += M
    return out


def percentile(S: np.ndarray, cand: np.ndarray) -> np.ndarray:
    """Each score as its percentile rank among the trait's candidates (ties share the mean rank)."""
    from scipy.stats import rankdata
    out = np.zeros(S.shape, np.float32)
    for j in range(S.shape[1]):
        c = cand[:, j]
        if c.any():
            out[c, j] = rankdata(S[c, j], method="average") / c.sum()
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    idx = {"A": 0, "B": 1, "C": 2, "D": 3}
    A = np.array([[1, 0], [0, 0], [1, 1], [0, 0]], bool)
    P = pd.DataFrame({"a": ["A", "B", "A", "Z"], "b": ["B", "C", "C", "A"]})
    X1 = partner_count(P, idx, A)
    ck("a partner's training link is counted, a pair listed either way round is one, and an unknown partner is dropped",
       X1.tolist() == [[1, 1], [2, 1], [1, 0], [0, 0]])
    M = pd.DataFrame({"target": ["A", "B", "C", "B", "D", "Q"], "grp": ["g1", "g1", "g1", "g2", "g2", "g2"]})
    X2 = group_share(M, idx, A)
    # g1 = {A,B,C}, size 3; g2 = {B,D,Q}, size 3 (Q is outside the index but counts in the size)
    # trait 0: linked in g1 = A,C. B gets 2/2 from g1 and 0 from g2. A gets (2-1)/2. C gets (2-1)/2. D gets 0.
    ck("a group votes with the share of its other members linked, and a target does not vote for itself",
       np.allclose(X2[:, 0], [0.5, 1.0, 0.5, 0.0]) and np.allclose(X2[:, 1], [0.5, 0.5, 0.0, 0.0]))
    ck("a group's size counts members outside the training set", abs(group_share(pd.DataFrame({"target": ["A", "B", "Q", "R"], "grp": ["g"] * 4}), idx, A)[1, 0] - 1 / 3) < 1e-6)
    L = pd.DataFrame({"source": ["s1", "s2", "s1"], "target": ["A", "A", "B"], "disease": ["MONDO:1", "MONDO:1", "MONDO:9"]})
    X5 = source_count(L, ["A", "B", "C", "D"], ["MONDO:1", "EFO:2"])
    ck("disease sources are counted on the trait that stands for the term, and nowhere else", X5.tolist() == [[2, 0], [0, 0], [0, 0], [0, 0]])
    S = np.array([[3.0], [1.0], [1.0], [9.0]])
    pc = percentile(S, np.array([[True], [True], [True], [False]]))
    ck("a percentile is among candidates only, and ties share a rank", np.allclose(pc[:, 0], [1.0, 0.5, 0.5, 0.0]))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("release")
    ap.add_argument("--cutoff", type=int, required=True)
    a = ap.parse_args(argv)
    M = baselines.matrices(split.load(a.cutoff))
    A, Y = M["train"], M["test"]
    idx = {t: i for i, t in enumerate(M["targets"])}
    C = paths.context(a.release)
    if not (C / "manifest.json").is_file():
        raise SystemExit(f"no context tables; run:  python -m pharos_gap_kg.ctxtables {a.release}")
    rd = lambda n: pd.read_parquet(C / f"{n}.parquet")
    terms = ontology.read(registry.require("efo_obo"))
    canon, dis = context.kinds(M["traits"], terms)
    R = controls.lineage_matrix(M["traits"], terms)
    tg = position.target_genes(a.release, position.genes(registry.require("ensembl_gtf_113")))
    pairs = position.neighbour_pairs(tg, idx)
    base = baselines.scorers(A)
    linked = (((A.astype(np.float32) @ R.astype(np.float32).T) > 0) | (position.neighbour_count(pairs, A) > 0)) & ~A
    S = {"B1": base["B1"], "B3": base["B3"],
         "X1": partner_count(rd("partners"), idx, A), "X2": group_share(rd("pathways"), idx, A),
         "X3": group_share(rd("go_terms"), idx, A), "X4": group_share(rd("compounds"), idx, A),
         "X5": source_count(rd("disease_links"), M["targets"], canon)}
    cand = ~A
    pct = {k: percentile(S[k], cand) for k in ("B3", "X1", "X2", "X3", "X4", "X5")}
    for k in NAMES:
        S[f"B3+{k}"] = (pct["B3"] + pct[k]) / 2
    S["B3+all"] = sum(pct.values()) / len(pct)
    scored = Y.any(axis=0)
    dcols = np.flatnonzero(dis & scored)
    new = Y & ~linked
    res, per = {}, {}
    cover = {k: {"candidates_scored_above_zero": int(((S[k] > 0) & cand).sum()), "positives_scored_above_zero": int(((S[k] > 0) & Y).sum())} for k in NAMES}
    for k in ["B1", "B3"] + list(NAMES) + [f"B3+{x}" for x in NAMES] + ["B3+all"]:
        t0 = time.time()
        full = score.by_axis(S[k], A, Y, axis=1)
        dz = score.by_axis(S[k][:, dcols], A[:, dcols], Y[:, dcols], axis=1)
        neu = score.by_axis(S[k], A, new, axis=1, exclude=linked)
        loc = controls.collapsed(S[k], A, Y, pairs)
        per[k] = {"primary": full["ap"], "disease_like": dz["ap"], "new_links": neu["ap"], "locus": loc["ap"]}
        res[k] = {v: score.interval(per[k][v]) for v in per[k]}
        r = res[k]
        print(f"  {k:8s} {NAMES.get(k, ''):16s} primary {r['primary'][0]:.4f}  disease-like {r['disease_like'][0]:.4f}  "
              f"new links {r['new_links'][0]:.4f}  by locus {r['locus'][0]:.4f}  ({time.time() - t0:.0f} s)", flush=True)
    paired = {}
    for k in [f"B3+{x}" for x in NAMES] + ["B3+all"]:
        paired[f"{k}-B3"] = {}
        for v in ("primary", "disease_like", "new_links", "locus"):
            m, lo, hi = score.interval(per[k][v] - per["B3"][v])
            paired[f"{k}-B3"][v] = {"difference": m, "low": lo, "high": hi}
        d = paired[f"{k}-B3"]
        print(f"  {k} minus B3: primary {d['primary']['difference']:+.4f} [{d['primary']['low']:+.4f}, {d['primary']['high']:+.4f}]  "
              f"by locus {d['locus']['difference']:+.4f} [{d['locus']['low']:+.4f}, {d['locus']['high']:+.4f}]  "
              f"new links {d['new_links']['difference']:+.4f} [{d['new_links']['low']:+.4f}, {d['new_links']['high']:+.4f}]")
    al = paired["B3+all-B3"]
    verdict = bool(al["primary"]["low"] > 0 and al["locus"]["low"] > 0)
    print(f"  context adds to the neighbours score (both intervals above zero): {verdict}")
    # An upper bound, not a result: the weight on B3 is chosen on the test links themselves. If no
    # weight below 1 beats B3 even so, no rank-averaging of this scorer with B3 can.
    bound = {}
    for k in list(NAMES):
        best = max(((w, float(np.nanmean(score.by_axis(w * pct["B3"] + (1 - w) * pct[k], A, Y, axis=1, seeds=(0,))["ap"]))) for w in WEIGHTS),
                   key=lambda x: x[1])
        bound[k] = {"best_weight_on_B3": best[0], "primary_at_best_weight": best[1]}
    alone = float(np.nanmean(score.by_axis(pct["B3"], A, Y, axis=1, seeds=(0,))["ap"]))
    print(f"  upper bound, weight on B3 tuned on the test links (B3 alone {alone:.4f}): "
          + "; ".join(f"{k} best at {v['best_weight_on_B3']:.2f} -> {v['primary_at_best_weight']:.4f}" for k, v in bound.items()))
    frozen = json.loads(split.frozen_path(a.cutoff).read_text())
    man = json.loads((C / "manifest.json").read_text())
    out = paths.DOCS_DATA / f"enrich_gene_trait_{a.cutoff}.json"
    out.write_text(json.dumps({"cutoff": a.cutoff, "release": a.release, "split": {k: frozen[k]["sha256"] for k in split.PARTS},
                               "context_tables": {k: v["sha256"] for k, v in man["tables"].items()}, "names": NAMES, "coverage": cover,
                               "traits_scored": int(scored.sum()), "disease_like_traits": int(len(dcols)),
                               "results": res, "paired": paired, "context_adds_to_neighbours": verdict,
                               "upper_bound_weight_tuned_on_test": {"B3_alone_one_seed": alone, "weights_tried": list(WEIGHTS), "by_scorer": bound}},
                              indent=1) + "\n")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
