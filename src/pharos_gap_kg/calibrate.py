"""How often a candidate ranked this high in the 2018 backtest was reported by 2022.

The front end ranks candidates by the neighbours score. A rank means nothing to a reader without the
rate at which candidates at that rank turned out to be links, so this computes that rate once, on the
frozen split, and the front end quotes it. Defined in docs/PROTOCOL.md before it was run.

Run:  python -m pharos_gap_kg.calibrate pharos319 --cutoff 2018
      python -m pharos_gap_kg.calibrate --selftest
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from pharos_gap_kg import baselines, context, controls, ontology, paths, position, registry, score, split

BINS = {"by_trait": [(1, 10), (11, 50), (51, 100), (101, 500)], "by_target": [(1, 5), (6, 20), (21, 50)]}


def tally(S: np.ndarray, train: np.ndarray, test: np.ndarray, groups: dict[str, np.ndarray], bins: list[tuple[int, int]],
          seeds: tuple[int, ...] = score.SEEDS) -> dict:
    """For each rank bin and each named group of candidates: candidates and positives, pooled over
    the scored columns and averaged over seeds. A group is a boolean matrix; "all" is every candidate."""
    units = np.flatnonzero(test.any(axis=0))
    n = np.zeros((len(seeds), len(bins), len(groups)))
    k = np.zeros_like(n)
    names = list(groups)
    for si, seed in enumerate(seeds):
        rng = np.random.default_rng(seed)
        for j in units:
            cand = np.flatnonzero(~train[:, j])
            order = cand[np.lexsort((rng.random(len(cand)), -S[cand, j]))]
            y = test[order, j]
            for bi, (lo, hi) in enumerate(bins):
                rows = order[lo - 1:hi]
                yy = y[lo - 1:hi]
                for gi, g in enumerate(names):
                    m = groups[g][rows, j]
                    n[si, bi, gi] += m.sum()
                    k[si, bi, gi] += (yy & m).sum()
    n, k = n.mean(axis=0), k.mean(axis=0)
    return {f"{lo}-{hi}": {g: {"candidates": float(n[bi, gi]), "positives": float(k[bi, gi]),
                                "share": (float(k[bi, gi] / n[bi, gi]) if n[bi, gi] else None)}
                            for gi, g in enumerate(names)} for bi, (lo, hi) in enumerate(bins)}


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    train = np.zeros((6, 2), bool)
    train[5, 0] = True
    test = np.zeros((6, 2), bool)
    test[[0, 2], 0] = True
    S = np.array([[9, 0], [8, 0], [7, 0], [6, 0], [5, 0], [4, 0]], float)
    grp = {"all": np.ones((6, 2), bool), "even": (np.arange(6) % 2 == 0)[:, None].repeat(2, axis=1)}
    r = tally(S, train, test, grp, [(1, 2), (3, 5)], seeds=(0,))
    ck("a bin counts its candidates and its positives", r["1-2"]["all"] == {"candidates": 2.0, "positives": 1.0, "share": 0.5})
    ck("a training link is never a candidate, so the second bin holds three", r["3-5"]["all"]["candidates"] == 3.0 and r["3-5"]["all"]["positives"] == 1.0)
    ck("a group counts only its own candidates", r["1-2"]["even"] == {"candidates": 1.0, "positives": 1.0, "share": 1.0})
    ck("only columns with a positive are pooled", sum(v["all"]["candidates"] for v in r.values()) == 5.0)
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
    terms = ontology.read(registry.require("efo_obo"))
    R = controls.lineage_matrix(M["traits"], terms)
    tg = position.target_genes(a.release, position.genes(registry.require("ensembl_gtf_113")))
    pairs = position.neighbour_pairs(tg, {t: i for i, t in enumerate(M["targets"])})
    S = baselines.scorers(A)
    linked = (((A.astype(np.float32) @ R.astype(np.float32).T) > 0) | (position.neighbour_count(pairs, A) > 0)) & ~A
    _, dis = context.kinds(M["traits"], terms)
    disease = np.repeat(dis[None, :], A.shape[0], axis=0)
    groups = {"all": np.ones(A.shape, bool), "linked": linked, "not_linked": ~linked,
              "disease_like": disease, "measurement_or_other": ~disease,
              "disease_like_not_linked": disease & ~linked}
    res = {"by_trait": tally(S["B3"], A, Y, groups, BINS["by_trait"]),
           "by_target": tally(S["B3"].T, A.T, Y.T, {k: v.T for k, v in groups.items()}, BINS["by_target"]),
           "degree_by_trait": tally(S["B1"], A, Y, {"all": groups["all"]}, BINS["by_trait"])}
    base = float(Y.sum() / (~A).sum())
    for axis in ("by_trait", "by_target"):
        print(f"{axis}: share of candidates at each rank that were reported by 2022 (base rate {base:.4f})")
        print(f"  {'ranks':>8s} " + " ".join(f"{g:>24s}" for g in groups))
        for b, v in res[axis].items():
            print(f"  {b:>8s} " + " ".join(f"{(v[g]['share'] if v[g]['share'] is not None else float('nan')):>15.3f} of {v[g]['candidates']:>6.0f}" for g in groups))
    print("target degree, by trait: " + ", ".join(f"{b} {v['all']['share']:.3f}" for b, v in res["degree_by_trait"].items()))
    frozen = json.loads(split.frozen_path(a.cutoff).read_text())
    out = paths.DOCS_DATA / f"calibration_gene_trait_{a.cutoff}.json"
    out.write_text(json.dumps({"cutoff": a.cutoff, "release": a.release, "scorer": "B3 neighbours", "base_rate": base,
                               "data_ends": int(M["year"].max()), "split": {k: frozen[k]["sha256"] for k in split.PARTS},
                               "seeds": list(score.SEEDS), **res}, indent=1) + "\n")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
