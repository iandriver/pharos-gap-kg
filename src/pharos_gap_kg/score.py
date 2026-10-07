"""The measures of docs/PROTOCOL.md, computed over every candidate, with nothing sampled.

A scorer is a matrix S[target, trait]. `train` marks the links it was allowed to see, `test` the links
first reported later. Candidates are the pairs with no training link. Ties are broken at random under
a seed, so a scorer that says nothing (every candidate equal) scores as chance and not as whatever
order the rows happened to be in.

Run:  python -m pharos_gap_kg.score --selftest
"""
from __future__ import annotations

import sys

import numpy as np

SEEDS = (0, 1, 2, 3, 4)
TOP_K = 100
N_BOOT = 2000


def ranked(scores: np.ndarray, labels: np.ndarray, rng: np.random.Generator, k: int = TOP_K) -> tuple[float, float, float]:
    """Average precision, ROC area and recall in the top k, for one ranking with random tie-breaks."""
    n, n_pos = len(labels), int(labels.sum())
    if n_pos == 0 or n_pos == n:
        return float("nan"), float("nan"), float("nan")
    order = np.lexsort((rng.random(n), -scores))
    y = labels[order]
    hits = np.cumsum(y)
    at = np.flatnonzero(y)
    ap = float((hits[at] / (at + 1)).mean())
    negatives_below = (n - 1 - at) - (n_pos - hits[at])
    auc = float(negatives_below.sum() / (n_pos * (n - n_pos)))
    return ap, auc, float(y[:k].sum() / n_pos)


def by_axis(S: np.ndarray, train: np.ndarray, test: np.ndarray, axis: int, exclude: np.ndarray | None = None,
            seeds: tuple[int, ...] = SEEDS) -> dict[str, np.ndarray]:
    """Per trait (axis=1) or per target (axis=0): the three measures, meaned over seeds.

    Returns arrays over the units that have at least one positive, plus `index` (which units) and
    `by_seed` (the overall mean of average precision under each seed). `exclude` removes further
    pairs from the candidates."""
    if axis == 0:
        S, train, test = S.T, train.T, test.T
        exclude = None if exclude is None else exclude.T
    units = np.flatnonzero(test.any(axis=0))
    out = np.full((len(seeds), len(units), 3), np.nan)
    for si, seed in enumerate(seeds):
        rng = np.random.default_rng(seed)
        for ui, j in enumerate(units):
            cand = ~train[:, j] if exclude is None else ~train[:, j] & ~exclude[:, j]
            out[si, ui] = ranked(S[cand, j], test[cand, j], rng)
    m = np.nanmean(out, axis=0)
    return {"index": units, "ap": m[:, 0], "auc": m[:, 1], "recall": m[:, 2], "by_seed": np.nanmean(out[:, :, 0], axis=1)}


def overall(S: np.ndarray, train: np.ndarray, test: np.ndarray, seeds: tuple[int, ...] = SEEDS) -> list[float]:
    """Average precision over every candidate at once, one value per seed."""
    cand = ~train
    s, y = S[cand], test[cand]
    return [ranked(s, y, np.random.default_rng(seed))[0] for seed in seeds]


def interval(values: np.ndarray, seed: int = 0, n: int = N_BOOT) -> tuple[float, float, float]:
    """Mean and 95% percentile interval, resampling units."""
    v = values[~np.isnan(values)]
    rng = np.random.default_rng(seed)
    means = v[rng.integers(0, len(v), (n, len(v)))].mean(axis=1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return float(v.mean()), float(lo), float(hi)


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    rng = np.random.default_rng(0)
    ap, auc, rec = ranked(np.array([3.0, 2.0, 1.0]), np.array([True, False, True]), rng, k=1)
    ck("average precision of hits at ranks 1 and 3 is (1 + 2/3) / 2", abs(ap - (1 + 2 / 3) / 2) < 1e-12)
    ck("its ROC area is one of two pairs in order", abs(auc - 0.5) < 1e-12)
    ck("recall in the top 1 is half", rec == 0.5)
    ap, auc, _ = ranked(np.array([9.0, 8.0, 1.0, 0.0]), np.array([True, True, False, False]), rng)
    ck("a perfect ranking scores 1 on both", ap == 1.0 and auc == 1.0)
    ck("no positive, or nothing but positives, has no value", all(np.isnan(ranked(np.zeros(3), np.zeros(3, bool), rng))) and all(np.isnan(ranked(np.zeros(3), np.ones(3, bool), rng))))
    y = np.zeros(2000, bool)
    y[:200] = True
    aps = [ranked(np.zeros(2000), y, np.random.default_rng(s))[0] for s in range(20)]
    aucs = [ranked(np.zeros(2000), y, np.random.default_rng(s))[1] for s in range(20)]
    ck("a scorer that ties everything scores as chance: precision near prevalence, ROC area near a half",
       abs(np.mean(aps) - 0.1) < 0.01 and abs(np.mean(aucs) - 0.5) < 0.02)
    few = np.arange(50) < 5
    ck("the same seed gives the same tie-break, and another seed a different one",
       ranked(np.zeros(50), few, np.random.default_rng(3)) == ranked(np.zeros(50), few, np.random.default_rng(3))
       and ranked(np.zeros(50), few, np.random.default_rng(3)) != ranked(np.zeros(50), few, np.random.default_rng(4)))
    train = np.array([[1, 0, 0], [0, 0, 1], [0, 0, 0], [1, 1, 0]], bool)
    test = np.array([[0, 1, 0], [0, 0, 0], [0, 1, 0], [0, 0, 0]], bool)
    S = np.array([[0, 5, 0], [0, 1, 0], [0, 4, 0], [0, 9, 0]], float)
    r = by_axis(S, train, test, axis=1, seeds=(0,))
    ck("only traits with a positive are scored, over candidates without a training link", r["index"].tolist() == [1] and r["ap"][0] == 1.0)
    r0 = by_axis(S, train, test, axis=0, seeds=(0,))
    ck("by target, each target with a positive ranks its own candidate traits", r0["index"].tolist() == [0, 2] and np.allclose(r0["ap"], 1.0))
    ex = np.zeros_like(train)
    ex[0, 1] = True
    r2 = by_axis(S, train, test & ~ex, axis=1, exclude=ex, seeds=(0,))
    ck("an excluded pair leaves the candidates", r2["ap"][0] == 1.0 and r2["index"].tolist() == [1])
    ck("overall precision pools every candidate", abs(overall(S, train, test, seeds=(0,))[0] - 1.0) < 1e-12)
    m, lo, hi = interval(np.array([0.2, 0.4, 0.6, np.nan]))
    ck("an interval ignores units with no value and brackets the mean", abs(m - 0.4) < 1e-12 and lo <= m <= hi)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else 0)
