"""The four baselines of docs/PROTOCOL.md on the frozen split. None has a fitted parameter.

    B0 chance          every candidate ties
    B1 target degree   how many training traits the target has
    B2 degree product  target degree times trait degree
    B3 neighbours      the target's training traits, weighted by how similar each is to the candidate
                       trait (cosine over the training targets two traits share)

Run:  python -m pharos_gaps.baselines --cutoff 2018
      python -m pharos_gaps.baselines --selftest
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np

from pharos_gaps import paths, score, split

NAMES = {"B0": "chance", "B1": "target degree", "B2": "degree product", "B3": "neighbours"}


def matrices(parts: dict) -> dict:
    tr, te = parts["train"], parts["test"]
    targets = sorted(tr.target.unique())
    traits = sorted(tr.trait.unique())
    ti = {t: i for i, t in enumerate(targets)}
    di = {d: i for i, d in enumerate(traits)}
    A = np.zeros((len(targets), len(traits)), bool)
    A[tr.target.map(ti).to_numpy(), tr.trait.map(di).to_numpy()] = True
    Y = np.zeros_like(A)
    year = np.zeros(A.shape, np.int16)
    r, c = te.target.map(ti).to_numpy(), te.trait.map(di).to_numpy()
    Y[r, c] = True
    year[r, c] = te.year.to_numpy()
    if (A & Y).any() or int(A.sum()) != len(tr) or int(Y.sum()) != len(te):
        raise SystemExit("REFUSED: the training and test links overlap, or a link was lost building the matrices")
    return {"train": A, "test": Y, "year": year, "targets": targets, "traits": traits}


def scorers(A: np.ndarray) -> dict[str, np.ndarray]:
    a = A.astype(np.float32)
    deg_t, deg_d = a.sum(axis=1), a.sum(axis=0)
    norm = np.sqrt(deg_d)
    sim = (a.T @ a) / np.outer(norm, norm)
    np.fill_diagonal(sim, 0.0)
    return {"B0": np.zeros_like(a), "B1": np.repeat(deg_t[:, None], a.shape[1], axis=1),
            "B2": np.outer(deg_t, deg_d), "B3": a @ sim}


def evaluate(S: np.ndarray, M: dict, years: list[int]) -> dict:
    A, Y = M["train"], M["test"]
    by_trait = score.by_axis(S, A, Y, axis=1)
    by_target = score.by_axis(S, A, Y, axis=0)
    out = {"per_trait_ap": by_trait["ap"], "units": by_trait["index"],
           "primary": score.interval(by_trait["ap"]), "primary_by_seed": by_trait["by_seed"].tolist(),
           "trait_auc": score.interval(by_trait["auc"]), "trait_recall_at_100": score.interval(by_trait["recall"]),
           "target_ap": score.interval(by_target["ap"]), "target_auc": score.interval(by_target["auc"]),
           "overall_ap_by_seed": score.overall(S, A, Y), "by_year": {}}
    for y in years:
        this = Y & (M["year"] == y)
        r = score.by_axis(S, A, this, axis=1, exclude=Y & ~this)
        out["by_year"][y] = {"traits": int(len(r["index"])), "links": int(this.sum()), "primary": score.interval(r["ap"])}
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    import pandas as pd
    tr = pd.DataFrame({"target": ["A", "A", "B", "C", "C"], "trait": ["x", "y", "x", "y", "z"], "year": 2010})
    te = pd.DataFrame({"target": ["B", "A"], "trait": ["y", "z"], "year": [2019, 2020]})
    M = matrices({"train": tr, "test": te})
    ck("the matrices hold each link once, in sorted order", M["train"].sum() == 5 and M["test"].sum() == 2 and M["targets"] == ["A", "B", "C"] and M["traits"] == ["x", "y", "z"])
    ck("a test link keeps its year", M["year"][1, 1] == 2019 and M["year"][0, 2] == 2020 and M["year"].sum() == 4039)
    S = scorers(M["train"])
    ck("chance ties every candidate", not S["B0"].any())
    ck("target degree is the count of a target's training traits", S["B1"][:, 0].tolist() == [2, 1, 2])
    ck("the degree product ranks as target degree within a trait", (np.argsort(-S["B2"][:, 1], kind="stable") == np.argsort(-S["B1"][:, 1], kind="stable")).all())
    # traits x and y share target A; y and z share target C. B has only x, so its score for y is sim(x,y) = 1/sqrt(2*2) = 0.5
    ck("a neighbour score is the summed similarity of the target's own traits to the candidate", abs(S["B3"][1, 1] - 0.5) < 1e-6 and S["B3"][1, 2] == 0)
    ck("a trait does not vote for itself", abs(S["B3"][0, 0] - 0.5) < 1e-6)   # A has x,y: score for x is sim(y,x) only
    try:
        matrices({"train": tr, "test": tr.iloc[:1]})
        ck("a test link that is also a training link is refused", False)
    except SystemExit:
        ck("a test link that is also a training link is refused", True)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cutoff", type=int, required=True)
    a = ap.parse_args(argv)
    parts = split.load(a.cutoff)
    M = matrices(parts)
    years = sorted(int(y) for y in parts["test"].year.unique() if (parts["test"].year == y).sum() >= 100)
    S = scorers(M["train"])
    res, per_trait = {}, {}
    for k in NAMES:
        t0 = time.time()
        r = evaluate(S[k], M, years)
        per_trait[k] = r.pop("per_trait_ap")
        r.pop("units")
        res[k] = r
        print(f"  {k} {NAMES[k]:15s} primary {r['primary'][0]:.4f} [{r['primary'][1]:.4f}, {r['primary'][2]:.4f}]  "
              f"trait ROC {r['trait_auc'][0]:.3f}  recall@100 {r['trait_recall_at_100'][0]:.3f}  "
              f"by target {r['target_ap'][0]:.4f}  overall {np.mean(r['overall_ap_by_seed']):.4f}  ({time.time() - t0:.0f} s)", flush=True)
    pairs = {}
    for x, y in (("B1", "B0"), ("B3", "B1"), ("B3", "B2"), ("B2", "B1")):
        m, lo, hi = score.interval(per_trait[x] - per_trait[y])
        pairs[f"{x}-{y}"] = {"difference": m, "low": lo, "high": hi, "half_width": (hi - lo) / 2}
        print(f"  {x} minus {y}: {m:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    A = M["train"]
    out = paths.DOCS_DATA / f"baselines_gene_trait_{a.cutoff}.json"
    frozen = json.loads(split.frozen_path(a.cutoff).read_text())
    out.write_text(json.dumps({
        "cutoff": a.cutoff, "split": {k: frozen[k]["sha256"] for k in split.PARTS},
        "targets": len(M["targets"]), "traits": len(M["traits"]), "traits_scored": int(M["test"].any(axis=0).sum()),
        "targets_scored": int(M["test"].any(axis=1).sum()), "candidates": int((~A).sum()), "positives": int(M["test"].sum()),
        "seeds": list(score.SEEDS), "bootstrap": score.N_BOOT, "names": NAMES, "results": res, "paired": pairs}, indent=1) + "\n")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
