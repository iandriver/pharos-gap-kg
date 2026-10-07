"""The lineage, position and locus controls of docs/PROTOCOL.md, for every parameter-free scorer.

A later link can be predicted for reasons that are not biology. The trait may be an earlier trait
under a broader, narrower or second name (lineage). The gene may sit beside a gene already linked to
the trait (position). And one GWAS signal is assigned to several neighbouring genes, so positives
come in clumps (locus). This measures how much of each scorer's result those account for.

Run:  python -m pharos_gap_kg.controls pharos319 --cutoff 2018
      python -m pharos_gap_kg.controls --selftest
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np

from pharos_gap_kg import baselines, ontology, paths, position, registry, score, split

NAMES = {**baselines.NAMES, "B4": "lineage", "B5": "position"}


def lineage_matrix(traits: list[str], terms: dict) -> np.ndarray:
    lin = ontology.lineage(traits, terms)
    at = {t: i for i, t in enumerate(traits)}
    R = np.zeros((len(traits), len(traits)), bool)
    for a, others in lin.items():
        for b in others:
            R[at[a], at[b]] = True
    return R


def collapsed(S: np.ndarray, train: np.ndarray, test: np.ndarray, pairs: np.ndarray,
              seeds: tuple[int, ...] = score.SEEDS) -> dict:
    """Per trait, average precision with each locus counted once at its best-ranked member."""
    units = np.flatnonzero(test.any(axis=0))
    out = np.full((len(seeds), len(units)), np.nan)
    n_loci = 0
    loci = {}
    for j in units:
        members = np.flatnonzero(test[:, j])
        loci[j] = (members, position.components(pairs, members))
        n_loci += len(set(loci[j][1]))
    for si, seed in enumerate(seeds):
        rng = np.random.default_rng(seed)
        for ui, j in enumerate(units):
            cand = np.flatnonzero(~train[:, j])
            order = cand[np.lexsort((rng.random(len(cand)), -S[cand, j]))]
            members, lab = loci[j]
            label_of = dict(zip(members.tolist(), lab.tolist()))
            keep = np.ones(len(order), bool)
            y = np.zeros(len(order), bool)
            seen: set[int] = set()
            for pos in np.flatnonzero(test[order, j]):
                L = label_of[int(order[pos])]
                if L in seen:
                    keep[pos] = False        # a later member of a locus already found leaves the list
                else:
                    seen.add(L)
                    y[pos] = True
            y = y[keep]
            hits = np.cumsum(y)
            at = np.flatnonzero(y)
            out[si, ui] = (hits[at] / (at + 1)).mean()
    return {"ap": np.nanmean(out, axis=0), "loci": n_loci, "index": units}


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
id: MONDO:4
name: lung disease
is_a: MONDO:1
""".splitlines())
    R = lineage_matrix(["MONDO:1", "MONDO:2", "MONDO:4"], terms)
    ck("the lineage matrix is symmetric with an empty diagonal, and siblings are apart", (R == R.T).all() and not R.diagonal().any() and R[0, 1] and not R[1, 2])
    # one trait, six candidate targets 0..5; positives 0,1 (one locus) and 4 (alone)
    train = np.zeros((6, 1), bool)
    test = np.zeros((6, 1), bool)
    test[[0, 1, 4], 0] = True
    pairs = np.array([[0, 1], [1, 0]])
    S = np.array([[6.0], [5.0], [4.0], [3.0], [2.0], [1.0]])
    r = collapsed(S, train, test, pairs, seeds=(0,))
    # ranked 0,1,2,3,4,5: locus {0,1} found at rank 1, its second member dropped -> list 0,2,3,4,5 ; hits at ranks 1 and 4
    ck("a locus counts once and its later members leave the list", r["loci"] == 2 and abs(r["ap"][0] - (1 / 1 + 2 / 4) / 2) < 1e-12)
    S2 = np.array([[1.0], [6.0], [4.0], [3.0], [2.0], [5.0]])
    r2 = collapsed(S2, train, test, pairs, seeds=(0,))
    # ranked 1,5,2,3,4,0: locus found at rank 1 via target 1; target 0 dropped -> 1,5,2,3,4 ; hits at ranks 1 and 5
    ck("the locus is found at its best-ranked member, whichever that is", abs(r2["ap"][0] - (1 / 1 + 2 / 5) / 2) < 1e-12)
    r3 = collapsed(S, train, test, np.zeros((0, 2), np.int64), seeds=(0,))
    full = score.by_axis(S, train, test, axis=1, seeds=(0,))
    ck("with no neighbours the collapsed measure is the plain one", r3["loci"] == 3 and abs(r3["ap"][0] - full["ap"][0]) < 1e-12)
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
    R = lineage_matrix(M["traits"], terms)
    g = position.genes(registry.require("ensembl_gtf_113"))
    tg = position.target_genes(a.release, g)
    index = {t: i for i, t in enumerate(M["targets"])}
    pairs = position.neighbour_pairs(tg, index)
    S = baselines.scorers(A)
    S["B4"] = A.astype(np.float32) @ R.astype(np.float32).T
    S["B5"] = position.neighbour_count(pairs, A)
    lin, pos = (S["B4"] > 0) & ~A, (S["B5"] > 0) & ~A
    linked = lin | pos
    new = Y & ~linked
    canon = [ontology.canonical(t, terms) for t in M["traits"]]
    scored = Y.any(axis=0)
    cover = {
        "traits": len(M["traits"]), "traits_in_ontology": sum(c is not None for c in canon),
        "scored_traits_in_ontology": int(sum(c is not None for c, s in zip(canon, scored) if s)),
        "traits_with_a_lineage_partner": int(R.any(axis=1).sum()), "lineage_pairs": int(R.sum() // 2),
        "targets": len(M["targets"]), "targets_with_coordinates": int(tg.target[tg.target.isin(index)].nunique()),
        "targets_with_a_neighbour": int(len(set(pairs[:, 0].tolist()))), "neighbour_pairs": int(len(pairs) // 2),
        "window": position.WINDOW,
    }
    shares = {
        "positives": int(Y.sum()), "lineage_linked": int((Y & lin).sum()), "position_linked": int((Y & pos).sum()),
        "both": int((Y & lin & pos).sum()), "either": int((Y & linked).sum()), "neither": int(new.sum()),
        "candidates": int((~A).sum()), "candidates_lineage_linked": int(lin.sum()), "candidates_position_linked": int(pos.sum()),
        "candidates_either": int(linked.sum()),
    }
    print(f"coverage: {cover}")
    print(f"positives {shares['positives']:,}: lineage-linked {shares['lineage_linked']:,}, position-linked {shares['position_linked']:,}, "
          f"both {shares['both']:,}, either {shares['either']:,}, neither {shares['neither']:,}")
    res, per = {}, {}
    for k in NAMES:
        t0 = time.time()
        full = score.by_axis(S[k], A, Y, axis=1)
        low = score.by_axis(S[k], A, new, axis=1, exclude=Y & linked)          # linked positives removed
        neu = score.by_axis(S[k], A, new, axis=1, exclude=linked)               # every linked candidate removed
        loc = collapsed(S[k], A, Y, pairs)
        per[k] = {"full": full["ap"], "low": low["ap"], "neutral": neu["ap"], "locus": loc["ap"]}
        res[k] = {"primary": score.interval(full["ap"]), "new_links_lower_bound": score.interval(low["ap"]),
                  "new_links_neutral": score.interval(neu["ap"]), "locus_collapsed": score.interval(loc["ap"]),
                  "traits_new": int(len(low["index"])), "loci": loc["loci"]}
        r = res[k]
        print(f"  {k} {NAMES[k]:15s} all {r['primary'][0]:.4f}  new, linked positives out {r['new_links_lower_bound'][0]:.4f}  "
              f"new, all linked out {r['new_links_neutral'][0]:.4f}  locus-collapsed {r['locus_collapsed'][0]:.4f}  ({time.time() - t0:.0f} s)", flush=True)
    paired = {}
    for x, y in (("B3", "B1"), ("B4", "B0"), ("B5", "B0"), ("B5", "B1"), ("B4", "B1"), ("B3", "B5"), ("B3", "B4")):
        paired[f"{x}-{y}"] = {}
        for view in ("full", "low", "neutral", "locus"):
            m, lo, hi = score.interval(per[x][view] - per[y][view])
            paired[f"{x}-{y}"][view] = {"difference": m, "low": lo, "high": hi}
        d = paired[f"{x}-{y}"]
        print(f"  {x} minus {y}: all {d['full']['difference']:+.4f} [{d['full']['low']:+.4f}, {d['full']['high']:+.4f}]  "
              f"new/lower bound {d['low']['difference']:+.4f} [{d['low']['low']:+.4f}, {d['low']['high']:+.4f}]  "
              f"new/neutral {d['neutral']['difference']:+.4f} [{d['neutral']['low']:+.4f}, {d['neutral']['high']:+.4f}]  "
              f"locus {d['locus']['difference']:+.4f} [{d['locus']['low']:+.4f}, {d['locus']['high']:+.4f}]")
    side = lambda k: json.loads(registry.SOURCES[k].path.with_name(registry.SOURCES[k].filename + ".lineage.json").read_text())["sha256"]
    frozen = json.loads(split.frozen_path(a.cutoff).read_text())
    out = paths.DOCS_DATA / f"controls_gene_trait_{a.cutoff}.json"
    out.write_text(json.dumps({"cutoff": a.cutoff, "release": a.release, "split": {k: frozen[k]["sha256"] for k in split.PARTS},
                               "efo_sha256": side("efo_obo"), "gtf_sha256": side("ensembl_gtf_113"),
                               "coverage": cover, "shares": shares, "names": NAMES, "results": res, "paired": paired}, indent=1) + "\n")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
