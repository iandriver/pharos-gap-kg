"""Ranking targets for a trait with no link, and new targets for a trait (docs/PROTOCOL.md, "Cold start").

A third of the links first reported after the cutoff were never scored, because their trait or their
target has no training link and the neighbours baseline has nothing to say about either. Trait
lineage reaches a new trait through the training traits above and below it in the ontology. Gene
position reaches a new target through the training targets beside it on the chromosome.

    for a new trait, the score of a training target
    N1 target degree       how many training traits the target is linked to
    N2 lineage transfer    how many of the new trait's lineage relatives the target is linked to
    N3 lineage neighbours  over those relatives, 1 for a link plus the target's neighbours score for the relative
    N4, N5                 N2 and N3 with ties broken by N1

    for a training trait, the score of a new target
    P1 locus degree        how many training links the target's neighbours hold
    P2 position transfer   how many of the target's neighbours are linked to the trait
    P3                     P2 with ties broken by P1

None has a fitted parameter. A second period inside the training years (cutoff 2016, links of 2017 and
2018) is computed first and reported beside the one that counts.

Run:  python -m pharos_gap_kg.coldstart pharos319 --cutoff 2018
      python -m pharos_gap_kg.coldstart --selftest
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd

from pharos_gap_kg import baselines, controls, ontology, paths, position, registry, score, split

EARLIER = 2016
TRAIT_SCORERS = {"N0": "chance", "N1": "target degree", "N2": "lineage transfer", "N3": "lineage neighbours", "N4": "lineage transfer, ties by degree", "N5": "lineage neighbours, ties by degree"}
TARGET_SCORERS = {"P0": "chance", "P1": "locus degree", "P2": "position transfer", "P3": "position transfer, ties by locus degree"}
BY_TARGET = {"Q0": "chance", "Q1": "trait degree", "Q2": "position transfer", "Q3": "position transfer, ties by trait degree"}


def tie_break(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """A score that orders each column by `first` and, where `first` ties, by `second`."""
    out = np.zeros(first.shape, np.float64)
    span = float(second.max() - second.min()) + 1.0
    for j in range(first.shape[1]):
        rank = np.unique(first[:, j], return_inverse=True)[1]
        out[:, j] = rank + (second[:, j] - second.min()) / (span + 1.0)
    return out


def parts_of(train: pd.DataFrame, later: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """The later links no warm scorer reaches, by what is new in each."""
    T, D = set(train.target), set(train.trait)
    nt, nd = ~later.target.isin(T), ~later.trait.isin(D)
    return {"new_trait": later[nd & ~nt], "new_target": later[~nd & nt], "both_new": later[nd & nt]}


def relation(new_traits: list[str], traits: list[str], terms: dict) -> np.ndarray:
    """[new trait, training trait]: whether the two are in one lineage."""
    lin = ontology.lineage(traits + new_traits, terms)
    at = {t: i for i, t in enumerate(traits)}
    R = np.zeros((len(new_traits), len(traits)), bool)
    for i, t in enumerate(new_traits):
        for b in lin[t]:
            if b in at:
                R[i, at[b]] = True
    return R


def trait_scores(A: np.ndarray, R: np.ndarray) -> dict[str, np.ndarray]:
    """The scores of every training target for every new trait. A [target, trait], R [new trait, trait]."""
    a, r = A.astype(np.float32), R.astype(np.float32)
    n1 = np.repeat(a.sum(axis=1)[:, None], R.shape[0], axis=1)
    n2 = a @ r.T
    n3 = n2 + baselines.scorers(A)["B3"] @ r.T
    return {"N0": np.zeros_like(n1), "N1": n1, "N2": n2, "N3": n3, "N4": tie_break(n2, n1), "N5": tie_break(n3, n1)}


def beside(pairs: np.ndarray, n_train: int, n_all: int) -> np.ndarray:
    """[candidate target, training target]: whether the two are neighbours. Targets are numbered with
    the training targets first."""
    N = np.zeros((n_all - n_train, n_train), bool)
    keep = (pairs[:, 0] >= n_train) & (pairs[:, 1] < n_train)
    N[pairs[keep, 0] - n_train, pairs[keep, 1]] = True
    return N


def target_scores(A: np.ndarray, N: np.ndarray) -> dict[str, np.ndarray]:
    """The scores of every candidate target for every training trait. A [target, trait], N [candidate, target]."""
    p2 = N.astype(np.float32) @ A.astype(np.float32)
    p1 = np.repeat(p2.sum(axis=1)[:, None], A.shape[1], axis=1)
    return {"P0": np.zeros_like(p2), "P1": p1, "P2": p2, "P3": tie_break(p2, p1)}


def by_target_scores(A: np.ndarray, N: np.ndarray) -> dict[str, np.ndarray]:
    """The same pairs read the other way, [trait, candidate], for ranking the training traits of one new target."""
    p2 = (N.astype(np.float32) @ A.astype(np.float32)).T
    q1 = np.repeat(A.sum(axis=0).astype(np.float32)[:, None], N.shape[0], axis=1)
    return {"Q0": np.zeros_like(p2), "Q1": q1, "Q2": p2, "Q3": tie_break(p2, q1)}


def summary(ap: np.ndarray, groups: dict[str, np.ndarray]) -> dict:
    out = {}
    for g, m in groups.items():
        v = ap[m]
        out[g] = {"units": int(m.sum()), "value": score.interval(v) if len(v) else None}
    return out


def against(a: np.ndarray, b: np.ndarray, groups: dict[str, np.ndarray]) -> dict:
    out = {}
    for g, m in groups.items():
        if m.sum():
            mean, lo, hi = score.interval(a[m] - b[m])
            out[g] = {"difference": mean, "low": lo, "high": hi, "half_width": (hi - lo) / 2, "units_higher": int((a[m] > b[m]).sum()), "units": int(m.sum())}
    return out


def evaluate(release: str, train: pd.DataFrame, later: pd.DataFrame, terms: dict, tg: pd.DataFrame, all_targets: list[str], log=print) -> dict:
    """Both tasks for one period: `train` the links a scorer may see, `later` the later links it could not be scored on."""
    P = parts_of(train, later)
    targets, traits = sorted(train.target.unique()), sorted(train.trait.unique())
    ti, di = {t: i for i, t in enumerate(targets)}, {d: i for i, d in enumerate(traits)}
    A = np.zeros((len(targets), len(traits)), bool)
    A[train.target.map(ti).to_numpy(), train.trait.map(di).to_numpy()] = True
    out = {"links": {k: int(len(v)) for k, v in P.items()}, "training": {"links": int(len(train)), "targets": len(targets), "traits": len(traits)}}

    # ---- new traits ----
    nt = P["new_trait"]
    new_traits = sorted(nt.trait.unique())
    ni = {d: i for i, d in enumerate(new_traits)}
    Y = np.zeros((len(targets), len(new_traits)), bool)
    Y[nt.target.map(ti).to_numpy(), nt.trait.map(ni).to_numpy()] = True
    R = relation(new_traits, traits, terms)
    S = trait_scores(A, R)
    none = np.zeros_like(Y)
    pairs = position.neighbour_pairs(tg, ti)
    reach = R.any(axis=1)
    groups = {"all": np.ones(len(new_traits), bool), "with_a_relative": reach, "without": ~reach}
    in_onto = sum(ontology.canonical(d, terms) is not None for d in new_traits)
    linked = (S["N2"] > 0) & Y
    out["new_traits"] = {"units": len(new_traits), "candidates_per_unit": len(targets), "positives": int(Y.sum()), "in_the_ontology": int(in_onto),
                         "with_a_relative": int(reach.sum()), "positives_of_those": int(Y[:, reach].sum()), "positives_on_a_target_linked_to_a_relative": int(linked.sum()),
                         "relatives_per_reached_trait_median": float(np.median(R.sum(axis=1)[reach])) if reach.any() else None, "results": {}, "by_locus": {}, "paired": {}}
    per, loc = {}, {}
    for k in TRAIT_SCORERS:
        t0 = time.time()
        per[k] = score.by_axis(S[k], none, Y, axis=1)["ap"]
        loc[k] = controls.collapsed(S[k], none, Y, pairs)["ap"]
        out["new_traits"]["results"][k] = summary(per[k], groups)
        out["new_traits"]["by_locus"][k] = summary(loc[k], groups)
        r = out["new_traits"]["results"][k]
        log(f"    new traits  {k} {TRAIT_SCORERS[k]:36s} all {r['all']['value'][0]:.4f}  with a relative {r['with_a_relative']['value'][0]:.4f}  without {r['without']['value'][0]:.4f}  "
            f"by locus {out['new_traits']['by_locus'][k]['all']['value'][0]:.4f}  ({time.time() - t0:.0f} s)")
    for x, y in (("N1", "N0"), ("N2", "N1"), ("N4", "N1"), ("N5", "N4"), ("N5", "N1")):
        out["new_traits"]["paired"][f"{x}-{y}"] = {"primary": against(per[x], per[y], groups), "by_locus": against(loc[x], loc[y], groups)}
    d = out["new_traits"]["paired"]["N4-N1"]["primary"]["all"]
    out["new_traits"]["lineage_transfer_works"] = bool(d["low"] > 0)
    log(f"    new traits  N4 minus N1 over all {d['units']} traits: {d['difference']:+.4f} [{d['low']:+.4f}, {d['high']:+.4f}]  works: {out['new_traits']['lineage_transfer_works']}")

    # ---- new targets ----
    cand = sorted(set(all_targets) - set(targets))
    ci = {t: i for i, t in enumerate(cand)}
    ng = P["new_target"]
    if not set(ng.target) <= set(cand):
        raise SystemExit("REFUSED: a later link's target is not a target of the release")
    Y2 = np.zeros((len(cand), len(traits)), bool)
    Y2[ng.target.map(ci).to_numpy(), ng.trait.map(di).to_numpy()] = True
    index = {**ti, **{t: len(targets) + i for t, i in ci.items()}}
    N = beside(position.neighbour_pairs(tg, index), len(targets), len(index))
    S2 = target_scores(A, N)
    none2 = np.zeros_like(Y2)
    units = Y2.any(axis=0)
    g2 = {"all": np.ones(int(units.sum()), bool)}
    placed = set(tg.target)
    out["new_targets"] = {"units": int(units.sum()), "candidates_per_unit": len(cand), "positives": int(Y2.sum()), "new_targets_with_a_positive": int(Y2.any(axis=1).sum()),
                          "candidates_placed_on_a_chromosome": int(sum(c in placed for c in cand)), "candidates_with_a_training_neighbour": int(N.any(axis=1).sum()),
                          "positives_on_a_target_beside_a_linked_one": int(((S2["P2"] > 0) & Y2).sum()), "results": {}, "paired": {}}
    per2 = {}
    for k in TARGET_SCORERS:
        per2[k] = score.by_axis(S2[k], none2, Y2, axis=1)["ap"]
        out["new_targets"]["results"][k] = summary(per2[k], g2)
        log(f"    new targets {k} {TARGET_SCORERS[k]:36s} all {out['new_targets']['results'][k]['all']['value'][0]:.4f}")
    for x, y in (("P1", "P0"), ("P2", "P1"), ("P3", "P1")):
        out["new_targets"]["paired"][f"{x}-{y}"] = against(per2[x], per2[y], g2)
    d = out["new_targets"]["paired"]["P3-P1"]["all"]
    out["new_targets"]["position_transfer_works"] = bool(d["low"] > 0)
    log(f"    new targets P3 minus P1 over all {d['units']} traits: {d['difference']:+.4f} [{d['low']:+.4f}, {d['high']:+.4f}]  works: {out['new_targets']['position_transfer_works']}")

    # the same links read by new target: which training traits will it be linked to
    S3 = by_target_scores(A, N)
    Y3 = Y2.T
    none3 = np.zeros_like(Y3)
    tu = np.flatnonzero(Y3.any(axis=0))
    has = N.any(axis=1)[tu]
    g3 = {"all": np.ones(len(tu), bool), "with_a_neighbour": has, "without": ~has}
    out["by_new_target"] = {"units": int(len(tu)), "candidates_per_unit": len(traits), "with_a_neighbour": int(has.sum()), "results": {}, "paired": {}}
    per3 = {}
    for k in BY_TARGET:
        per3[k] = score.by_axis(S3[k], none3, Y3, axis=1)["ap"]
        out["by_new_target"]["results"][k] = summary(per3[k], g3)
        r = out["by_new_target"]["results"][k]
        log(f"    by new target {k} {BY_TARGET[k]:36s} all {r['all']['value'][0]:.4f}  with a neighbour {r['with_a_neighbour']['value'][0]:.4f}")
    for x, y in (("Q1", "Q0"), ("Q2", "Q1"), ("Q3", "Q1")):
        out["by_new_target"]["paired"][f"{x}-{y}"] = against(per3[x], per3[y], g3)
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    first = np.array([[2.0], [0.0], [2.0], [0.0], [1.0]])
    second = np.array([[1.0], [9.0], [5.0], [3.0], [100.0]])
    order = np.argsort(-tie_break(first, second)[:, 0]).tolist()
    ck("a tie in the first score is broken by the second, and the second never overturns the first", order == [2, 0, 4, 1, 3])
    f2 = np.array([[0.3], [0.3000001], [0.1]])
    ck("first scores that differ by very little are still told apart", np.argsort(-tie_break(f2, np.array([[9.0], [0.0], [5.0]]))[:, 0]).tolist() == [1, 0, 2])

    train = pd.DataFrame([("t1", "a"), ("t2", "a"), ("t2", "b"), ("t3", "b")], columns=["target", "trait"])
    later = pd.DataFrame([("t1", "z"), ("t9", "a"), ("t9", "z"), ("t3", "a")], columns=["target", "trait"])
    P = parts_of(train, later[:3])
    ck("a later link is filed by what is new in it: its trait, its target, or both",
       P["new_trait"].values.tolist() == [["t1", "z"]] and P["new_target"].values.tolist() == [["t9", "a"]] and P["both_new"].values.tolist() == [["t9", "z"]])

    terms = ontology.load("""[Term]
id: EFO:1
name: top

[Term]
id: EFO:2
name: child
is_a: EFO:1

[Term]
id: EFO:3
name: grandchild
is_a: EFO:2

[Term]
id: EFO:4
name: elsewhere
""".splitlines())
    R = relation(["EFO:3", "EFO:4"], ["EFO:1", "EFO:2"], terms)
    ck("a new trait is related to the training traits above it, and an unrelated one to none", R.tolist() == [[True, True], [False, False]])
    # 3 targets, 2 training traits; target 0 linked to both, target 1 to the first, target 2 to neither
    A = np.array([[1, 1], [1, 0], [0, 0]], bool)
    S = trait_scores(A, R)
    ck("degree is the same for every new trait, and lineage transfer counts the relatives a target is linked to",
       S["N1"][:, 0].tolist() == [2, 1, 0] and S["N2"][:, 0].tolist() == [2, 1, 0] and S["N2"][:, 1].tolist() == [0, 0, 0])
    b3 = baselines.scorers(A)["B3"]
    ck("lineage neighbours adds the target's neighbours score for each relative to its links", np.allclose(S["N3"][:, 0], A.sum(axis=1) + b3.sum(axis=1)))
    ck("for a trait with no relative, lineage with ties by degree is degree", np.argsort(-S["N4"][:, 1], kind="stable").tolist() == np.argsort(-S["N1"][:, 1], kind="stable").tolist())
    # targets 0,1,2 are training targets; candidates 3 and 4; 3 is beside 0 and 1, 4 beside nothing
    pairs = np.array([[0, 3], [3, 0], [1, 3], [3, 1], [0, 1], [1, 0]])
    N = beside(pairs, 3, 5)
    ck("a candidate is beside the training targets it neighbours, and pairs among training targets are left out", N.tolist() == [[True, True, False], [False, False, False]])
    S2 = target_scores(A, N)
    ck("position transfer counts the candidate's neighbours linked to the trait, and locus degree all their links",
       S2["P2"].tolist() == [[2, 1], [0, 0]] and S2["P1"].tolist() == [[3, 3], [0, 0]])
    S3 = by_target_scores(A, N)
    ck("read by new target, the candidates are traits and the plain baseline is the trait's degree", S3["Q2"].tolist() == [[2, 0], [1, 0]] and S3["Q1"][:, 0].tolist() == [2, 1])
    g = {"all": np.array([True, True, True]), "some": np.array([True, False, False]), "none": np.array([False, False, False])}
    s = summary(np.array([0.5, 0.1, 0.3]), g)
    ck("a summary gives each group its units, and no value where there are none", s["all"]["units"] == 3 and abs(s["some"]["value"][0] - 0.5) < 1e-12 and s["none"]["value"] is None)
    d = against(np.array([0.5, 0.1, 0.3]), np.array([0.2, 0.1, 0.4]), g)
    ck("a paired difference is taken unit by unit within a group", abs(d["all"]["difference"] - (0.3 + 0 - 0.1) / 3) < 1e-12 and d["all"]["units_higher"] == 1 and "none" not in d)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    import pyarrow.parquet as pq
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("release")
    ap.add_argument("--cutoff", type=int, required=True)
    a = ap.parse_args(argv)
    parts = split.load(a.cutoff)
    terms = ontology.read(registry.require("efo_obo"))
    tg = position.target_genes(a.release, position.genes(registry.require("ensembl_gtf_113")))
    prot = pq.read_table(paths.tables(a.release) / "protein.parquet", columns=["uniprot"]).to_pandas()
    all_targets = sorted(set("UniProtKB:" + prot.uniprot.dropna()))
    earlier = split.split(parts["train"], EARLIER)
    out = {"cutoff": a.cutoff, "release": a.release, "targets_in_the_release": len(all_targets), "scorers": {**TRAIT_SCORERS, **TARGET_SCORERS, **BY_TARGET}, "periods": {}}
    print(f"  cutoff {EARLIER}, inside the training years: the links of {EARLIER + 1} to {a.cutoff} whose trait or target was new in {EARLIER}", flush=True)
    out["periods"][str(EARLIER)] = evaluate(a.release, earlier["train"], earlier["unscored"], terms, tg, all_targets, log=lambda s: print(s, flush=True))
    print(f"  cutoff {a.cutoff}: the later links the test left out", flush=True)
    out["periods"][str(a.cutoff)] = evaluate(a.release, parts["train"], parts["unscored"], terms, tg, all_targets, log=lambda s: print(s, flush=True))
    side = lambda k: json.loads(registry.SOURCES[k].path.with_name(registry.SOURCES[k].filename + ".lineage.json").read_text())["sha256"]
    frozen = json.loads(split.frozen_path(a.cutoff).read_text())
    out.update({"split": {k: frozen[k]["sha256"] for k in split.PARTS}, "efo_sha256": side("efo_obo"), "gtf_sha256": side("ensembl_gtf_113"), "window": position.WINDOW})
    f = paths.DOCS_DATA / f"coldstart_gene_trait_{a.cutoff}.json"
    f.write_text(json.dumps(out, indent=1) + "\n")
    print(f"  -> {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
