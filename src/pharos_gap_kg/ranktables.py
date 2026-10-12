"""The rankings the front end serves, as tables, computed by matrix code that does not call the question layer.

`ask` computes a ranking for one trait or one target at a time. This computes all of them at once from
the bundle's tables, a different way: so the usability test can read its expected answers from here,
give the same tables to the SQL condition, and check that the question layer agrees. Where the two
disagree, one of them is wrong.

    candidates        the highest-ranked candidate targets for every trait, with the parts of each score
    rank_rates        how often candidates at each rank were reported in the test on later links, by group
    no_link_targets   for every trait, the targets with no link that lie beside a gene with the link
    no_link_traits    for every target with no link, the traits its neighbouring genes are linked to

Run:  python -m pharos_gap_kg.ranktables --selftest
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

GROUPS = {"linked": "explained", "disease_like_not_linked": "not explained, trait is a disease or phenotype", "not_linked": "not explained, other trait"}


def matrices(nodes: pd.DataFrame, links: pd.DataFrame, lineage: pd.DataFrame, neighbours: pd.DataFrame, weights: dict) -> dict:
    """The link matrix and everything the two rankings are made of. Targets and traits keep the order of `nodes`."""
    from scipy import sparse
    t = nodes[nodes.kind == "target"].reset_index(drop=True)
    d = nodes[nodes.kind == "trait"].reset_index(drop=True)
    ti, di = {v: i for i, v in enumerate(t.id)}, {v: i for i, v in enumerate(d.id)}
    A = np.zeros((len(t), len(d)), np.float64)
    A[links.target.map(ti).to_numpy(), links.trait.map(di).to_numpy()] = 1.0
    deg_t, deg_d = A.sum(axis=1), A.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        sim = np.where(np.outer(deg_d, deg_d) > 0, (A.T @ A) / np.sqrt(np.outer(deg_d, deg_d)), 0.0)
    np.fill_diagonal(sim, 0.0)
    nb = A @ sim
    N = sparse.csr_matrix((np.ones(len(neighbours)), (neighbours.a.map(ti).to_numpy(), neighbours.b.map(ti).to_numpy())), shape=(len(t), len(t)))
    pos = np.asarray(N @ A)
    score = weights["B3"] * np.log1p(nb) + weights["B1"] * np.log1p(deg_t)[:, None] + weights["B5"] * np.log1p(pos)
    score[deg_t == 0] = 0.0
    R = np.zeros((len(d), len(d)))
    R[lineage.a.map(di).to_numpy(), lineage.b.map(di).to_numpy()] = 1.0
    wording = (A @ R.T) > 0                       # the target is linked to a trait in the query trait's lineage
    return {"t": t, "d": d, "A": A > 0, "deg_t": deg_t, "deg_d": deg_d, "neighbours": nb, "position": pos, "score": score, "wording": wording,
            "locus": np.asarray(N @ deg_t).ravel()}


def candidates(M: dict, top: int = 200) -> pd.DataFrame:
    """For every trait with a link, its `top` highest candidates: by score, ties by the target's links, then by the order of the node table."""
    rows = []
    idx = np.arange(len(M["t"]))
    sym, tid, did = M["t"].name.to_numpy(), M["t"].id.to_numpy(), M["d"].id.to_numpy()
    for j in np.flatnonzero(M["deg_d"] > 0):
        s = M["score"][:, j]
        order = np.lexsort((idx, -M["deg_t"], -s))
        order = order[~M["A"][order, j] & (s[order] > 0)][:top]
        for rank, i in enumerate(order, 1):
            rows.append((did[j], rank, tid[i], sym[i], float(s[i]), float(M["neighbours"][i, j]), int(M["deg_t"][i]), int(M["position"][i, j]),
                         int(M["position"][i, j] > 0 or M["wording"][i, j]), int((s[~M["A"][:, j]] > s[i]).sum()) + 1))
    return pd.DataFrame(rows, columns=["trait_id", "rank", "target_id", "symbol", "score", "shared_trait_neighbours", "traits_already_linked",
                                       "genes_nearby_with_this_link", "explained", "targets_scoring_higher_plus_one"])


def rank_rates(calibration: dict) -> pd.DataFrame:
    """The rate the front end quotes beside a candidate for a trait: by its ranks and its group."""
    return pd.DataFrame([(b, GROUPS[g], round(100 * v[g]["share"], 1)) for b, v in calibration["by_trait"].items() for g in GROUPS],
                        columns=["ranks", "candidate_group", "percent"])


def no_link(M: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The two position lists for targets with no link: by trait and by target."""
    U = np.flatnonzero(M["deg_t"] == 0)
    P = M["position"][U]
    sym, tid, did, dname = M["t"].name.to_numpy(), M["t"].id.to_numpy(), M["d"].id.to_numpy(), M["d"].name.to_numpy()
    by_trait, by_target = [], []
    for j in np.flatnonzero((P > 0).any(axis=0)):
        listed = np.flatnonzero(P[:, j] > 0)
        order = sorted(listed, key=lambda u: (-P[u, j], -M["locus"][U[u]], U[u]))
        by_trait += [(did[j], r, tid[U[u]], sym[U[u]], int(P[u, j]), int(M["locus"][U[u]])) for r, u in enumerate(order, 1)]
    for u in np.flatnonzero((P > 0).any(axis=1)):
        listed = np.flatnonzero(P[u] > 0)
        order = sorted(listed, key=lambda j: (-P[u, j], -M["deg_d"][j], j))
        by_target += [(tid[U[u]], sym[U[u]], r, did[j], dname[j], int(P[u, j]), int(M["deg_d"][j])) for r, j in enumerate(order, 1)]
    return (pd.DataFrame(by_trait, columns=["trait_id", "rank", "target_id", "symbol", "genes_nearby_with_this_link", "links_held_by_genes_nearby"]),
            pd.DataFrame(by_target, columns=["target_id", "symbol", "rank", "trait_id", "trait_name", "genes_nearby_with_this_link", "targets_linked_to_this_trait"]))


def build(bundle: Path, top: int = 200) -> dict[str, pd.DataFrame]:
    bundle = Path(bundle)
    bt = json.loads((bundle / "backtest.json").read_text())
    if "scorer" not in bt or "no_link" not in bt:
        raise SystemExit(f"the bundle at {bundle} names no scorer or does not record the position rankings; refresh it with:  python -m pharos_gap_kg.bundle <release> --backtest-only")
    rd = lambda n: pd.read_parquet(bundle / f"{n}.parquet")
    M = matrices(rd("nodes"), rd("links"), rd("lineage"), rd("neighbours"), bt["scorer"]["weights"])
    a, b = no_link(M)
    return {"candidates": candidates(M, top), "rank_rates": rank_rates(bt["calibration"]), "no_link_targets": a, "no_link_traits": b,
            "targets_with_no_link": M["t"][M["deg_t"] == 0][["id", "name"]].reset_index(drop=True)}


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    T = lambda i, s: {"id": f"t{i}", "kind": "target", "name": s}
    D = lambda i, s: {"id": f"d{i}", "kind": "trait", "name": s}
    nodes = pd.DataFrame([T(0, "AAA"), T(1, "BBB"), T(2, "CCC"), T(3, "DDD"), T(4, "EEE"), D(0, "x"), D(1, "y"), D(2, "z")])
    # AAA: x, y. BBB: x, y, z. CCC: y. DDD and EEE: no link. DDD lies beside CCC; EEE beside nothing.
    links = pd.DataFrame([("t0", "d0"), ("t0", "d1"), ("t1", "d0"), ("t1", "d1"), ("t1", "d2"), ("t2", "d1")], columns=["target", "trait"])
    lineage = pd.DataFrame([("d0", "d2", "narrower"), ("d2", "d0", "broader")], columns=["a", "b", "relation"])
    nb = pd.DataFrame([("t2", "t3", 100), ("t3", "t2", 100), ("t0", "t2", 9), ("t2", "t0", 9)], columns=["a", "b", "gap"])
    w = {"B1": 1 / 32, "B3": 1.0, "B5": 0.25}
    M = matrices(nodes, links, lineage, nb, w)
    c = candidates(M)
    z = c[c.trait_id == "d2"]
    # for z (held by BBB): AAA shares x and y with BBB; CCC shares y. AAA is linked to x, a wording of z: explained
    ck("a trait's candidates are the linked targets without it, best first", z.symbol.tolist() == ["AAA", "CCC"] and z["rank"].tolist() == [1, 2])
    sim_xz = 1 / np.sqrt(2 * 1)
    sim_yz = 1 / np.sqrt(3 * 1)
    want = np.log1p(sim_xz + sim_yz) + np.log1p(2) / 32
    ck("a score is log(1 + neighbours) plus the two weighted parts", abs(z.score.iloc[0] - want) < 1e-12 and abs(z.shared_trait_neighbours.iloc[0] - (sim_xz + sim_yz)) < 1e-12)
    ck("a candidate linked to another wording of the trait is explained, and one that is not is not", z.explained.tolist() == [1, 0])
    x = c[c.trait_id == "d0"]
    ck("a gene beside one that has the link counts it, and is explained by it", x.symbol.tolist() == ["CCC"] and x.genes_nearby_with_this_link.iloc[0] == 1 and x.explained.iloc[0] == 1)
    ck("a target with no link is not a candidate", not c.symbol.isin(["DDD", "EEE"]).any())
    a, b = no_link(M)
    ck("a target with no link is listed for the traits of the gene beside it", a.values.tolist() == [["d1", 1, "t3", "DDD", 1, 1]] and b.values.tolist() == [["t3", "DDD", 1, "d1", "y", 1, 3]])
    ck("a target with no link and nothing beside it is in neither list", "EEE" not in set(a.symbol) | set(b.symbol))
    cal = {"by_trait": {"1-10": {g: {"share": s} for g, s in zip(GROUPS, (0.111, 0.051, 0.04))}}}
    r = rank_rates(cal)
    ck("a rate is given for each group of candidates at each rank", r.values.tolist() == [["1-10", "explained", 11.1], ["1-10", "not explained, trait is a disease or phenotype", 5.1],
                                                                                           ["1-10", "not explained, other trait", 4.0]])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else 0)
