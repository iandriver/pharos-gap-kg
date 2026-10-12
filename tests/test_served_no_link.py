"""The position rankings served for targets with no link are the rankings that were scored.

Two artifacts answer the same question: the bundle, read by `ask`, and the prospective split's
training matrix, on which the cold-start scorers were run. For a trait the prospective split did not
merge with another, the number of linked genes beside each target with no link must be the same, and
so must the order `ask` lists them in.

Skipped where the bundle or the prospective split has not been built."""
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest

from pharos_gap_kg import ask, coldstart, paths, position, prospect, registry, split

BUNDLE = paths.bundle("pharos319")
TRAIN = split.out_dir(prospect.CUTOFF) / "train.parquet"
pytestmark = pytest.mark.skipif(not ((BUNDLE / "manifest.json").is_file() and TRAIN.is_file()), reason="the bundle or the prospective split is not built")


def test_served_position_counts_are_the_scored_ones():
    G = ask.Graph(BUNDLE)
    assert "no_link" in (G.backtest or {}), "the bundle does not record the position rankings; refresh it with:  python -m pharos_gap_kg.bundle pharos319 --backtest-only"
    M = prospect.train_matrix(prospect.training())
    A, targets, traits = M["train"], M["targets"], M["traits"]
    prot = pq.read_table(paths.tables("pharos319") / "protein.parquet", columns=["uniprot"]).to_pandas()
    cand = sorted(set("UniProtKB:" + prot.uniprot.dropna()) - set(targets))
    tg = position.target_genes("pharos319", position.genes(registry.require("ensembl_gtf_113")))
    index = {**{t: i for i, t in enumerate(targets)}, **{t: len(targets) + i for i, t in enumerate(cand)}}
    N = coldstart.beside(position.neighbour_pairs(tg, index), len(targets), len(index))
    S = coldstart.target_scores(A, N)
    key = np.array([c if c else d for d, c in zip(G.d.id, G.d.live_term)])
    counts = pd.Series(key).value_counts()
    col = {t: j for j, t in enumerate(traits)}
    kept = [i for i in range(len(G.d)) if counts[key[i]] == 1 and key[i] in col]
    mask = G._mask(G.latest)
    at = np.array([G.ti[c] for c in cand])
    assert not np.isin(at, G.lt[mask]).any() and len(cand) == 2407           # the bundle's targets with no link are the split's candidates
    # the tie-breaker: the links held by the genes beside a target. The bundle counts a merged trait under each
    # of its ids, so it can only be higher than the split's count, and is equal where no gene nearby holds one.
    merged = np.array([counts[k] > 1 for k in key])
    holds_merged = np.zeros(len(G.t), bool)
    holds_merged[G.lt[mask & merged[G.ld]]] = True
    deg = np.bincount(G.lt[mask], minlength=len(G.t))
    served = np.array([sum(deg[o] for o, _ in G.neigh.get(int(i), [])) for i in at])
    clean = np.array([not any(holds_merged[o] for o, _ in G.neigh.get(int(i), [])) for i in at])
    tested = S["P1"][:, 0]
    assert (served >= tested).all() and np.array_equal(served[clean], tested[clean]) and clean.sum() > 500
    rng = np.random.default_rng(0)
    listed = same_order = 0
    for i in rng.choice(kept, size=60, replace=False):
        order, pos, locus = G._no_link_for_trait(int(i), mask)
        j = col[key[i]]
        assert np.array_equal(pos[at], S["P2"][:, j]), G.d.id[i]
        assert all(locus[int(at[c])] == served[c] for c in np.flatnonzero(S["P2"][:, j] > 0))
        want = sorted(np.flatnonzero(S["P2"][:, j] > 0).tolist(), key=lambda c: (-S["P2"][c, j], -S["P1"][c, j], at[c]))
        if all(clean[c] for c in want):                    # every listed target has the tested tie-breaker: the order must be the tested order
            assert [int(at[c]) for c in want] == order, G.d.id[i]
            same_order += 1
        listed += len(order)
    assert listed > 1000 and same_order >= 5
