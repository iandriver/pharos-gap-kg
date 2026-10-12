"""The ranking the question layer serves is the ranking that was tested.

Two artifacts answer the same question: the bundle, read by `ask`, and the prospective split's
training matrix, on which the combiner was scored. The prospective split keys a trait by its live
ontology term, so two of the bundle's traits can be one trait there. For a trait that was not merged
with another, and a target that holds no merged trait, the two must give the same score.

Skipped where the bundle or the prospective split has not been built."""
import json

import numpy as np
import pandas as pd
import pytest

from pharos_gap_kg import ask, combine, model, paths, prospect, split

BUNDLE = paths.bundle("pharos319")
TRAIN = split.out_dir(prospect.CUTOFF) / "train.parquet"
pytestmark = pytest.mark.skipif(not ((BUNDLE / "manifest.json").is_file() and TRAIN.is_file()), reason="the bundle or the prospective split is not built")


def test_served_scores_are_the_tested_scores():
    G = ask.Graph(BUNDLE)
    assert G.weights is not None, "the bundle names no scorer; refresh it with:  python -m pharos_gap_kg.bundle pharos319 --backtest-only"
    rec = json.loads(prospect.weights_path().read_text())["combiners"][prospect.SERVED]
    assert G.weights == dict(zip(("B1", "B3", "B5"), rec["weights"])) and rec["scores"] == ["B1", "B3", "B5"] and rec["form"] == "F1"
    M = prospect.train_matrix(prospect.training())
    A = M["train"]
    pairs, R, _ = model.load_structure("pharos319", M["targets"], M["traits"], ("position", "lineage"))
    S = combine.combined(combine.raw_scores(A, pairs, R), A, tuple(rec["scores"]), rec["form"], rec["weights"])
    # the prospective split keys a trait by its live ontology term; a bundle trait is merged if another shares its key
    key = np.array([c if c else d for d, c in zip(G.d.id, G.d.live_term)])
    counts = pd.Series(key).value_counts()
    merged = np.array([counts[k] > 1 for k in key])
    col = {t: j for j, t in enumerate(M["traits"])}
    kept = [i for i in range(len(G.d)) if not merged[i] and key[i] in col]
    assert len(kept) > 2000 and 0 < merged.sum() < 100
    mask = G._mask(G.latest)
    holds_merged = np.zeros(len(G.t), bool)
    holds_merged[G.lt[mask & merged[G.ld]]] = True
    at = np.array([G.ti[t] for t in M["targets"]])
    clean = ~holds_merged[at]
    assert clean.mean() > 0.5
    rng = np.random.default_rng(0)
    checked = 0
    for i in rng.choice(kept, size=40, replace=False):
        served = G._score_trait(int(i), mask)[0][at]
        tested = S[:, col[key[i]]]
        linked = A[:, col[key[i]]]
        ok = clean & ~linked
        assert np.allclose(served[ok], tested[ok], rtol=1e-5, atol=1e-5), G.d.id[i]
        checked += int(ok.sum())
    assert checked > 100_000
