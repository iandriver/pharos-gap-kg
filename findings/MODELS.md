# Fitted models: tuning inside the training years

8 Oct 2026. `python -m pharos_gap_kg.model tune pharos319 --cutoff 2018 --model <name>` writes
`docs/data/model_<name>_2018.json`, and the `validation` stage writes
`docs/data/model_<name>_2018_validation.json`. The models, the grid of settings and the rule for
choosing among them were committed to `docs/PROTOCOL.md` before anything was trained.

**No test link has been read.** Everything below is on the inner split of the training links: fitted
on the 57,411 links first reported in or before 2016, validated on the 22,610 links of 2017 and 2018,
over 284 traits. On those validation links the baselines score 0.0084 (chance), 0.0394 (target degree)
and 0.0712 (neighbours).

Two of the four declared models are tuned, both on a laptop CPU. The two that join the context tables
on are not.

| model | settings tried | chosen setting | validation | against neighbours |
|---|---|---|---|---|
| `mf` | 8 | vectors of 64, rate 0.001, weight 1e-4, 300 passes | 0.0546 | -0.0166 [-0.0273, -0.0082] |
| `lightgcn` | 24 | vectors of 128, 3 steps, rate 0.005, weight 1e-6, 60 passes | 0.0743 | +0.0031 [-0.0058, +0.0120] |

The last column is the paired difference over the 284 traits with a 95% interval. It flatters the
model, because the setting was chosen on these same links.

- **Matrix factorisation is worse than the neighbours baseline.** All eight settings score between
  0.0518 and 0.0546. The chosen one peaks at the cap of 300 passes, with a flat curve behind it:
  0.0539 at 250 passes, 0.0530 at 280, 0.0546 at 300. A longer fit would not close a gap of 0.0166.
- **LightGCN's best setting cannot be told from the neighbours baseline.** It is higher on 157 of the
  284 traits. The interval's half-width, 0.0089, is the smallest difference this comparison resolves.
- **More smoothing scores higher.** With one step the settings score 0.0644 to 0.0670, with two 0.0691
  to 0.0726, with three 0.0703 to 0.0743. The chosen setting sits at the edge of the declared grid in
  both vector length and steps, so a longer vector or a fourth step was not tried.
- **Both beat target degree**: LightGCN by +0.0349 [+0.0172, +0.0553] and matrix factorisation by
  +0.0152 [+0.0001, +0.0319].

Fitting all 8 and all 24 settings took 7 and 43 minutes.

## What follows

- The two settings are committed. Under the protocol each model is next fitted on all 99,435 training
  links under five seeds and scored on the test links once.
- The grid's edge is the open question for LightGCN. Widening the grid is an amendment to make before
  any model is scored on the test links, not after.
- Validation gives no sign that either model beats the parameter-free neighbours baseline.
