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

## Added 8 Oct 2026: the settings widened for LightGCN

The chosen LightGCN setting above sat at the edge of the declared settings, so 20 more were added by a
protocol amendment before any model was scored on the test links: vectors of 256 with 3, 4 or 5 steps,
and vectors of 128 with 4 or 5 steps. The first 24 results were kept and not fitted again. The table
and bullets above describe the first 24 and are left as written.

| | settings | chosen setting | validation | against neighbours |
|---|---|---|---|---|
| first settings | 24 | vectors of 128, 3 steps, rate 0.005, weight 1e-6, 60 passes | 0.0743 | +0.0031 [-0.0058, +0.0120] |
| all settings | 44 | vectors of 256, 5 steps, rate 0.005, weight 1e-4, 60 passes | 0.0766 | +0.0054 [-0.0034, +0.0145] |

- **The wider settings raise validation by 0.0023, and the model still cannot be told from the
  neighbours baseline.** The interval spans zero, and the setting is higher on 152 of the 284 traits.
  The comparison flatters the model, as before.
- **The chosen setting is again at the edge, and the gains have flattened.** The best of each added
  combination, by vector length and steps:

  | | 3 steps | 4 steps | 5 steps |
  |---|---|---|---|
  | vectors of 128 | 0.0743 | 0.0759 | 0.0764 |
  | vectors of 256 | 0.0760 | 0.0765 | 0.0766 |

  The last step adds 0.0001 at 256 and the five best settings lie within 0.0003 of each other. As the
  amendment said, nothing more is added.
- Against target degree the chosen setting is +0.0372 [+0.0189, +0.0580].

The 20 added settings took 95 minutes of fitting, run as three processes.

**Two checks were added to the code before any scoring**, after a review of it by another session
(`docs/AUDIT_20261008.md`):

- A kept result is tied to the code that made it and, for the context models, to the content of the
  graph and not its counts. The fitting functions at the three commits that produced the 52 kept
  results are identical to the present ones, so those results were carried over, and the two records
  are unchanged.
- `final` now refuses a record whose settings are not the ones the protocol declares, one made by
  other fitting code, and one whose file is not the committed one a `freeze` step recorded. A record
  of the first 24 settings can no longer authorise a scoring.

Both settings are frozen. Neither model has been scored on the test links.
