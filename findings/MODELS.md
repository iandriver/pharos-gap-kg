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

## Added 8 Oct 2026: both models scored on the test links, once

`python -m pharos_gap_kg.model final pharos319 --cutoff 2018 --model <name>`, then `score`, which
writes `docs/data/model_<name>_2018_scored.json`. Each model was fitted on all 99,435 training links
under seeds 0 to 4 with its frozen setting, and its score for a pair is the mean of the five. Before
the scoring, the controls were run again and reproduced the committed baseline file exactly. Each
scoring wrote a marker before it read a test link, and a second start is refused.

The primary measure, mean average precision over the 540 scored traits, on the four views of the
protocol:

| | all links | new links, lower bound | new links, neutral | by locus |
|---|---|---|---|---|
| chance | 0.0178 | | | |
| target degree | 0.0733 | 0.0496 | 0.0556 | 0.0255 |
| neighbours | 0.0969 | 0.0539 | 0.0662 | 0.0512 |
| `mf` | 0.0906 | 0.0465 | 0.0564 | 0.0662 |
| `lightgcn` | 0.0825 | 0.0433 | 0.0530 | 0.0527 |

Each model against the neighbours baseline, paired over traits, with 95% intervals:

| | all links | new links, lower bound | new links, neutral | by locus |
|---|---|---|---|---|
| `mf` | -0.0063 [-0.0155, +0.0025] | -0.0075 [-0.0143, -0.0012] | -0.0098 [-0.0185, -0.0023] | +0.0150 [+0.0079, +0.0222] |
| `lightgcn` | -0.0144 [-0.0232, -0.0064] | -0.0106 [-0.0180, -0.0043] | -0.0133 [-0.0220, -0.0056] | +0.0015 [-0.0055, +0.0078] |

- **Neither fitted model beats the neighbours baseline.** LightGCN is worse than it on all links.
  Matrix factorisation cannot be told from it: the interval spans zero, and its half-width, 0.0090, is
  the smallest difference this test resolves. By the protocol's rule neither model helps.
- **On new links both are worse than neighbours, both ways.**
- **By locus, matrix factorisation is better than neighbours**, 0.0662 against 0.0512, and LightGCN
  ties it. This is the one view of four in which a fitted model is ahead. The protocol asks for the
  primary measure, new links and loci together before a claim, so it is reported and not claimed.
- **Against target degree** matrix factorisation is ahead on all links (+0.0173 [+0.0046, +0.0310])
  and by locus (+0.0407 [+0.0294, +0.0535]) and not on new links (-0.0032 [-0.0119, +0.0043] and
  +0.0008 [-0.0085, +0.0090]). LightGCN is ahead by locus only (+0.0272 [+0.0174, +0.0372]); on all
  links it is +0.0092 [-0.0025, +0.0215] and on new links -0.0064 [-0.0152, +0.0011] and -0.0026
  [-0.0120, +0.0057].

**Validation ordered the two models the other way, and did not predict the comparison with the
baseline.**

| | validation, one seed | test, each seed alone | test, mean of five seeds |
|---|---|---|---|
| neighbours | 0.0712 | | 0.0969 |
| `mf` | 0.0546 | 0.0730 to 0.0799 | 0.0906 |
| `lightgcn` | 0.0766 | 0.0791 to 0.0831 | 0.0825 |

One seed at a time, LightGCN is above matrix factorisation on the test links as it was on validation.
Averaging five fits raises matrix factorisation by 0.0150 over the mean of its single seeds and
LightGCN by 0.0015, and that reverses their order. Against neighbours, LightGCN was +0.0054
[-0.0034, +0.0145] on validation and is -0.0144 [-0.0232, -0.0064] on the test links.

Secondary measures of the protocol, as values. No paired interval was computed for these.

| | trait ROC area | recall in the top 100 | average precision by target | ROC area by target | all candidates at once |
|---|---|---|---|---|---|
| target degree | 0.6856 | 0.0613 | 0.0117 | 0.4994 | 0.0145 |
| neighbours | 0.6944 | 0.1007 | 0.0875 | 0.7930 | 0.0391 |
| `mf` | 0.6744 | 0.1031 | 0.0436 | 0.5041 | 0.0225 |
| `lightgcn` | 0.6597 | 0.1049 | 0.1065 | 0.6854 | 0.0721 |

Matrix factorisation orders the traits of a target no better than chance (ROC area 0.5041). Both models
are trained only to order targets within a trait. LightGCN is the highest of the four by target and
over all candidates at once, measures on which it was not tuned and for which no interval is given.

By year of the later link (2019, 2020, 2021, 2022): neighbours 0.0696, 0.0916, 0.0712, 0.0703; `mf`
0.0702, 0.0789, 0.0631, 0.0584; `lightgcn` 0.0608, 0.0762, 0.0590, 0.0493.

### What follows

- These two models are scored and are not scored again.
- The best scorer on all links and on new links is still the neighbours baseline, which has no fitted
  parameter. A graph model on Pharos alone does not beat it.
- The two models that join the context on (`lightgcn_ctx`, `lightgcn_ctx_pruned`) are the remaining
  declared ones. LightGCN's 0.0825 is what they are compared with, and neighbours' 0.0969 is what a
  claim of improvement has to beat.
- Averaging over seeds mattered more than the choice of model. The protocol fixed five seeds before
  either model was scored and that stands for the remaining two.

## Added 9 Oct 2026: models with lineage and position as inputs

Six further models, added after the two above were scored, use trait lineage and gene position as
inputs. They are fitted and frozen on validation and none is scored on the test links yet. They are
in `findings/INPUTS.md`.
