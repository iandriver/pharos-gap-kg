# Trait lineage and gene position as model inputs

9 Oct 2026. The models, the forms, the fitting and what counts were committed to `docs/PROTOCOL.md`
("Lineage and position as inputs") before anything was fitted. Every model here was **added after
`mf` and `lightgcn` had been scored on the test links**, and the decision to use these two inputs was
taken after seeing the controls, which are computed on the test links.

**No test link has been read for any model in this file.** Everything below is on the inner split of
the training links: fitted on the 57,411 links first reported in or before 2016, validated on the
22,610 links of 2017 and 2018 over 284 traits.

(Added 10 Oct 2026: the six models have since been scored on the test links, once each. That is the
last section. The sections before it are as written on 9 Oct, when the statement above held.)

## The two inputs on the validation links

Two traits are in one lineage if one is an ancestor of the other in the Experimental Factor Ontology
v3.94.0. Two targets are neighbours if genes of theirs lie within 500,000 bases on one chromosome
(Ensembl 113). On the inner split there are 855 lineage pairs among the 850 traits and 76,011
neighbour pairs among the 10,648 targets.

| | candidates | validation positives among them | share of the 22,610 positives |
|---|---|---|---|
| every candidate of the 284 traits | 2,994,372 | 22,610 (0.76%) | |
| lineage-linked (B4 above zero) | 35,854 | 1,270 (3.5%) | 5.6% |
| position-linked (B5 above zero) | 58,546 | 2,367 (4.0%) | 10.5% |
| either | | 3,477 | 15.4% |

A linked candidate is about five times as likely to be a later link as a candidate at random. Alone,
the two scores reach 0.0331 (B4 lineage) and 0.0372 (B5 position) on validation, beside 0.0394 for
target degree and 0.0712 for the neighbours baseline B3.

## The combiners as first declared

A combiner scores a candidate by a weighted sum of parameter-free scores. C0 combines target degree
(B1) and neighbours (B3) and uses neither input; C1 adds position (B5), C2 adds lineage (B4), C3 adds
both. As first declared the weights maximised a likelihood: the probability of a trait's validation
positives under a softmax over its candidates (code at commit `86b57d5`,
`docs/data/combine_gene_trait_2018_first.json`). Each figure is cross-validated: a trait is scored by
weights fitted on the other four fifths of the traits.

The form was chosen for C3 among three: F1 0.0730, F2 0.0774, F3 0.0715. Under F2, each score enters
as its percentile rank among the trait's candidates.

| | scores | cross-validated | against neighbours, 0.0712 |
|---|---|---|---|
| C0 `combine` | B1, B3 | 0.0511 | -0.0201 [-0.0345, -0.0082] |
| C1 `combine_pos` | B1, B3, B5 | 0.0717 | +0.0005 [-0.0139, +0.0133] |
| C2 `combine_lin` | B1, B3, B4 | 0.0556 | -0.0156 [-0.0318, -0.0016] |
| C3 `combine_lin_pos` | B1, B3, B4, B5 | 0.0774 | +0.0062 [-0.0086, +0.0207] |

- **Fitted this way, a combiner that holds the neighbours score is worse than the neighbours score.**
  A weight of zero on degree would make C0 rank exactly as B3 does. The fitted weights put it 0.0201
  below. The likelihood counts a positive the same wherever it ranks, and average precision counts the
  top of the list most.
- **The comparison this was built for would have misled.** C3 is 0.0264 above C0 [+0.0144, +0.0405],
  which reads as a large gain from the two inputs, while C3 is 0.0062 above B3 with an interval that
  spans zero.

The protocol was amended on this evidence, before any model of this section was scored on the test
links. This fit stays on record as a validation result and is not scored on the test links.

## The combiners as amended

The weights now maximise the primary measure itself on the validation positives. The weight on B3 is
fixed at 1 and the others start at 0, so the search starts at B3 alone. It moves one weight at a time
over 23 fixed values (0, and plus or minus each power of two from 1/64 to 16), and nothing in it is
random (code at commit `74fe460`, `docs/data/combine_gene_trait_2018.json`, frozen at `4e5859a`).

The form chosen for C3: F1 0.0757, F2 0.0751, **F3 0.0810**. Under F3 a score enters as log(1 + score),
and again multiplied by log(1 + the number of links the trait has), so a weight can differ between
small and large traits.

| | scores | cross-validated | against neighbours, 0.0712 | against C0 | fitted on all 284 traits |
|---|---|---|---|---|---|
| C0 `combine` | B1, B3 | 0.0732 | +0.0020 [-0.0030, +0.0069] | | 0.0759 |
| C1 `combine_pos` | B1, B3, B5 | 0.0813 | +0.0101 [+0.0036, +0.0173] | +0.0081 [+0.0006, +0.0160] | 0.0853 |
| C2 `combine_lin` | B1, B3, B4 | 0.0757 | +0.0045 [-0.0023, +0.0116] | +0.0025 [-0.0029, +0.0092] | 0.0824 |
| C3 `combine_lin_pos` | B1, B3, B4, B5 | 0.0810 | +0.0098 [+0.0019, +0.0177] | +0.0077 [-0.0005, +0.0165] | 0.0901 |

- **Position adds to the neighbours score on validation.** C1 is 0.0101 above B3 and 0.0081 above C0,
  both intervals above zero. Its weight on log(1 + B5) is 0.5, beside 1 on log(1 + B3).
- **Lineage is not shown to add.** C2 is 0.0025 above C0 with an interval that spans zero, and C3
  (0.0810) is level with C1 (0.0813). The half-width of C2 against C0, 0.0060, is the smallest gain
  from lineage this comparison would have resolved.
- **Fitting weights to degree and neighbours alone does not beat neighbours.** C0 is +0.0020
  [-0.0030, +0.0069], so the gains of C1 and C3 come from position and not from the fitting.
- **Each fit flatters itself, more with more weights.** Scored on the traits it was fitted on, a
  combiner reads 0.0027 (C0, 3 free weights), 0.0039 (C1, 5), 0.0067 (C2, 5) and 0.0092 (C3, 7) above
  its cross-validated figure. For C3 that gap is as large as its margin over B3.
- **The chosen form's figure is a little high too**, because F3 was picked on C3's cross-validated
  score. F1 and F2 are also above B3's 0.0712; no interval was computed for them.

The fitted weights. A term not listed has weight 0. "size" is log(1 + the trait's links).

| | weights |
|---|---|
| C0 | B3 1; B1 1/64; B1 x size -1/64; B3 x size -1/8 |
| C1 | B3 1; B1 1/64; B5 1/2; B3 x size 1/8; B5 x size 1/64 |
| C2 | B3 1; B1 1/32; B4 1/16; B1 x size -1/64; B3 x size -1/8; B4 x size 1/32 |
| C3 | B3 1; B1 1/32; B4 -1/2; B5 1/8; B1 x size -1/64; B3 x size -1/8; B4 x size 1/4 |

The sign of the weight on B3 x size differs between C1 and the other three, and C3 sets B5 x size to
0 where C1 uses it. Weights that move this much between neighbouring models are not individually
meaningful; the cross-validated score is what the fit is judged on.

## The graph models

`lightgcn_pos` is `lightgcn`'s graph with an edge between every two targets that are neighbours, and
`lightgcn_lin_pos` adds an edge between every two traits in one lineage. Both were tuned over the same
44 settings as `lightgcn`, seed 0, on the inner split (`docs/data/model_<name>_2018.json` and
`_validation.json`). The fitting code is the code that fitted `lightgcn`.

| | edges joined on | chosen setting | validation | against neighbours, 0.0712 |
|---|---|---|---|---|
| `lightgcn`, for reference | none | vectors of 256, 5 steps, rate 0.005, weight 1e-4, 60 passes | 0.0766 | +0.0054 [-0.0034, +0.0145] |
| `lightgcn_pos` | 76,011 neighbour edges | vectors of 128, 4 steps, rate 0.001, weight 1e-6, 250 passes | 0.0777 | +0.0065 [-0.0049, +0.0183] |
| `lightgcn_lin_pos` | 76,011 neighbour and 855 lineage edges | vectors of 64, 2 steps, rate 0.005, weight 1e-4, 110 passes | 0.0792 | +0.0080 [-0.0050, +0.0199] |

The last column is the paired difference over the 284 traits. It flatters each model, because its
setting was chosen on these links.

- **Neither graph model can be told from the neighbours baseline on validation.** Both intervals span
  zero. Their half-widths, 0.0116 and 0.0124, are the smallest differences these comparisons resolve.
- **The edges raise the best validation score by 0.0011 and 0.0026 over `lightgcn`.** No paired
  interval between the graph models was computed on validation.
- **With the edges, smaller settings do best.** For `lightgcn` the score rose with vector length and
  steps up to the edge of the settings (0.0750 to 0.0766 at vectors of 256 and 5 steps). With neighbour
  edges that corner is the lowest of the settings with three or more steps (0.0707 to 0.0748), and
  neither chosen setting is at an edge.
- **One step of smoothing gains most.** With one step `lightgcn` scores 0.0644 to 0.0670,
  `lightgcn_pos` 0.0689 to 0.0737 and `lightgcn_lin_pos` 0.0699 to 0.0764.
- **On validation the combiner with position is above both graph models**: 0.0813 cross-validated,
  against 0.0777 and 0.0792 for the best of 44 settings chosen on the same links. That is two figures
  side by side and not a paired test. The protocol's comparison is on the test links.

Fitting the 44 settings took 3.3 and 3.4 hours of processor time for the two models, on a laptop CPU.

## What is waiting

The six models are frozen: the four combiners' weights (`docs/data/combine_gene_trait_2018_frozen.json`)
and the two graph models' settings (`docs/data/model_lightgcn_pos_2018_frozen.json`,
`docs/data/model_lightgcn_lin_pos_2018_frozen.json`). None has been scored on the test links. Under
the protocol the four combiners are scored together once, and each graph model once under five seeds.

What the test can and cannot show is fixed in the protocol:

- A model helps if it beats the neighbours baseline (0.0969 on the test links) with the paired
  interval excluding zero.
- A gain from a combiner that uses lineage or position is a gain on all links and on loci. With every
  linked candidate set aside the two inputs are zero for all that remain, so such a gain says nothing
  about new links.
- A graph model is kept over a combiner only if it beats the combiner with the same inputs.
- Validation has not predicted the test comparison before: `lightgcn` was +0.0054 against neighbours
  on validation and -0.0144 on the test links (`findings/MODELS.md`).

## Added 10 Oct 2026: all six scored on the test links, once

Scored on Ian's word. `python -m pharos_gap_kg.combine final pharos319 --cutoff 2018` scored the four
combiners together from the weights frozen at `4e5859a`, and `python -m pharos_gap_kg.model final`
scored each graph model from the setting frozen at `2f48edb`, fitted on all 99,435 training links
under seeds 0 to 4, its score for a pair being the mean of the five. Each scoring wrote a marker before
it read a test link, and a second start is refused. The `score` and `versus` stages write
`docs/data/combine_gene_trait_2018_scored.json`, `docs/data/model_<name>_2018_scored.json` and
`docs/data/versus_gene_trait_2018.json`. Each recomputes the two baselines and refuses unless they
equal the committed controls.

The primary measure, mean average precision over the 540 scored traits, on the four views of the
protocol:

| | all links | new links, lower bound | new links, neutral | by locus |
|---|---|---|---|---|
| target degree | 0.0733 | 0.0496 | 0.0556 | 0.0255 |
| neighbours | 0.0969 | 0.0539 | 0.0662 | 0.0512 |
| `mf`, scored 8 Oct | 0.0906 | 0.0465 | 0.0564 | 0.0662 |
| `lightgcn`, scored 8 Oct | 0.0825 | 0.0433 | 0.0530 | 0.0527 |
| `combine` | 0.0931 | 0.0504 | 0.0621 | 0.0524 |
| `combine_pos` | 0.1037 | 0.0440 | 0.0687 | 0.0553 |
| `combine_lin` | 0.1010 | 0.0494 | 0.0665 | 0.0554 |
| `combine_lin_pos` | 0.0995 | 0.0390 | 0.0665 | 0.0547 |
| `lightgcn_pos` | 0.0884 | 0.0373 | 0.0549 | 0.0557 |
| `lightgcn_lin_pos` | 0.0962 | 0.0379 | 0.0537 | 0.0661 |

Each model against the neighbours baseline, paired over traits, with 95% intervals:

| | all links | new links, lower bound | new links, neutral | by locus |
|---|---|---|---|---|
| `combine` | -0.0039 [-0.0070, -0.0010] | -0.0036 [-0.0057, -0.0014] | -0.0041 [-0.0068, -0.0015] | +0.0012 [-0.0016, +0.0035] |
| `combine_pos` | +0.0067 [-0.0009, +0.0151] | -0.0100 [-0.0138, -0.0066] | +0.0025 [+0.0009, +0.0046] | +0.0041 [-0.0030, +0.0118] |
| `combine_lin` | +0.0041 [-0.0000, +0.0093] | -0.0045 [-0.0071, -0.0019] | +0.0002 [-0.0022, +0.0027] | +0.0042 [+0.0006, +0.0092] |
| `combine_lin_pos` | +0.0026 [-0.0056, +0.0114] | -0.0149 [-0.0196, -0.0108] | +0.0002 [-0.0022, +0.0027] | +0.0035 [-0.0034, +0.0114] |
| `lightgcn_pos` | -0.0085 [-0.0174, -0.0002] | -0.0166 [-0.0231, -0.0110] | -0.0114 [-0.0190, -0.0047] | +0.0045 [-0.0022, +0.0111] |
| `lightgcn_lin_pos` | -0.0008 [-0.0108, +0.0092] | -0.0160 [-0.0231, -0.0098] | -0.0126 [-0.0212, -0.0048] | +0.0149 [+0.0069, +0.0238] |

- **None of the six is shown to beat the neighbours baseline.** The protocol's rule is a paired
  interval that excludes zero on the primary measure. The highest model, `combine_pos`, scores 0.1037
  against 0.0969, and its interval reaches below zero. The interval's half-width, 0.0080, is the
  smallest gain this test would have resolved for it. `combine_lin`'s interval ends just below zero.
- **Two are worse than neighbours on all links**: the reference combiner `combine`, which uses
  neither input, and `lightgcn_pos`. `lightgcn_lin_pos` cannot be told from neighbours.

What each input adds to the model it is put into:

| what is added | comparison | all links | by locus |
|---|---|---|---|
| position, in the combiner | `combine_pos` against `combine` | +0.0106 [+0.0031, +0.0189] | +0.0030 [-0.0034, +0.0105] |
| lineage, in the combiner | `combine_lin` against `combine` | +0.0080 [+0.0042, +0.0129] | +0.0030 [-0.0006, +0.0081] |
| both, in the combiner | `combine_lin_pos` against `combine` | +0.0065 [-0.0009, +0.0151] | +0.0023 [-0.0037, +0.0096] |
| position, in the graph model | `lightgcn_pos` against `lightgcn` | +0.0059 [+0.0006, +0.0116] | +0.0031 [-0.0020, +0.0088] |
| lineage, in the graph model with position | `lightgcn_lin_pos` against `lightgcn_pos` | +0.0078 [+0.0021, +0.0144] | +0.0104 [+0.0051, +0.0160] |

- **Each input adds on all links, in both kinds of model.** Position adds 0.0106 in the combiner and
  0.0059 in the graph model. Lineage adds 0.0080 in the combiner and 0.0078 in the graph model that
  already has position. All four intervals exclude zero.
- **The combiner with both inputs is not above the reference**, and it is below each combiner with
  one. It has seven free weights, and on validation it read 0.0092 higher on the traits it was fitted
  on than in cross-validation.
- **Lineage was not shown to add on validation and does add on the test links.** Two things differ.
  The test split has 1,630 lineage pairs among 1,404 traits where the inner split had 855 among 850.
  And the ontology is the 2026 release, so a relation between two traits may have been recorded after
  2018. The protocol's caveat for any result that uses lineage applies to these two rows.
- **The gain stops short of the baseline because fitting costs about what the inputs give.** The
  reference combiner is 0.0039 below neighbours, and `lightgcn` is 0.0144 below.
- **By locus, the graph model with both inputs is ahead of neighbours**, 0.0661 against 0.0512, and
  the lineage edges are worth 0.0104 of that over `lightgcn_pos`. `combine_lin` is ahead by locus too.
  Matrix factorisation scored 0.0662 on that view on 8 Oct. The protocol asks for the primary measure,
  new links and loci together before a claim, so this is reported and not claimed.
- **Nothing here is a result about new links.** With the linked positives removed, all five models
  that use an input are below neighbours, which is what the protocol said of the combiners. With every
  linked candidate removed, both graph models are below neighbours. `combine_lin` and
  `combine_lin_pos` score the same, 0.0665, because their weights on degree and neighbours are the
  same and both inputs are zero there. `combine_pos` is 0.0025 above neighbours on that view. Position
  is zero for every candidate in it, so that comes from its weights on degree and neighbours, and its
  lower-bound view is 0.0100 below, so the two-way requirement is not met.
- **Neither graph model is kept over a combiner.** `lightgcn_pos` against `combine_pos` is
  -0.0153 [-0.0229, -0.0087]. `lightgcn_lin_pos` against `combine_lin_pos` is
  -0.0033 [-0.0128, +0.0061], a tie, which the protocol gives to the combiner. By locus
  `lightgcn_lin_pos` is ahead of `combine_lin_pos`, +0.0115 [+0.0033, +0.0192].

**Validation against the test links.** Each model against neighbours on the primary measure:

| | validation | test links |
|---|---|---|
| `combine` | +0.0020 [-0.0030, +0.0069] | -0.0039 [-0.0070, -0.0010] |
| `combine_pos` | +0.0101 [+0.0036, +0.0173] | +0.0067 [-0.0009, +0.0151] |
| `combine_lin` | +0.0045 [-0.0023, +0.0116] | +0.0041 [-0.0000, +0.0093] |
| `combine_lin_pos` | +0.0098 [+0.0019, +0.0177] | +0.0026 [-0.0056, +0.0114] |
| `lightgcn_pos` | +0.0065 [-0.0049, +0.0183] | -0.0085 [-0.0174, -0.0002] |
| `lightgcn_lin_pos` | +0.0080 [-0.0050, +0.0199] | -0.0008 [-0.0108, +0.0092] |

- **Every model stands lower against neighbours on the test links than on validation**: the combiners
  by 0.0004 to 0.0072, the graph models by 0.0150 and 0.0088.
- **Validation ordered the models that share inputs correctly.** It put `combine_pos` above both
  graph models and `lightgcn_lin_pos` above `lightgcn_pos`, and both held. Among the combiners it
  did not: on validation the combiner with both inputs was level with `combine_pos` (0.0810 and
  0.0813) and above `combine_lin` (0.0757), and on the test links it is below both.
- **Averaging seeds mattered for one graph model.** `lightgcn_pos` scores 0.0867 to 0.0884 one seed
  at a time and 0.0884 for the mean of five. `lightgcn_lin_pos` scores 0.0893 to 0.0911 one seed at a
  time and 0.0962 for the mean.

Secondary measures of the protocol, as values. No paired interval was computed for these.

| | trait ROC area | recall in the top 100 | average precision by target | ROC area by target | all candidates at once |
|---|---|---|---|---|---|
| neighbours | 0.6944 | 0.1007 | 0.0875 | 0.7930 | 0.0391 |
| `combine` | 0.6446 | 0.1062 | 0.0116 | 0.3568 | 0.0082 |
| `combine_pos` | 0.7215 | 0.1172 | 0.1002 | 0.8034 | 0.0529 |
| `combine_lin` | 0.6972 | 0.1078 | 0.0212 | 0.3764 | 0.0103 |
| `combine_lin_pos` | 0.6856 | 0.1141 | 0.0463 | 0.4430 | 0.0188 |
| `lightgcn_pos` | 0.6775 | 0.1164 | 0.0990 | 0.6784 | 0.0675 |
| `lightgcn_lin_pos` | 0.6914 | 0.1153 | 0.0855 | 0.6405 | 0.0512 |

The three combiners with a negative weight on neighbours times trait size order the traits of one
target worse than chance (ROC area by target 0.3568, 0.3764 and 0.4430): their scores are fitted to
order targets within a trait and are not comparable across traits. `combine_pos`, whose weight on that
term is positive, is the highest of all seven rows on trait ROC area, recall in the top 100 and ROC
area by target.

By year of the later link (2019, 2020, 2021, 2022): neighbours 0.0696, 0.0916, 0.0712, 0.0703; `combine` 0.0683, 0.0843, 0.0699, 0.0630; `combine_pos` 0.0712, 0.0967, 0.0666, 0.0949; `combine_lin` 0.0739, 0.0892, 0.0804, 0.0689; `combine_lin_pos` 0.0690, 0.0905, 0.0731, 0.0834; `lightgcn_pos` 0.0680, 0.0808, 0.0605, 0.0693; `lightgcn_lin_pos` 0.0767, 0.0792, 0.0641, 0.0888.

### What follows

- The neighbours baseline stays the scorer behind the typed questions. No fitted model has been shown
  to beat it on this split.
- Gene position and trait lineage each carry signal a link-prediction score counts: four paired
  comparisons exclude zero. All of it is on links that follow from an earlier link at the same locus
  or under a related trait name.
- The combiners are closed. The amendment of 9 Oct said they would be fitted once more and nothing
  further tried. A simpler way of adding position to neighbours, such as one fixed weight with no
  trait-size terms, would now be chosen knowing these test results and would have to be reported so.
- The two models that join the context tables on (`lightgcn_ctx`, `lightgcn_ctx_pruned`) are still
  not tuned.
