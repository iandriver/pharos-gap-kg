# The prospective test: links first reported from 2023

10 Oct 2026. The test, its link rule and what is scored on it were committed to `docs/PROTOCOL.md`
("A prospective test", and three dated additions of 10 Oct) before the steps they govern.

## The split

**The source.** The GWAS Catalog's release of 5 Oct 2026, registered and pinned by digest:
1,198,394 associations from 76,336 studies in 6,905 papers, published 2005 to 2026
(`python -I -m pharos_gap_kg.gwascat extract`). Pharos 3.19 holds GWAS links to 2022.

**The rule declared first failed.** Four rules were declared, built on the Catalog's own gene columns
and on trait ids as written. On the 12,260 studies that Pharos and the Catalog both hold, all published
by 2022, against Pharos's 720,702 (gene, trait, study) triples (`gwascat first-rule`,
`docs/data/gwascat_rule_first_pharos319.json`):

| rule | triples produced | of them Pharos's | of Pharos's, recovered |
|---|---|---|---|
| genes the variant lies in, p at most 5e-8 | 119,752 | 68.1% | 11.3% |
| genes the variant lies in, any p | 135,726 | 60.3% | 11.4% |
| those and the nearest gene each side, p at most 5e-8 | 175,035 | 67.1% | 16.3% |
| those and the nearest gene each side, any p | 199,731 | 59.0% | 16.4% |

The protocol had said that under 90% the test would go ahead and no further rule would be tried. That
sentence was set aside in a dated amendment, with this evidence (`gwascat look`,
`docs/data/gwascat_look_pharos319.json`), all of it from the same studies:

- **The Catalog has replaced half of its trait terms since Pharos's snapshot.** The terms differ for
  6,433 of the 12,260 studies, mostly an EFO term replaced by an OBA or MONDO term. As written, 50.0%
  of the Catalog's (trait, study) pairs are Pharos's. Through the ontology's live term, an obsolete
  term standing for its replacement, 97.1% are, and 98.7% of Pharos's are recovered.
- **Pharos's source took every gene near the variant, not the nearest.** Pharos lists 3 times as
  many genes for a study, at the median, as the Catalog's columns give. Taking every gene within a
  distance of a variant significant at 5e-8:

| distance from the variant | of the rule's (gene, study) pairs, Pharos's | of Pharos's, recovered |
|---|---|---|
| 0 bases | 98.5% | 16.6% |
| 10,000 bases | 98.4% | 24.8% |
| 25,000 bases | 98.4% | 35.2% |
| 50,000 bases | 98.3% | 50.7% |
| 100,000 bases | 98.1% | 79.3% |
| 250,000 bases | 50.9% | 81.7% |

  What the rule produces stays Pharos's up to 100,000 bases and stops being so beyond. Without the
  p-value filter, 87.7% of the pairs at 100,000 bases are Pharos's.

**The amended rule.** A link joins a trait, as the live term of its ontology id, to every target with
a gene within 100,000 bases of a variant associated with it at 5e-8 or smaller. At the level a link is
made, (target, live trait, study), on the shared studies (`gwascat rule`,
`docs/data/gwascat_rule_pharos319.json`): 95.1% of the rule's 636,104 triples are Pharos's, and 92.8% of
Pharos's 651,671 are recovered. For pairs it is 95.3% and 92.0%. Both are over the 90% asked for.
The amendment had expected about 77% recovered, a figure taken gene by gene that left out targets
Pharos ties only to gene ids on alternate assemblies; a dated note corrects it.

**The split** (`gwascat build`, `docs/data/split_gene_trait_2022.json`), frozen by code that reports
counts and runs no scorer:

| | pairs | |
|---|---|---|
| training | 294,954 | Pharos 3.19 by live trait term: 18,005 targets, 2,225 traits. As written they were 295,327 pairs over 2,247 traits |
| prospective test | 194,017 | first reported in 2023 or later, between a training target and a training trait: 1,110 traits, 17,456 targets |
| cold start | 394,383 | 377,156 to a new trait (8,425 new traits), 6,729 to a new target (1,529), 10,498 to both |
| set aside | 228,538 | pairs Pharos lacks that were first reported by 2022; neither positives nor candidates. 139,482 are between a training target and trait |

The test links by first year: 2023 34,067, 2024 58,594, 2025 85,222, 2026 16,134. The Catalog gives 93.3% of Pharos's
training pairs back under the rule.

- **The set-aside pairs are many.** They are pairs from studies published by 2022 that Pharos does not
  hold: the Catalog curated the study after Pharos's snapshot, or Pharos's source left the study out.
  Which, this protocol cannot see, and that is why they are counted neither way.
- **The last year is partial.** 2026 holds what had been curated by 5 Oct.

## What was scored

Fixed in a dated addition to the protocol before any scorer read a prospective link: the
parameter-free scorers; the 2018 split's four combiners and two of its models exactly as frozen; and
three candidates from the rehearsal, tuned or fitted with a cutoff of 2020 on training links only
(`prospect tune`, `prospect weights`, `prospect fit`; `docs/data/prospect_candidate_2022.json`,
`docs/data/prospect_weights_2022.json`). The fitted models are the mean of five seeds on all 294,954
training links. The addition also wrote down what the rehearsal led one to expect: no graph model
above neighbours, `combine[B1+B3+B5]` a little above it, and lineage transfer well above degree for
new traits.

`python -m pharos_gap_kg.prospect score pharos319` writes `docs/data/prospect_gene_trait_2022.json`. The
candidates of a trait are the 18,005 training targets less its training links and less the pairs
set aside. 68,846 of the 194,017 test links (35.5%) are lineage- or position-linked to a training
link; on the 2018 split's test it was 20.2%.

## The test links

Mean average precision over the 1,110 scored traits:

| | | all links | new links, lower bound | new links, neutral | by locus |
|---|---|---|---|---|---|
| `B1` | target degree | 0.0851 | 0.0412 | 0.0548 | 0.0261 |
| `B3` | neighbours | 0.0947 | 0.0447 | 0.0582 | 0.0475 |
| `B4` | lineage | 0.0316 | 0.0071 | 0.0087 | 0.0238 |
| `B5` | position | 0.0439 | 0.0066 | 0.0087 | 0.0315 |
| `combine` | 2018 combiner, frozen | 0.0870 | 0.0406 | 0.0524 | 0.0475 |
| `combine_pos` | 2018 combiner, frozen | 0.0997 | 0.0344 | 0.0587 | 0.0506 |
| `combine_lin` | 2018 combiner, frozen | 0.0930 | 0.0381 | 0.0547 | 0.0480 |
| `combine_lin_pos` | 2018 combiner, frozen | 0.0919 | 0.0321 | 0.0547 | 0.0508 |
| `mf` | 2018 recipe, frozen | 0.0620 | 0.0271 | 0.0368 | 0.0486 |
| `lightgcn_lin_pos` | 2018 recipe, frozen | 0.0676 | 0.0241 | 0.0364 | 0.0519 |
| `combine[B1+B3+B5]:F1` | rehearsal candidate, carried | 0.1012 | 0.0387 | 0.0587 | 0.0517 |
| `combine[B1+B3+B4+B5]:F1` | rehearsal candidate, did not carry | 0.1031 | 0.0309 | 0.0590 | 0.0523 |
| `lightgcn_lin_pos:L05:steps` | rehearsal candidate, did not carry | 0.0913 | 0.0356 | 0.0523 | 0.0569 |

Each against the neighbours baseline, paired over traits, with 95% intervals:

| | all links | new links, lower bound | new links, neutral | by locus |
|---|---|---|---|---|
| `combine` | -0.0077 [-0.0102, -0.0049] | -0.0041 [-0.0058, -0.0026] | -0.0058 [-0.0080, -0.0039] | -0.0000 [-0.0018, +0.0024] |
| `combine_pos` | +0.0051 [+0.0001, +0.0098] | -0.0104 [-0.0132, -0.0081] | +0.0005 [+0.0002, +0.0007] | +0.0031 [-0.0014, +0.0079] |
| `combine_lin` | -0.0016 [-0.0042, +0.0010] | -0.0066 [-0.0084, -0.0050] | -0.0035 [-0.0054, -0.0018] | +0.0005 [-0.0010, +0.0020] |
| `combine_lin_pos` | -0.0028 [-0.0063, +0.0011] | -0.0126 [-0.0158, -0.0099] | -0.0035 [-0.0054, -0.0018] | +0.0033 [-0.0001, +0.0067] |
| `mf` | -0.0327 [-0.0390, -0.0266] | -0.0176 [-0.0225, -0.0130] | -0.0214 [-0.0272, -0.0159] | +0.0011 [-0.0040, +0.0059] |
| `lightgcn_lin_pos` | -0.0271 [-0.0339, -0.0205] | -0.0206 [-0.0257, -0.0160] | -0.0218 [-0.0277, -0.0161] | +0.0044 [-0.0012, +0.0096] |
| `combine[B1+B3+B5]:F1` | +0.0066 [+0.0031, +0.0102] | -0.0060 [-0.0081, -0.0043] | +0.0005 [-0.0006, +0.0012] | +0.0042 [+0.0008, +0.0078] |
| `combine[B1+B3+B4+B5]:F1` | +0.0084 [+0.0034, +0.0135] | -0.0138 [-0.0170, -0.0111] | +0.0008 [-0.0007, +0.0020] | +0.0047 [-0.0000, +0.0097] |
| `lightgcn_lin_pos:L05:steps` | -0.0033 [-0.0092, +0.0024] | -0.0091 [-0.0130, -0.0053] | -0.0059 [-0.0112, -0.0006] | +0.0094 [+0.0050, +0.0139] |

- **Three combiners that use position beat the neighbours baseline.** `combine[B1+B3+B5]:F1`, the arm
  that carried in rehearsal, is +0.0066 [+0.0031, +0.0102], higher on 628 of the 1,110 traits.
  `combine[B1+B3+B4+B5]:F1` is +0.0084 [+0.0034, +0.0135]. The 2018 split's `combine_pos`, with weights
  untouched since 9 Oct, is +0.0051 [+0.0001, +0.0098]. By the protocol's rule each helps. These are the
  first models here to beat the baseline, on links no model had been scored on, and the first of
  them was named in advance.
- **The gain is 0.0066 to 0.0084 on a baseline of 0.0947**, seven to nine parts in a hundred. The weights are
  few: B1 1/32; B3 1; B5 1/4 for the first, and B1 1/16; B3 1; B4 1/4; B5 1/2 for the second.
- **It is not a gain on new links.** With the linked positives removed, every model that uses an input
  is below neighbours. With every linked candidate removed the two candidates are level with it
  (+0.0005 [-0.0006, +0.0012] and +0.0008 [-0.0007, +0.0020]). The gain is on links at a locus, or under
  a trait name, that training already holds, which is what the protocol said such a combiner can show.
- **Counted by locus `combine[B1+B3+B5]:F1` is still ahead**, +0.0042 [+0.0008, +0.0078].
- **No graph model beats neighbours.** The candidate tuned at 2020, `lightgcn_lin_pos:L05:steps`, is
  -0.0033 [-0.0092, +0.0024]: not told apart. It is below on new links and above by locus
  (+0.0094 [+0.0050, +0.0139]), as `lightgcn_lin_pos` was on the 2018 split.
- **The two recipes frozen on the 2018 split fail at three times the links.** `mf` is -0.0327 [-0.0390, -0.0266]
  and `lightgcn_lin_pos` -0.0271 [-0.0339, -0.0205]. Their settings and passes were chosen on 57,411 links,
  scored on the 2018 split after a fit on 99,435, and are applied here to 294,954. The same graph model retuned at 2020 is +0.0237 [+0.0190, +0.0286] above the
  frozen recipe.

**Do the 2018 split's results repeat?**

| | the 2018 split's test | this split |
|---|---|---|
| `combine` against neighbours | -0.0039 [-0.0070, -0.0010] | -0.0077 [-0.0102, -0.0049] |
| `combine_pos` against neighbours | +0.0067 [-0.0009, +0.0151] | +0.0051 [+0.0001, +0.0098] |
| `combine_lin` against neighbours | +0.0041 [-0.0000, +0.0093] | -0.0016 [-0.0042, +0.0010] |
| `combine_lin_pos` against neighbours | +0.0026 [-0.0056, +0.0114] | -0.0028 [-0.0063, +0.0011] |
| `combine_pos` against `combine` | +0.0106 [+0.0031, +0.0189] | +0.0128 [+0.0080, +0.0175] |
| `combine_lin` against `combine` | +0.0080 [+0.0042, +0.0129] | +0.0061 [+0.0032, +0.0086] |
| `combine_lin_pos` against `combine` | +0.0065 [-0.0009, +0.0151] | +0.0049 [+0.0015, +0.0081] |
| `mf` against neighbours | -0.0063 [-0.0155, +0.0025] | -0.0327 [-0.0390, -0.0266] |
| `lightgcn_lin_pos` against neighbours | -0.0008 [-0.0108, +0.0092] | -0.0271 [-0.0339, -0.0205] |

- **For the combiners, yes.** The reference combiner is below neighbours on both. Position adds to it
  on both, and so does lineage. The one with position is above neighbours on both, with the interval
  reaching zero on the first and clearing it here. The two with lineage are not told apart from
  neighbours on either.
- **For the two fitted models, no**: both were near neighbours on the 2018 split and are far below
  here, for the reason above.

**Did the rehearsal say what would happen?**

| | in rehearsal | on the prospective links |
|---|---|---|
| `combine[B1+B3+B5]:F1` | +0.0038 [+0.0008, +0.0077] | +0.0066 [+0.0031, +0.0102] |
| `combine[B1+B3+B4+B5]:F1` | +0.0095 [-0.0049, +0.0233] | +0.0084 [+0.0034, +0.0135] |
| `lightgcn_lin_pos:L05:steps` | +0.0079 [-0.0044, +0.0195] | -0.0033 [-0.0092, +0.0024] |

- The arm that carried in rehearsal beats neighbours here. The combiner that did not carry had the
  same sign and size in rehearsal with an interval nearly three times as wide, over 284 traits against 1,110
  here. The graph arm cannot be told from neighbours in either.

Each fitted model one seed at a time and as the mean of five: `mf` 0.0504 to 0.0519 and 0.0620;
`lightgcn_lin_pos` 0.0616 to 0.0622 and 0.0676; `lightgcn_lin_pos:L05:steps` 0.0844 to 0.0892 and
0.0913.

Two secondary measures, as values with no interval:

| | trait ROC area | recall in the top 100 |
|---|---|---|
| `B1` | 0.7331 | 0.1054 |
| `B3` | 0.7502 | 0.1156 |
| `B4` | 0.5859 | 0.0495 |
| `B5` | 0.5473 | 0.0664 |
| `combine` | 0.7127 | 0.1106 |
| `combine_pos` | 0.7669 | 0.1310 |
| `combine_lin` | 0.7366 | 0.1180 |
| `combine_lin_pos` | 0.7368 | 0.1222 |
| `mf` | 0.6949 | 0.0831 |
| `lightgcn_lin_pos` | 0.7147 | 0.0927 |
| `combine[B1+B3+B5]:F1` | 0.7666 | 0.1275 |
| `combine[B1+B3+B4+B5]:F1` | 0.7723 | 0.1334 |
| `lightgcn_lin_pos:L05:steps` | 0.7336 | 0.1235 |

By first year of the link, as values under one seed with no interval:

| first year | links | traits | `B3` | `combine_pos` | `combine[B1+B3+B5]:F1` | `combine[B1+B3+B4+B5]:F1` | `lightgcn_lin_pos:L05:steps` |
|---|---|---|---|---|---|---|---|
| 2023 | 34,067 | 436 | 0.0581 | 0.0655 | 0.0661 | 0.0668 | 0.0595 |
| 2024 | 58,594 | 504 | 0.0729 | 0.0656 | 0.0704 | 0.0676 | 0.0711 |
| 2025 | 85,222 | 674 | 0.0878 | 0.1000 | 0.0983 | 0.1025 | 0.0820 |
| 2026 | 16,134 | 226 | 0.0408 | 0.0385 | 0.0395 | 0.0387 | 0.0366 |

The combiners are ahead of neighbours in 2023 and 2025 and behind in 2024 and in the partial 2026.
The result over all years does not hold in each year.

## Cold start on this split

`python -m pharos_gap_kg.prospect coldstart pharos319` writes `docs/data/prospect_coldstart_2022.json`:
the scorers of `findings/COLDSTART.md`, with the pairs set aside left out of the candidates.

| | cutoff 2016 | cutoff 2018 | this split |
|---|---|---|---|
| new traits scored, and those with a relative | 497, 346 | 829, 543 | 8,404, 7,879 |
| N1 target degree | 0.0497 | 0.0610 | 0.0639 |
| N4 lineage transfer, ties by degree | 0.0663 | 0.1019 | 0.0613 |
| N5 lineage neighbours, ties by degree | 0.0672 | 0.0853 | 0.0729 |
| N4 against N1, the declared test | +0.0166 [+0.0078, +0.0268] | +0.0408 [+0.0300, +0.0519] | -0.0026 [-0.0040, -0.0011] |
| N5 against N1 | +0.0176 [+0.0094, +0.0277] | +0.0242 [+0.0163, +0.0329] | +0.0090 [+0.0068, +0.0114] |
| N5 against N4 | +0.0009 [-0.0037, +0.0053] | -0.0166 [-0.0233, -0.0104] | +0.0116 [+0.0095, +0.0137] |
| traits with a link to a new target, and candidate targets | 210, 9,764 | 322, 5,928 | 714, 2,407 |
| P1 locus degree | 0.0049 | 0.0102 | 0.1177 |
| P3 position transfer, ties by locus degree | 0.0236 | 0.0250 | 0.1427 |
| P3 against P1, the declared test | +0.0187 [+0.0072, +0.0336] | +0.0148 [+0.0074, +0.0243] | +0.0250 [+0.0087, +0.0418] |
| by new target, Q3 against Q1 | +0.0261 [+0.0219, +0.0305] | +0.0162 [+0.0130, +0.0197] | +0.0249 [+0.0182, +0.0318] |

- **Position transfer for new targets repeats.** P3 against P1 is +0.0250 [+0.0087, +0.0418] over 714 traits,
  the third period of three in which it is above zero. Read by new target it repeats too.
- **Lineage transfer for new traits does not repeat.** The declared test, N4 against N1, is
  -0.0026 [-0.0040, -0.0011] over 8,404 new traits: a little below degree, where the 2018 split gave
  +0.0408 [+0.0300, +0.0519]. The expectation written down before scoring was wrong.
- **The new traits here are of another kind.** 7,879 of the 8,404 have a relative in training, 94% against
  66% on the 2018 split, and 167,952 of their 353,045 links (48%) are to a target linked to a relative,
  against 24%. A relative that nearly every new trait has, and that is linked to half of what the new
  trait will be linked to, is a broad trait with very many links, and counting links to it sorts
  targets less well than degree does.
- **The scorer that passes on the relatives' neighbours scores is above degree in all three periods.**
  N5 against N1 is +0.0176 [+0.0094, +0.0277], +0.0242 [+0.0163, +0.0329] and +0.0090 [+0.0068, +0.0114]. N5 was a
  declared scorer and not the declared test, and on the 2018 split it was the weaker of the two, so
  this is a pattern across periods and not a passed test.

## What follows

- A combiner of neighbours, degree and position, with two fitted weights, is the first model here
  shown to beat the neighbours baseline on untouched links. Its gain is small and is on links at loci
  that training already holds.
- Fitted graph models remain unshown: level with neighbours when retuned, far below when a recipe is
  carried to a training set three times the size.
- For a target with no link, position transfer has held in three periods. For a trait with no link,
  lineage transfer depends on the kind of trait: it held where a new trait had a specific relative
  and failed where nearly every new trait hangs under one broad term.
- The question layer still ranks by neighbours alone. Whether to serve the combiner, and the
  cold-start rankings, is a decision not taken here.

## Added 11 Oct 2026: the question layer ranks by the combiner

On Ian's decision the scorer behind `predict` and `why` changed from the neighbours score to
`combine[B1+B3+B5]:F1`, the model that beat it above. The change, the calibration and a rule for the direction
the combiner had not been tested in were written into `docs/PROTOCOL.md` first.

**What is served.** For a trait, a target's score is log(1 + neighbours) plus 1/32 of log(1 + the
traits it is already linked to) plus 1/4 of log(1 + the genes within 500,000 bases that already have
the link), from the links as of the year asked. The weights are the ones scored above. They travel in
the bundle's `backtest.json`, and a bundle that names no scorer is ranked by neighbours alone as
before. Each candidate shows the three parts of its score.

**How often a rank held** (`python -m pharos_gap_kg.prospect calibrate pharos319`,
`docs/data/calibration_prospect_2022.json`). The share of candidates at each rank that were reported in
2023 to 2026, pooled over the 1,110 scored traits, against a base rate of 0.49%. The stage
refuses unless the scorer it recomputes scores what the test recorded.

| ranks | all candidates | an earlier link explains | none does | a disease, and none does | neighbours alone, all candidates |
|---|---|---|---|---|---|
| 1-10 | 17.8% | 20.5% | 9.7% | 9.1% | 17.3% |
| 11-50 | 13.8% | 17.5% | 9.0% | 7.0% | 13.2% |
| 51-100 | 11.1% | 14.6% | 7.9% | 6.4% | 10.5% |
| 101-500 | 5.5% | 9.2% | 3.6% | 2.9% | 5.1% |

Ranking the traits for one target, pooled over 17,456 targets:

| ranks | all candidates | an earlier link explains | none does | neighbours alone, all candidates |
|---|---|---|---|---|
| 1-5 | 7.8% | 8.0% | 6.0% | 8.1% |
| 6-20 | 5.6% | 6.1% | 4.4% | 5.4% |
| 21-50 | 3.7% | 5.2% | 2.9% | 3.4% |

- **A top-ten candidate for a trait was reported within four years 18 times in a hundred**, 36 times the
  base rate. Where no earlier link explains the candidate it is 10 in a hundred.
- **Rank by rank the two rankings are close.** In the top ten for a trait the combiner's share is 17.8%
  and that of neighbours alone 17.3%. Ranking the traits for one target, neighbours alone is a little
  higher in the top five, 8.1% against 7.8%, and lower below. These are values with no interval. The
  difference that was tested is the one in mean average precision above.
- **Every share is a lower bound.** A candidate not reported by 2026 is not known to be absent, and
  the answers say so wherever they quote one.
- **The front end quotes the share for the candidate's own group**: explained by an earlier link, not
  explained, or not explained and a disease.

**The direction that had not been tested.** The combiner was fitted and scored for ranking targets
within a trait. Ranking the traits for one target, it scores 0.0629 against 0.0615 for neighbours
alone, +0.0014 [+0.0006, +0.0022] over 17,456 targets. By the rule written beforehand it is served
that way too. That comparison was made on the prospective links after the test, so it is not a
declared result, and the answers say so.

**The served scores are the tested scores** (`tests/test_served_scorer.py`). The prospective split
keys a trait by its live ontology term and the bundle keeps Pharos's ids, so 38 of the bundle's 2,247
traits are one of two ids for the same trait there. For a trait that was not merged and a target that
holds no merged trait, the score `ask` computes from the bundle equals the score in the tested matrix;
the test checks 40 traits. For the 3,739 of 18,005 linked targets that hold a merged trait the served
score can differ a little from the tested one.

**What a reader will notice.**

- The highest candidates for a trait are now mostly genes beside a gene that already has the link,
  because that is where the gain over neighbours was. Every candidate says what explains it, and
  `--unexplained-only` lists the others, which the combiner orders about as neighbours did.
- A target with no link as of the year is not ranked, as the scorer was tested only on targets that
  have one. `why` still names a linked gene beside it.
- A ranking as of an earlier year is an illustration. The two weights were fitted on links to 2022.

**Not done.**

- No usability test has been run on the ranking, before or after this change. The usability tests of
  7 and 9 Oct asked 105 questions and none asks for a rank or a candidate, so their results stand for
  the questions they asked.
- Rankings for a trait or a target with no link (`findings/COLDSTART.md`) are not served.
