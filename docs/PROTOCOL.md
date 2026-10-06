# Scoring protocol: target and trait links first reported after 2018

Written and committed on 6 Oct 2026, before any scorer was run on the test links. The commit that
adds this file is the evidence of that order.

## What is scored

The split is the one frozen in `docs/data/split_gene_trait_2018.json`, identified by the SHA-256 of
each part. A link is a Pharos target and a GWAS trait. Its year is the PubMed publication year of the
earliest study cited for it.

- **Training links:** 99,435, first reported in or before 2018, over 14,484 targets and 1,404 traits.
- **Candidates:** the 20,236,101 pairs of those targets and traits with no training link.
- **Positives:** the 130,369 candidates first reported after 2018. They span 540 traits.
- **Negatives:** every other candidate.

A negative is not a known absence. It is a pair no study in this source had reported by the time the
data ends in 2022, so every precision figure is a lower bound.

65,523 later links involve a target or trait with no training link. They are counted and not scored.

## The measures

**Primary.** For each of the 540 traits with at least one positive, rank that trait's candidate
targets and take the average precision. The primary measure is the mean over traits. It answers:
given what was known about this trait in 2018, which targets came next.

**Secondary.**
- the mean over traits of the area under the ROC curve;
- the mean over traits of recall in the top 100 candidates;
- average precision and ROC area by target (for each target with a positive, rank its candidate
  traits), meaned over targets;
- average precision over all 20,236,101 candidates at once;
- the primary measure with positives restricted to one year at a time (2019, 2020, 2021, 2022),
  the other years' positives removed from the candidates.

**Ties** are broken at random, with five fixed seeds. A trait's value is its mean over the seeds, and
the spread of the overall mean across seeds is reported.

## Uncertainty

The unit is the trait. Intervals are 95% percentile intervals from 2,000 bootstrap resamples of
traits. Two scorers are compared by the paired bootstrap of their per-trait difference in the primary
measure, and a difference counts only if its interval excludes zero. The half-width of that interval
is reported with every comparison: it is the smallest difference this test can resolve.

## Baselines

Every later result is shown beside these, and a claim that a graph or a model helps requires beating
the best of them on the primary measure.

| | score of a candidate (target, trait) |
|---|---|
| B0 chance | none: every candidate ties |
| B1 target degree | the number of training traits linked to the target |
| B2 degree product | target degree times trait degree; within a trait it ranks as B1 does |
| B3 neighbours | the sum, over the target's training traits, of the cosine similarity between that trait and the candidate trait, where two traits are similar by the training targets they share |

None has a fitted parameter. B1 knows only how well-studied a target is.

## What a scorer may see

The graph as of 2018 (`data/graph/pharos319_asof2018`) and anything else dated 2018 or earlier. A
table with no dates may be used as context only if the result names it and is reported with and
without it.

## Not yet defined

These must be defined in this file, and committed, before a fitted model is scored. The baselines
above may be scored first, since they have nothing to tune.

- near-duplicate traits: a later link to a trait that is another wording of one the target already has;
- chromosomal neighbours: a later link that a neighbouring gene's training link predicts by position.

## Not claimed

That performance on GWAS links carries to any other kind of link. That a model which ranks later
links well has learned biology and not the order in which cohorts grew.

## Amendments

A change to this file after a fitted model has been scored is added below with its date and reason.
Results under the earlier text are kept.
