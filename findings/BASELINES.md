# Baselines on the frozen 2018 split

6 Oct 2026. `python -m pharos_gap_kg.baselines --cutoff 2018` writes
`docs/data/baselines_gene_trait_2018.json`. The protocol (`docs/PROTOCOL.md`) was committed before this
was run. None of the four scorers has a fitted parameter.

14,484 targets, 1,404 traits, 20,236,101 candidates, 130,369 positives over 540 traits.

| scorer | primary: mean average precision over traits | 95% interval | ROC area by trait | recall in top 100 | average precision by target | over all candidates |
|---|---|---|---|---|---|---|
| B0 chance | 0.0178 | 0.0149 to 0.0209 | 0.503 | 0.007 | 0.0117 | 0.0064 |
| B1 target degree | 0.0733 | 0.0621 to 0.0852 | 0.686 | 0.061 | 0.0117 | 0.0145 |
| B2 degree product | 0.0733 | 0.0621 to 0.0852 | 0.686 | 0.061 | 0.0754 | 0.0614 |
| B3 neighbours | 0.0969 | 0.0838 to 0.1117 | 0.694 | 0.101 | 0.0875 | 0.0391 |

Paired over traits, on the primary measure:

| comparison | difference | 95% interval |
|---|---|---|
| target degree minus chance | +0.0555 | +0.0455 to +0.0663 |
| neighbours minus target degree | +0.0236 | +0.0130 to +0.0359 |

## What it shows

- **How well-studied a target was in 2018 ranks later links at four times chance** (0.0733 against
  0.0178). Nothing about the trait enters that score.
- **Shared-trait neighbours add 0.0236 to that.** The ROC area hardly moves (0.686 to 0.694); the gain
  is at the top of each ranking, where recall in the first 100 candidates goes from 0.061 to 0.101.
- **The number a fitted model has to beat is 0.0969.** The paired interval for neighbours against
  degree is 0.0115 wide on each side, which is the scale of difference this test resolves.
- **The pooled ranking tells a different story.** Over all candidates at once the degree product
  (0.0614) beats neighbours (0.0391). A few traits hold thousands of links each and dominate a pooled
  ranking; the mean over traits weights each trait once. The primary measure was fixed as the latter
  before either was computed.
- By target, target degree equals chance (0.0117 both), as it must: a target's degree is the same for
  every trait it could be linked to.

By year of the later link, primary measure:

| year | links | traits | chance | target degree | neighbours |
|---|---|---|---|---|---|
| 2019 | 23,972 | 222 | 0.0083 | 0.0361 | 0.0696 |
| 2020 | 48,018 | 222 | 0.0169 | 0.0750 | 0.0916 |
| 2021 | 42,639 | 265 | 0.0122 | 0.0598 | 0.0712 |
| 2022 | 15,739 | 130 | 0.0094 | 0.0584 | 0.0703 |

Random tie-breaking moves the primary measure by at most 0.0003 across the five seeds.

## What it does not show

The neighbours gain has not been separated from two things the protocol lists as undefined controls.
A trait reported under two wordings gives a target a later link that its earlier link to the other
wording predicts. And genes that sit together on a chromosome share GWAS links by position. Either
could account for some or all of the 0.0236, and neither has been measured.

## Since: the two controls were measured (7 Oct 2026)

`findings/CONTROLS.md`. A fifth of the positives are linked to a training link by trait lineage or by
chromosomal position. On the positives that neither explains, the neighbours gain over target degree
is +0.0121 with every linked candidate removed, and +0.0060 with an interval that includes zero when
only the linked positives are removed. The 0.0236 above is therefore not a gain on new links.
