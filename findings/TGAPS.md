# Target-level gaps: what a target lacks for the next level, and whether that predicted what came

9 Oct 2026. `python -m pharos_gap_kg.tgaps pharos319` writes `docs/data/target_gaps_pharos319.json`
and two tables of gaps. The definitions, the scores and how each would be read were committed to
`docs/PROTOCOL.md` ("Target-level gaps") before any target's level in Pharos 4.0.0 was joined to them.
Before that, the moves between levels and the number of ligand gains had been seen as counts only.

Pharos grades each target by explicit evidence. A gap here is what a target lacks for the next level
in the 3.19 database, with the nearest way to close it. The back-test asks Pharos 4.0.0 (May 2026)
what happened.

## The rule is Pharos's own

A target with no qualifying compound is Tdark if it meets at least two of three criteria (PubMed
score under 5, at most 3 GeneRIFs, at most 50 antibodies) and has no experimental GO leaf annotation.
Applied to each release's own evidence, the rule gives the release's grade for 99.90% of the 17,793
Tdark and Tbio targets of 3.19 (18 differ) and for all 16,689 of 4.0.0.

## The dark gap

Of 5,516 Tdark targets, 5,499 have a gap (for 17 the evidence and the grade disagree). 2,425 need one
criterion to clear and 3,074 need two. The gap names the nearest criteria and gives a route cost: the
shortfall on each, as a share of its threshold, added up. 5,314 are in both releases, and 788 of them
(14.8%) are no longer Tdark in 4.0.0.

| ranking the 5,314 by | ROC area | average precision |
|---|---|---|
| route cost | 0.859 [0.845, 0.872] | 0.567 [0.529, 0.605] |
| PubMed score | 0.793 [0.777, 0.809] | 0.318 [0.293, 0.342] |
| recent attention | 0.774 [0.758, 0.790] | 0.303 [0.279, 0.328] |
| GeneRIF count | 0.751 [0.733, 0.769] | 0.327 [0.298, 0.356] |
| antibody count | 0.729 [0.711, 0.745] | 0.299 [0.270, 0.329] |
| flips needed | 0.718 [0.702, 0.733] | 0.267 [0.238, 0.283] |

- **The distance to the next level predicts who gets there.** Route cost ranks a promoted target above
  an unpromoted one 85.9% of the time.
- **It is better than how much a target is studied.** Against the best single count, the PubMed score,
  the difference in area is +0.066 [+0.055, +0.077]. Both readings stated in advance are met.

## What a user is shown

The share of each band that was promoted, and beside it the mirror image: Tbio targets that are bio by
their counts alone, the same distance above the thresholds, and the share of them that fell to Tdark.
If levels moved only as knowledge grew, none would fall, so the second column is what the pipeline
alone moves at that distance. (The mirror was added after the first run.)

| band | promoted | fell, in the mirror |
|---|---|---|
| 1 flip, nearest third | 53.6% [50.1%, 57.0%] of 806 | 5.9% of 698 |
| 1 flip, middle third | 19.8% [17.2%, 22.7%] of 804 | 3.0% of 700 |
| 1 flip, farthest third | 6.6% [5.1%, 8.7%] of 737 | 5.2% of 688 |
| 2 flips, nearest third | 11.4% [9.6%, 13.5%] of 1,021 | 0.1% of 1,589 |
| 2 flips, middle third | 2.6% [1.8%, 3.8%] of 994 | 0.0% of 1,588 |
| 2 flips, farthest third | 0.6% [0.3%, 1.4%] of 952 | 0.0% of 1,587 |

- **A target one criterion short and close to it was promoted more often than not**, nine times the
  rate at which the pipeline alone moves targets at that distance.
- **The farthest third of one-flip targets is not told apart from the pipeline**: 6.6% promoted,
  5.2% falling in the mirror.
- The mirror does not depend on distance (5.9%, 3.0%, 5.2%), so what moves those targets is not being
  near a threshold. 100 of 6,850 fell in all.

## The proposal names the route that was taken

640 of the promoted targets needed one criterion. 498 of them left through the counts, and for 85.3%
[82.0%, 88.2%] the criterion the gap had named is one that cleared. Naming one of the two at random
would be right 52.1% of the time. The difference is +0.332 [+0.300, +0.364].

The other 142 left another way, by a GO annotation or a compound, for which the gap has no distance.
Over all 788 promoted targets, 334 hold a GO term in 4.0.0 and 46 a compound. Among the criteria,
the PubMed score cleared for 457, antibodies for 130 and GeneRIFs for 48.

## The counts are not one instrument across the two releases

For the 4,526 targets still Tdark: the PubMed score fell for 15.2% and rose for 62.5% (rank correlation
0.82); GeneRIFs fell for 8.0% and rose for 2.8% (0.85); antibodies fell for 11.3% and rose for 51.7%
(0.95). A count that falls for a tenth of targets is being recomputed as well as growing. The mirror
column above is the size of that effect on the level.

## The ligand gap

17,517 targets have no activity row in 3.19, and 1,359 (7.8%) hold a ligand grade in 4.0.0. What each
lacks is a compound at its family's cutoff: 1 uM for 16,653 of them, 100 nM for 560, 30 nM for 179 and
10 uM for 125. The release keeps only compounds that pass, so this gap has no distance.

| ranking the 17,517 by | ROC area | average precision |
|---|---|---|
| partners with a ligand | 0.746 [0.734, 0.759] | 0.168 [0.156, 0.182] |
| antibody count | 0.729 [0.717, 0.741] | 0.157 [0.146, 0.169] |
| GeneRIF count | 0.722 [0.710, 0.734] | 0.149 [0.138, 0.162] |
| PubMed score | 0.718 [0.707, 0.729] | 0.138 [0.129, 0.149] |
| recent attention | 0.714 [0.703, 0.725] | 0.139 [0.129, 0.150] |

The number of interaction partners that already hold a ligand is ahead of the best single count by
+0.017 [+0.004, +0.030]. All five lie between 0.714 and 0.746: how much a target is studied predicts a
first ligand about as well as anything read here.

## What this is and is not

- It is Pharos's grading read forward. The gap says what the grade lacks. A PubMed score crossing 5 is
  attention to a protein, not knowledge of its function.
- The 3.19 evidence does not reach one date. PubMed scores are complete to 2021 and GeneRIFs run to
  early 2024, so "what came" covers two to four years.
- It is one back-test on one pair of releases.

## What follows

- A typed question for the front end: what a target lacks, the route proposed, and the rate at which
  targets in its band were promoted, beside the mirror.
- A distance for the ligand gap needs compounds that fall short of the cutoff, which ChEMBL holds and
  Pharos does not.
- The GO route, by which 334 targets left, has no distance here.
