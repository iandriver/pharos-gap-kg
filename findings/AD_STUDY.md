# What the graph said of Alzheimer disease, against a study published after it

28 Oct 2026. `python -I -m pharos_gap_kg.adstudy audit pharos319` writes
`docs/data/adstudy_pharos319.json`. The four tests, their readings and their limits were committed to
`docs/PROTOCOL.md` ("A study published after the graph: Alzheimer disease", 036e506), and the code
(c687327), before any gene of the study was matched to the graph. One check was added after the
result, as a dated amendment (test 3).

The study is Uffelmann et al., Nature, 28 Oct 2026 (doi:10.1038/s41586-026-11206-y): 127 loci, 48 of
them new, at an effective size of 416,493. Its predicted effector genes are 191 symbols, of which 181
name a graph gene, GBA1 through its Ensembl symbol. The other 10 are RNA genes, HLA, KIAA0125, and two
earlier names of protein-coding genes, FAM105B and C6orf57, which the matching rule does not reach.
The graph's side is Alzheimer disease (MONDO:0004975) as of 2023, as `gaps "Alzheimer disease"`
serves it: 605 genes with a GWAS link, 74 tied by sources or a drug with no GWAS link, and 90 loci that
carry the link on more than one gene. The sets were computed in full and agree with the served
answer's counts and rows.

## 1. Gaps that close

A gene's gap closes if the study names it.

| the graph's genes with no GWAS link to the trait | genes | named, any locus | named, the 48 new loci |
|---|---:|---|---|
| tied without a link (two or more sources, or a drug's mechanism) | 74 | 3, 4.0% [1.4%, 11.2%] | 0, 0.0% [0.0%, 4.9%] |
| one source | 57 | 4, 7.0% [2.8%, 16.7%] | 1, 1.8% [0.3%, 9.3%] |
| withheld rows only | 1 | 0 | 0 |
| an expression observation only | 251 | 0 | 0 |
| nothing | 19,266 | 133, 0.69% [0.58%, 0.82%] | 84, 0.44% [0.35%, 0.54%] |

Tied minus nothing is +0.034 [+0.007, +0.106] over all loci. By the reading stated in advance the
interval lies above 0: the genes the graph ties to Alzheimer disease with no GWAS link gained one more
often than the genes it does not tie. Three genes carry the figure, ABI3, APP and PLCG2, all at loci
the study calls known. At the 48 new loci none of the 74 is named, and the difference is -0.004
[-0.005, +0.045]; a rise smaller than about 10 points (0.099) had less than 80% power to be seen there.
No difference is detected between tied genes and genes one source asserts (tied minus one source
-0.030 [-0.130, +0.054]). The study names 41 of the 605 genes the graph already links.

## 2. The order of a locus's genes

28 of the graph's 90 loci hold a gene the study names, and at 20 of them the graph's other evidence
points to a first gene. The first gene is one the study names at 17 of the 20, where chance gives 5.7
and the gene with the most studies behind its link gives 7.7. Graph minus chance is +11.3 [+7.2,
+14.8] and graph minus that baseline +9.3 [+5.2, +12.9], so both intervals lie above 0. Where other
evidence singles out one gene, the study names it at 12 of 13 loci. The three loci where the graph's
first gene is not the study's are MPO (the study names TSPOAP1), NDUFAF6 (TP53INP1) and ACTR2
(SPRED2).

The first genes are mostly the genes the literature had already tied to the disease (APOE, CD33, MAPT,
BIN1, TREM2, CLU), and at known loci the study prefers exome and FLAMES evidence and otherwise keeps the gene earlier
studies reported. The data-only
reading, which keeps an effector gene only where FLAMES, the exome analysis or colocalisation also
names it, gives 17 against a chance of 5.5 at 19 loci (+11.5 [+7.9, +14.5]).

## 3. What the largest study could have found

The largest study behind the trait's links is GCST007320 (2019), held as 24,087 cases and 383,378
controls, for which the graph states an odds ratio of 1.077 at a variant carried by 20% of chromosomes.
At each new locus's own frequency, none of the 48 should have been seen by it: their odds ratios run
from 1.029 to 1.114, each below what that study could find at its frequency. For an odds ratio of 1.05
the graph calls for an effective size of 207,947; the study has 416,493.

The study-size table holds GCST007320's 47,793 proxy cases as individuals, which the floor's formula
for cases and controls leaves out. Counted as cases, the floor at 20% is 1.046, and 2 of the 48 new
loci could have been found: DOK3 (odds ratio 1.066 against 1.056 at its frequency) and CEP89/CEBPA
(1.048 against 1.039). A proxy case carries less information than a case, so the two figures bound
the study's reach, and the floor the graph serves for Alzheimer disease, 1.077, is the upper end: the
study could have found somewhat smaller effects than it says.

## 4. The ranking

The served `predict` ranks 17,285 candidate genes for the trait. The study names 85 of them at its
new loci. Their ROC area is 0.682 under the served score, 0.675 under the neighbour score and 0.679
under target degree. Served minus neighbours is +0.008 [-0.014, +0.038] and served minus degree
+0.004 [-0.024, +0.038], so no difference is detected, and the interval excludes a gain of more than 0.038. Over
all 127 loci (139 positives) served minus neighbours is +0.006 [-0.008, +0.023]. All three scores rank a
gene the study names above a randomly chosen candidate about two times in three.

## Limits

- One disease and one study.
- The effector genes are the study's predictions; any of several given for a locus counts.
- The graph's GWAS links stop in early 2022, and some loci the study calls known were first reported
  after that, so the 48 new loci are the cleaner test; there the tied genes closed at 0 of 74.
- The order in test 2 and the study's choice at known loci both draw on the literature; the data-only
  reading narrows this and does not remove it.
- Nothing in the graph changed on this result. Adding the study as a dated source is a separate
  protocol entry.

## Added 29 Oct 2026: the study as a dated source

The rule is `docs/PROTOCOL.md`, "A study added after the release" (43d943d), with its amendment of the
same day (16254ac). `python -I -m pharos_gap_kg.later build pharos319` writes the table to
`data/later/pharos319/`, beside the bundle, and its record to `docs/data/later_pharos319.json`.

The table holds 181 links, one for each predicted effector gene that names a graph gene (181 of the
191 symbols), at 122 of the 127 loci; 88 of the genes sit at loci the study calls new. 41 of the 181
already carry a GWAS link to Alzheimer disease in the release and 140 do not. The study's FLAMES
column also names the gene in 91 rows, Coloc_MR in 140, its literature column in 83 and the exome
analysis in 15.

`gaps` for Alzheimer disease now marks 3 of the 74 genes tied with no GWAS link (ABI3, APP and PLCG2)
and 28 of the 90 loci, the figures of test 1 and test 2 above; no list or count changes. `why`,
`known` and `whats_new` say what the study reported, apart from what the release holds.

The release's files and the release archive are unchanged. Answers to questions that do not involve
the trait are the same with and without the table (`tests/test_later.py` checks a sample), and both
frozen person-test sets pass their check. The hosted server serves the release alone until it is
redeployed with the table.
