# Second pruning round: the literature, and the sources' own fields

7 Oct 2026. `python -m pharos_gap_kg.prune2 <stage> pharos319` writes `docs/data/pruning2_pharos319.json`,
`pruning2_reader_pharos319.json` and `pruning2_sample_pharos319.json`. The planting, the sample, the
checks, the measures and the bar were committed to `docs/PROTOCOL.md` before any pair was planted or
searched. Three amendments made during the round are in that file with their dates and are marked
below where they bear on a figure.

Round one (`PRUNING.md`) found that nothing inside the graph tells a real pair from a planted one at an
acceptable cost. This round asks two things from outside it.

## What could be tested, and on what

A planted pair is a row no source wrote, so it has no confidence, z-score or p-value. A rule that
reads a source's own field cannot be scored on planted pairs. So:

- a **literature check** is scored on planted pairs, as the in-graph rule was;
- each **source's own field** is scored on that source's real pairs, against the literature check.

The sample is 200 planted and 200 real one-source pairs for each of the eleven sources, 4,400 pairs.
Planted pairs are crossed among a source's one-source pairs only. 22 pairs are unsearchable (PubMed
did not find the disease name) and are left out of every share. No planted pair is one that a source
asserts under another identifier with the same name (0 of 2,200).

The bar is round one's: a rule passes on a source if it withholds at least 80% of planted pairs and at
most 50% of real one-source pairs.

## The code-only check

One PubMed search per pair: the gene symbol in the title or abstract, and the disease by PubMed's own
mapping. 8,576 searches in all, counting the symbols and the diseases alone.

**L0, withhold a pair no paper mentions.** Passes on two sources.

| source | planted withheld | real withheld | passes |
|---|---|---|---|
| UniProt Disease | 93.8% | 24.5% | yes |
| Monarch | 88.2% | 36.4% | yes |
| JensenLab Knowledge AmyCo | 78.0% | 61.5% | no |
| JensenLab Knowledge MedlinePlus | 77.0% | 9.5% | no |
| DrugCentral Indication | 73.5% | 53.0% | no |
| Expression Atlas | 72.5% | 61.0% | no |
| DisGeNET | 65.0% | 26.0% | no |
| eRAM | 61.5% | 27.0% | no |
| CTD | 48.0% | 23.5% | no |
| JensenLab Text Mining | 43.7% | 8.7% | no |
| JensenLab Knowledge UniProtKB-KW | 43.5% | 6.0% | no |
| all sources | 67.6% | 30.7% | no |

A third of pairs made at random are mentioned together in some paper. For a much-studied gene and a
common disease, one paper with both is not evidence.

**L0b, withhold a pair whose lift is at or below a cut.** The lift is the number of papers that
mention both, over the number expected if the gene and the disease were mentioned independently. This
rule was added after the counts above were seen (protocol amendment). The cut for each source is set
on 100 of its planted pairs and scored on its other 100 planted and 100 real pairs. Passes on five
sources.

| source | planted withheld | real withheld | passes |
|---|---|---|---|
| UniProt Disease | 94.9% | 21.4% | yes |
| JensenLab Knowledge UniProtKB-KW | 91.0% | 20.0% | yes |
| Monarch | 86.5% | 45.0% | yes |
| JensenLab Knowledge MedlinePlus | 84.0% | 14.0% | yes |
| CTD | 83.0% | 35.0% | yes |
| DisGeNET | 79.0% | 36.0% | no |
| JensenLab Text Mining | 78.0% | 11.2% | no |
| Expression Atlas | 78.0% | 64.0% | no |
| eRAM | 77.0% | 32.0% | no |
| JensenLab Knowledge AmyCo | 76.0% | 61.0% | no |
| DrugCentral Indication | 69.0% | 53.0% | no |
| all sources | 81.5% | 35.8% | yes |

Three more sources miss on the planted side by one to three points at a cost of 11% to 36% (DisGeNET,
JensenLab Text Mining, eRAM). With 100 planted pairs a share of 80% has a 95% interval of about 71% to
87%, so those three are neither passes nor clear failures.

The lift ranks a real pair above a planted one 75.9% of the time over all sources. Round one's best
in-graph score, the neighbours score, does so 66.3% of the time on this planting.

## The reader

For a pair some paper mentions, a model reads the first five abstracts and answers `asserts`, with a
PMID and a passage copied from that record, or `mentions_only`. Code checks that the PMID was one it
was given and that the passage is in that record. The reader can remove a pair L0 kept. It cannot
restore one.

**Choosing it.** Three hosted models read the development part (20 planted and 20 real pairs per
source, 225 of which some paper mentions). As first read, 84, 213 and 208 of 225 answers passed the
check, and none reached the 95% the protocol required. Nearly all the failures were empty answers: the
model had used its 1,500 output tokens before answering. The limit was raised to 4,000 and those pairs
read again (protocol amendment, made before any model read the test part). Then:

| model | answers that pass the check | planted withheld | real withheld | difference |
|---|---|---|---|---|
| deepseek/deepseek-v4-flash | 98.7% | 84.9% | 41.4% | +0.436 |
| z-ai/glm-5.3-flash | 97.8% | 82.6% | 38.6% | +0.440 |
| inception/mercury-2.5 | 94.7% | 81.7% | 40.9% | +0.408 |

Two qualify. Their differences are within 0.02, so the cheaper was chosen: `deepseek/deepseek-v4-flash`.

**L1 on the test part** (80 planted and 80 real pairs per source; 879 pairs read, 9 with no verdict).
Passes on three sources.

| source | planted withheld | real withheld | passes |
|---|---|---|---|
| UniProt Disease | 94.9% | 28.2% | yes |
| DisGeNET | 90.0% | 43.8% | yes |
| JensenLab Knowledge MedlinePlus | 87.5% | 17.5% | yes |
| Monarch | 94.8% | 60.0% | no |
| JensenLab Knowledge AmyCo | 86.2% | 72.5% | no |
| Expression Atlas | 85.0% | 71.2% | no |
| DrugCentral Indication | 82.5% | 66.2% | no |
| eRAM | 76.2% | 43.8% | no |
| CTD | 71.2% | 37.5% | no |
| JensenLab Knowledge UniProtKB-KW | 62.5% | 10.0% | no |
| JensenLab Text Mining | 61.3% | 29.5% | no |
| all sources | 81.1% | 43.7% | yes |

**The reader does not add to the lift, and does worse than it.** Over all sources on the test part, L1 withholds
81.1% of planted pairs and 43.7% of real ones. A cut on the lift that withholds the same 43.7% of real
pairs withholds 87.8% of the planted ones. The difference is -0.067 [-0.103, -0.037]. The stated test
was whether the reader beats that cut, and it does not.

**Why.** Of the 283 planted pairs the reader read, it kept 165 with a passage that code found in the
cited abstract. Eighteen of those, drawn at random, were read by hand. In 14 the passage does state a relation between the gene
and the disease (PLOD1 and gastric cancer, PMID 41186886; CCR7 and hepatocellular carcinoma, PMID
40685403). In 4 it does not: a measurement with no stated association, a study's aim, or a symbol that
means another gene there. That is a hand reading of 18, not a measurement. It says the reader is
mostly doing what it was asked, and that what it was asked is too easy to satisfy: when a paper
mentions a gene and a disease together, it usually states some relation between them, whether or not
any database holds the pair. "Mentioned together more often than chance" separates real from planted
pairs better than "some abstract states a relation".

It follows that a planted pair is false only in the sense that no source asserts it. The shares
withheld are lower bounds on what each rule does with pairs that are false.

## The sources' own fields

For each field, on its source's 200 sampled real one-source pairs: the area under the ROC curve for
"at least one paper mentions the pair". The protocol said a field grades its source at 0.60 or more
with the interval above 0.50.

| source | field | area [95% interval] | grades |
|---|---|---|---|
| DisGeNET | score | 0.642 [0.588, 0.692] | yes, but not independent |
| Expression Atlas | size of the fold change | 0.592 [0.511, 0.673] | no |
| eRAM | number of upstream sources | 0.548 [0.502, 0.593] | no |
| CTD | therapeutic | 0.528 [0.489, 0.562] | no |
| Expression Atlas | minus log10 p | 0.524 [0.439, 0.610] | no |
| DrugCentral Indication | number of drugs | 0.524 [0.469, 0.585] | no |
| DisGeNET | PubMed ids | 0.514 [0.435, 0.589] | no |
| DisGeNET | number of upstream sources | 0.507 [0.482, 0.530] | no |
| Monarch | S2O | 0.467 [0.385, 0.549] | no |
| JensenLab Text Mining | z-score | not measurable | |

**No field that is independent of the label grades its source.** The one that does, the DisGeNET
score, counts papers and upstream sources, so it is partly the label. For six of the other eight an
area of 0.60 is outside the interval. It is inside for the two Expression Atlas fields, whose
intervals reach 0.673 and 0.610, so a weak grading by fold change is not excluded. The JensenLab
z-score could not be measured: 17 of its 195 sampled pairs have no paper, under the 20 the protocol
required, which is expected of a source built from co-mention.

On all of a source's one-source pairs, against round one's lineage criterion, no field's area is
above 0.60 either (0.389 to 0.591).

The sources' fields therefore add nothing to pruning in this round. UniProt Disease and the three
JensenLab Knowledge sources carry no field that varies.

## The in-graph rule again, and a correction to round one

Measured again on the new planting (12,431 planted pairs, 5% of each source's one-source pairs), the
round-one rule withholds 32.7% of planted pairs and 16.4% of real ones, and passes on Monarch only
(88.2% and 29.9%). The 508 Monarch pairs withheld in round one stay withheld.

**Round one's account of its design fault was incomplete.** It said planted pairs were drawn from all
of a source's pairs and compared with its one-source pairs. Planting among one-source pairs did not
remove the fault: the stated check, the product of a pair's two degrees less the pair itself, still
separates real from planted in four sources (0.132 UniProt Disease, 0.240 UniProtKB-KW, 0.323
MedlinePlus, 0.349 Monarch; 0.40 to 0.60 was the range for a matched design). On degrees that count
the pair itself those four are 0.48 to 0.55, so the planted pairs are as common as real ones. The
cause is that a planted pair is added to a target, so its target has one more other pair than a real
pair's target does. In a source where the typical target has one or two diseases (the median is 1 in
UniProt Disease and UniProtKB-KW, 2 in MedlinePlus, 4 in Monarch) that alone separates them.

The fault runs against the rule: a planted pair has more around it to be corroborated by. Weighting
the planted pairs to the real pairs' counts of other pairs raises Monarch's 88.2% to 92.0%, so
Monarch's pass does not come from it. The same weighting takes MedlinePlus from 73.5% to 82.7% on 49
planted pairs. That figure was computed after the run and the rule is not applied to MedlinePlus on it.

The literature checks do not read the table and are not affected.

## Two kinds of evidence agree on the real pairs

Added as description after the stated measures. Among the sampled real one-source pairs, those the
in-graph rule leaves uncorroborated have no paper 39.9% of the time, against 27.9% for those it
corroborates. For Monarch it is 52.8% against 30.3%, and for MedlinePlus 24.5% against 4.6%. The two
checks share no evidence. 203 of the 2,189 searchable sampled real pairs are flagged by both.

## What each source can be pruned by

| source | in-graph | no paper (L0) | lift cut (L0b) | reader (L1) | own field |
|---|---|---|---|---|---|
| UniProt Disease | no | yes | yes | yes | none |
| Monarch | yes | yes | yes | no | no |
| JensenLab Knowledge MedlinePlus | no | no | yes | yes | none |
| CTD | no | no | yes | no | no |
| JensenLab Knowledge UniProtKB-KW | no | no | yes | no | none |
| DisGeNET | no | no | no (79%) | yes | not independent |
| JensenLab Text Mining | no | no | no (78%) | no | not measurable |
| eRAM | no | no | no (77%) | no | no |
| Expression Atlas | no | no | no | no | no |
| DrugCentral Indication | no | no | no | no | no |
| JensenLab Knowledge AmyCo | no | no | no | no | none |

- **The lift cut reaches the most sources**, five, with three more within three points on the
  planted side. It needs no model.
- **The reader reaches DisGeNET**, which the lift cut misses by one point, at a cost of 43.8% of its
  real pairs. It reaches no source the lift cut clearly fails on.
- **Three sources are out of reach of every rule tried**: Expression Atlas, DrugCentral Indication and
  JensenLab Knowledge AmyCo. More than half of their real one-source pairs have no paper that mentions
  the gene with the disease (61.0%, 53.0%, 61.5%), so the literature cannot tell them from random
  pairs. For two of them the reason is what the source states. An Expression Atlas row is a change in
  expression in a disease sample and a DrugCentral row is a drug's target and that drug's indication.
  Neither is a claim an abstract makes about a gene and a disease. The graph already types them as
  `biolink:correlated_with` and `biolink:target_for` and not as gene-disease associations.

## Not done in this round

Nothing new was withheld. The protocol set this round to measure on samples. Applying the lift cut to
the five sources it passes on means searching each of their 17,254 one-source pairs, about 2.7 hours
at NCBI's limit without a key and no model calls.

## Cost

8,576 PubMed searches. Model calls cost 1.22 dollars of the 2.00 cap, about 0.08 cents for each pair
the chosen reader read. No local model was run: the machine had no free memory for one.

## What follows

- The pruner for the disease sources should be the lift cut, applied to the five sources it passes on,
  with each withheld pair written as a gap that names the count of papers, the lift and the cut.
- A larger sample would settle DisGeNET, JensenLab Text Mining and eRAM, which sit one to three points
  under the bar on 100 planted pairs each.
- A reader that only confirms a stated relation adds nothing. If a model is used at all, it should be
  on a harder question than the one asked here, such as whether the paper's own data support the pair.
- Expression Atlas and DrugCentral pairs need a check of their own kind: the expression study behind
  the row, and the drug's activity on the target with the drug's indication.
- Planting by adding pairs gives a planted pair's target one more neighbour than a real pair's. A
  future in-graph test should plant by swapping the diseases of two real pairs.
