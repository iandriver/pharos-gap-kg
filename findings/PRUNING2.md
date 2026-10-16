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

## Added 8 Oct 2026: the lift cut applied

`python -m pharos_gap_kg.applycut search pharos319` and `write` produce
`docs/data/pruning2_applied_pharos319.json`. How the cut is applied was added to the protocol before
any pair outside the sample was searched. Every one-source pair of the five sources was searched,
25,303 searches in all, and a pair at or below its source's cut is withheld.

| source | cut | one-source pairs | withheld by the cut | share | share on the 100 pairs the cut was scored on |
|---|---|---|---|---|---|
| CTD | 2.626 | 12,860 | 5,982 | 46.5% | 35.0% [26.4%, 44.7%] |
| Monarch | 0 | 1,692 | 547 | 32.8% | 45.0% [35.6%, 54.8%] |
| JensenLab Knowledge UniProtKB-KW | 2.963 | 1,329 | 273 | 20.5% | 20.0% [13.3%, 28.9%] |
| JensenLab Knowledge MedlinePlus | 0.383 | 972 | 103 | 10.6% | 14.0% [8.5%, 22.1%] |
| UniProt Disease | 0 | 401 | 94 | 23.8% | 21.4% [14.5%, 30.5%] |

Shares are of searchable pairs. 39 pairs were unsearchable and are kept. All 1,000 sampled real pairs
of these sources got the counts they had in the sample.

**Two sources fall outside the interval the check stated**, CTD above it and Monarch below. On all 200
of a source's sampled real pairs, and not only the half the cut was scored on, the share at the cut is
44.5% [37.8%, 51.4%] for CTD and 36.4% [30.0%, 43.3%] for Monarch, and both whole-source figures are
inside those. (This second comparison was added after the stated check was read.) The halves the cut
was scored on were unrepresentative, in opposite directions.

**CTD's cost is 46.5% of its one-source pairs, not the 35.0% in the table further up.** That is still
under the bar's 50%, by 3.5 points. The share of planted pairs withheld on CTD, 83.0% [74.5%, 89.1%],
does not depend on which real pairs were sampled, and rests on 100 planted pairs.

**Monarch has two rules.** Its in-graph list was recomputed on the real table alone: round one had
decided it with planted pairs still in the table, and a planted pair can corroborate a real one. It
is 509 pairs, where round one wrote 508. On the sampled Monarch pairs that carry both verdicts (83
planted, since 2 of the 85 are unsearchable, and 198 real):

| | planted withheld | real withheld |
|---|---|---|
| in-graph rule | 89.2% | 26.8% |
| lift cut | 86.7% | 36.4% |
| either | 95.2% | 49.0% |
| both | 80.7% | 14.1% |

Either rule together stays inside the bar's cost, by one point, so both are applied, as the protocol
said. On the whole source 796 of Monarch's 1,692 one-source pairs are withheld (47.0%): 287 by the cut
alone, 249 by the in-graph rule alone, and 260 by both. The 260 are flagged by two kinds of evidence
that share nothing, and each gap says which rules made it.

**The graph.** 7,248 of 324,882 rows are withheld. The export is 317,634 edges and 7,248 gaps and
passes all ten conformance checks. A gap from the cut carries the number of papers that mention both,
the lift and the cut. Of the 6,999 pairs the cut withholds, 3,672 have no paper at all and 3,327 are
mentioned together, but no more often than the cut for their source.

**What this rests on, and what it does not show.**

- 5,982 of the 7,248 withheld pairs are CTD's (82.5%). CTD's pass rests on 100 planted pairs, and its
  cost in real pairs is within 3.5 points of the bar. A larger planted sample for CTD is the first
  check to make before relying on these gaps.
- Withholding is not deletion. Every pair is in the gaps file with what would restore it.
- No link-prediction result uses the pruned table yet. Whether pruning changes prediction is what the
  fourth fitted model (`lightgcn_ctx_pruned`, `docs/PROTOCOL.md`) is for.

## Added 9 Oct 2026: full text for the sources abstracts could not reach, first day

`python -m pharos_gap_kg.fulltext run pharos319`, then `report`, which writes
`docs/data/pruning2_fulltext_pharos319.json`. The rule, the pairs and the order they are asked in were
added to the protocol after a pilot on five hand-picked pairs and before any sample pair was searched
in this form. The service (Paperclip) allows 100 requests a day, so the 300 pairs take four days.
This is the first: 32 Expression Atlas pairs, 16 planted and 16 real. It is not complete and is not a
verdict.

Rule F withholds a pair if no full text holds the gene symbol within 200 characters of one of the
disease's names.

| Expression Atlas, first 16 of each | planted withheld | real withheld |
|---|---|---|
| full text, rule F | 56.2% [33.2%, 76.9%] | 87.5% [64.0%, 96.5%] |
| abstracts, rule L0, same pairs | 75.0% | 81.2% |

- **So far full text does not find Expression Atlas's real pairs.** 14 of the 16 real pairs have no full
  text that mentions the gene near the disease. The interval's lower end, 64.0%, is already above
  the 50% the bar allows.
- **Planted pairs are found at least as often as real ones**, 7 of 16 against 2 of 16. At this size
  that difference is not established.
- Of the 13 real pairs that no abstract mentions, full text finds 2. Of the 12 planted pairs that no
  abstract mentions, it finds 4.

DrugCentral Indication and JensenLab Knowledge AmyCo have not been asked about yet.

**Expression Atlas was stopped at these 16 of each**, by a protocol amendment made the same day after
these answers were read, on Ian's decision: the interval for real pairs withheld lies above both the
bar and the margin that would have earned an extension. It is reported as stopped, not as having
completed the test. The requests saved go to the other two sources, 200 pairs, asked at up to 90 a day.

## Added 15 Oct 2026: CTD's cut on 400 more planted pairs

`python -m pharos_gap_kg.prune2 extend pharos319 CTD`. The pairs, the measure and what follows were set
in `docs/PROTOCOL.md` ("15 Oct 2026, CTD's lift cut on more planted pairs") before any of them was
searched. The record is `docs/data/pruning2_extend_CTD_pharos319.json`, and the answers are committed
under `recorded/pharos319/`.

The pairs are CTD's planted pairs 200 to 599 in the drawn order of the frozen planting. They were drawn
on 9 Oct with the others, left out of the sample, and never searched. All 400 were searchable. The 400
pair counts were asked on 15 Oct (400 new searches); the gene and disease counts and PubMed's size
(41,252,707 records) are the answers kept from 8 Oct.

| planted pairs | withheld at the cut (2.626) |
|---|---|
| 200 to 599, new | 82.2% [78.2%, 85.7%], 329 of 400 |
| 0 to 99, the pairs the cut was first scored on | 83.0% [74.5%, 89.1%], 83 of 100 |
| both | 82.4% [78.8%, 85.5%], 412 of 500 |

**CTD's cut stands.** On the new pairs it withholds 82.2% of planted pairs, against the bar's 80%, at a
cost of 46.5% of CTD's searchable one-source pairs, against the bar's 50%. The 5,982 CTD rows stay
withheld. The share agrees with the first 100 pairs within a point, and the interval's lower end
moves from 74.5% to 78.2%. It still reaches below 80%: the pass is on the point share, as the bar is
stated, and is not provisional.

A planted pair is a crossing of two of CTD's own pairs that no source asserts, not a pair known to be
false. Of the 71 the cut keeps, 48 are mentioned together in five papers or more, and the most
mentioned, CDKN2A with lung adenocarcinoma, in 214. If some planted pairs are true associations that no
source records, the planted share understates how much of the false pairs the cut withholds, by an
amount this does not measure.

## Added 16 Oct 2026: full text, complete

`python -m pharos_gap_kg.fulltext run pharos319`, then `report`; the record is
`docs/data/pruning2_fulltext_pharos319.json`. Every pair the protocol asks about is answered: 332 pairs,
394 requests over five days, the pilot and the session openings included. DrugCentral Indication was
extended to its first 100 of each by the rule of 9 Oct, which the module applied from 14 Oct, since
both its shares at 50 were within 10 points of the bar. AmyCo's were not, and it stops at 50.
Expression Atlas was stopped at 16 of each on 9 Oct.

Rule F withholds a pair if no full text holds the gene symbol within 200 characters of one of the
disease's names. L0 withholds a pair if no PubMed abstract mentions the gene with the disease.

| source | pairs of each | full text, planted withheld | full text, real withheld | abstracts, same pairs |
|---|---|---|---|---|
| DrugCentral Indication | 100 | 76.0% [66.8%, 83.3%] | 53.0% [43.3%, 62.5%] | 68.0% and 53.0% |
| DrugCentral Indication, first 50 | 50 | 76.0% [62.6%, 85.7%] | 56.0% [42.3%, 68.8%] | |
| JensenLab Knowledge AmyCo | 50 | 80.0% [67.0%, 88.8%] | 70.0% [56.2%, 80.9%] | 80.0% and 70.0% |
| Expression Atlas, stopped | 16 | 56.2% [33.2%, 76.9%] | 87.5% [64.0%, 96.5%] | 75.0% and 81.2% |

**Full text does not bring any of the three sources inside the bar.** DrugCentral Indication misses on
both sides, withholding 76.0% of planted pairs against the 80% wanted and 53.0% of real pairs against
the 50% allowed. AmyCo withholds 70.0% of its real pairs, 20 points over. Expression Atlas withholds
more real pairs than planted ones.

**Full text withholds as many real pairs as abstracts do, on the same pairs:** 53 of 100 for
DrugCentral Indication under both, 35 of 50 for AmyCo under both, and 14 of 16 against 13 for
Expression Atlas. On planted pairs it withholds more for DrugCentral Indication (76 against 68), the
same for AmyCo (40), and fewer for Expression Atlas (9 against 12).

**The two disagree on single pairs, in both directions.** Of DrugCentral Indication's real pairs, full
text finds 10 of the 53 that no abstract mentions, and finds no mention for 10 of the 47 that an
abstract does mention. For AmyCo it finds 10 of 35 and misses 10 of 15; for Expression Atlas, 2 of 13
and 3 of 3. Rule F asks for the gene and the disease within 200 characters on one line, in the papers
the service holds in full text; an abstract counts a mention anywhere in it. Neither contains the
other.

**What follows.** No rule tried in either round cleans these three sources, and their rows are not
withheld. Since 13 Oct DrugCentral Indication's rows are exported as paths through a drug and
Expression Atlas's as observations, which is a change of form and not a pruning result. AmyCo's rows
stay as the source asserts them.
