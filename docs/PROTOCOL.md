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

## Controls

Defined on 7 Oct 2026, after the baselines were scored and before any control was computed. No fitted
model had been scored. Each control removes positives that something other than biology could
predict, so each can only make the test harder.

**Trait lineage.** From the Experimental Factor Ontology v3.94.0, `is_a` relations only, with an
obsolete term standing for the term that replaced it. Two traits are in one lineage if they are the
same term or one is an ancestor of the other. A positive (target, trait) is *lineage-linked* if the
target has a training link to another trait in that trait's lineage: the later link is the earlier
one under a broader, narrower or second name.

**Position.** From Ensembl release 113, GRCh38, primary chromosomes. Two targets are neighbours if
genes of theirs lie on one chromosome with at most 500,000 bases between them. A positive is
*position-linked* if a neighbour of the target has a training link to the same trait. The window is
500,000 bases and is not tuned.

**Locus.** Within one trait, positives whose targets chain together as neighbours form one locus. In
the locus-collapsed measure a locus counts once, at the rank of its best-ranked member under the
scorer being measured, and its other members leave the candidates. GWAS assigns one signal to several
neighbouring genes, so positives are not independent; this counts signals, not genes.

Reported for every scorer:

- the share of positives that are lineage-linked, position-linked, and either;
- the primary measure on positives that are neither, with linked positives removed from the candidates;
- the same with every linked candidate removed, positive or not;
- the locus-collapsed primary measure;
- the paired comparisons, repeated on the positives that are neither, both ways.

The two ways of setting linked pairs aside differ in what they ask. Removing only the linked
positives leaves linked non-links among the candidates, so a scorer that ranks linked pairs high now
ranks mostly negatives high: that reading is a lower bound for any scorer that uses lineage or
position. Removing every linked candidate asks only how well the scorer orders the pairs nothing else
explains. A claim that one scorer beats another on new links needs both. (Second way added on 7 Oct
before either was computed.)

Two further baselines, also without a fitted parameter:

| | score of a candidate (target, trait) |
|---|---|
| B4 lineage | the number of the target's training traits that are in the candidate trait's lineage |
| B5 position | the number of the target's neighbours with a training link to the candidate trait |

The ontology is the 2026 release. As a control that is harmless. As a scorer, B4 may use a relation
the 2018 ontology did not have, and is reported with that caveat. A target with no gene on a primary
chromosome has no neighbours, and a trait the ontology does not hold has no lineage; both counts are
reported.

## Disease traits, and the undated disease table

Defined on 7 Oct 2026 before either was computed.

**By kind of trait.** A trait is *disease-like* when the live ontology term it stands for is a MONDO,
Orphanet or HP term, and a *measurement or other* trait otherwise. The primary measure is reported
for each kind as well as overall. What was already known when this was written: 14,227 of the
130,369 positives are on disease-like traits.

**The undated disease table.** Pharos holds target and disease associations from twelve sources with
no dates, as of 2022. A candidate (target, trait) is *in a source* when that source has a row for the
target whose MONDO id is the trait's live term. For each source, over the scored traits whose live
term is a MONDO id:

- how many candidates are in it, and how many of those are positives;
- how many training links are in it;
- the primary measure of the source used as a scorer (1 if the candidate is in it, else 0).

A source that ranks later links well is either prior knowledge or the later links themselves under
another name. Without dates the two cannot be told apart, which is why a result that uses such a
table must also be reported without it.

## Calibration by rank, for the question-answering front end

Defined on 7 Oct 2026 before it was computed. The front end offers ranked candidates for a trait or a
target, scored by neighbours (B3). Beside each it states how often a candidate ranked that high in
the 2018 backtest was reported by 2022. Those shares are computed once, here.

For each scored trait, rank its candidates by B3. Pool over traits the candidates at ranks 1 to 10,
11 to 50, 51 to 100 and 101 to 500, and report the share that are positives: overall, by kind of
trait, and separately for candidates that are lineage- or position-linked and those that are not. The
same by target, for ranks 1 to 5, 6 to 20 and 21 to 50 of each scored target's candidate traits.
Ties are broken at random with the five seeds and the shares averaged.

Every share is a lower bound for the reason given at the top: a candidate not reported by 2022 is not
known to be absent. The front end must say so wherever it quotes one.

## Usability: can a small model answer through the typed questions

Defined on 7 Oct 2026 before any model was run on any question.

**What is tested.** Whether a model given only what the server itself says (its instructions and its
tool descriptions) answers questions a person would put to the graph, and whether it does so better
than with direct query access to the same data.

**Questions.** 35, written by code from the bundle with a fixed seed, in ten kinds: the year a link
was first reported; the studies behind a link as of a year; a link that did not yet exist as of a
year; a pair with no link at all; a gene that is not in the graph; a name that means several things;
the best-supported target of a trait and how many traits that target has; the targets two traits
share; the links first reported after a year; and, for an absent link, a broader or narrower form of
the trait the target already had. The expected answers are computed from the bundle's tables by code
that does not call the question layer. The set is frozen by digest in
`docs/data/usability_questions.json`.

**Conditions.** The same rules of reading are given in each (what "as of" means, that an absence is
not a negative, that an unknown or ambiguous name is not to be resolved by guessing).

- *typed*: the seven tools, with the instructions and descriptions exactly as the MCP server serves them;
- *sql*: one tool that runs read-only SQL over the same tables, with the schema described;
- *none*: no tool.

**Models.** Chosen by price before any was run: three hosted through OpenRouter (`qwen/qwen3.7-flash`,
`deepseek/deepseek-v4-flash`, `z-ai/glm-5.3-flash`) and one local (LFM2.5 2.6B through Ollama). The
local 35B model is not run: this machine is short of memory. Temperature 0 where the endpoint honours
it, at most 8 tool calls a question, one run per cell, total spend capped at 1.50 US dollars.

**Grading.** By code, from a last line `FINAL: {...}` that every question asks for. No model grades.
The status must match, and then the fields of that kind of question: a year, a count, a symbol, or a
trait that is one of the acceptable ones.

**Measures.** Correct of 35 for each model and condition, and by kind. On the three not-in-graph
questions, how often a model answered as if the gene were known. Typed minus sql, paired by question,
for each model (exact sign test) and pooled over the hosted models (bootstrap over questions).

**How it will be read.** The typed questions are usable by a small model if at least two of the three
hosted models answer 28 or more of 35 under *typed* and none of them answers a not-in-graph question
as if the gene were known. Typed is better than SQL if the pooled paired interval excludes zero.
Both verdicts are reported whichever way they fall.

**What it does not test.** A person. And with one run per cell, a difference of one or two questions
between two cells is noise.

## Context: does the rest of Pharos add to the neighbours score

Defined on 7 Oct 2026 before any context scorer was computed.

Context is what the 3.19 database holds besides the dated GWAS links: interactions, pathways, GO
annotations, compound activities and the eleven disease sources that are not the GWAS data. None is
dated except compound activities, and all are as of 2022 at the latest, so a result that uses them is
an upper bound on what context known in 2018 could offer. That bound is the point of running them:
if undated context adds nothing, dated context will not.

Five scorers, none with a fitted parameter. For a candidate (target, trait), "linked" means a
training link.

| | score |
|---|---|
| X1 partners | the number of the target's interaction partners linked to the trait. A partner is a BioPlex or Reactome interaction, or a STRING interaction at combined score 700 or more. Rows that pair a protein with itself are not interactions. |
| X2 pathways | over the pathways that contain the target, the sum of the share of each pathway's other members linked to the trait |
| X3 GO terms | the same over GO terms |
| X4 shared compounds | the same over compounds active on the target and on at least one other target |
| X5 disease sources | for a trait whose live term is a MONDO id, the number of the eleven sources that assert the target with that term |

Each is also combined with the neighbours score B3 by averaging the two percentile ranks among a
trait's candidates, and all five are combined with B3 the same way.

Reported for each: the primary measure, the same on disease-like traits, on new links with every
linked candidate removed, and with each locus counted once; and the paired difference from B3.

**How it will be read.** Context adds to the neighbours score if "B3 with all five" minus B3 has a
95% interval above zero on the primary measure and on the locus-collapsed measure. Each scorer alone
is reported against chance and against target degree.

## Pruning: planted false pairs

Defined on 7 Oct 2026 before any pair was planted.

**The table.** The (target, disease) pairs of the eleven undated disease sources, keeping rows whose
disease has an identifier. 276,658 pairs, of which 248,550 are asserted by one source only.

**Planting.** For each source, false pairs numbering 5% of its pairs, each made by drawing two of the
source's pairs at random and crossing them (the target of one with the disease of the other), kept
only if no source asserts the crossed pair. A planted pair is therefore asserted by one source, and
its target and its disease are as common in that source as a real pair's. Fixed seed.

**The rule under test.** A pair is kept if it is corroborated and withheld if it is not.
Corroborated means any of:

- (a) a second source asserts the same pair;
- (b) some source asserts the same target with a disease in this disease's lineage (the ontology's
  `is_a`, as for traits);
- (c) a GWAS link, of any year, joins the target to the disease or to one in its lineage;
- (d) an interaction partner of the target (as in X1) is asserted with the same disease.

A withheld pair is not deleted. It is written as a gap that names the rule and what would restore it.

**Measures**, for each source: the share of planted pairs withheld; the share of the source's real
single-source pairs withheld, which is the rule's cost (those pairs are not known to be false); and
each of (b), (c), (d) alone. For comparison with a score rather than a rule: the area under the ROC
curve, planted against real single-source pairs, of the neighbours score computed on this table with
the pair itself left out, and of the product of the two degrees, which by construction should be
near a half.

**How it will be read.** The rule is applied to a source if, on that source, it withholds at least
80% of planted pairs and at most 50% of real single-source pairs. A source where it does not is
reported as such and is left unpruned.

## Pruning, second round: evidence from outside the graph

Defined on 7 Oct 2026, after the first round (`findings/PRUNING.md`) and before any second-round pair
was planted, any field was read against a label, or any paper was searched.

**What can and cannot be scored on planted pairs.** A planted pair is a row no source wrote, so it has
no confidence, z-score or p-value. A rule that reads a source's own field cannot be scored on planted
pairs. A literature check can be asked of any pair. This round therefore scores the literature check
on planted pairs, as round one scored the in-graph rule, and scores each source's field on that
source's real pairs against the literature check, with the planted pairs giving the floor.

**Planting, matched.** Seed 20261009. For each source, two of its *one-source* pairs are drawn and
crossed, kept only if no source asserts the crossed pair. The number is 5% of the source's one-source
pairs, and at least 200 where crossing can supply them. The first 5% in drawn order are the in-graph
set; the first 200 in drawn order are the literature sample. Check on the design: the product of the
two degrees should separate real one-source from planted pairs about half the time; a source outside
0.40 to 0.60 is reported as unmatched.

**In-graph rule, repeated.** The round-one rule and the neighbours score are measured again on the
in-graph set, unchanged. The 508 Monarch pairs withheld in round one stay withheld only if the rule
clears the round-one bar on Monarch under matched planting. If it does not, they are restored and the
restoration recorded.

**The literature sample.** For each source, 200 planted pairs and 200 real one-source pairs drawn with
the seed (all of them where fewer). In each 200, the first 20 are the development part and the next 80
the test part for the reader; all 200 are used by the code-only check. The sample is written and its
digest committed before any search.

**The code-only check.** One PubMed search per pair (NCBI E-utilities `esearch`, relevance order):
`"<gene symbol>"[tiab] AND (<disease name>)`, the disease name being the release's normalised name
for that identifier with query punctuation removed, left untagged so that PubMed applies its own
mapping to subject headings. Also one search for the symbol alone, one for the disease alone and one
for the size of PubMed. From these: the number of papers that mention both, and the lift, which is
that number over the number expected if the gene and the disease were mentioned independently.

- *Rule L0:* withhold a pair that no paper mentions.
- A pair whose target has no symbol, or whose disease name PubMed reports as not found, is
  unsearchable. It is counted, left out of the shares, and never withheld on that ground.
- Papers of any year are read. A table pruned this way is not an input to the 2018 split.

**The reader.** For a pair that at least one paper mentions, a model is given the gene symbol, the
disease name and the first five records in relevance order (title and abstract, each cut at 1,800
characters), without being told whether the pair is planted. Its brief is fixed in
`pharos_gap_kg.litcheck.BRIEF`: answer `asserts` with the PMID and a passage copied from that record,
or `mentions_only`. Code checks the answer: the PMID must be one of those given and the passage, of
at least 20 characters, must occur in that record's text. An answer that fails is asked for once more
and then counted as no verdict.

- *Rule L1:* withhold a pair that no paper mentions, and a pair the reader answers `mentions_only`.
  A pair with no verdict is kept. The reader can only remove a pair that L0 kept; it cannot restore one.
- *Choosing the reader.* The three hosted models of the usability test read the development part. A
  model qualifies if at least 95% of its answers pass the code check. Among those, the one with the
  largest difference between planted and real pairs withheld under L1 is chosen; within 0.02, the
  cheaper. Only the chosen model reads the test part. No local model is run: the machine has no free
  memory for one this week.
- Spend on model calls is capped at 2.00 dollars for this round.

**Measures.** For each source and over all sources: the share of planted pairs and of real one-source
pairs withheld, under L0 on all 200 and under L1 on the test part, each with a Wilson 95% interval;
the area under the ROC curve, planted against real, of the number of papers and of the lift; and what
a cut on the lift would cost in real pairs if set to withhold 80% of planted ones.

**How it will be read.**

- *The bar is round one's.* A rule passes on a source if it withholds at least 80% of planted pairs
  and at most 50% of real one-source pairs. A pass whose interval for planted pairs reaches below 70%
  is called provisional.
- *The reader against the simple score.* On the test part, over all sources, a cut on the lift is set
  to withhold the same share of real pairs as L1 does. The reader adds something if L1 withholds more
  planted pairs than that cut, with the 95% interval of the difference (2,000 resamples of pairs)
  above zero.
- *What "false" means here.* A planted pair is one no source asserts. Some will be real associations
  that no source holds, so the share withheld is a lower bound on what the check does with pairs that
  are false. Planted pairs the reader keeps are listed with the passage it copied.
- This round measures on samples. Applying a rule that passes means checking every one-source pair of
  the source, which is a separate run.

**The sources' own fields.** Read for each pair as the strongest value among the source's rows:

| source | field |
|---|---|
| DisGeNET | `score`; the number of PubMed ids in `evidence`; the number of upstream sources in `source` |
| Expression Atlas | minus log10 of `pvalue`; the size of `log2foldchange` |
| JensenLab Text Mining | `zscore` |
| Monarch | `S2O` |
| eRAM | the number of upstream sources in `source` |
| DrugCentral Indication | the number of drugs |
| CTD | whether `evidence` includes "therapeutic" |

UniProt Disease and the three JensenLab Knowledge sources carry no field that varies.

For each field, on the source's 200 sampled real one-source pairs: the area under the ROC curve of the
field for "at least one paper mentions the pair", with a 95% interval (2,000 resamples of pairs), and
the share no paper mentions in the lower and upper half of the field, beside the planted pairs'
share. On all the source's one-source pairs: the same area for "some source asserts this target with
a disease in this one's lineage" (round one's criterion b), resampling targets.

- A field *grades* its source if the first area is at least 0.60 with the interval above 0.50. Where
  fewer than 20 sampled pairs fall on either side of the label, the field is reported as not
  measurable on this sample.
- A field is not used alone to withhold a pair, because its recall on false pairs cannot be measured.
  A field that grades its source orders which pairs the literature check reads first.
- Two fields are not independent of the label and are marked so: the JensenLab z-score is itself a
  co-mention statistic, and the DisGeNET score counts papers and upstream sources.

## Fitted models

Defined on 7 Oct 2026, before any fitted model was trained or scored on anything.

**Tuning sees training links only.** The training links are split again by year: links first reported
in or before 2016 are fitted on (57,411, over 10,648 targets and 850 traits), and links of 2017 and
2018 between a target and a trait already seen are the validation positives (22,610, over 284 traits;
8,993,389 candidates). On that inner split the baselines score 0.0084 (chance), 0.0393 (target degree)
and 0.0712 (neighbours) on the primary measure. Every setting of a model, including when to stop, is
chosen by the primary measure on the validation positives, ties broken with one seed.

**The models**, each scoring a (target, trait) pair by the dot product of two learned vectors and
trained to rank a linked target above a target drawn at random for the same trait:

| | what it is fitted on |
|---|---|
| M1 `mf` | the training links, as two tables of vectors and a number per target and per trait |
| M2 `lightgcn` | the training links as a graph, vectors smoothed over one to three steps |
| M3 `lightgcn_ctx` | the same with the undated context joined on: interaction partners, pathway, GO and compound memberships, and the disease sources as released |
| M4 `lightgcn_ctx_pruned` | M3 with the disease pairs the pruning rules withheld left out |

The disease source that repeats the GWAS data ("JensenLab Experiment TIGA") is in none of them.

**Settings tried**, all of them, for each model: vector length 64 or 128; steps 1, 2 or 3 (M2 to M4);
learning rate 0.001 or 0.005; weight on the squared length of the vectors 1e-6 or 1e-4. Batches of
4,096 links, one random target per link, at most 300 passes over the links, validation every 10
passes, stopping after 50 passes without improvement. Tuning uses seed 0.

**Scoring on the test links, once per model.** The chosen setting and its number of passes are
written to `docs/data` and committed. The model is then fitted on all 99,435 training links under
seeds 0 to 4, with that setting and that number of passes and nothing else changed. Its score for a
pair is the mean over the five. That one matrix is scored with every measure and control above. The
primary measure under each seed alone is reported beside it.

**What counts.** As stated under Baselines and Controls: a model helps if it beats the neighbours
baseline on the primary measure with the paired interval excluding zero, and the claim that it
predicts new links needs the same on the positives neither lineage nor position explains, both ways,
and on loci. M3 against M2 is the reading with and without undated context that "What a scorer may
see" requires. M4 against M3 says whether pruning the disease sources changes prediction.

**What is not claimed in advance.** That any of the four beats the neighbours baseline. A model added
after one of these has been scored on the test links is marked as such.

## Target-level gaps

Defined on 9 Oct 2026, before any target's level in Pharos 4.0.0 was joined to anything computed
here. What had been seen of 4.0.0 before: the table of moves between levels (of the 5,331 Tdark
targets in both releases, 759 are Tbio and 46 Tchem in 4.0.0), and the number of targets that gained a
ligand, in all and by family (1,359 of 17,517; `findings/SNAPSHOTS.md`). No score had been computed
against either.

**What a gap is here.** Pharos grades each target by what is known about it. A target's gap is what
it lacks for the next level, read from the release's own evidence, with the nearest way to close it.

**The dark gap.** In the 3.19 database a target with no qualifying compound is Tdark if it meets at
least two of three criteria (PubMed score under 5, at most 3 GeneRIFs, at most 50 antibodies) and has
no experimental Gene Ontology function or process leaf annotation. Applied to the release's own
tables, that rule gives the release's grade for 99.90% of its 17,793 Tdark and Tbio targets. For a
Tdark target:

- *flips needed*: the number of dark criteria it meets, less one (1 or 2);
- *shortfall* on each dark criterion it meets: 5 minus the score, 4 minus the GeneRIFs, 51 minus the
  antibodies, each divided by its threshold (5, 4, 51);
- *route cost*: the sum of its smallest shortfalls, as many as the flips needed. The criteria in that
  sum are the cheapest route, and they are what the gap proposes;
- the other way out, an experimental GO annotation, has no distance.

**The ligand gap.** A target with no activity row in the 3.19 ligand table lacks a compound at its
family's potency cutoff (30 nM for kinases, 100 nM for GPCRs and nuclear receptors, 10 uM for ion
channels, 1 uM for the rest). The release holds only activities that pass the cutoff, so how close a
target has come is not in it, and this gap is typed without a distance.

**Back-test.** Labels are from Pharos 4.0.0 (May 2026), targets matched by UniProt accession. A 3.19
Tdark target is *promoted* if its 4.0.0 level is not Tdark. A ligand-gap target *gained* if its 4.0.0
level is Tchem or Tclin. The 3.19 evidence does not all reach one date: PubMed scores are complete to
2021 and GeneRIFs run to early 2024.

**Scores**, all from the 3.19 database and none with a fitted parameter.

- Dark gap: PubMed score; GeneRIF count; antibody count; recent attention (the PubMed score added in
  2019, 2020 and 2021); flips needed (1 above 2); route cost (lower above higher).
- Ligand gap: PubMed score; GeneRIF count; antibody count; recent attention; the number of
  interaction partners, as in the context tables, that hold a ligand in 3.19.

Measured by the area under the ROC curve and by average precision, ties broken at random under the
five seeds, with 95% intervals from 2,000 resamples of targets, and paired differences the same way.

**How it will be read.**

- *The distance is worth showing* if route cost ranks promoted targets above the rest with an interval
  that excludes a half. It is *better than how much a target is studied* only if its area exceeds the
  best of the three single counts with the paired interval excluding zero. If it does not, a ranking
  shown to a user is by that single count, and says so.
- *The proposal.* For promoted targets that needed one flip and left through the counts, the
  criterion the gap named is compared with the criteria 4.0.0 shows cleared, by its own evidence
  columns. Naming one of two at random is right half the time. The proposal is right if the interval
  for its share lies above a half.
- *What a user is shown*: the share promoted in each band, by flips needed and by thirds of route cost
  within each, with its interval and the releases it spans.
- *A level can move because the pipeline moved.* 118 of 12,242 Tbio targets are Tdark in 4.0.0. How
  the three counts compare between the releases is reported for targets that stayed Tdark, and
  4.0.0's level is checked against its own evidence columns by the rule in its release notes.

## Usability of the gaps question

Defined on 9 Oct 2026, before any model was put a question about gaps.

**What is tested.** Whether a small model, given only what the server says, answers questions about
what a target lacks (the eighth question, `gaps`), alone and together with the other questions, and
whether it does so better than with SQL over the same tables.

**Questions.** 35, written by code with a fixed seed, in ten kinds: how many criteria a Tdark target
must clear (4); which criterion is nearest for a target one criterion short (4); how many more
GeneRIFs or antibodies such a target needs (4); the share of its band that was promoted in the
back-test, as a percentage (4); the potency cutoff a target with no compound lacks, across families
(4); a Tclin target, which has no gap (3); a gene that is not in the graph (3); a symbol prefix that
means several targets (3); the Tdark target of a family nearest to leaving (3); and, of the three
nearest in a family, how many have any GWAS link, which needs two kinds of question (3). The expected
answers are computed from the bundle's tables by code that does not call the question layer.

**Conditions.** As in the first usability test: *typed* (the eight tools as the server serves them),
*sql* (the same tables, with the gap table and the band rates added and described), and *none*. One
rules text is given in every condition. It includes what a reader needs to answer at all: the levels,
the three criteria and where each clears, what "criteria to clear" means, and the cutoffs by family.

**Models, limits, grading.** The three hosted models of the first test, temperature 0, at most 8 tool
calls a question, one run per cell, graded by code from the `FINAL` line. A percentage is right within
one point. Spend on this test is capped at 0.75 dollars, which keeps all usability testing under the
1.50 stated above.

**Two sets.** A development set (seed 20261010) is run first, under *typed* and *sql*, to find faults
in the test itself. Its results are reported and are not read against the verdicts. The evaluation set
(seed 20261011, with other invented symbols and prefixes) is generated and frozen after the
development round and before any model sees it. Whatever is changed between the rounds is recorded.

**How it will be read.** On the evaluation set, by the two verdicts of the first test, unchanged. The
question is usable by a small model if at least two of the three hosted models answer 28 or more of 35
under *typed* and none answers a not-in-graph question as if the gene were known. Typed is better than
SQL if the pooled paired interval excludes zero.

**What it does not test.** A person. Nor whether these are the questions a user would ask: they were
written by the author of the tool.

## Usability of the rankings

Defined on 12 Oct 2026, before any model was put a question about a ranking.

**What is tested.** Whether a small model, given only what the server says, answers questions about
the rankings the front end serves: the candidates for a trait, where a target ranks, what explains a
candidate, how often a rank has held, and the position lists for targets with no link. And whether it
does so better than with SQL over tables that hold the same rankings. The two earlier usability tests
asked no question of this kind.

**Questions.** 35, written by code with a fixed seed, in eleven kinds:

- the highest-ranked candidate target for a trait (4);
- the rank of a named target that has no link to a named trait (4);
- how often candidates ranked where that target ranks were reported in the test on later links, as a
  percentage, for the candidate's own group (4);
- whether the highest-ranked candidate for a trait is explained by an earlier link, two where it is and
  two where it is not (4);
- for the highest-ranked candidate, how many genes within 500,000 bases already have the link (3);
- the highest-ranked candidate once those an earlier link explains are left out (3);
- among targets with no link at all, the one ranked highest for a trait (3);
- for a target with no link, the trait ranked highest for it (3);
- a target with no link and no linked gene nearby, for which there is nothing to rank (2);
- a gene that is not in the graph (3);
- a trait word that means several traits (2).

The expected answers are computed from the bundle's tables by matrix code that does not call the
question layer. Before a set is frozen, the question layer is asked each question directly and must
give the expected answer; a disagreement means one of the two is wrong and nothing is written.

**Conditions.** As before. *typed*: the eight tools as the server serves them. *sql*: the tables of the
first test, and four more that hold the rankings themselves, described to the model: the 200 highest
candidates for every trait with the three parts of each score and what explains each; the rate at
which each rank held, by group; and the two position lists for targets with no link. *none*: no tool.
One rules text is given in every condition. It says what a reader needs to answer at all: what the
score is made of, what "explained by an earlier link" means, that a target with no link is ranked
apart by position, and what the rate beside a candidate is.

**Models, limits, grading.** The three hosted models of the earlier tests, temperature 0, at most 8
tool calls a question, one run per cell, graded by code from the `FINAL` line. A percentage is right
within one point. A trait may be given by its name or its id. A question is put up to three times,
and one the provider fails each time counts as wrong. Spend on this test is capped at 0.60 dollars,
which keeps all usability testing under the 1.50 stated above: 0.72 has been spent.

**Two sets.** A development set (seed 20261012) is run first, under *typed* and *sql*, to find faults
in the test itself. Its results are reported and are not read against the verdicts. The evaluation set
(seed 20261013, with other invented symbols and words) is generated and frozen after the development
round and before any model sees it, and is run under all three conditions. Whatever is changed between
the rounds is recorded.

**How it will be read.** On the evaluation set, by the two verdicts of the earlier tests, unchanged.
The rankings are usable by a small model if at least two of the three hosted models answer 28 or more
of 35 under *typed* and none answers a not-in-graph question as if the gene were known. Typed is
better than SQL if the pooled paired interval excludes zero.

**What it does not test.** A person. Whether a ranking is right: that is the prospective test. Nor
whether these are the questions a user would ask: they were written by the author of the tool.

## Lineage and position as inputs

Defined on 9 Oct 2026, before any model here was fitted. Every model in this section is **added after
M1 and M2 were scored on the test links** and is marked so wherever it is reported.

**What was known when this was written.** The test scores of every parameter-free scorer and of M1 and
M2 (neighbours 0.0969, `mf` 0.0906, `lightgcn` 0.0825), and the controls: 26,346 of the 130,369 test
positives (20.2%) are lineage-linked or position-linked, and the two scorers built on those relations
alone reach 0.0459 (B4 lineage) and 0.0468 (B5 position). The decision to use lineage and position as
inputs was taken after seeing those figures, which are computed on the test links. Nothing else here
is: every weight and setting below is fitted or chosen inside the training years.

**The inputs** are the two relations of the Controls section, unchanged. Two traits are in one lineage
if one is an ancestor of the other in the Experimental Factor Ontology v3.94.0. Two targets are
neighbours if genes of theirs lie within 500,000 bases on one chromosome (Ensembl 113). A gene's
position has no date. The ontology is the 2026 release and may hold a relation the 2018 ontology did
not, so by "What a scorer may see" every result that uses lineage is reported beside the same model
without it.

**Used in two ways.**

*As scores, combined with fitted weights.* A combiner scores a candidate by a weighted sum of
parameter-free scores of that candidate:

| | scores combined |
|---|---|
| C0 `combine` | B1 target degree and B3 neighbours |
| C1 `combine_pos` | C0 and B5 position |
| C2 `combine_lin` | C0 and B4 lineage |
| C3 `combine_lin_pos` | C0, B4 and B5 |

C0 uses neither input. It is there so that what the inputs add can be told from what fitting two
weights adds.

The weights are fitted on the inner split of "Fitted models". The scores are computed from the links
of 2016 and earlier. For each of the 284 traits with a validation positive, a candidate's probability
is the softmax of its weighted sum over that trait's candidates, and the weights maximise the mean log
probability of the trait's validation positives, averaged over traits so that each trait counts once.
There is no penalty on the weights, no sampling and no seed, so one fit gives one answer.

Three forms of the sum are declared:

| | each score enters as |
|---|---|
| F1 | log(1 + score) |
| F2 | its percentile rank among the trait's candidates, ties sharing the mean rank |
| F3 | F1, and again multiplied by log(1 + the number of fitting links the trait has) |

The form is chosen once, for C3, by five-fold cross-validation over the 284 traits (folds from a
permutation under seed 0): the weights are fitted on four fifths of the traits and the primary measure
is taken on the rest, with ties broken under seed 0. The form with the highest mean over the 284
traits is used for all four combiners, the earlier form on a tie. The weights are then fitted on all
284 traits. The form, the weights and the cross-validated scores are written to `docs/data` and
committed.

Scored on the test links once, the four together: the scores are computed again from all 99,435
training links, combined with the committed weights, and each matrix is scored with every measure and
control. Nothing in a combiner is random, so there are no seeds to average.

*As edges, in the graph model.*

| | what it is fitted on |
|---|---|
| M5 `lightgcn_pos` | M2's graph with an edge between every two training targets that are neighbours |
| M6 `lightgcn_lin_pos` | M5 with an edge between every two training traits in one lineage |

The same 44 settings as M2, tuned, chosen and scored on the test links as "Fitted models" says: seed 0
for tuning, the chosen setting committed, five seeds on all training links, one scoring per model.

**What counts.**

- A model helps if it beats the neighbours baseline B3 on the primary measure with the paired interval
  excluding zero, as before.
- What the inputs add is read from paired differences on the primary measure and on the
  locus-collapsed measure: C1, C2 and C3 each against C0; M5 against M2; M6 against M5.
- A graph model is kept over a combiner only if it beats the combiner with the same inputs (M5 against
  C1, M6 against C3) on the primary measure with the interval excluding zero. A combiner has at most
  eight fitted numbers, so on a tie the combiner is the model.

**What these inputs cannot show.** A candidate is lineage-linked exactly when B4 is above zero and
position-linked exactly when B5 is above zero. With every linked candidate removed, B4 and B5 are zero
for all that remain, so C1 to C3 order them by degree and neighbours alone. With only the linked
positives removed, the linked candidates that remain are all negatives and a combiner that ranks them
high is expected to score below C0. A gain from C1 to C3 is therefore a gain on all links and on loci,
and says nothing about new links. It measures how much of later GWAS reporting is an earlier link
under a related trait name or at the same locus, which a link-prediction score counts as correct. The
graph models can pass a linked trait's or a neighbour's vectors on to pairs that are not linked, so
for M5 and M6 the claim about new links is open and needs what it always needed: both ways of setting
linked pairs aside, and loci.

**What is not claimed in advance.** That any of these beats the neighbours baseline. That a rank
average of B3 with B4 or B5 is a fair simple baseline: it is not computed, because the same average
lost to B3 at every weight for each of the five context scorers (`findings/ENRICH.md`).

## Cold start: new traits and new targets

Defined on 10 Oct 2026, before any scorer was run on the links scored here.

**The links.** 65,523 links first reported after 2018 were left out of the test because their trait
or their target has no training link. The neighbours baseline has no score for a trait it has never
seen and none for a target it has never seen. Trait lineage and gene position, as defined under
Controls, reach both.

**What was known when this was written.** Every result on the test links. And these counts, which say
which traits and targets the 65,523 links involve and nothing about how any scorer ranks them: 52,601
join a new trait to a training target (843 new traits; 550 have a lineage relative among the training
traits, and those hold 37,644 of the links); 10,008 join a training trait to a new target (3,521 new
targets in all, 3,088 with a neighbour among the training targets, holding 9,039 of the links); 2,914
join a new trait to a new target.

**Two tasks.**

- *New traits.* The unit is a trait with no training link. Its candidates are the 14,484 training
  targets, and its positives are its later links to those targets.
- *New targets.* The unit is a training trait. Its candidates are every target of the release with no
  training link, and its positives are its later links to those targets. Read the other way, as a
  secondary measure: the unit is a new target, its candidates are the 1,404 training traits.

Links whose trait and target are both new are scored by neither.

**Scorers.** None has a fitted parameter.

| | for a new trait, the score of a training target |
|---|---|
| N0 chance | none |
| N1 target degree | the number of training traits linked to the target |
| N2 lineage transfer | the number of the new trait's lineage relatives, among the training traits, that the target is linked to |
| N3 lineage neighbours | over those relatives, the sum of 1 where the target is linked to the relative and of its neighbours score B3 for the relative |
| N4 | N2, ties broken by N1 |
| N5 | N3, ties broken by N1 |

| | for a training trait, the score of a new target |
|---|---|
| P0 chance | none |
| P1 locus degree | the number of training links held by the target's neighbours |
| P2 position transfer | the number of the target's neighbours with a training link to the trait |
| P3 | P2, ties broken by P1 |

For the secondary reading by new target, the scorers are the trait's degree (the number of training
targets linked to it), P2, and P2 with ties broken by the trait's degree.

**Measures.** The primary measure of this protocol: average precision per unit over all its
candidates, ties broken at random under five seeds, the mean over units, with a bootstrap interval
over units. For new traits also the locus-collapsed measure. Each is reported over all units, over
the units the input reaches (a new trait with a lineage relative in training; a new target with a
neighbour in training), and over the rest, where a transfer scorer ties every candidate.

**How it will be read.** Lineage transfer works for new traits if N4 beats N1 on the primary measure
over all scored new traits, with the paired interval excluding zero. Position transfer works for new
targets if P3 beats P1 the same way. N5 against N4 says whether passing the relatives' neighbours
scores on adds to counting their links.

**A second period, inside the training years.** The same scorers with a cutoff of 2016, on the links
of 2017 and 2018 whose trait or target was new in 2016. It is computed first and reported beside the
other. It decides nothing.

**What a pass would mean.** A ranked list for a trait or a target that has no link, with a measured
precision. It is not a claim about biology. A new trait with a lineage relative is often a narrower or
broader form of a trait already studied, and a new target beside a linked gene is often the same GWAS
signal assigned to one more gene. The ontology is the 2026 release: a new trait's term and its place
in the ontology were fixed after the cutoff, in some cases because the trait had been studied, though
the placement does not use gene links.

**Not done here.** A fitted model for either task, and relations wider than lineage (two traits under
one parent).

## Rehearsal inside the training years

Defined on 10 Oct 2026, before any of it was run.

**Why.** Every fitted model stood lower against the neighbours baseline on the test links than on the
validation that chose it, the graph models by 0.009 to 0.020. The procedure is therefore run one step
back in time, where its outcome can be read without reading a test link: tune with a cutoff of 2014,
refit on the links to 2016, and score on the links of 2017 and 2018.

**The splits.** Tuning fits on the 26,503 links first reported in or before 2014 (7,047 targets, 562
traits) and validates on the 10,025 links of 2015 and 2016 between a target and a trait already seen
(156 traits). The refit is on the 57,411 links to 2016, and the rehearsal's target is the 22,610 links
of 2017 and 2018 over 284 traits. Those are the fitting links and the validation positives of "Fitted
models". The refit has 2.17 times the links tuning had; in the scored procedure the ratio was 1.73.

**What was known when this was written.** The links of 2017 and 2018 are not untouched. Six models
were tuned on them and the combiners' weights fitted on them. On 10 Oct they were also looked at by
trait size, outside the committed code: the largest third of the 284 traits holds 87% of the fitting
links and the smallest third 2%, and `lightgcn` at its chosen setting is 0.0025 below neighbours on
the smallest third and about 0.010 above on the other two. The tuning curves show that the number of
passes chosen on 57,411 links was applied to 99,435, which is the number of steps of 1.73 times as
many passes; for `lightgcn` that is past the peak of its validation curve. The arms below were chosen
knowing all of this and every result on the test links. So a rehearsal result is evidence about a
procedure, and a claim that a model beats the baseline needs links no model has been scored on.

**Arms.** Each is a whole procedure: tune at 2014, refit at 2016 under seeds 0 to 4, score the mean
of the five matrices on 2017 and 2018.

- *Models:* `mf`, `lightgcn`, `lightgcn_lin_pos`.
- *How a link is drawn in fitting:* each link once a pass, as before (L0); or as many draws a pass as
  there are links, each link drawn with probability in proportion to one over the square root of its
  trait's number of links (L½), or to one over that number, so that every trait is drawn equally
  often (L1). The primary measure counts each trait once, and L1 is the weighting that matches it.
- *How the length of fitting is carried to the refit:* the number of passes that was best in tuning
  (by passes, as before); or that number multiplied by the tuning links over the refit links and
  rounded, at least one, which keeps the number of steps (by steps).
- *Settings:* for the graph models, vectors of 64 or 128, 2 or 4 steps, learning rate 0.001 or 0.005,
  weight 1e-4: eight settings. For `mf`, the four without steps. Batches, the cap of 300 passes,
  validation every 10 and the stopping rule are those of "Fitted models". Tuning uses seed 0.

That is 3 models, 3 ways of drawing and 2 ways of carrying: 18 arms.

- *Combiners:* the weights fitted to the primary measure as amended on 9 Oct, on scores computed from
  the links to 2014 and the positives of 2015 and 2016, all 156 traits, then applied to scores
  computed from the links to 2016. Score sets: B3 with B5; B1 with B3; B1, B3 and B5; B1, B3 and B4;
  all four. Forms: F1, and F3 as scored. Ten arms. Nothing is cross-validated, because the target is
  later in time than the fit.

**Measures.** For each arm, on the rehearsal target: the primary measure, and its paired difference
from the neighbours baseline computed from the links to 2016, with a 95% interval. For description:
the same difference on each third of the traits by size, each seed alone, and the difference from
neighbours the arm showed in tuning at 2014 beside the one it shows on the target.

**How it will be read.**

- A procedure carries if its difference from neighbours on the rehearsal target has an interval above
  zero.
- One way of drawing, or of carrying, is better than another if the paired difference between the two
  arms that differ only in that has an interval that excludes zero, for both graph models. Otherwise
  the earlier way stays.
- With 28 arms read against one baseline, one interval above zero is expected by chance alone even if
  nothing carries. The rehearsal names candidates. It does not confirm them.

**Not rehearsed.** Choosing a setting by the mean of several seeds, which multiplies the cost of
tuning by the number of seeds. An earlier origin still (2012), which has 18,010 fitting links.

## A prospective test: links first reported from 2023

Defined on 10 Oct 2026, before the file it rests on was fetched.

**Why.** The test links of the 2018 split have now been used to score eight fitted models, and a
further model chosen with those results in view could not be tested on them. Pharos 3.19 holds GWAS
links to 2022. The GWAS Catalog, which Pharos's GWAS source aggregates, has gone on being curated.

**The file.** The GWAS Catalog's release of 5 Oct 2026, associations with ontology annotations
(`gwas-catalog-associations_ontology-annotated-full.zip`, 72,826,184 bytes, from
`ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/10/05/`), registered and pinned by digest.

**From an association to a link.** An association gives a variant, the Ensembl genes the Catalog maps
it to (the genes it lies in, and the nearest gene on each side when it lies in none), a p-value, the
ontology terms of the trait, and the study with its publication date. A gene becomes a target through
the release's own cross-references, as for position. A trait is the ontology term, written as Pharos
writes it. Pharos's source applied filters this protocol does not know in full, so the rule is
chosen where it can be checked, on the studies Pharos already holds, and before any later study is
turned into links:

- Four rules are declared: genes the variant lies in only, or those and the nearest gene on each side;
  with a p-value of 5e-8 or smaller, or any p-value the Catalog lists.
- For each rule, the (target, trait, study) triples it produces for the studies that Pharos 3.19 also
  cites are compared with Pharos's own triples for those studies. The rule with the highest F1 is
  used, the earlier in the order above on a tie.
- If the chosen rule recovers less than 90% of Pharos's triples or less than 90% of its triples are
  Pharos's, the test goes ahead and is reported with that mismatch stated. No further rule is tried.

**Dating and the split.** A pair's first year is the earliest publication year among the Catalog's
studies that give it under the rule. Training is every pair of Pharos 3.19, 295,327 in all. A
prospective positive is a pair that Pharos 3.19 does not hold and whose first year is 2023 or later.
A pair that Pharos does not hold and whose first year is 2022 or earlier is set aside: it is neither
a positive nor a candidate, since it was reported before the cutoff and is missing from Pharos for a
reason this protocol cannot see. Positives between a target and a trait that both have a training link
are the test. The others are the cold-start links of this split.

The split is built and frozen by code that reports counts and reads no score. Its digests are
committed before any scorer is run on it.

**What is scored on it.** Fixed in a dated addition below, written after the rehearsal has been read
and before any scorer reads a prospective link. Decided now:

- the parameter-free scorers B1, B3, B4 and B5, and the cold-start scorers of the section above on the
  cold-start links;
- the four combiners with the weights frozen on 9 Oct, and `mf` and `lightgcn_lin_pos` with their
  frozen settings and passes, refitted on all training links under five seeds. These ask whether the
  results on the 2018 split repeat on a second period;
- the best graph arm and the best combiner arm of the rehearsal, each tuned with a cutoff of 2020 on
  the links of 2021 and 2022, refitted on all training links and marked with whether it carried in
  rehearsal. Scoring both whether or not they carried also tests the rehearsal.

Each is scored once, with a marker written before a prospective link is read.

**What counts.** As throughout: a model helps if it beats the neighbours baseline on the primary
measure with the paired interval excluding zero; a claim about new links needs both ways of setting
lineage- and position-linked pairs aside, and loci, with lineage and position taken against all
training links.

**Limits known in advance.** The Catalog's studies of 2023 on differ from those of 2019 to 2022 in
kind, with more protein-level and biobank-wide studies, so the mix of traits differs. Curation lags
publication, so the last year or two are incomplete. And the link rule is fitted to Pharos's, not
taken from it.

## Not claimed

That performance on GWAS links carries to any other kind of link. That a model which ranks later
links well has learned biology and not the order in which cohorts grew.

## Amendments

A change to this file after a fitted model has been scored is added below with its date and reason.
Results under the earlier text are kept.

**7 Oct 2026, usability models, before any question reached a model.** `qwen/qwen3.7-flash` has no
endpoint that this account's zero-data-retention setting admits; the provider refused both smoke-test
requests before they reached a model. It is replaced by the next cheapest model the setting admits
that can call a tool, `inception/mercury-2.5`. The other two hosted models stand.

**7 Oct 2026, usability, second round. Written after the first round was graded and before the second
set of questions was generated.**

The first round, graded as stated: typed 27, 27 and 26 of 35; sql 27, 30 and 30; no tool 0 for all
three. Neither stated verdict is met. Reading the wrong answers showed two faults in the test and one
in the tools, none of which the stated grading can separate from what it meant to measure.

- *The conditions were not given the same rules.* The sql and no-tool prompts spelled out which
  status word goes with which situation. The typed prompt was the server's instructions and one line
  about the format, without that mapping. On the two kinds of question where the status word was the
  only error, typed answers said "answered" 19 times of 21 and sql answers 13 of 21.
- *Two kinds of question had an expected status a correct answer does not naturally give.* "As of
  2019, was X linked to Y? If not, when was a link reported?" answered "no, in 2020" was marked wrong
  unless it also said the status was "no_link_reported", and every one of the 18 answers to the
  related-wording questions named an acceptable trait under the status "answered".
- *`find` invites a choice.* Asked about "alzheimer", it lists four traits and calls the first "the
  closest". All three models then answered for Alzheimer disease, where the server's own rule is to
  ask which is meant.

The first set of questions is now the development set. The second round changes three things:

1. The same rules text, in full, is given in every condition. In the typed condition it follows the
   server's instructions.
2. The answer format no longer asks for a status. Every answer has a `problem` field: null when the
   model found what was asked (including when what it found is that there is no link), or
   "not_in_graph", "ambiguous" or "unknown". A question about whether a link existed asks for
   `linked: true or false`.
3. `find` says whether any name matches exactly, and when none does it says to ask which is meant.

A second set of 35 questions of the same ten kinds is generated with a different seed, different
invented symbols and different ambiguous words, and frozen before any model sees it. Models, limits,
grading by code, measures and the two verdict rules are those stated above, applied to the second set.
The spend cap of 1.50 dollars covers both rounds together. The first round's results stand as recorded.

**7 Oct 2026, second pruning round. Written after the code-only counts were seen and before any reading
was scored.**

The stated measures include what a cut on the lift would cost in real pairs if set to withhold 80% of
planted ones. Seen on all 200 pairs of each source, that cost is under 50% for eight sources, which
the count of papers alone (rule L0) does not achieve. But the cut is set on the pairs it is scored
on, so the figure is not a test of a rule. Added, and marked as added after those figures were seen:

- *Rule L0b.* For each source, the cut is the lift at or below which 80% of the planted pairs in the
  part of the sample the reader does not use (the last 100 in drawn order) fall. A pair is withheld if
  its lift is at or below the cut. The rule is scored on the other 100 planted and 100 real pairs of
  the source, against the same bar.

The comparison of the reader with a lift cut of the same cost stands as stated.

**7 Oct 2026, second pruning round, the reader's output limit. Written after the development part was
read and before any model read the test part.**

The development part as first read, 225 pairs each: answers that passed the code check were 84 for
`inception/mercury-2.5`, 213 for `deepseek/deepseek-v4-flash` and 208 for `z-ai/glm-5.3-flash`. None
reaches 95%. Of the answers that did not pass, 141 of 141, 10 of 12 and 17 of 17 were empty: the model
spent the whole output limit of 1,500 tokens before writing an answer, on both attempts. That measures
the limit, not the reading. Two changes to the harness, made on the development part, which is what
that part is for:

1. The output limit is 4,000 tokens. A pair whose stored answers were all empty under the old limit is
   read again; every other stored answer stands. What the first reading cost stays in the total.
2. A PMID written with the letters "PMID" before the number is read as the number.

The brief, the code check of the passage, the qualifying rule, the choice rule and the spend cap are
unchanged, and are applied to the development part after the re-reading.

**7 Oct 2026, applying the lift cut. Written after the second round's results (`findings/PRUNING2.md`)
and before any pair outside the sample was searched.**

- *Where.* The five sources on which rule L0b passed, each with the cut that was scored for it: CTD
  2.626, JensenLab Knowledge MedlinePlus 0.383, JensenLab Knowledge UniProtKB-KW 2.963, Monarch 0,
  UniProt Disease 0. A cut of 0 withholds a pair no paper mentions, which is rule L0.
- *What.* Every one-source pair of those sources is searched as in the code-only check. A pair whose
  lift is at or below its source's cut is withheld. An unsearchable pair is kept and counted. A pair
  with a second source is not examined.
- *How it is written.* A withheld pair is a gap that names the rule, the number of papers that mention
  both, the lift, the cut, and what would restore it: a lift above the cut, or a second source.
- *Checks.* Each sampled real pair must get the verdict it had in the sample. The share of a source's
  one-source pairs withheld is expected inside the 95% interval of the share on the 100 real pairs the
  cut was scored on; a source outside it is reported.
- *The in-graph list for Monarch is recomputed.* Round one decided which real pairs were
  uncorroborated while its planted pairs were still in the table, and a planted pair can corroborate a
  real one. The list is now computed on the real table alone. The change in its size is reported.
- *Monarch has two passing rules.* A pair is withheld if either withholds it, and the gap names which.
  The union and the overlap are measured on the sampled pairs that carry both verdicts (85 planted,
  198 real). If the union withholds more than 50% of those real pairs it is over the bar's cost, and
  only the rule that withholds fewer of Monarch's real pairs is applied. From figures already
  reported, the union is 97 of the 198, so both are expected to apply.

**8 Oct 2026, fitted models, the settings tried for the graph models. Written after M1 and M2 were
tuned on validation (`findings/MODELS.md`) and before any model was scored on the test links.**

M2's best setting on validation sits at the edge of the declared settings in both directions: the
longest vector (128) and the most steps (3). Validation rises with each step added: at most 0.0670
with one, 0.0726 with two and 0.0743 with three. The settings stopped where the score was still
rising, so the best the model can do on validation was not found.

Added for M2, M3 and M4, to the 24 settings already declared: vectors of 256 with 3, 4 or 5 steps, and
vectors of 128 with 4 or 5 steps, each with both learning rates and both weights. That is 20 more
settings, 44 in all, the same for the three models so that M3 and M4 are compared with M2 on equal
terms. Everything else stands: the batches, the cap of 300 passes, validation every 10, the stopping
rule, seed 0, and the choice of the best on validation, the first of equals in the order tried.

- The 24 settings already tried for M2 are not fitted again. Their recorded values are used.
- The settings are widened once. If the chosen setting is again at the edge, that is reported and
  nothing more is added.
- M1's setting stands as committed. Its eight settings differ by 0.003 and its curve is flat at the cap.
- No test link has been read, and tuning reads none.

**9 Oct 2026, full text for the three sources no rule reached. Written after a pilot on five hand-picked
Expression Atlas pairs outside the sample, and before any sample pair was searched in the form fixed
here.**

More than half of the real one-source pairs of Expression Atlas, DrugCentral Indication and JensenLab
Knowledge AmyCo have no PubMed abstract that mentions the gene with the disease (61.0%, 53.0%, 61.5%),
so no literature rule built on abstracts cleans them. This asks whether full text does.

- *The service and its allowance.* Paperclip (paperclip.gxl.ai) searches the full text of PubMed
  Central and preprints. Its free allowance is 100 requests a day and 25 a minute. This check sends
  at most 90 a day and one request for each pair. Every request is written to a log before it is sent,
  and an answer is kept and not asked for again. The check is optional: nothing else in the package
  needs the service or its key.
- *Rule F.* A pair is withheld if no full text holds the gene symbol within 200 characters, on one
  line, of one of the disease's names, in either order and any letter case, by an exhaustive scan. The
  names are the release's own, the ontology's label and its exact synonyms (EFO v3.94.0), cleaned as
  for PubMed, none shorter than 5 characters, at most 12.
- *Pairs.* From the frozen sample: for each of the three sources, the first 50 planted and the first
  50 real one-source pairs in drawn order. They are asked source by source in the order above, and
  within a source planted and real in turn, so that wherever a day stops the two groups are level. 300
  requests, four days.
- *Measures.* For each source, the shares of planted and of real pairs withheld with Wilson intervals,
  against the round-two bar (80% and 50%), beside the abstract rule L0 on the same pairs; and how many
  of the real and of the planted pairs that no abstract mentions are found in full text.
- *How it will be read.* With 50 pairs a share of 80% has an interval of about 67% to 89%, so a pass at
  this size is provisional whatever the code's flag says. A source whose two shares are each inside
  the bar or within 10 points of it is extended to the first 100 of each; another source stops at 50.
- *What it cannot do.* At 100 requests a day the rule cannot be applied to Expression Atlas's 158,479
  one-source pairs, and DrugCentral's 3,145 would take five weeks. This measures whether full text
  carries the signal.
- *What was seen before this was written.* The pilot's five pairs, none in the sample: under the form
  that asks for the symbol and every word of one name anywhere in a paper, the answers stop counting
  at 20 and a pair named "malignant pancreatic neoplasm" finds nothing, while Rule F with synonyms
  finds 10 papers for it; a crossing of two real pairs (SLC20A1 with breast cancer) finds 45. One
  sample pair, CENPV with lung carcinoma (Expression Atlas, real, first in drawn order), was searched
  in the whole-paper form while the tool was being explored and returned 2 papers. It stays in the sample.

**9 Oct 2026, full text: Expression Atlas is stopped. Written after the first day's 32 answers were
read, on Ian's decision, so this is a change made after seeing data.**

After 16 planted and 16 real Expression Atlas pairs, Rule F withholds 14 of the 16 real ones. The
interval for that share, 64.0% to 96.5%, lies above the 50% the bar allows and above the 60% that
would have earned an extension. The remaining 34 pairs of each kind are not asked about, and the
requests they would have used go to the other two sources. Expression Atlas is reported on 16 of each,
as stopped, and is not said to have completed the test. DrugCentral Indication and JensenLab Knowledge
AmyCo are asked about as stated above: the first 50 planted and the first 50 real pairs of each.

**9 Oct 2026, usability of the gaps question, between the two rounds. Written after the development
round and before the evaluation questions were generated.**

The development round, 35 questions: under *typed* the three models answered 35, 35 and 35; under
*sql* 35, 32 and 32. No answer was invented. Reading the wrong answers found no fault in the
questions, the expected answers or the tool, so none of those is changed. One thing is added:

- *A provider that does not answer.* For one model the provider timed out on one question in each of
  several attempts before it answered. From here on a question is put up to three times, and one that
  the provider fails each time is recorded as not answered and counts as wrong.

The evaluation set is generated and frozen next, and no model has seen it.

**9 Oct 2026, how the combiners' weights are fitted. Written after the four combiners were fitted on
validation as first declared and before any model of "Lineage and position as inputs" was scored on
the test links.**

Fitted as declared (`docs/data/combine_gene_trait_2018_first.json`, code at commit `86b57d5`), the form
chosen is F2, and the cross-validated scores on the 284 validation traits are C0 0.0511, C1 0.0717,
C2 0.0556 and C3 0.0774, beside 0.0712 for the neighbours baseline B3 alone. C0 combines B3 with
target degree and scores 0.0201 below B3 [-0.0345, -0.0082], although a weight of zero on degree
would rank every trait's candidates exactly as B3 does. The likelihood the weights maximise counts a
positive the same wherever it ranks, and the primary measure counts the top of the list most, so the
fit gives up the top of B3's ranking for a better average position. That defeats what C0 is for: it
would make the inputs look as if they added 0.026 when C3 is 0.006 above B3. It would also make C3
an easy bar for the graph models.

Changed, once:

1. *What the weights maximise.* The primary measure itself on the validation positives: the mean over
   traits of average precision, ties broken under seed 0.
2. *Where the search starts.* The weight on B3's own term is fixed at 1, since multiplying every weight
   by one positive number changes no ranking, and every other weight starts at 0. The search starts
   at B3 alone and cannot end below it on the links it is fitted on.
3. *The search.* One weight at a time, in the order the terms are listed, each is set to the value
   that scores highest with the others held, among 0 and plus or minus 2 to the power -6, -5, up to 4,
   keeping its current value on a tie. Passes over the weights repeat until one changes nothing, at
   most ten. Nothing is random.

Unchanged: the four combiners, the three forms (in F3, B3's product with the trait's size has a free
weight), the five folds, the choice of form by C3's cross-validated score, what is committed, one
scoring of the four on the test links, and what counts.

The combiners are fitted this once more. If C3 again cannot be told from B3 in cross-validation, that
is what is reported, and nothing further is tried on the combiners. The first fit stays on record as
a validation result. It is not scored on the test links.

**10 Oct 2026, the link rule of the prospective test. Written after the four declared rules were scored
on the studies Pharos holds, and before any study published after 2022 was turned into a link.**

The four rules, on the 12,260 studies that Pharos 3.19 and the Catalog both hold (all published by
2022, 720,702 triples in Pharos; `docs/data/gwascat_rule_first_pharos319.json`): the best recovers 16.3%
of Pharos's triples, and 67.1% of what it produces is Pharos's. The protocol said that below 90% the
test would go ahead and no further rule would be tried. That sentence is set aside. The four rules
rested on two wrong assumptions about how Pharos's source makes a link, and a test built on them
would not test links of Pharos's kind. Looking further at those same studies, outside the committed
code, showed both:

- *Traits.* The Catalog has re-curated its trait terms since Pharos's snapshot. For 6,407 of the 12,260
  studies the terms now differ, mostly an EFO term replaced by an OBA or MONDO term. Compared as
  written, 50% of the Catalog's (trait, study) pairs are Pharos's and 51% of Pharos's are recovered.
  Compared through the ontology's live term, an obsolete term standing for the term that replaced it
  as under Controls, it is 97.1% and 98.8%.
- *Genes.* Pharos lists three times as many genes for a study, at the median, as the Catalog's own
  gene columns give. Taking instead every gene within a fixed distance of the variant, with a p-value
  of 5e-8 or smaller: 98% of the (gene, study) pairs produced are Pharos's at every distance up to
  100,000 bases, and 51% at 250,000; the share of Pharos's pairs recovered is 50% at 50,000 bases, 79%
  at 100,000 and 81% at 250,000. Without the p-value filter, 88% of what is produced is Pharos's.
  Pharos's source took the genes within 100,000 bases of a variant significant at 5e-8.

Changed:

1. *A trait is the live term of its ontology id*, in Pharos's pairs and in the Catalog's alike. Pharos's
   pairs are re-keyed before the split is made, so two Pharos traits with one live term become one
   trait. A term the ontology does not hold is kept as written.
2. *The rule.* A link joins a trait to every target with a gene within 100,000 bases of a variant
   associated with it at a p-value of 5e-8 or smaller. A target is placed on the chromosome as for
   position. The Catalog's own gene columns are not used. This is the rule that scored best in the
   look above. It is fixed here and is not chosen again.
3. *Its agreement with Pharos* is computed by the committed code at the level a link is made, (target,
   live trait, study), on the shared studies, and reported. From the look it will recover about 77% of
   Pharos's triples, which is under the 90% asked for, with about 97% of its own being Pharos's. The
   test goes ahead with that stated, and nothing further is tried on the rule.

Unchanged: the dating of a pair, what is training, what is a prospective positive, what is set aside,
that the split is frozen by code that reports counts and runs no scorer, and what is scored on it.

Known limit after the change: about a quarter of Pharos's triples for the shared studies are not
reproduced, spread over some 5,300 studies. The cause is not established. The rule is narrower than
Pharos's by that much, so some true links of Pharos's kind will be counted as negatives for every
scorer alike.

Not read: no study published after 2022 has been turned into a link, and nothing about those studies
has been counted beyond the totals of the file.

**10 Oct 2026, what is scored on the prospective split. Written after the rehearsal was read
(`findings/REHEARSAL.md`) and the split was frozen (`docs/data/split_gene_trait_2022.json`), and before
any scorer read a prospective link.**

The split: 294,954 training pairs (Pharos 3.19 by live trait term, 18,005 targets, 2,225 traits);
194,017 prospective links between a training target and a training trait, over 1,110 traits, first
reported in 2023 to 2026; 394,383 cold-start links; 228,538 pairs set aside, which are neither
positives nor candidates anywhere below.

On the 194,017 test links, each scored once:

1. *The parameter-free scorers* B1 target degree, B3 neighbours, B4 lineage and B5 position, from all
   training links.
2. *The 2018 split's models, as frozen.* The four combiners with the weights frozen on 9 Oct, applied
   to scores computed from all training links. `mf` and `lightgcn_lin_pos` with their frozen settings
   and passes, fitted on all training links under seeds 0 to 4, the mean of the five matrices scored.
   These ask whether the results of the 2018 split repeat on a second period.
3. *The rehearsal's candidates.* The protocol said the best graph arm and the best combiner arm. The
   combiner arm with the highest value in rehearsal and the one arm that carried are different arms,
   so three are scored, each marked with what the rehearsal said of it:
   - `lightgcn_lin_pos:L05:steps`, the best graph arm, which did not carry. Its eight settings are
     tuned with a cutoff of 2020 (197,553 fitting links, validated on the 72,234 links of 2021 and
     2022 over 430 traits, seed 0), and the chosen setting is refitted on all training links under
     five seeds with the length of fitting carried by steps.
   - `combine[B1+B3+B5]:F1`, the arm that carried, and `combine[B1+B3+B4+B5]:F1`, the highest in
     value, which did not. Their weights are fitted on scores computed from the links to 2020 against
     the links of 2021 and 2022, and applied to scores computed from all training links.

On the cold-start links, the scorers of "Cold start" as defined there: for a new trait the training
targets are ranked, and for a training trait the targets of the release with no training link.

**Measures and what counts** are those of the protocol throughout: the primary measure with its
paired difference from B3 and a 95% interval; the two ways of setting lineage- and position-linked
pairs aside, with lineage and position taken against all training links; and the locus-collapsed
measure. A model helps if it beats B3 on the primary measure with the interval excluding zero. For
cold start, lineage transfer works if N4 beats N1 and position transfer if P3 beats P1, as before.

**What the rehearsal leads one to expect**, stated so that the outcome can be read against it: that
no fitted graph model beats B3; that `combine[B1+B3+B5]` is a little above B3; and, from the 2018
split, that lineage transfer with ties by degree is well above degree for new traits.

**Known about this split before scoring.** Its new traits are of another kind than the 2018 split's:
8,425 of them, against 829, from four years in which the Catalog took in many studies of protein and
metabolite levels. A trait of that kind is often linked to the gene that encodes the measured
protein, which neither lineage nor degree knows.

**10 Oct 2026, the figures behind the link rule, regenerated. Written after the split was frozen and
before any scorer read a prospective link.**

The figures in the amendment on the link rule came from a look outside the committed code.
`python -m pharos_gap_kg.gwascat look pharos319` now regenerates them
(`docs/data/gwascat_look_pharos319.json`), and four differ by a little: the trait terms differ for
6,433 of the 12,260 studies, counting a difference either way, where the text says 6,407; by live term 98.7% of
Pharos's (trait, study) pairs are recovered, where it says 98.8%; and the share of Pharos's (gene,
study) pairs recovered is 51% at 50,000 bases and 82% at 250,000, where it says 50% and 81%. None of
these bears on the rule.

One statement in that amendment was wrong, and it understated the rule. It expected the rule to recover about
77% of Pharos's triples and so to miss the 90% asked for. That figure was taken gene by gene, and it
left out the targets that Pharos ties only to gene ids on alternate assemblies of a chromosome, which
`position` places by their symbol. At the level a link is made, (target, live trait, study), the
committed check (`docs/data/gwascat_rule_pharos319.json`) gives 92.8% of Pharos's triples recovered and
95.1% of the rule's being Pharos's, and for pairs 92.0% and 95.3%. The rule passes the bar. The triples
it does not reproduce are 7.2% of Pharos's, in 4,709 studies, and not the quarter stated there.

**11 Oct 2026, the scorer behind the question layer. Written on Ian's decision to switch it, after the
prospective test was scored and before anything below was computed.**

- *What changes.* The front end ranked candidates by neighbours (B3). It now ranks the targets for a
  trait by `combine[B1+B3+B5]:F1` with the weights that were scored on the prospective test: log(1 + B3),
  plus 1/32 of log(1 + B1), plus 1/4 of log(1 + B5), each computed from the links as of the year asked.
  This is the model named in advance and shown to beat B3 on links no model had been scored on.
- *Calibration.* "Calibration by rank" is computed again for this scorer, on the prospective split, as
  defined in that section: the same rank bins and groups and the five seeds. A candidate is a pair with
  no training link that is not set aside, and a positive is a link first reported in 2023 to 2026. The
  front end quotes these shares and says which test they come from. This reads the prospective links
  again, for a scorer already scored on them.
- *Ranking the traits for one target.* The combiner was fitted and tested for ranking targets within a
  trait. Read the other way, the degree term is the same for every trait of one target, so the order
  is by log(1 + B3) plus 1/4 of log(1 + B5). Before it is served that way, its mean average precision
  over targets is compared with B3's on the prospective split, paired over the targets with a
  prospective link. It is served by target unless the interval of that difference lies below zero, in
  which case B3 stays for that direction. This choice is made on the prospective links, so for that
  direction what the front end quotes is not from a clean test, and its documentation says so.
- *As of an earlier year.* The two weights were fitted on links to 2022. A ranking as of an earlier
  year uses weights fitted on later links, so what it shows as reported later is an illustration and
  not a backtest. The clean test is the prospective one.
- *Trait ids.* The prospective test keyed traits by live ontology term, 2,225 of them. The bundle
  keeps Pharos's ids, 2,247. Where two ids were merged the served score differs from the tested one. A
  check compares the served scores with the tested matrix for the traits that were not merged.
- *What does not change.* Every other question, the grouping of candidates into leads, the cautions,
  and the two statements every answer carries: a candidate is a lead, and one not reported is not
  known to be absent.
- *Usability.* The tests of 7 and 9 Oct were run on the layer as it was. None of their 105 frozen
  questions asks for a rank or a candidate, so their expected answers are unchanged. No usability test
  has been run on the ranking itself, before or after this change.

**12 Oct 2026, position rankings for targets with no link, and a correction to calibration by rank.
Written on Ian's decision to serve the rankings, before either was computed.**

*What is served.* Since 11 Oct a target with no link as of the year asked is not ranked. It is now
ranked by position alone, with no fitted number, in the two readings scored under "Cold start":

- *For a trait*, the targets with no link that have a gene within 500,000 bases linked to the trait,
  ordered by the number of such genes, ties broken by the number of links those genes hold in all
  (P3). It is a list of its own, given when asked for. It is not merged into the combiner's list,
  because the two were never scored on one list.
- *For a target with no link*, the traits that a gene within 500,000 bases of it is linked to, ordered
  by the number of such genes, ties broken by the number of targets linked to the trait (Q3).

P3 beat P1, and Q3 beat Q1, with the interval above zero in each of three periods: the cutoffs of 2016
and 2018 and the prospective split. A target with no link and no linked gene nearby is not ranked.
Targets with no link that lie beside each other are shown as one lead. Nothing is fitted, so a ranking
as of an earlier year is a backtest in the plain sense.

*Calibration for it.* On the cold-start links of the prospective split. For every training trait that
has a listed candidate (a target with no training link, not set aside, with a gene nearby linked to
the trait), the listed candidates are ranked by P3, and over all those traits the share of candidates
at ranks 1 to 5, 6 to 20 and 21 to 50 that are positives is taken. The same for every target with no
training link that has a listed trait, ranked by Q3. Ties are broken under the five seeds.

*The correction.* "Calibration by rank" pooled over scored units only, that is, over the traits (or
targets) that had at least one positive. Someone asking about a trait cannot know it will be one of
those, and the sentence the front end quotes did not say the share was conditional. On the
prospective split 1,110 of the 2,225 traits are scored, so a share quoted for the top of a trait's list
was about twice what holds for a trait asked about without that knowledge. From here on a share is
pooled over every unit the front end ranks: every trait with a training link, and every target with
one. The shares over scored units are kept in the record under their own name. The front end quotes
the others, and says what they are pooled over. The 17.8% quoted on 11 Oct for the ten highest
candidates of a trait is a share over scored traits.

**12 Oct 2026, usability of the rankings, between the two rounds. Written after the development round
and before the evaluation questions were generated.**

The development round, 35 questions (`docs/data/usability_results_rank_dev.json`): under *typed* the
three models answered 32, 30 and 30; under *sql* 34, 34 and 35. No answer was invented. All 13 typed
misses are of two kinds, the rank of one named target and the rate for that rank. Those questions are
answered by `why`, which gives a target with no link its own rank among the trait's candidates and
the rate that goes with it. In none of the 13 misses was `why` called: the models went to `predict`
and looked for the target in the list. That list shows the top of the ranking only, folds neighbouring
genes into one lead under the rank of its first gene, and was cut short by the harness's limit of
8,000 characters on a tool result. Of the 11 answers of those two kinds that were right, `why` had
been called in 8. Reading the misses found no fault in the questions, the expected answers or the
grading. The fault is in what the server tells a model about its own tools.

Changed, in the tools and not in the test:

1. *The server's instructions* send "where does X rank for Y" and "how often has a rank like that
   held" to `why`, and say that `predict` shows only the top of a list.
2. *The descriptions of `why` and `predict`* say the same, each in its own words.
3. *A list of candidates for a trait* carries one sentence on how to read it: a row is a lead under the
   rank of its first gene, and a gene listed with it has its own rank, which `why` gives.

Unchanged: the eleven kinds of question and how they are written, the expected answers and how they
are computed, the rules text, the tables of the *sql* condition, the grading, the models, the limits
and the two verdicts. The typed tools returned about three times as many tokens a question as the SQL
tool in this round (15,691 against 4,816 on average); that is reported and not changed.

The evaluation set is generated and frozen next, and no model has seen it. The development round cost
0.23 of the 0.60 dollars.

**13 Oct 2026, one graph, and each source exported as what it asserts. Written before the merged
graph was built or any question read it.**

Why. Three exports exist (the dated links as of 2018, the disease sources, the target-level gaps),
with separate node files, and the questions read only the dated links and the target gaps. The dated
links after 2018 are in no KGX file. And two of the three sources that no pruning rule reaches are
not exported as what the release says they are. Looked up in the release's own tables on 13 Oct:

- *DrugCentral Indication.* All 13,919 rows name a drug. Each is a row of the release's drug-activity
  table (a drug with a measured or stated action on a target) joined to an indication of that drug:
  10,327 through a target the release marks as the drug's mechanism of action, 3,592 through another
  activity. 7,738 rows name the indication without an identifier, and they are every disease row
  `ctxtables` leaves out for having none.
- *Expression Atlas.* All 159,846 rows carry a log2 fold change of at least 1.0 in size and a p-value
  below 0.05, over 142 conditions. The export carries neither.

These counts are a look. The build recomputes each and writes it to its record; a difference is
reported here as a dated note.

What changes.

1. *One export*, `data/graph/<release>/`, with one node file: every dated GWAS link of the release
   with its first year, studies and papers; the disease-source assertions; the expression
   observations; the drug edges; every withheld row as a gap; every target-level gap. It is written
   only if it passes all ten conformance checks, and every input row is an edge or a gap. The three
   earlier exports stay as built. The graph as of 2018 is the training graph of the 2018 test and
   stays pinned to that split's digest.
2. *DrugCentral Indication* is no longer a direct edge from target to disease. It is two edges
   through a drug node (`DrugCentral:<id>`, category `biolink:Drug`): drug to target
   (`biolink:directly_physically_interacts_with`, carrying the release's action type and whether the
   target is the drug's mechanism of action) and drug to indication (`biolink:treats`). An
   indication with no identifier is a node under the local prefix, flagged as not grounded. The
   target and disease pair is the path, and an answer that reports it says whether it runs through
   the drug's mechanism target.
3. *Expression Atlas* stays `biolink:correlated_with` and gains what the release measured: the
   subject aspect `expression`, a direction from the sign, the fold change as the effect size, the
   p-value, and the source's own name for the condition as context. Its floor is 0.999: the release
   holds no row below 1.0 at three decimals, so a pair that is absent may have a fold change up to
   that size. This is the release's inclusion rule and not a statement of any experiment's power.
4. *The other nine sources* keep their mapping. A kept edge also carries the source's own evidence
   field where the release gives one.
5. *The questions* read the merged tables. `known` and `why` report, beside the dated GWAS links,
   which sources assert a pair, which rows are withheld and by which rule with what was counted, the
   drugs that act on a target with their indications, and expression observations in a list of
   their own, apart from the assertions. `connect` may pass through a drug. `predict` is unchanged:
   the tested scorer reads the dated links only.

What does not change. The pruning tests and their results stand as measured on the pair table of
the eleven sources, and so does "three sources are out of reach of every rule tried". Nothing here
prunes them. No withheld row belongs to either remodelled source, and the build refuses if one
does. A pair reached through a target that is not the drug's mechanism is flagged and not withheld:
no planted-pair test of that flag has been run, and a planted pair has no drug. The context arm's
result, the fitted models and every scored test split are untouched. The usability results on
record are of the question layer as it was before this change.
