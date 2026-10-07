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
