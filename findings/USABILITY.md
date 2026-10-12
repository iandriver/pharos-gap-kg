# Can a small model answer through the typed questions

7 Oct 2026. `python -m pharos_gap_kg.usability run` puts frozen questions to a model and
`python -m pharos_gap_kg.usability report [--set eval]` grades the stored answers by code. Results:
`docs/data/usability_results.json` (first round) and `docs/data/usability_results_eval.json` (second).
The protocol, and the amendment that separates the two rounds, were each committed before the runs
they govern (`docs/PROTOCOL.md`).

Three hosted models chosen for price (0.04 to 0.15 dollars per million input tokens), 35 questions a
round, three conditions: the seven typed tools as the MCP server serves them, read-only SQL over the
same tables, and no tool. One run per cell.

## Second round: the one the verdicts rest on

A new set of 35 questions, generated after the first round's faults were fixed and frozen before any
model saw it.

| model | typed tools | SQL | no tool |
|---|---|---|---|
| mercury-2.5 | 35 | 31 | 1 |
| deepseek-v4-flash | 35 | 33 | 4 |
| glm-5.3-flash | 35 | 35 | 1 |

**Usable by a small model: met.** The stated bar was 28 of 35 for two of the three models with no
not-in-graph question answered as if the gene were known. All three answered 35 of 35.

**Better than SQL: met, narrowly.** Pooled over the three models and paired by question, typed minus
SQL is +0.057 of a question, 95% interval +0.010 to +0.105. That is six answers of 105: typed right
and SQL wrong six times, and never the reverse (sign test on the six, p = 0.031). For no single model
is the difference significant (p = 0.125, 0.5 and 1.0), and one model made no error either way.

The six answers SQL got wrong:

- "How many targets is arthritis linked to" was answered "not in the graph" by two models. They
  looked for a trait named exactly "arthritis", found none, and stopped. Several traits contain the
  word.
- Two questions of the form "as of 2010, was X linked to Y" were answered "yes" by one model for
  links first reported in 2013 and 2019. It found the link and did not apply the year.
- One count of shared targets came back 2 where it is 17.
- One related trait was not found.

**No model invented an answer** for a gene that is not in the graph, in any condition of either round.

**Without a tool** the models answered 1, 4 and 1 of 35, each time by declaring a problem that
happened to be the expected one.

Typed questions took 2.4 tool calls a question against 4.0 for SQL, and 0.0007 dollars a question
against 0.0010.

## First round: now the development set

Graded as stated in advance:

| model | typed tools | SQL | no tool |
|---|---|---|---|
| mercury-2.5 | 27 | 30 | 0 |
| deepseek-v4-flash | 27 | 27 | 0 |
| glm-5.3-flash | 26 | 30 | 0 |

Neither verdict was met: no model reached 28 under typed, and typed minus SQL was −0.067 of a
question (−0.143 to +0.010). These results stand.

The wrong answers showed two faults in the test and one in the tools.

- **The conditions had not been given the same rules.** The SQL and no-tool prompts said which status
  word goes with which situation; the typed prompt did not. On the two kinds of question where the
  status word was the only error, typed answers used "answered" 19 times of 21 and SQL answers 13 of 21.
- **Two kinds of question expected a status that a correct answer does not naturally give.** With the
  status word set aside, the facts were right in 34, 34 and 32 of 35 typed answers and 33, 34 and 33 of
  35 SQL answers. That regrading was done after the fact and is not the stated grading.
- **`find` invited a choice.** Asked about "alzheimer" it listed four traits and called the first "the
  closest", and all three models answered for Alzheimer disease. `find` now says whether any name
  matches exactly, and when none does it says to ask which is meant.

## What this shows and what it does not

- Three cheap hosted models answered every question of a fresh set through the typed questions,
  including two-step questions, dated questions, absences, an unknown gene and an ambiguous name.
- The advantage over SQL is small and comes from the cases the typed questions were built for: a
  name that means several things, and a date that has to be applied.
- Between the rounds typed went from 27, 27, 26 to 35, 35, 35 and SQL from 30, 27, 30 to 31, 33, 35.
  Three things changed at once (the rules given, the answer format, `find`), so the gain cannot be
  divided among them.
- The second set is new, but it is the same ten kinds of question, written by the author of the tools.
  Nothing here tests what only the typed questions offer, a ranked prediction with its track record,
  because SQL has no counterpart to compare.
- No person has used it. This is a test of whether a model can, which is necessary and not sufficient.
- One run per cell. The pooled interval clears zero by a hundredth of a question.

## Cost

630 hosted conversations over the two rounds cost 0.46 dollars as recorded, against a cap of 1.50.
Requests that returned an API error were retried and are not in that total. One model named in the
protocol could not be run: the account's zero-data-retention setting excludes its only endpoint, and
the next cheapest admitted model replaced it before any question reached a model.

## Added 9 Oct 2026: the gaps question

`python -m pharos_gap_kg.usability run pharos319 --set gaps ...` and `report --set gaps` write
`docs/data/usability_results_gaps.json`; the development set is `gaps_dev`. The test was written into
`docs/PROTOCOL.md` ("Usability of the gaps question") before any model was asked about gaps: 35
questions in ten kinds about what a target lacks, three conditions, the same three hosted models, and
the two verdicts of the first test. The expected answers come from the bundle's tables, and the
question layer, asked separately, gives the same answer to all 35 of each set.

**Development round** (typed and SQL only; reported, not read against the verdicts):

| model | typed | SQL |
|---|---|---|
| deepseek-v4-flash | 35 | 32 |
| glm-5.3-flash | 35 | 35 |
| mercury-2.5 | 35 | 32 |

Reading the wrong answers found no fault in the questions, the expected answers or the tool, so
nothing was changed but one rule: a question the provider fails three times is recorded as not
answered. None was in either round.

**Evaluation round**, on a set generated and frozen after the development round:

| model | typed | SQL | no tool |
|---|---|---|---|
| deepseek-v4-flash | 35 | 33 | 11 |
| glm-5.3-flash | 35 | 34 | 9 |
| mercury-2.5 | 34 | 33 | 3 |

- **Usable by a small model: met.** All three answer 28 or more of 35 through the typed questions, and
  no model in any condition answered a not-in-graph question as if the gene were known.
- **Better than SQL: not met.** Pooled over the three models, typed minus SQL is +0.038
  [+0.000, +0.086] of a question, and the interval reaches zero. Typed was right where SQL was wrong
  on five answers, and SQL where typed was wrong on one.
- **Where SQL went wrong.** Two models reported a target named by its accession (Q9BRP9) as not in the
  graph, though it is in the table. The other errors were on the nearest target of a family and on
  the two-step question, which one model also failed under typed by answering "unknown".
- **The two-step question works through the tools.** Asked how many of a family's three nearest
  targets have a GWAS link, the models called `gaps` and then `known` and were right in 8 of 9 answers.
- **Without a tool the models score 3 to 11**, well above the 1 to 4 of the first test. A quarter of
  these questions can be answered from the rules every condition is given and from general knowledge:
  the level of a target with a well-known drug, or the cutoff for a kinase.

The typed condition used 1.62 tool calls a question (`gaps` 89 times, `find` 57, `known` 24) and SQL
2.25.

**What this shows and what it does not.**

- Three small models can use the eighth question, alone and with the others.
- SQL over the same tables does nearly as well. The gap table already holds each answer in a column,
  so little is left for the tool to add beyond resolving a name.
- Two of the 35 evaluation questions are also in the development set, because the family questions
  draw on ten families. No model carries anything from one round to the next.
- The questions were written by the author of the tool, and no person has used it.

**Cost.** 525 conversations over the two rounds cost 0.26 dollars against this test's cap of 0.75.
All usability testing to date has cost 0.72 dollars of the 1.50 stated.

## Added 12 Oct 2026: the rankings

`python -m pharos_gap_kg.usability run pharos319 --set rank ...` and `report --set rank` write
`docs/data/usability_results_rank.json`; the development set is `rank_dev`. The test was written into
`docs/PROTOCOL.md` ("Usability of the rankings") before any model was asked about a ranking: 35
questions in eleven kinds, three conditions, the same three hosted models, and the two verdicts of the
first test. The questions ask for the highest candidate for a trait, one named target's rank and the
rate for that rank, whether and by how many genes the highest candidate is explained, the highest
candidate no earlier link explains, the two position lists for targets with no link, a target with
nothing to rank, an invented gene and an ambiguous trait word.

The expected answers come from `ranktables`, which computes every served ranking from the bundle's
tables by matrix code that does not call the question layer. Asked separately, the question layer
gives the same answer to all 35 questions of each set. The SQL condition is given those same tables.

**Development round** (typed and SQL only; reported, not read against the verdicts):

| model | typed | SQL |
|---|---|---|
| deepseek-v4-flash | 32 | 34 |
| glm-5.3-flash | 30 | 35 |
| mercury-2.5 | 30 | 34 |

Pooled, typed minus SQL was -0.105 [-0.210, -0.010] of a question: typed was worse.

- **All 13 typed misses were of two kinds**: the rank of one named target, and the rate for that rank.
  The models answered 11 of those 24 questions.
- **The tool that answers both was not called.** `why` gives a target with no link its own rank among
  the trait's candidates and the rate that goes with it. It was called in none of the 13 misses, and in
  8 of the 11 right answers. The models went to `predict` and looked for the target in its list, which shows
  the top of the ranking only, folds neighbouring genes into one row under the rank of its first gene,
  and was cut short by the harness's limit of 8,000 characters on a tool result.
- **The fault was in what the server says about its own tools**, not in the questions, the expected
  answers or the grading. Between the rounds the server's instructions and the descriptions of `why`
  and `predict` were changed to say where one target's rank is found, and a list of candidates now says
  that a row is a lead. Nothing in the test was changed. This is recorded in the protocol.

**Evaluation round**, on a set generated and frozen after those changes:

| model | typed | SQL | no tool |
|---|---|---|---|
| deepseek-v4-flash | 35 | 34 | 0 |
| glm-5.3-flash | 34 | 35 | 0 |
| mercury-2.5 | 35 | 32 | 0 |

- **Usable by a small model: met.** All three answer 28 or more of 35 through the typed questions, and
  no model in any condition answered a not-in-graph question as if the gene were known.
- **Better than SQL: not met.** Pooled over the three models, typed minus SQL is +0.029
  [-0.019, +0.086] of a question. Typed was right where SQL was wrong on 4 answers, and SQL where
  typed was wrong on 1.
- **The rank questions are now answered.** All 24 answers about one target's rank or its rate are right
  under typed, and `why` was called for every one of them.
- **Where typed went wrong.** One answer of 105: asked for the highest candidate that no earlier link
  explains, one model called `predict` eight times, said each result was cut short, and reported that
  it could not tell.
- **Where SQL went wrong.** 3 of its 4 misses are the ambiguous trait word: given "cardiomyopathy" or "colitis",
  two models found no trait of exactly that name and said the trait is not in the graph, where the
  graph holds several. The typed questions say a name is ambiguous and list what it could mean.
- **Without a tool the models score 0.** They answered "unknown" to 102 of the 105 questions. Nothing about a
  ranking can be answered from the rules text or from general knowledge.

The typed condition used 1.96 tool calls a question (`find` 106 times, `predict` 76 times, `why` 24 times) and SQL 2.99.
A typed conversation read 9,561 tokens on average against 5,268 for SQL; in the development round,
before the change, it was 15,691 against 4,816. A list of candidates repeats a sentence on
how often the rank held beside every row.

**What this shows and what it does not.**

- Three small models can use the rankings through the typed questions, once the server says which
  question places a single target. Before it said so they missed 13 of 24 questions about one target's
  rank or its rate. The evaluation result follows a change to the tool made on the development round,
  which is what that round is for.
- SQL over tables that already hold every rank does about as well. The typed questions' margin is in
  resolving names: every SQL miss but one is an ambiguous name.
- It does not show that a ranking is right. That is the prospective test.
- The questions were written by the author of the tool, and no person has used it.

**Cost.** 525 conversations over the two rounds cost 0.40 dollars against this test's cap of 0.60.
All usability testing to date has cost 1.12 dollars of the 1.50 stated.
