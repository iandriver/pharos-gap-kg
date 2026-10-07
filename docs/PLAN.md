# Plan for the Bio-Agent KG Construction Challenge

Challenge window 6 Oct to 17 Nov 2026. Progress report by email on 27 Oct. Judging to 8 Dec, results
15 Dec. Open-source track. Solo entry.

## What is asked

An agent that takes public Pharos datasets and builds an enriched knowledge graph for link prediction,
pulling context from other public sources. The organisers' key question: does a graph built this way
predict held-out links better than an existing public graph such as Pharos, using graph neural
networks or embeddings. The agent is expected to run a loop: ingest, filter, clean, limit, maintain.

| criterion | weight | what it asks |
|---|---|---|
| Usability | 40% | can someone who is not a KG expert use it, and do multi-hop, contextual work |
| Accuracy and pruning | 30% | how well the agent removes wrong or contextually irrelevant links |
| Explainability | 20% | can it trace the lineage of one edge: why does A connect to B |
| Standards | 10% | nodes and properties map to existing ontology classes; specific; few classes |

## The entry

Given a table of measured relationships and a Pharos snapshot, the agent builds a KGX graph in which:

1. each edge id is a digest of the claim, and resolves to source row, published file and file digest;
2. each removed link stays as a typed gap naming the rule that removed it and what would restore it;
3. a small set of typed questions answers from the graph, and a language model only routes and
   renders. A gene or disease the graph has not seen returns "not in the graph".

Agents propose (a column mapping, a suspect edge, a literature record). A code check decides what
enters. Two measurements made before the challenge, in private projects and so not reproducible from
this repository, are why: a structured-decision model named a gap's type 16 of 50 times, and when the
right ontology term was absent from its candidates it still chose one in 12 of 29 cases. Both are to
be repeated here on public data before the entry leans on them.

## The loop, and what checks each step

| step | the agent does | the check that can refuse it |
|---|---|---|
| Ingest | proposes a mapping from a table's columns to BioLink categories and predicates | `kgx-gaps` conformance suite; ungrounded ids are refused, not forced |
| Filter | turns a question into a typed query and a subgraph | the query resolves its arguments against the graph |
| Clean | flags suspect edges; applies declared pruning rules | planted false edges: recall and the cost in true edges |
| Limit | nothing by default | |
| Maintain | re-derives what a changed input affects | digests: an output built from a replaced input is stale |

## The measurement

Three controls are fixed before any model is trained.

- **A baseline that knows only how well-studied each node is.** On graphs of this kind node degree
  alone ranks held-out links well. Every result is reported beside it.
- **Hold out the pair, not the row.** One target and disease pair reaches Pharos through several
  sources. Holding out one copy leaves the answer in the graph.
- **Score on what appeared later.** Links are held out by date, not at random.

Pruning is scored by planting known-false edges and counting how many the agent removes, and at what
cost in true ones. A model's score may flag an edge. Only a declared rule or a cited contradiction
removes one.

### The prospective test, as measured on 6 Oct

The archive's dumps cannot be split by publication date: the April 2023 and April 2024 databases hold
byte-identical link tables (`findings/SNAPSHOTS.md`). What the files do support:

1. **Gene and trait links, split by GWAS study year.** 295,327 pairs, every one dated. At a 2018 cutoff
   99,435 train, and 130,369 later links join a gene and a trait already in the graph. This is the
   pair-level test. The cutoff year is not frozen yet.
2. **Ligand gain per target, from the 3.19 database to Pharos 4.0.0.** 17,517 targets with no activity
   row; 1,359 hold a ligand grade in May 2026. Features are strictly earlier than labels.

Dropped: target and ligand pairs (after any cutoff nearly every new pair involves a compound the graph
has not seen), and 3.16 against 3.19.

To settle before a model is scored:

- GWAS assigns a signal to a gene by position, so neighbouring genes share links. The split, or the
  negatives, must not let a neighbour's link stand in for the held-out one.
- Associations per year follow cohort size (4,101 in 2010, 63,978 in 2020).
- Only these links carry dates. Every other table is as of 2022, so as context it can hold what the
  held-out link later made known. Either date the context or report the test with and without it.
  Measured 7 Oct for the disease table (`findings/CONTEXT.md`): the source "JensenLab Experiment TIGA"
  is the held-out data itself and is excluded from any context; the other eleven are undated.
- Whether each of the 1,359 ligand gains rests on evidence published after 3.19's ends, from ChEMBL.

## Schedule

| dates | work |
|---|---|
| 6 to 12 Oct | done 6 Oct: snapshots read and compared, split frozen at 2018, protocol committed, graph as of 2018 built, four baselines scored. 7 Oct: the lineage, position and locus controls defined and measured. Next: the context tables, each with a pruning rule and its planted-edge test |
| 13 to 26 Oct | baselines and one graph model on Pharos alone; first pruning rules |
| 27 Oct | progress report; the public repository `iandriver/pharos-gap-kg` is created with it |
| 27 Oct to 9 Nov | the typed questions as an MCP server; lineage; one enrichment arm; planted-edge test |
| 10 to 17 Nov | usability test, standards pass, reproduction from a clean clone, submission |

Another deadline-bearing project shares this machine and these weeks (22 Oct, 5 Nov). The enrichment
arm is the first thing cut if time runs short; the baselines, the pruning test and the typed
questions are not.

## Not built

A live Amazon Neptune service, GPU graph analytics, and retrieval-augmented generation over the graph.
None is scored and no compute is provided. The graph ships as KGX files with a loader script for
Neptune.

## Open questions for the organisers

- Which link type is held out, and who runs the scoring models?
- Is any compute provided during the challenge?
- How should a stored absence be written in BioLink? (For the Monarch mentors.)
- Is the Pharos 4.0.0 database or graph itself available to teams? The public release is a table of
  levels and counts; with its activities and disease links, 3.19 to 4.0.0 could be scored link by link.
