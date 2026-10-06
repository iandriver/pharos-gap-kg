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

### Open: which dates the prospective test uses

The archive holds database dumps dated July 2022, April 2023, August 2023, October 2023 and April
2024, and a May 2026 table of development levels. The publication date of a dump is not the date of
its contents: the April 2024 dump's own `dataset` table lists its newest source loads as September
2022. Before a split is frozen, `findings/SNAPSHOTS.md` must say how many links differ between each
pair of dumps. Candidate designs, in order of preference once that is known:

1. two dumps whose contents differ, scored on links present only in the later one;
2. within one dump, a split by the publication year of the evidence behind each link;
3. target level: from the April 2024 graph, which targets gained a qualifying ligand or moved up a
   development level by May 2026 (the release notes report 2,137 such moves).

## Schedule

| dates | work |
|---|---|
| 6 to 12 Oct | fetch and read the snapshots; measure how they differ; freeze the split and the pass rules |
| 13 to 26 Oct | baselines and one graph model on Pharos alone; first pruning rules |
| 27 Oct | progress report |
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
