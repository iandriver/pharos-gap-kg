# pharos-gap-kg

A knowledge graph built by an agent from dated [Pharos](https://pharos.nih.gov) snapshots, for the
Bio-Agent Knowledge Graph Construction Challenge (6 Oct to 17 Nov 2026).

Three properties are the point of it:

- **Every edge carries its lineage.** An edge id is a digest of what the edge claims, and each edge
  resolves to the source row, the published file and that file's SHA-256.
- **Every pruned edge stays, as a typed gap.** Removing a link is a claim too. The graph records what
  was removed, by which rule, and what evidence would put it back. The rules and their conformance
  suite are [`kgx-gaps`](https://github.com/iandriver/kgx-gaps).
- **Predictions are scored on links that became known later.** A link is dated by the first paper it
  cites, and every score is reported beside a baseline that knows only how well-studied each node is.

## Status

Started 6 Oct 2026. No model has been fitted. What exists:

- **The data, read without a database server.** A registry with pinned digests, a resumable fetcher,
  and a streaming reader that turns a published MySQL dump into typed tables.
- **What the snapshots can date** ([`findings/SNAPSHOTS.md`](findings/SNAPSHOTS.md)). The April 2023
  and April 2024 databases hold identical link tables, so the archive cannot be split by publication
  date. GWAS target and trait links carry their own dates (295,327 pairs, each dated by its first
  study).
- **A frozen split and a protocol** ([`docs/PROTOCOL.md`](docs/PROTOCOL.md)). Links first reported by
  2018 train; 130,369 later links are scored. Each part is pinned by digest, and the protocol was
  committed before any scorer ran.
- **The graph as of 2018.** 20,412 targets, 1,404 traits, 99,435 edges carrying only the studies dated
  2018 or earlier. It passes the `kgx-gaps` conformance suite, and `graph lineage` answers why one
  edge exists, down to the papers and the digest of the file they were read from.
- **Baselines** ([`findings/BASELINES.md`](findings/BASELINES.md)). How well-studied a target was in
  2018 ranks later links at four times chance; a parameter-free neighbour score adds to that. A fitted
  model has to beat 0.0969.
- **Controls** ([`findings/CONTROLS.md`](findings/CONTROLS.md)). A fifth of the later links follow from
  an earlier link by trait lineage or by chromosomal position, and the ontology alone or the
  chromosome alone ranks later links at two and a half times chance. On links neither explains, the
  neighbour score is not shown to beat target degree. The 130,369 later links are 31,885 loci, and
  counted by locus target degree is at twice chance, not four times.
- **Context audit** ([`findings/CONTEXT.md`](findings/CONTEXT.md)). Nine tenths of the later links are
  to measurements, not diseases. Of Pharos's twelve undated disease sources, one is the held-out GWAS
  data under another name: used as context it would score 0.60 by reading the answers.

- **Seven questions, and an MCP server** ([`docs/ASK.md`](docs/ASK.md)). `find`, `known`, `why`,
  `connect`, `predict`, `whats_new`, `status`, answered from an 11 MB bundle whose files are checked
  against a manifest of digests. A name that matches nothing is "not in the graph"; an absent link is
  never reported as a negative; a prediction comes with how often its rank held in the backtest.

- **A usability test** ([`findings/USABILITY.md`](findings/USABILITY.md)). On a fresh set of 35
  code-graded questions, three cheap hosted models each answered 35 through the typed questions, and
  31, 33 and 35 with raw SQL over the same tables. None invented an answer for a gene that is not in
  the graph. A first round failed both of its stated verdicts and exposed two faults in the test and
  one in the tools; both rounds are recorded. No person has tested it yet.

- **Context** ([`findings/ENRICH.md`](findings/ENRICH.md)). Interactions, pathways, GO terms, compounds
  and eleven disease sources, each scored alone and blended with the neighbour score. Each carries
  some signal alone; every blend is worse than the neighbour score, and with the blending weight
  tuned on the test links the best weight is still none.
- **Pruning** ([`findings/PRUNING.md`](findings/PRUNING.md)). 16,250 false pairs were planted in the
  disease sources. The corroboration rule stated in advance withheld 39% of them and cleared its bar
  on one source of eleven; the criterion that fails is "an interaction partner has the same disease",
  which is true of most random pairs. The 508 pairs it did withhold are in the graph as typed gaps.
- **Pruning, second round** ([`findings/PRUNING2.md`](findings/PRUNING2.md)). Evidence from outside
  the graph, on 2,200 planted and 2,200 real pairs. A cut on how much more often than chance PubMed
  mentions a gene with a disease clears the bar on five sources of eleven and needs no model. A model
  that reads the abstracts and must copy the passage it relies on clears it on three, and at the same
  cost in real pairs it withholds fewer planted ones than the cut does. No source's own confidence
  field orders its pairs the way the literature does. Three sources are out of reach of every rule.
  The cut is applied to the five sources: 7,248 of 324,882 rows are withheld, each a typed gap that
  carries the papers counted, the lift and the cut. 5,982 of them are CTD's, where the cut costs 46.5%
  of one-source pairs against a bar of 50%.
- **Fitted models** ([`findings/MODELS.md`](findings/MODELS.md)) are defined and not yet scored on the
  test links: four declared models, tuned inside the training years, each scored once
  ([`docs/PROTOCOL.md`](docs/PROTOCOL.md), "Fitted models"). Two are tuned. On validation the best
  LightGCN setting cannot be told from the neighbour baseline (+0.0031 [-0.0058, +0.0120]) and matrix
  factorisation is worse than it. [`docs/BREV.md`](docs/BREV.md) is the route to a rented GPU for the
  two models that join the context on.

The plan is [`docs/PLAN.md`](docs/PLAN.md).

## Inputs

Public files only, each listed in `src/pharos_gap_kg/registry.py` with its size, publication date and
fetch command. A missing input names the command that fetches it.

| key | published | what |
|---|---|---|
| `pharos316` | 2023-04-12 | Pharos 3.16 database; same links as 3.19 |
| `pharos319` | 2024-04-29 | Pharos 3.19, the last full legacy database; the graph is built from it |
| `pharos400_canonicals` | 2026-05-28 | development levels for 20,654 canonical proteins |
| `efo_obo` | 2026-09-15 | trait ontology, for the lineage control |
| `ensembl_gtf_113` | 2024-08-18 | gene coordinates, for the position and locus controls |

## Outputs

`data/tables/<release>/<table>.parquet` (typed tables), `data/graph/<release>/` (KGX nodes, edges and
gaps), `data/eval/` (frozen splits and scores), `data/bundle/<release>/` (what the questions read). `data/` is not committed; small summaries that the
documents quote are committed under `docs/data/`.

## Run

```bash
uv venv && uv pip install -e ".[test]"
python -m pharos_gap_kg.registry            # what is needed and what is on disk
python -m pharos_gap_kg.fetch pharos319     # ranged, resumable, refuses a file that changed
python -m pharos_gap_kg.sqldump extract pharos319 --tables tiga,tiga_provenance,protein
NCBI_EMAIL=you@example.org python -m pharos_gap_kg.pubyear pharos319
python -m pharos_gap_kg.split pharos319 --cutoff 2018          # rebuilds the split and checks its digests
python -m pharos_gap_kg.graph build pharos319 --cutoff 2018
python -m pharos_gap_kg.graph lineage pharos319 --cutoff 2018 --target APOE --trait MONDO:0004975
python -m pharos_gap_kg.baselines --cutoff 2018
python -m pharos_gap_kg.fetch efo_obo ensembl_gtf_113
python -m pharos_gap_kg.sqldump extract pharos319 --tables xref
python -m pharos_gap_kg.controls pharos319 --cutoff 2018
python -m pharos_gap_kg.sqldump extract pharos319 --tables disease
python -m pharos_gap_kg.context pharos319 --cutoff 2018
python -m pharos_gap_kg.sqldump extract pharos319 --tables ncats_ppi,pathway,goa,ncats_ligand_activity
python -m pharos_gap_kg.ctxtables pharos319                  # the rest of the release as link tables
python -m pharos_gap_kg.enrich pharos319 --cutoff 2018       # needs the [eval] extra
python -m pharos_gap_kg.planted pharos319
python -m pharos_gap_kg.evidence pharos319                   # each source's own evidence field, per pair
python -m pharos_gap_kg.prune2 freeze pharos319              # second pruning round; stages in its docstring
python -m pharos_gap_kg.applycut search pharos319            # the lift cut on the sources it passed on (hours; resumes)
python -m pharos_gap_kg.applycut write pharos319
python -m pharos_gap_kg.model tune pharos319 --cutoff 2018 --model mf    # fitted models, [model] extra; docs/BREV.md for a GPU
python -m pharos_gap_kg.ctxgraph pharos319                   # the disease sources as KGX, withheld pairs as gaps
python -m pharos_gap_kg.calibrate pharos319 --cutoff 2018
python -m pharos_gap_kg.sqldump extract pharos319 --tables target,t2tc
python -m pharos_gap_kg.bundle pharos319                     # what the questions read
python -m pharos_gap_kg.ask why APOE "Alzheimer disease" --as-of 2010
python -m pharos_gap_kg.ask predict "type 2 diabetes mellitus" --as-of 2018 --unexplained-only
python -m pharos_gap_kg.server                               # the same questions over MCP (needs the [server] extra)
make check
```

## Layout

`src/pharos_gap_kg/paths.py` is the only module that names a location. `findings/` is an append-only
record of what was measured. `docs/` says how things work.

Apache-2.0.
