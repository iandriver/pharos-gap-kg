# pharos-gap-kg

A knowledge graph built by an agent from the [Pharos](https://pharos.nih.gov) 3.19 database, with each
GWAS link dated by its first study, for the Bio-Agent Knowledge Graph Construction Challenge (6 Oct to
17 Nov 2026).

Three properties are the point of it:

- **Every edge carries its lineage.** An edge id is a digest of what the edge claims; each edge names
  its source, and the graph records the published file and its SHA-256. A GWAS link names its studies
  and papers, an expression observation its row in the release, and a disease source's assertion the
  source's own evidence where it gives any.
- **Every pruned edge stays, as a typed gap.** Removing a link is a claim too. The graph records what
  was removed, by which rule, and what evidence would put it back. The rules and their conformance
  suite are [`kgx-gaps`](https://github.com/iandriver/kgx-gaps).
- **Predictions are scored on links that became known later.** A link is dated by the first paper it
  cites, and every score is reported beside a baseline that knows only how well-studied each node is.

## Try it

The questions read a bundle of files that each release carries as a 16 MB archive (20 MB
unpacked), so trying them needs neither the 2.5 GB of published files nor the build:

```bash
git clone https://github.com/iandriver/pharos-gap-kg && cd pharos-gap-kg
uv venv && source .venv/bin/activate
uv pip install -e ".[server]" -c constraints.txt
python -m pharos_gap_kg.release get pharos319          # downloads the bundle, checks every digest, unpacks it
python -m pharos_gap_kg.ask gaps "Parkinson disease"   # the genes tied to a disease, what each lacks, the study that would close it
python -m pharos_gap_kg.ask why TREM2 "Alzheimer disease" --as-of 2012
python -m pharos_gap_kg.server                         # the same eight questions over MCP, the Model Context Protocol
python -m pharos_gap_kg.desktop install                # or put them in Claude Desktop on macOS; quit it first (docs/ASK.md)
```

`get` refuses an archive, or a file in it, whose digest is not the one committed in
`docs/data/bundle_archive_pharos319.json`. Without `uv`, `python -m venv .venv` and `pip install` take
the same arguments. [`docs/ASK.md`](docs/ASK.md) describes the eight questions and how to add the
server to an MCP client. To build the bundle from the published files instead, see [Run](#run).

## Status

Started 6 Oct 2026. What exists, by the challenge's criteria. The progress report, [`docs/PROGRESS_20261027.md`](docs/PROGRESS_20261027.md), gives the detail, and `make check` compares each of its figures with the records:

- **The graph** ([`findings/GRAPH.md`](findings/GRAPH.md)). One KGX export of 640,321 edges and 30,484
  gaps over 59,071 nodes, with every input row an edge or a gap: 295,327 GWAS links between targets
  and traits, each dated by its first study; 154,363 assertions of a gene with a condition from nine
  disease sources, each naming the protein it was for, and the 7,248 rows a pruning rule withheld, as
  gaps; 20,153 edges that make each protein the product of its gene; 159,846 expression
  observations, kept apart from the assertions; 1,782 drugs with their targets and indications; and
  23,236 target-level gaps. Pharos's disease table has twelve sources: one, TIGA, repeats the GWAS
  data and is left out; eleven were tested for pruning; and two of those, DrugCentral Indication and
  Expression Atlas, are exported as paths through a drug and as observations. Interactions, pathways,
  GO terms and compound memberships are read as tables by the scorers and are not exported.
- **Gaps for a disease** ([`findings/DGAPS.md`](findings/DGAPS.md)). `gaps` lists the genes the graph
  ties to a disease and what each lacks, with what would close each gap: for a missing genetic
  association, the size of study that would find one. Over 513 disease traits, 92.8% of GWAS-linked
  genes share their locus with another gene carrying the same link, 77.3% of tied genes have no drug
  and no compound at their family's potency cutoff, and for 270 of 457 traits the largest study
  lacked the cases to find a common variant with an odds ratio of 1.2. For a Tdark target, Pharos's
  least-studied level, the distance to the next level ranks the targets Pharos 4.0.0 promoted with a
  ROC area of 0.859 ([`findings/TGAPS.md`](findings/TGAPS.md)).
- **Usability** ([`findings/USABILITY.md`](findings/USABILITY.md)). Eight typed questions, from the
  command line or an MCP client. On four frozen sets of 35 code-graded tasks, three small hosted
  models (mercury-2.5, deepseek-v4-flash and glm-5.3-flash, in that order) answered 35, 35 and 35;
  34, 35 and 35; 35, 35 and 34; and 34, 35 and 35 through the typed questions, against 31, 33 and 35;
  33, 33 and 34; 32, 34 and 35; and 28, 34 and 35 with SQL over the same tables. Typed beat SQL on two
  of the four sets. A first frozen set, before the faults it showed were fixed, gave 27, 27 and 26
  typed against 30, 27 and 30 with SQL. No model answered about an invented gene, drug or disease as
  if it were known. No researcher outside the project has used it yet; the test is written and frozen
  ([`docs/PERSON_TEST.md`](docs/PERSON_TEST.md)), and the builder takes it first on a separate set,
  reported apart.
- **Pruning** ([`findings/PRUNING.md`](findings/PRUNING.md), [`findings/PRUNING2.md`](findings/PRUNING2.md)).
  Each rule is scored on planted false pairs against a bar stated in advance: at least 80% of planted
  pairs withheld, at most 50% of one-source pairs. A cut on how much more often than chance PubMed
  mentions a gene with a disease passes on five of the eleven sources and is applied: 7,248 rows are
  withheld, 6,999 by the cut and 249 by in-graph rules on Monarch. On CTD, which holds 5,982 of them,
  it withholds 82.2% [78.2%, 85.7%] of 400 further planted pairs. A model that reads the abstracts
  does worse than the cut. Three sources are out of reach of every rule tried, full text included.
  Against ClinGen's expert curations ([`findings/CLINGEN.md`](findings/CLINGEN.md)), where they judge
  a kept row (6,562 rows) they support 86.7% and are against 4.6%, and they support the 26 withheld
  rows they judge less often, though not on UniProt Disease (10 of 11 supported). The rule counts a
  second source by identifier: 4,081 withheld rows have a kept row from another source, 4,059 of them
  DisGeNET's, on a condition of the same name.
- **Explainability** ([`docs/ASK.md`](docs/ASK.md)). An edge id is a digest of the claim, and each edge
  names its source. `why` gives the studies behind a
  link; where there is none, what lies nearby, what the trait's largest study could have found and
  the study that would close the gap; and for a withheld row, the rule, the count and what would
  restore it. `subgraph` returns the nodes and edges behind any answer as KGX, with the graph's own
  ids, and every row it can give equals the export's.
- **Standards** ([`findings/STANDARDS.md`](findings/STANDARDS.md)). BioLink categories and predicates,
  with a knowledge level and agent type on every edge. The export passes this entry's own `kgx-gaps`
  conformance suite (10 of 10). The reference KGX validator (kgx 2.7.0, Biolink 4.4.3) finds no
  unknown predicate and no missing property, and flags the identifiers of 1,044 of 59,071 nodes for
  four stated reasons. It does not check Biolink's own constraints on classes, domains and
  qualifiers; against those all 640,321 edges pass, after a remodel on 19 Oct that put the disease
  sources' assertions on gene nodes (614,400 of 620,168 edges failed before it).
- **Ingest by an agent** ([`findings/INGEST.md`](findings/INGEST.md)). A model proposes how a Pharos
  table the graph does not hold becomes Biolink edges, and a validator in code accepts or refuses it.
  Claude Haiku 5.5, Sonnet 5.5 and Opus 5.5 qualify on the nine tables the graph already maps; a local
  qwen3.6 35B does not. All three map three of four new tables (orthologs' diseases, OMIM and JAX/MGI
  phenotypes), which are exported apart and not yet joined; none maps the Gene Ontology annotations,
  which Biolink 4.4.3's domains refuse for genes and proteins.
- **The key question: link prediction** ([`findings/BASELINES.md`](findings/BASELINES.md),
  [`findings/MODELS.md`](findings/MODELS.md), [`findings/INPUTS.md`](findings/INPUTS.md),
  [`findings/PROSPECTIVE.md`](findings/PROSPECTIVE.md)). Scored as mean average precision over traits,
  on links reported after the training data. A neighbour score from Pharos links alone, which ranks a
  target by how similar a trait is to the traits the target is already linked to and has no fitted
  parameter, scores 0.0969 on links after 2018. No graph model beats it: matrix factorisation 0.0906,
  LightGCN 0.0825, LightGCN with interactions, pathways, GO terms, compounds and disease sources
  joined on 0.0866, and the same without the pruned rows 0.0864. The context adds to LightGCN
  (+0.0041 [+0.0006, +0.0085]), also on links neither trait lineage nor gene position explains, but
  it is undated, so that is an upper bound; pruning changes nothing (-0.0002 [-0.0007, +0.0002]).
  Blended into the neighbour score, each kind of context lowers it. Three small combiners of
  neighbours, degree and gene position beat it on links first reported from 2023; the one named in
  advance, which the questions rank by, scores 0.1012 against 0.0947 (+0.0066 [+0.0031, +0.0102]), at
  loci the training links already hold. On links nothing earlier explains it is level with neighbours
  (+0.0005 [-0.0006, +0.0012]), or behind on the reading that is a lower bound for it (-0.0060
  [-0.0081, -0.0043]). Cold
  start, rehearsal inside the training years and the controls are in
  [`findings/COLDSTART.md`](findings/COLDSTART.md), [`findings/REHEARSAL.md`](findings/REHEARSAL.md)
  and [`findings/CONTROLS.md`](findings/CONTROLS.md).
- **Keeping it current** ([`findings/MAINTENANCE.md`](findings/MAINTENANCE.md)). A traced map of what
  each build stage reads and writes. When the recorded PubMed answers for 20 pairs changed, 6 of the 20
  stages were stale and rebuilt, in 72 seconds, and every file came out byte for byte as a full build
  gives.
- **Reproducible** ([`findings/REPRODUCE.md`](findings/REPRODUCE.md)). A clean clone of commit
  e4e581f following the README rebuilt 111 of the 112 files the build touches byte for byte as on
  record, the other being this checkout's longer extraction manifest, and rewrote every committed
  record unchanged.
- **For other teams** ([`findings/SNAPSHOTS.md`](findings/SNAPSHOTS.md),
  [`findings/CONTEXT.md`](findings/CONTEXT.md)). The April 2023 and April 2024 Pharos databases hold
  identical link tables, so no later links come from comparing them. One Pharos disease source,
  JensenLab Experiment TIGA, is the GWAS data under another name, and used as context it reads the
  answers to a GWAS test.

The plan is [`docs/PLAN.md`](docs/PLAN.md).

## Inputs

Public files only, each listed in `src/pharos_gap_kg/registry.py` with its size, publication date and
fetch command. A missing input names the command that fetches it. The five the build reads come to
2.5 GB. They go in `data/raw`, which may be a symlink to a larger volume, or in the directory that
`PHAROS_GAP_KG_RAW` names.

| key | published | size | what |
|---|---|---|---|
| `pharos319` | 2024-04-29 | 2.26 GB | Pharos 3.19, the last full legacy database; the graph is built from it |
| `pharos400_canonicals` | 2026-05-28 | 45 MB | Pharos 4.0.0 development levels, for the target-level gaps' back-test |
| `efo_obo` | 2026-09-15 | 88 MB | trait ontology, for lineage and for the names a disease is searched under |
| `ensembl_gtf_113` | 2024-08-18 | 64 MB | gene coordinates, for position and loci |
| `gwas_catalog_20261005` | 2026-10-05 | 73 MB | the GWAS Catalog, for study sizes and the prospective test |
| `pharos316` | 2023-04-12 | 2.20 GB | optional: Pharos 3.16, only to show its link tables equal 3.19's |

Three stages ask a live service: NCBI for the year of each cited paper, and PubMed for co-mention
counts (33,879 searches). PubMed's answers move with time, so the answers the graph was built from are
committed under `recorded/` and restored in place of asking again.

## Outputs

`data/tables/<release>/<table>.parquet` (typed tables), `data/graph/<release>/` (the merged graph: KGX
nodes, edges and gaps; the three earlier exports sit beside it), `data/eval/` (frozen splits and scores), `data/bundle/<release>/` (what the questions read). `data/` is not committed; small summaries that the
documents quote are committed under `docs/data/`.

## Run

Install, with the library versions the committed records were written with (a parquet file carries
the version of the library that wrote it, so other versions give other digests):

```bash
uv venv && source .venv/bin/activate
uv pip install -e ".[test,server]" -c constraints.txt
```

Fetch, extract every table the build reads in one pass, and restore the recorded answers:

```bash
python -m pharos_gap_kg.registry            # what is needed and what is on disk
python -m pharos_gap_kg.fetch pharos319 pharos400_canonicals efo_obo ensembl_gtf_113 gwas_catalog_20261005
python -m pharos_gap_kg.sqldump extract pharos319 --tables tiga,tiga_provenance,gwas,protein,target,t2tc,xref,disease,ncats_ppi,pathway,goa,ncats_ligand_activity,drug_activity,tdl_info,generif,pmscore
python -m pharos_gap_kg.recorded restore pharos319     # PubMed and NCBI answers as they were asked
python -m pharos_gap_kg.pubyear pharos319              # the year of each cited paper; nothing to ask once restored
```

Build the graph and the bundle the questions read (about four minutes):

```bash
python -m pharos_gap_kg.split pharos319 --cutoff 2018           # rebuilds the split and checks its digests
python -m pharos_gap_kg.graph build pharos319 --cutoff 2018
python -m pharos_gap_kg.ctxtables pharos319                     # the rest of the release as link tables
python -m pharos_gap_kg.bundle pharos319                        # the dated links, which the next step reads
python -m pharos_gap_kg.planted pharos319                       # planted false pairs, for scoring pruning rules
python -m pharos_gap_kg.evidence pharos319                      # each source's own evidence field, per pair
python -m pharos_gap_kg.prune2 freeze pharos319                 # refuses unless it draws the committed sample
python -m pharos_gap_kg.prune2 ingraph pharos319
python -m pharos_gap_kg.applycut write pharos319                # the lift cut, applied from the recorded searches
python -m pharos_gap_kg.tgaps pharos319                         # what each target lacks for the next level
python -m pharos_gap_kg.tgapgraph pharos319                     # those gaps as KGX gap edges
python -m pharos_gap_kg.ctxgraph pharos319                      # the disease sources as KGX, withheld pairs as gaps
python -m pharos_gap_kg.kg build pharos319                      # ONE graph: dated links, sources, drugs, observations, every gap
python -I -m pharos_gap_kg.gwascat extract                      # the GWAS Catalog's columns, from the registered zip
python -m pharos_gap_kg.gwascat rule pharos319                  # how well the link rule reproduces Pharos
python -m pharos_gap_kg.gwascat build pharos319                 # links first reported from 2023; refuses unless it is the frozen split
python -I -m pharos_gap_kg.studysize build pharos319            # the size of each study, read from the GWAS Catalog file
python -m pharos_gap_kg.bundle pharos319                        # again, now with the gaps, the other parts and the study sizes
python -m pharos_gap_kg.diseasegaps census pharos319            # every disease trait, counted by kind of gap
make check
python -m pharos_gap_kg.release archive pharos319               # the bundle as the archive a release carries; its digest is committed
```

To keep it current ([`findings/MAINTENANCE.md`](findings/MAINTENANCE.md)): `python -m
pharos_gap_kg.maintain trace pharos319` runs the build block once more and maps what each stage reads
and writes; after an input changes, `maintain stale pharos319` says which stages that makes stale and
why, and `maintain rebuild pharos319` runs only those.

Ask:

```bash
python -m pharos_gap_kg.ask gaps "Parkinson disease"            # the genes tied to a disease, what each lacks, the study that would close it
python -m pharos_gap_kg.ask why ACHE "Alzheimer disease"        # no GWAS link; four sources and four drugs tie them
python -m pharos_gap_kg.ask why APOE "Alzheimer disease" --as-of 2010
python -m pharos_gap_kg.ask predict "type 2 diabetes mellitus" --as-of 2018 --unexplained-only
python -m pharos_gap_kg.ask predict "Parkinson disease" --unlinked          # targets with no link, by position
python -m pharos_gap_kg.server                                  # the same questions over MCP; docs/ASK.md
```

The evaluations. Each was scored once, and its figures are on record under `docs/data/`; running them
again reproduces those figures.

```bash
python -m pharos_gap_kg.graph lineage pharos319 --cutoff 2018 --target APOE --trait MONDO:0004975
python -m pharos_gap_kg.baselines --cutoff 2018
python -m pharos_gap_kg.controls pharos319 --cutoff 2018
python -m pharos_gap_kg.context pharos319 --cutoff 2018
python -m pharos_gap_kg.enrich pharos319 --cutoff 2018
python -m pharos_gap_kg.model tune pharos319 --cutoff 2018 --model mf    # fitted models, [model] extra; docs/BREV.md for a GPU
python -m pharos_gap_kg.combine fit pharos319 --cutoff 2018     # fitted weights over the baselines, lineage and position
python -m pharos_gap_kg.coldstart pharos319 --cutoff 2018       # a trait or a target with no link: lineage and position transfer
python -m pharos_gap_kg.rehearse report pharos319               # the procedure one step back in time; stages in its docstring
python -m pharos_gap_kg.prospect score pharos319                # the prospective test; stages in its docstring
python -m pharos_gap_kg.prospect calibrate pharos319            # how often a rank held, for the scorer the questions rank by
python -m pharos_gap_kg.calibrate pharos319 --cutoff 2018
python -m pharos_gap_kg.bundle pharos319 --backtest-only        # put the scorer and its calibration into a built bundle
```

The export against the Biolink model, with the reference KGX validator, offline, in an environment of
its own (`findings/STANDARDS.md`):

```bash
uv venv .venv-standards && uv pip install --python .venv-standards/bin/python kgx==2.7.0 biolink-model==4.4.3
.venv-standards/bin/python -I src/pharos_gap_kg/standards.py pharos319
```

The ingest stage (`findings/INGEST.md`). A hosted model needs `OPENROUTER_API_KEY`; the local one
needs `ollama serve`:

```bash
.venv-standards/bin/python -I src/pharos_gap_kg/biolinkslice.py    # the Biolink 4.4.3 slice the validator reads
python -m pharos_gap_kg.ingest cards pharos319                     # the cards, frozen before any model reads one
python -m pharos_gap_kg.ingest control pharos319                   # the graph's own mappings through the validator
python -m pharos_gap_kg.ingest propose pharos319 --model anthropic/claude-haiku-5.5 --set scored
python -m pharos_gap_kg.ingest report pharos319                    # who qualifies; then --set new for those that do
python -m pharos_gap_kg.ingest export pharos319                    # accepted new tables, under data/graph/pharos319/ingested
python -m pharos_gap_kg.ingest graphcheck pharos319                # Biolink's constraints on every kind of edge exported
```

To ask PubMed and NCBI again instead of restoring (about three hours at NCBI's rate; the counts will
have moved, so the withheld rows will differ from the record), set `NCBI_EMAIL` to your address and run
`python -m pharos_gap_kg.pubyear pharos319` after extracting, `python -m pharos_gap_kg.prune2 search
pharos319` after `prune2 freeze`, and `python -m pharos_gap_kg.applycut search pharos319` before
`applycut write`.
The pruning round's other stages (`read`, `choose`, `report`) call a model through OpenRouter and are
not needed to build the graph; they are listed in `prune2`'s docstring.

## Layout

`src/pharos_gap_kg/paths.py` is the only module that names a location. `findings/` is an append-only
record of what was measured. `docs/` says how things work.

Apache-2.0.
