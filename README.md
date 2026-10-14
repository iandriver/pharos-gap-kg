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
- **One graph** ([`findings/GRAPH.md`](findings/GRAPH.md)). Since 13 Oct 2026 everything is one KGX
  export with one node file: 620,168 edges and 30,484 gaps over 39,060 nodes, passing 10 of 10
  conformance checks with every input row an edge or a gap. It holds every dated GWAS link of the
  release, what nine disease sources assert (154,363 edges, and the 7,248 withheld rows as gaps),
  159,846 expression observations with direction, fold change and p-value, 1,782 drugs with their
  targets and indications, and the 23,236 target-level gaps. A DrugCentral indication row is a path
  through a drug, and 10,327 of the 13,919 run through the drug's mechanism target. 394 `same_as`
  edges tie a retired GWAS trait id to its current term, which joins 346 traits to the other sources.
- **Gaps for a disease** ([`findings/DGAPS.md`](findings/DGAPS.md)). `gaps` with a trait or a condition
  lists the genes the graph ties to it (a GWAS link, two or more disease sources, or the mechanism
  target of an indicated drug) and what they lack, in six lists, each with what would close the gap
  and a rate where one was measured. Over the 513 disease traits, 92.8% of GWAS-linked genes share
  their locus with another gene carrying the same link, and 77.3% of tied genes have no drug and no
  compound at their family's cutoff. An absent GWAS link comes with the smallest effect the largest
  study of the trait could have found (a median odds ratio of 1.26 where the study has cases and
  controls), read from 11,952 of the release's 12,335 studies. A gap is evidence the graph does not
  hold. It is not a finding that the gene matters in the disease.
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

- **Eight questions, and an MCP server** ([`docs/ASK.md`](docs/ASK.md)). `find`, `known`, `why`, `gaps`,
  `connect`, `predict`, `whats_new`, `status`, answered from a 20 MB bundle whose files are checked
  against a manifest of digests. `known`, `why` and `connect` read every part of the merged graph:
  beside the dated links they say which sources assert a pair, which rows are withheld and on what
  count, which drugs act on a target and treat the condition, and where expression differs, with
  observations kept apart from assertions. A name that matches nothing is "not in the graph"; an absent link is
  never reported as a negative; a prediction comes with the parts of its score and how often its rank
  held on later links. Since 11 Oct 2026 candidates are ranked by the position combiner of the
  prospective test below, and before that by the neighbour score alone.

- **A usability test** ([`findings/USABILITY.md`](findings/USABILITY.md)). On a fresh set of 35
  code-graded questions, three cheap hosted models each answered 35 through the typed questions, and
  31, 33 and 35 with raw SQL over the same tables. None invented an answer for a gene that is not in
  the graph. A first round failed both of its stated verdicts and exposed two faults in the test and
  one in the tools; both rounds are recorded. No person has tested it yet.
  The eighth question, `gaps`, was tested the same way on a fresh set: 35, 35 and 34 of 35 through the
  typed questions and 33, 34 and 33 with SQL, a difference that does not clear zero.
  The rankings were tested the same way on 12 Oct: 35, 35 and 34 of 35 through the typed questions,
  34, 32 and 35 with SQL (+0.029 [-0.019, +0.086] of a question, not shown better) and 0 with no tool.
  Its development round found that the server did not say which question gives one named target's
  rank; that was fixed before the evaluation questions were written.
  The other sources, the drugs and the gaps of a disease were tested the same way on 13 Oct: 34, 35
  and 35 of 35 through the typed questions, 28, 34 and 35 with SQL (+0.067 [+0.019, +0.124], six of
  the seven answers that differ being one model's) and 0, 0 and 1 with no tool. Every typed answer a
  model gave was right; the one miss is a question a provider refused. Its development round found
  that `find`, asked for a condition, did not look among diseases that are GWAS traits; that was
  fixed before the evaluation questions were written.

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
- **Target-level gaps** ([`findings/TGAPS.md`](findings/TGAPS.md)). Pharos grades targets by explicit
  evidence, so what a target lacks for the next level is on record, and so is how far it is from having
  it. For each of 5,499 Tdark targets the gap names the nearest criterion and a distance. That distance
  ranks the targets Pharos 4.0.0 promoted above the rest 85.9% of the time, better than how much a
  target is studied (+0.066 [+0.055, +0.077]). Of targets one criterion short and close to it, 53.6%
  were promoted, against 5.9% that the pipeline alone moves at that distance. Where a target was
  promoted through the counts, the gap had named the criterion that cleared 85.3% of the time.
  The gaps are in the graph as 23,236 gap edges, each with what would close it and what would retire it.
- **Fitted models** ([`findings/MODELS.md`](findings/MODELS.md)). Four declared models, tuned inside
  the training years and each scored on the test links once ([`docs/PROTOCOL.md`](docs/PROTOCOL.md),
  "Fitted models"). The two on Pharos links alone are scored, and neither beats the neighbour baseline's
  0.0969: matrix factorisation scores 0.0906 (-0.0063 [-0.0155, +0.0025]) and LightGCN 0.0825 (-0.0144
  [-0.0232, -0.0064]). On links that lineage and position do not explain both are worse than it. Counted
  by locus, matrix factorisation is ahead of it (0.0662 against 0.0512). The two models that join the
  context on are not fitted yet; [`docs/BREV.md`](docs/BREV.md) is the route to a rented GPU for them.
- **Trait lineage and gene position as inputs** ([`findings/INPUTS.md`](findings/INPUTS.md)). Added
  after the two models above were scored, and marked so. Four combiners of parameter-free scores with
  fitted weights and two graph models with neighbour and lineage edges, fitted inside the training
  years and each scored on the test links once. None is shown to beat the neighbours baseline: the
  highest, the combiner that adds gene position, scores 0.1037 (+0.0067 [-0.0009, +0.0151]). Each
  input does add to the model it is put into (position +0.0106 [+0.0031, +0.0189] and lineage +0.0080
  [+0.0042, +0.0129] in the combiner), on links that follow from an earlier link at the same locus or
  under a related trait name. Nothing here is a result about new links.
- **Cold start** ([`findings/COLDSTART.md`](findings/COLDSTART.md)). A third of the later links were
  never scored because their trait or target has no earlier link, and the neighbours baseline has no
  score for them. With no fitted parameter, ranking targets for a new trait by its ontology relatives'
  links, ties broken by degree, scores 0.1019 against 0.0610 for degree alone (+0.0408 [+0.0300,
  +0.0519] over 829 new traits). Ranking new targets for a trait by their chromosome neighbours scores
  0.0250 against 0.0102 (+0.0148 [+0.0074, +0.0243]). An earlier period agrees on both. The position
  ranking held again on the prospective split and is served by the question layer since 12 Oct, for
  the 2,407 targets with no link: a listed pair was linked within four years 2.8% of the time, 5.0% in
  the first five for a trait.
- **Rehearsal** ([`findings/REHEARSAL.md`](findings/REHEARSAL.md)). The tune, refit, score procedure
  run one step back in time, inside the training years, over 28 arms. One small combiner (neighbours,
  degree, position) carries, +0.0038 [+0.0008, +0.0077]; no graph or factorisation arm does, and
  neither carrying the length of fitting by steps nor drawing links by trait helps. Rehearsed, each
  scored model stands about where it stood on the test links.
- **The prospective test** ([`findings/PROSPECTIVE.md`](findings/PROSPECTIVE.md)). 194,017 links first
  reported from 2023 on, over 1,110 traits, built from the GWAS Catalog's release of 5 Oct 2026 by a
  rule checked against Pharos's own links (on shared studies 95.1% of the rule's triples are Pharos's
  and 92.8% of Pharos's are recovered), with all of Pharos 3.19 as training. What is scored was fixed
  before any scorer read a link. A combiner of neighbours, degree and position, with two fitted
  weights, beats the neighbours baseline: 0.1012 against 0.0947, +0.0066 [+0.0031, +0.0102]. It is the
  arm that carried in rehearsal. The gain is on links at loci that training already holds and not on
  new links. No graph model beats neighbours. For cold start, position transfer for new targets
  repeats (+0.0250 [+0.0087, +0.0418]) and lineage transfer for new traits does not (-0.0026 [-0.0040,
  -0.0011]). The question layer ranks by this combiner. Counted over every trait it ranks, 8.9% of
  its ten highest candidates for a trait were reported in 2023 to 2026, against a base rate of 0.49%;
  where no earlier link explains the candidate, 4.0%. (The 17.8% and 9.7% given here on 11 Oct were
  counted only over traits that had a later link.)

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

`data/tables/<release>/<table>.parquet` (typed tables), `data/graph/<release>/` (the merged graph: KGX
nodes, edges and gaps; the three earlier exports sit beside it), `data/eval/` (frozen splits and scores), `data/bundle/<release>/` (what the questions read). `data/` is not committed; small summaries that the
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
python -m pharos_gap_kg.combine fit pharos319 --cutoff 2018   # fitted weights over the baselines, lineage and position; [eval] extra
python -m pharos_gap_kg.coldstart pharos319 --cutoff 2018     # a trait or a target with no link: lineage and position transfer
python -m pharos_gap_kg.rehearse report pharos319             # the procedure one step back in time; stages in its docstring
python -m pharos_gap_kg.fetch gwas_catalog_20261005           # the GWAS Catalog as curated since, for the prospective test
python -I -m pharos_gap_kg.gwascat extract                    # then `rule` and `build`: links first reported from 2023
python -m pharos_gap_kg.prospect score pharos319              # the prospective test; stages in its docstring
python -m pharos_gap_kg.prospect calibrate pharos319          # how often a rank held, for the scorer the questions rank by
python -m pharos_gap_kg.bundle pharos319 --backtest-only      # put that scorer and its calibration into a built bundle
python -m pharos_gap_kg.tgaps pharos319                     # what each target lacks for the next level, against Pharos 4.0.0
python -m pharos_gap_kg.tgapgraph pharos319                 # those gaps as KGX gap edges
python -m pharos_gap_kg.bundle pharos319 --gaps-only        # and into the bundle the questions read
python -m pharos_gap_kg.ctxgraph pharos319                   # the disease sources as KGX, withheld pairs as gaps
python -m pharos_gap_kg.sqldump extract pharos319 --tables drug_activity
python -m pharos_gap_kg.kg build pharos319                   # ONE graph: dated links, sources, drugs, observations, every gap
python -m pharos_gap_kg.calibrate pharos319 --cutoff 2018
python -m pharos_gap_kg.sqldump extract pharos319 --tables target,t2tc
python -m pharos_gap_kg.bundle pharos319                     # what the questions read
python -m pharos_gap_kg.bundle pharos319 --context-only      # the merged graph's other parts, into a built bundle
python -I -m pharos_gap_kg.studysize build pharos319         # the size of each study, read from the GWAS Catalog file
python -m pharos_gap_kg.bundle pharos319 --sizes-only        # and into the bundle
python -m pharos_gap_kg.ask gaps "Parkinson disease"         # the genes tied to a disease and what each lacks
python -m pharos_gap_kg.diseasegaps census pharos319         # every disease trait, counted by kind of gap
python -m pharos_gap_kg.ask why ACHE "Alzheimer disease"     # no GWAS link; four sources and four drugs tie them
python -m pharos_gap_kg.ask why APOE "Alzheimer disease" --as-of 2010
python -m pharos_gap_kg.ask predict "type 2 diabetes mellitus" --as-of 2018 --unexplained-only
python -m pharos_gap_kg.ask predict "Parkinson disease" --unlinked          # targets with no link, by position
python -m pharos_gap_kg.server                               # the same questions over MCP (needs the [server] extra)
make check
```

## Layout

`src/pharos_gap_kg/paths.py` is the only module that names a location. `findings/` is an append-only
record of what was measured. `docs/` says how things work.

Apache-2.0.
