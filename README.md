# pharos-gaps

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
  neighbour score is not shown to beat target degree. The 130,369 later links are 34,149 loci.

The plan is [`docs/PLAN.md`](docs/PLAN.md).

## Inputs

Public files only, each listed in `src/pharos_gaps/registry.py` with its size, publication date and
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
gaps), `data/eval/` (frozen splits and scores). `data/` is not committed; small summaries that the
documents quote are committed under `docs/data/`.

## Run

```bash
uv venv && uv pip install -e ".[test]"
python -m pharos_gaps.registry            # what is needed and what is on disk
python -m pharos_gaps.fetch pharos319     # ranged, resumable, refuses a file that changed
python -m pharos_gaps.sqldump extract pharos319 --tables tiga,tiga_provenance,protein
NCBI_EMAIL=you@example.org python -m pharos_gaps.pubyear pharos319
python -m pharos_gaps.split pharos319 --cutoff 2018          # rebuilds the split and checks its digests
python -m pharos_gaps.graph build pharos319 --cutoff 2018
python -m pharos_gaps.graph lineage pharos319 --cutoff 2018 --target APOE --trait MONDO:0004975
python -m pharos_gaps.baselines --cutoff 2018
python -m pharos_gaps.fetch efo_obo ensembl_gtf_113
python -m pharos_gaps.sqldump extract pharos319 --tables xref
python -m pharos_gaps.controls pharos319 --cutoff 2018
make check
```

## Layout

`src/pharos_gaps/paths.py` is the only module that names a location. `findings/` is an append-only
record of what was measured. `docs/` says how things work.

Apache-2.0.
