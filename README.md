# pharos-gaps

A knowledge graph built by an agent from dated [Pharos](https://pharos.nih.gov) snapshots, for the
Bio-Agent Knowledge Graph Construction Challenge (6 Oct to 17 Nov 2026).

Three properties are the point of it:

- **Every edge carries its lineage.** An edge id is a digest of what the edge claims, and each edge
  resolves to the source row, the published file and that file's SHA-256.
- **Every pruned edge stays, as a typed gap.** Removing a link is a claim too. The graph records what
  was removed, by which rule, and what evidence would put it back. The rules and their conformance
  suite are [`kgx-gaps`](https://github.com/iandriver/kgx-gaps).
- **Predictions are scored on links that appeared later.** The graph is built from the April 2023
  snapshot and scored against April 2024, beside a baseline that knows only how well-studied each
  node is.

## Status

Started 6 Oct 2026. Built so far: the data registry, the resumable fetcher, and the streaming reader
that turns a published MySQL dump into typed tables. Nothing here is a result yet. The plan, with what
is decided and what is open, is [`docs/PLAN.md`](docs/PLAN.md).

## Inputs

Public files only, each listed in `src/pharos_gaps/registry.py` with its size, publication date and
fetch command. A missing input names the command that fetches it.

| key | published | what |
|---|---|---|
| `pharos316` | 2023-04-12 | the snapshot the graph is built from |
| `pharos319` | 2024-04-29 | the snapshot links are scored against |
| `pharos400_canonicals` | 2026-05-28 | development levels for 20,654 canonical proteins |

## Outputs

`data/tables/<release>/<table>.parquet` (typed tables), `data/graph/<release>/` (KGX nodes, edges and
gaps), `data/eval/` (frozen splits and scores). `data/` is not committed; small summaries that the
documents quote are committed under `docs/data/`.

## Run

```bash
uv venv && uv pip install -e ".[test]"
python -m pharos_gaps.registry            # what is needed and what is on disk
python -m pharos_gaps.fetch pharos319     # ranged, resumable, refuses a file that changed
make check
```

## Layout

`src/pharos_gaps/paths.py` is the only module that names a location. `findings/` is an append-only
record of what was measured. `docs/` says how things work.

Apache-2.0.
