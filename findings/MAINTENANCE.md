# Keeping the graph current

18 Oct 2026. `python -m pharos_gap_kg.maintain` (`trace`, then `demo`), run in two clean clones of
commit 8d413c3 that read the published files read-only. The definitions and what counts were written
in `docs/PROTOCOL.md` ("Keeping the graph current") before it was run. The records are
`docs/data/build_graph_pharos319.json` (the map), `docs/data/maintenance_demo_pharos319.json` and
`docs/data/maintenance_check_pharos319.json`.

**The map.** The README's build block was run once with every file each stage opens recorded,
including the parquet files pyarrow reads and writes itself. Its 20 stages read and write 112 files.
The two clones traced the same map.

**The change.** PubMed asked again: for the 20 CTD one-source pairs outside the sample with the
smallest lift above CTD's cut, one fewer paper mentions both, with the lift and the cut's verdict
recomputed. 19 of the 20 fall to the cut.

**What it made stale, and what ran.** 6 of the 20 stages, each for the reason recorded:

| stage | why it was stale | seconds |
|---|---|---|
| `applycut write` | its input, the recorded PubMed answers, changed | 2.4 |
| `ctxgraph` | reads what `applycut` writes | 4.3 |
| `kg build` | reads what `applycut` writes | 15.6 |
| `bundle` | reads what `kg` writes | 12.8 |
| `diseasegaps census` | reads what `bundle` writes | 34.5 |
| `release archive` | reads what `bundle` writes | 1.0 |

The six ran in 72 seconds, where the full build takes about 3.5 minutes. The other 14 stages did not
run. 19 of the 112 files changed: the changed input, and 18 outputs of the six stages that ran (the
withheld and kept disease rows, 2 files; the context graph, 4; the merged graph's edges, gaps and
summary, 3; two bundle files and the archive, 3; and six committed records). The merged graph's node
file did not change. After the rebuild nothing was stale.

**Against a full build.** In the second clone the same change was made before a full build. All 112
files are byte for byte the same in the two clones.

**What this does not show.**

- A change to an input the map does not see. The map records what a stage opened in this build; a
  stage that would read another file under other data is caught when a rebuild runs it, which
  refuses a stage that reads a file the map does not know or writes other files than it gives.
- A changed published release. That passes through the registry, which pins each file's digest, so a
  new release is a new registry entry before it is an input here.
- The evaluations. They are scored once and are outside the build block.

A defect found on the first run is fixed in commit 8d413c3: the build writes the bundle twice, and
the rebuild compared the early stage's outputs with files the final stage had overwritten, so it ran
that stage without need.
