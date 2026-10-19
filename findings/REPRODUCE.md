# Rebuilds from a clean clone

Whether someone with only the repository and the published files gets the same graph.

**18 Oct 2026, commit 8d413c3.** A clean clone, installed as the README says, read the published files
read-only and ran the README's extraction, `recorded restore` and build block, with the network
blocked after the install. `python -m pharos_gap_kg.maintain` traced the build and took the digest of
every file it reads or writes, 112 in all, and compared them with this checkout's
(`docs/data/reproduce_check_pharos319.json`). 111 are byte for byte the same. The one that differs is
the extraction manifest, `data/tables/pharos319/_manifest.json`: this checkout's lists 43 more tables,
extracted here while exploring, that no stage of the build reads. Every table the build reads is the
same.

**14 Oct 2026, commit d599dc6.** An agent given only the README, in a clean clone, followed it as
written with no file of the author's: install, extraction, `recorded restore`, build, the question
examples and the MCP server over stdio. Every step passed, and after the build and `make check`, git
showed no change: every record the build writes under `docs/data` was rewritten byte for byte.

**14 Oct 2026, commit 6fb3446,** the first such run, found what the two above rest on. Every rebuilt
table matched the author's frame for frame and the KGX files byte for byte, but the README ran stages
before their inputs, the base install lacked scipy, the withheld rows depended on 33,879 live PubMed
searches, and every parquet file's digest changed with the version of pyarrow that wrote it. The
answers are now committed under `recorded/`, the versions pinned in `constraints.txt`, and the README
reordered.

**What this does not show.** The two runs of 14 Oct are recorded in their logs and commit messages,
not in a committed record. The evaluations, which are scored once, and the stages that ask a live
service are outside the build block and were not rerun.
