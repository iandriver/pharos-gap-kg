# Rebuilds from a clean clone

Whether someone with only the repository and the published files gets the same graph.

**Final submission, revision 3 (November 17).** A fresh checkout of `87e9309` extracted the 16 source tables, restored the recorded service answers and completed all 20 build stages. The three graph files matched the members of the public archive byte for byte. The rebuilt serving archive matched its download, and the checkout had no tracked changes. KGX 2.7.0 reproduced the counts and error classes listed below. The project’s module self-tests and 180 pytest tests passed. All 31 public checks passed against the deployed revision, including the Parkinson and type 2 diabetes reports, the Parkinson subgraph export, MCP, and both submission downloads. [Local checks](../docs/data/submission_r3_local_checks.json) · [Public checks](../docs/data/submission_r3_public_checks.json)

**November 17 submission.** A clean clone of `0b70bcf` (tag `submission-2026-11-17-r2`) completed all 20 build stages using cached copies of the published inputs and recorded service answers. All three graph files and the serving archive matched the published bytes. Git reported no tracked changes after the build. The wheel exposed 12 MCP tools outside the checkout and passed a Parkinson query and mapping preview. All 180 tests and the module self-tests passed. This run did not repeat model calls or prediction scoring. [Build record, hashes and installation checks](../docs/data/submission_r2_local_checks.json)

On the released KGX archive, KGX 2.7.0 read 59,071 nodes and 640,321 edges. It flagged 1,042 nodes and 1,062 subject–object pairs under its identifier and category checks. Its scope excludes the separate gaps file. All 26 public checks of the r2 service, tag, documents and archive downloads passed. Both records identify the artifact they checked. [KGX result](../docs/data/submission_kgx_20261117.json) · [Public checks](../docs/data/submission_r2_public_checks.json)

**20 Oct 2026, commit e4e581f,** after the remodel for Biolink 4.4.3 (`findings/GRAPH.md`). The same
run as on 18 Oct: a clean clone, installed as the README says, the published files read-only, the
network blocked after the install, and the build traced by `python -m pharos_gap_kg.maintain`. Of the
112 files the build reads or writes, 111 are byte for byte this checkout's, and the one that differs is
again the extraction manifest (`docs/data/reproduce_check_pharos319_e4e581f.json`). After the build,
git showed no change in the clone: every record the build writes under `docs/data` was rewritten byte
for byte, among them the remodelled graph's records and the release archive's digest. Each stage read
and wrote the files the committed map gives it.

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
