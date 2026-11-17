# pharos-gap-kg: final submission, November 17, 2026

Bio-Agent Knowledge Graph Construction Challenge, open-source track. Ian Driver.

## Open the demo

- [Browser demo](https://pm4ie6tprh.us-east-1.awsapprunner.com/): choose Parkinson disease to inspect known evidence and missing evidence.
- [Graph-building walkthrough](https://pm4ie6tprh.us-east-1.awsapprunner.com/add?run=D): inspect an unseen study's sheet selection, proposed mapping, validator feedback and resulting links.
- [MCP endpoint](https://pm4ie6tprh.us-east-1.awsapprunner.com/mcp): add this URL in a compatible MCP client. A browser GET is not an MCP session.
- [Service health and deployed commit](https://pm4ie6tprh.us-east-1.awsapprunner.com/health).

The public demo is read-only and needs no model/API key. The questions are live graph queries; the import walkthrough displays recorded runs, not a public file-upload service. A local MCP installation additionally provides `add_card`, `add_check` and `add_join` for a scientist's file. Models propose, code validates, and a join requires the person's instruction. The new local workflow reports unresolved source rows and preserves immutable import receipts. See [the import guide](IMPORT_WORKFLOW.md).

## Download and reproduce

- [Public repository](https://github.com/iandriver/pharos-gap-kg).
- [Frozen submission source](https://github.com/iandriver/pharos-gap-kg/tree/submission-2026-11-17).
- [KGX graph archive](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-graph-pharos319.tar.gz), with `nodes.tsv`, `edges.tsv`, `gaps.tsv` and a file-digest manifest.
- [Graph archive checksum](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-graph-pharos319.tar.gz.sha256).
- [Question-serving bundle](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-bundle-pharos319.tar.gz), for the local questions/MCP server.
- [Detailed results and progress report](PROGRESS_20261027.md), [scoring protocol](PROTOCOL.md), and [reproduction findings](../findings/REPRODUCE.md).

Graph archive SHA-256: `4891f7c80671a5771348817af16a8c9a8bc2652e4792db14d5221e7474d8c72c`.
Bundle archive SHA-256: `a6b1554c9cac320836652f1f198920c3013491ccc8081fbfd38ba112d8b0a3d3`.

The graph archive packages the existing measured graph; it does not rerun scientific evaluations. Its manifest records the packaging checkout and the exact file bytes. The submission tag freezes the code and write-up. The deployed image reports its exact code commit at `/health`. The bundle contains no API secrets or licensed later-study file. Repository code is Apache-2.0; upstream data provenance and source-specific terms remain applicable.

## What the agent builds

The entry ingests public Pharos tables, normalises identifiers, maps them to Biolink, applies declared pruning rules, exports a graph and serves bounded contextual questions. The graph has 640,321 edges and 30,484 typed gaps over 59,071 nodes. Disease/target gaps state what evidence is missing, the observation that would close the gap, and historical closure rates where measured. A missing relationship is not reported as a biological negative.

Every removed association remains a gap naming the rule and restoration condition. Every ordinary edge identifies its source; GWAS edges carry their studies/papers, and build records identify the source files and digests. A traced dependency map rebuilds only outputs affected by changed inputs; a controlled demonstration matched a full rebuild byte for byte.

The table-mapping agent was evaluated on separate study files, including unseen studies outside the original disease example. Schema validation and identifier coverage are measured separately from comparison with Catalog gene assignments. Paper-named genes and genes assigned by genomic position can differ; a valid mapping does not prove a causal gene.

## Evidence by judging criterion

| Criterion | Evidence | Limitation |
|---|---|---|
| Usability, 40% | Browser questions, MCP tools, contextual disease gaps and recorded import walkthrough; multiple frozen task sets answered by three hosted models | No completed independent human session is claimed. The hosted walkthrough is not an arbitrary-file importer |
| Accuracy and pruning, 30% | Dated pair-level prediction tests, degree/neighbour baselines, lineage/position/locus controls, planted-edge tests and a CTD extension | No graph neural model beats neighbours; pruning has no resolved prediction benefit; planted crossings are not independently adjudicated real errors |
| Explainability, 20% | Source tracing, typed gaps/restoration conditions, import diagnostics and receipts, changed-input rebuild demonstration | A gap in Pharos's evidence grade is not a functional or disease-causal conclusion |
| Standards, 10% | KGX export, Biolink 4.4.3 semantic checks after the remodel, reference-validator report and recorded clean-clone reproduction | Known identifier warnings and the gap-file extension remain explicitly documented |

The prespecified prospective combiner scores mean average precision 0.1012 versus neighbours 0.0947, paired difference +0.0066 [0.0031, 0.0102]. The gain concerns links at previously known loci. On candidates unexplained by earlier lineage/position evidence it is level or worse under the two protocol views. This is evidence for bounded ranking utility, not discovery of new biological mechanisms or an official competition score.

The submission includes the negative results. No new model was fitted or rescored to prepare this release. Fresh code tests and deployed-link checks are recorded separately from the historical clean-clone scientific build. A demo video script is provided in [VIDEO.md](VIDEO.md); no completed video or independent-user result is claimed.
