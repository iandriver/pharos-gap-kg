# Pharos gap graph

Bio-Agent Knowledge Graph Construction Challenge · Open-source track · Ian Driver · November 17, 2026

## Try it with a scientist’s question

**Which genes does this refractive-error study add to Pharos?** Open the recorded Cheng study. A scientist supplied the workbook and a description of the paper. The model picked Supplementary Table 3 from 15 sheets. Open “2. The mapping” to see why the validator refused a proposed gene relationship whose subject was a protein. A second model supplied a mapping that passed. The scientist could then inspect the rows and decide whether to join the study. [Open the walkthrough](https://pm4ie6tprh.us-east-1.awsapprunner.com/add?run=D) · [Captioned demo video](https://github.com/iandriver/pharos-gap-kg/releases/download/submission-2026-11-17-r2/pharos-gap-kg-demo.mp4)

The accepted mapping accounted for all 539 rows. With the full Pharos reference tables it yielded **500 distinct genes**. An installation with the downloadable bundle yields **499** because only the full lookup resolves the symbol `C14orf39`. The study names a nearby protein-coding gene at each locus. Of the 500 genes, 327 also appear among the GWAS Catalog’s 737 positional gene assignments for this study. The scientist left the table unjoined. [Study D results](../findings/ADD.md)

**Why does ACHE appear in an Alzheimer disease query?** The live graph traces donepezil to ACHE, its mechanism target, and to Alzheimer disease, its indication. Three more drugs follow this path. Pharos 3.19 contains no GWAS link for this gene–disease pair. The answer names the source of each kind of evidence and lets you download the supporting subgraph. [Ask the live question](https://pm4ie6tprh.us-east-1.awsapprunner.com/?target=ACHE&trait=Alzheimer%20disease)

The browser serves recorded table imports and live graph queries without an API key. To check a file of your own, use the local MCP tools. They inspect the sheet, preview what would join, and write the table beside the release after you instruct the client to join it. Each import records its source file, mapping and validator identity. [Install and try a local client](ASK.md#from-an-mcp-client) · [Import guide](IMPORT_WORKFLOW.md) · [Hosted MCP](https://pm4ie6tprh.us-east-1.awsapprunner.com/mcp)

## What the evaluations showed

| Question | Result and scope |
|---|---|
| Can a scientist use the agent? | Three hosted models each answered 34 or 35 of 35 tasks on four frozen sets through the typed questions. The typed interface beat SQL on two sets. An earlier set exposed faults later fixed. Independent researcher testing remains open. The ingest pipeline also exported 57,873 JAX/MGI gene–phenotype edges separately, an example beyond the external GWAS workbooks. [Usability](../findings/USABILITY.md) · [Human test](PERSON_TEST.md) · [Ingest](../findings/INGEST.md) |
| Does context improve held-out link prediction? | With the same graph model, LightGCN, mean average precision on post-2018 links rose from **0.0825** on Pharos links alone to **0.0866** with context. The paired gain was 0.0041 [0.0006, 0.0085]. The context was undated and could include later information. A simpler neighbour method scored **0.0969**. In the separate 2023–2026 test, the prespecified ranking that used gene position scored 0.1012 versus 0.0947 for neighbours. Its gain concerned previously known loci. [Matched model test](../findings/MODELS.md) · [Later-link test](../findings/PROSPECTIVE.md) |
| Does pruning work? | The literature rule met its prespecified thresholds on five of eleven sources. The rules withheld 7,248 rows. ClinGen judged 26 withheld rows, supporting 15, rating six limited and judging five against. Ten of the eleven withheld UniProt Disease rows it judged had support. The small, selectively curated subset cannot estimate how many of all 7,248 were wrong. The prediction test found no clear overall benefit from pruning. [Pruning](../findings/PRUNING2.md) · [ClinGen audit](../findings/CLINGEN.md) |
| Can a reviewer trace or restore an edge? | An ordinary edge names its source. GWAS edges name their studies and papers. Each pruned row survives as a typed gap with its rule and restoration condition. Changing the recorded evidence for 20 pairs triggered a six-stage rebuild whose files matched a full rebuild. [Graph](../findings/GRAPH.md) · [Maintenance](../findings/MAINTENANCE.md) |
| Does the export follow a standard? | The graph contains 59,071 nodes, 640,321 edges and 30,484 typed gaps. The remodelled edges pass the project’s Biolink 4.4.3 semantic checks. The reference KGX validator flagged identifier and category issues, including the extension used for gaps. [Standards report](../findings/STANDARDS.md) |

## Download and reproduce

- [Corrected source and installation instructions](https://github.com/iandriver/pharos-gap-kg/tree/submission-2026-11-17-r2) · [Original freeze](https://github.com/iandriver/pharos-gap-kg/tree/submission-2026-11-17-final)
- [KGX graph archive](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-graph-pharos319.tar.gz) · [Graph checksum](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-graph-pharos319.tar.gz.sha256)
- [Question-serving bundle](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-bundle-pharos319.tar.gz) · [Deployed commit](https://pm4ie6tprh.us-east-1.awsapprunner.com/health) · [Public checks](https://github.com/iandriver/pharos-gap-kg/releases/download/submission-2026-11-17-r2/public-checks.json)
- [Protocol](PROTOCOL.md) · [Detailed results](PROGRESS_20261027.md) · [Reproduction history](../findings/REPRODUCE.md)

A clean build of the original November 17 freeze used the published source files and recorded service answers. All 20 stages completed. The three graph TSVs and the serving archive matched the release digests byte for byte. The graph archive SHA-256 is `4891f7c80671a5771348817af16a8c9a8bc2652e4792db14d5221e7474d8c72c`. The serving archive SHA-256 is `a6b1554c9cac320836652f1f198920c3013491ccc8081fbfd38ba112d8b0a3d3`. Repository code is Apache-2.0. Upstream data keep their source-specific terms.
