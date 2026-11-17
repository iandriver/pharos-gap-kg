# Pharos gap graph

Bio-Agent Knowledge Graph Construction Challenge · Open-source track · Ian Driver · November 17, 2026

## What evidence is missing for this disease?

Open the Parkinson disease report. Pharos 3.19 ties 345 genes to Parkinson disease. For 97, another source or a drug mechanism ties the gene to the disease, yet the release contains no genome-wide association study (GWAS) link. For 251, it records neither a drug nor a compound active enough for the protein’s family. Each **typed gap** names the missing evidence and what would close it. TMEM175, for example, has ten studies behind its genetic link. Its compound gap would close with measured activity of 1 micromolar or better. [Try the live Parkinson report](https://pm4ie6tprh.us-east-1.awsapprunner.com/?id=MONDO%3A0005180&name=Parkinson%20disease&kind=trait&view=gaps) · [Watch the captioned walkthrough](https://github.com/iandriver/pharos-gap-kg/releases/download/submission-2026-11-17-r3/pharos-gap-kg-gaps-demo.mp4)

The graph also puts a scale on missing genetic evidence. The largest Parkinson study in the release had 20,184 cases. For a variant carried by 20% of chromosomes, it could detect an odds ratio of about 1.08. At the same ratio of controls to cases, detecting 1.05 would take about 55,000 cases under the graph’s power assumptions. Select another covered disease, such as type 2 diabetes, to generate the same report from that disease’s studies, genes and compounds. The “The graph behind this answer (KGX)” link exports the supporting nodes and edges. [Disease gap methods](../findings/DGAPS.md) · [Export a disease graph](ASK.md#the-subgraph-behind-an-answer)

## How a later study updates the answer

On October 29 we added a dated local source from the Alzheimer’s GWAS that Uffelmann and colleagues had published the day before. We mapped 181 of its 191 predicted effector genes across 122 loci. Of the 74 genes that Pharos tied to Alzheimer disease through other evidence but lacked a GWAS link for, the study named **ABI3, APP and PLCG2**, all at previously known loci. The local answer shows the publication beside the earlier gap. The public demo serves the original release, and the study’s licensed table stays local. [Original paper](https://www.nature.com/articles/s41586-026-11206-y) · [Recorded update](../findings/AD_STUDY.md#added-29-oct-2026-the-study-as-a-dated-source)

For another paper, give the local agent its workbook, DOI, date and disease. It proposes a sheet and a mapping. Review the contributing rows and unresolved identifiers, then decide whether to add the dated source. On the Alzheimer workbook this workflow recovered the same 181 gene–disease pairs as the hand-built update. [Import instructions](IMPORT_WORKFLOW.md) · [Workbook results](../findings/ADD.md)

## What the agent does with a scientist’s file

**Which genes does this refractive-error study add to Pharos?** Open the recorded Cheng study. A scientist supplied the workbook and a description of the paper. The model picked Supplementary Table 3 from 15 sheets. Open “2. The mapping” to see why the validator refused a proposed gene relationship whose subject was a protein. A second model supplied a mapping that passed. [Open the recorded import](https://pm4ie6tprh.us-east-1.awsapprunner.com/add?run=D)

The accepted mapping accounted for all 539 rows. With the full Pharos reference tables it yielded **500 distinct genes**. An installation with the downloadable bundle yields **499** because only the full lookup resolves the symbol `C14orf39`. The study names a nearby protein-coding gene at each locus. Of the 500 genes, 327 also appear among the GWAS Catalog’s 737 positional gene assignments for this study. The scientist left the table unjoined. [Study D results](../findings/ADD.md)

**Can it handle a different type of study?** Higginbotham and colleagues reported 271 candidate cerebrospinal-fluid protein biomarkers for Alzheimer disease. The same tools accepted a mapping that preserved the measured protein differences. The compact bundle matched 264 protein symbols. Two expanded to multiple targets, yielding 285 proposed links. The table stayed at preview for review of those matches. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7577712/) · [Protein panel test](../findings/ADD.md#november-17-test-with-a-protein-biomarker-panel)

The browser serves recorded table imports and live graph queries without an API key. To check a file of your own, use the local MCP tools. They inspect the sheet, preview what would join, and write the table beside the release after you instruct the client to join it. Each import records its source file, mapping and validator identity. [Install and try a local client](ASK.md#from-an-mcp-client) · [Import guide](IMPORT_WORKFLOW.md) · [Hosted MCP](https://pm4ie6tprh.us-east-1.awsapprunner.com/mcp)

## What the evaluations showed

| Question | Result and scope |
|---|---|
| Can a scientist use the agent? | Three hosted models each answered 34 or 35 of 35 tasks on four frozen sets through the typed questions. The typed interface beat SQL on two sets. An earlier set exposed faults later fixed. Independent researcher testing remains open. The CSF protein panel tested a different evidence type through the local import tools. [Usability](../findings/USABILITY.md) · [Human test](PERSON_TEST.md) · [Protein panel](../docs/data/csf_panel_preview_pharos319.json) |
| Does context improve held-out link prediction? | With the same graph model, LightGCN, mean average precision on post-2018 links rose from **0.0825** on Pharos links alone to **0.0866** with context. The paired gain was 0.0041 [0.0006, 0.0085]. The context was undated and could include later information. A simpler neighbour method scored **0.0969**. In the separate 2023–2026 test, the prespecified ranking that used gene position scored 0.1012 versus 0.0947 for neighbours. Its gain concerned previously known loci. [Matched model test](../findings/MODELS.md) · [Later-link test](../findings/PROSPECTIVE.md) |
| Does pruning work? | The literature rule met its prespecified thresholds on five of eleven sources. The rules withheld 7,248 rows. ClinGen judged 26 withheld rows, supporting 15, rating six limited and judging five against. Ten of the eleven withheld UniProt Disease rows it judged had support. The small, selectively curated subset cannot estimate how many of all 7,248 were wrong. The prediction test found no clear overall benefit from pruning. [Pruning](../findings/PRUNING2.md) · [ClinGen audit](../findings/CLINGEN.md) |
| Can a reviewer trace or restore an edge? | An ordinary edge names its source. GWAS edges name their studies and papers. Each pruned row survives as a typed gap with its rule and restoration condition. Changing the recorded evidence for 20 pairs triggered a six-stage rebuild whose files matched a full rebuild. [Graph](../findings/GRAPH.md) · [Maintenance](../findings/MAINTENANCE.md) |
| Does the export follow a standard? | The graph contains 59,071 nodes, 640,321 edges and 30,484 typed gaps. The remodelled edges pass the project’s Biolink 4.4.3 semantic checks. The reference KGX validator flagged identifier and category issues, including the extension used for gaps. [Standards report](../findings/STANDARDS.md) |

## Download and reproduce

- [Corrected source and installation instructions](https://github.com/iandriver/pharos-gap-kg/tree/submission-2026-11-17-r3) · [Original freeze](https://github.com/iandriver/pharos-gap-kg/tree/submission-2026-11-17-final)
- [KGX graph archive](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-graph-pharos319.tar.gz) · [Graph checksum](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-graph-pharos319.tar.gz.sha256)
- [Question-serving bundle](https://github.com/iandriver/pharos-gap-kg/releases/download/pharos319-bundle/pharos-gap-kg-bundle-pharos319.tar.gz) · [Deployed commit](https://pm4ie6tprh.us-east-1.awsapprunner.com/health) · [Public checks](https://github.com/iandriver/pharos-gap-kg/releases/download/submission-2026-11-17-r3/public-checks.json)
- [Protocol](PROTOCOL.md) · [Detailed results](PROGRESS_20261027.md) · [Reproduction history](../findings/REPRODUCE.md)

A clean build of the original November 17 freeze used the published source files and recorded service answers. All 20 stages completed. The three graph TSVs and the serving archive matched the release digests byte for byte. The graph archive SHA-256 is `4891f7c80671a5771348817af16a8c9a8bc2652e4792db14d5221e7474d8c72c`. The serving archive SHA-256 is `a6b1554c9cac320836652f1f198920c3013491ccc8081fbfd38ba112d8b0a3d3`. Repository code is Apache-2.0. Upstream data keep their source-specific terms.
