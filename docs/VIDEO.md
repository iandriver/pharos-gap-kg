# Demo walkthrough

[Watch the captioned screen walkthrough](https://github.com/iandriver/pharos-gap-kg/releases/download/submission-2026-11-17-r2/pharos-gap-kg-demo.mp4). For a live presentation, use the steps below. The import page displays the recorded agent decisions. The query page reads the serving bundle.

| Step | Show | Explain |
|---|---|---|
| Choose the table | Open [the refractive-error study](https://pm4ie6tprh.us-east-1.awsapprunner.com/add?run=D). Expand the 15-sheet list under “1. The sheet” and point to Supplementary Table 3. | The scientist supplied the paper’s workbook. The model had to find the table of loci before it could map a gene. |
| Read a refusal | Under “2. The mapping”, find the rejected protein subject for a predicate that requires a gene. Continue to the accepted proposal. | Identify the exact disagreement a scientist can inspect. The accepted mapping reads the nearest protein-coding gene column. |
| Decide whether to join | Under “3. The person”, read the preview, then the Catalog comparison. | The full reference tables yield 500 genes, the downloadable bundle 499. Only the full lookup resolves C14orf39. Of 500 genes, 327 overlap the Catalog’s 737 positional assignments. Scientific agreement is a separate check. The researcher left this study unjoined. |
| Follow a drug through the graph | Open [ACHE and Alzheimer disease](https://pm4ie6tprh.us-east-1.awsapprunner.com/?target=ACHE&trait=Alzheimer%20disease). Point to donepezil, its mechanism target and its indication. | Four indicated drugs act on ACHE. This release contains no GWAS link for the pair. The graph keeps those evidence types distinct. For the Parkinson example, the corresponding wording is “Pharos 3.19 contains no GWAS link for PARK7 and Parkinson disease.” |
| Answer the prediction question | Open [the matched model comparison](../findings/MODELS.md). | Adding context raised LightGCN’s mean average precision from 0.0825 to 0.0866. The undated context could include later information. Neighbours scored 0.0969. The separate later-link ranking gained at known loci. |

For a question about other data types, show the separate JAX/MGI phenotype export. For a question about human usability, show the frozen 12-task protocol and distinguish model scores from a researcher’s experience. [Ingest results](../findings/INGEST.md) · [Researcher session](PERSON_TEST.md)
