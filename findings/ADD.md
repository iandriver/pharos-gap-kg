# Adding a table: a model proposes, code checks, the person joins

3 Nov 2026. The stages, the two studies and their checks are in `docs/PROTOCOL.md`, "Adding a table"
(be71e37). Its amendment (00acb2a) was written after the code and the control ran and before any model
read a card. `python -I -m pharos_gap_kg.add` writes the records
`docs/data/add_uffelmann2026_supp_pharos319.json`, `docs/data/add_kim2023_pd_supp_pharos319.json`,
`docs/data/add_uffelmann2026_supp_control_pharos319.json` and `docs/data/add_compare_pharos319.json`.
Hosted spend was $0.036 of the $0.50 cap.

## The control

The hand mapping of study A, written in the proposal's form, is accepted on all nine checks. Its join
holds the same 181 pairs of a graph gene and a trait as the table `later.py` built by hand.

## Study A: reproduced

The model was given the whole file, 23 sheets, and the description. Claude Haiku 5.5's first proposal
named SD2 and its header row, row 2. Code refused both of Haiku's mappings: the first used
`biolink:genetic_association`, which is abstract, and the second `biolink:gene_associated_with_condition`,
whose domain is a gene and not the protein a graph gene is. Claude Sonnet 5.5's first mapping was
accepted. It maps the predicted effector gene column, split on "/" and read by gene symbol, to
Alzheimer disease, with `biolink:genetically_associated_with` in a `GeneToDiseaseAssociation`, and
carries nine of the table's columns. Its join holds 181 links from 181 graph genes, the same pairs as
the table built by hand. 5 of the 127 rows give no graph gene: FAM105B, HLA and three gene clusters.

## Study B: not added

Haiku named "Supplementary Table 5 - Genomic", one of 14 sheets. Its columns are a locus number, the
lead variant, the chromosome and the locus's coordinates; none names a gene. Code accepted the choice,
because the stage checks that a sheet and its header row can be read, not what they hold. The model's
own note on that proposal says the genes are in Table 3's "Nearest Protein Coding Gene" column. Given
Table 5, its mapping was that the table should not become edges, because "the table has no gene
column", and nothing joined.

The comparison with the GWAS Catalog cannot be computed. The Catalog release of 5 Oct 2026 gives no
reported gene for any of the study's 121 associations. It maps them by position to 138 genes, which is
a different field from the one the protocol named.

## The tools

`add_card`, `add_check` and `add_join` are served by the MCP server over stdio, on the person's own
machine; the server over the web serves none of them (`tests/test_server_http.py`). For each of the
six proposals the command made, `add_check` gives the command's verdict, word for word
(`tests/test_add.py`).

## What this shows

From the whole file, with no code written for that table, the stages reproduce a table built by hand.
One model chose the sheet and a second the mapping. On the unseen file the stages stopped at the
sheet: the model chose a sheet from which the description could not be met, and the stages as fixed
give no way back from the mapping to the sheet.

## Limits

- Two studies, one model choice of sheet each.
- Study A's table had been adapted by hand in this project before its description was written.
- Code checks a sheet for whether it can be read, not for whether it holds what the description asks.
- `find` matched no trait for a disease named with a possessive ("Parkinson's disease") until the
  amendment; the questions now match names with a possessive dropped.
