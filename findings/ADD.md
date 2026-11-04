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

## Added 4 Nov 2026: a way back to the sheet, study B again, and study C

The amendment of 4 Nov (98be7b8) was written after the result above, before the rerun and before
any model read study C's file. The records are `docs/data/add_kim2023_pd_supp_rerun_pharos319.json`,
`docs/data/add_mcgrail2026_t1d_supp_pharos319.json` and the parts `B_rerun` and `C` of
`docs/data/add_compare_pharos319.json`. Hosted spend for the whole entry is $0.139 of the $0.50 cap.

**The Catalog's genes.** The GWAS Catalog release of 5 Oct 2026 gives a reported gene for no
association published from 2022 on. Its positional mapping is there for study B. For study C it is
empty too: the Catalog gives C's 101 associations as positions, with no mapped gene, so neither
Catalog check can be made for C.

**Study B again, as development.** Given the mapping stage's reason, Haiku named "Supplementary
Table 3 - Lead SN" (the name as the file holds it, cut at 31 characters), whose rows are lead
variants with their nearest protein-coding gene. Code refused both of Haiku's mappings for it (a
predicate whose domain is a gene; then the Gene category, whose identifiers are not UniProtKB) and
Sonnet's first (279 of 303 rows resolve, 92.1%, under the 95% that check 8 requires). Sonnet's
second was accepted. It splits the rows by the method column, keeping the lead variants of the two
methods (random effects and MR-MEGA) and leaving out the legend of the columns written below them,
and it joins 191 links from 131 graph genes to Parkinson disease. Against the Catalog's positional
mapping (107 graph genes), 93 genes are shared: 71% of the table's genes [63%, 78%] and 87% of the
Catalog's [79%, 92%], partial by the stated reading. The sheet's 303 rows lie outside 39 to 156,
half to twice the abstract's 78 loci, so by the stated check it is not a table of the study's loci.
It was not joined.

**Study C, unseen.** Haiku named Supplementary Table 2, whose 102 rows are the study's loci, each
named by a gene in its "Signal" column. 102 lies within 80 to 320, so by the stated check it is a
table of the study's loci. Code refused both of Haiku's mappings (a predicate whose domain is a gene;
then an abstract predicate and a mixin category) and Sonnet's first, which read Signal by gene symbol
and was refused on check 8 alone: 95 of 102 rows (93.1%) resolve. Sonnet's second read Signal with
the resolver `name:gene`, which makes an ungrounded local node from the text and so resolves every
row. Code accepted it, because check 8 counts such a node as resolved. The mapping links 102 local
nodes to type 1 diabetes and no graph gene, so it joins nothing, and `add join` refuses it.

**What this shows.** The way back to the sheet turned study B's wrong sheet into one that names
genes. On study C the model met the check it had failed by choosing identifiers that cannot fail to
resolve: the refusal named the 95% bar, and the menu held a resolver that clears it without naming a
graph gene. The mapping code refused was the one that would have joined the study's genes.

**Limits.** Study B's rerun is development. Study C ran once. The two checks against the Catalog
could not be made for C, and the check from the abstract says only whether the row count fits.

## Added 4 Nov 2026, second: check 8 for an added table, study C again, and study D

The second amendment of 4 Nov (3f98ae2) was written after the results above, before the code for the new check
(70726c0), before the control ran under it and before any model read study D's file. The records are
`docs/data/add_mcgrail2026_t1d_supp_rerun_pharos319.json`, `docs/data/add_cheng2026_re_supp_pharos319.json`,
`docs/data/add_candidates_pharos319.json` and the parts `C_rerun`, `D` and `rechecked` of
`docs/data/add_compare_pharos319.json`. Hosted spend for the whole entry is $0.165 of the $0.50 cap.

**The new check.** A row of the table counts when it gives an edge between a target and a trait of the graph, and 80% of
each rule's rows must. The control passes under it: 122 of 127 rows count, and its join equals `later.py`'s 181 pairs.
Checked again under both rules, calling no model, 11 of the 15 mapping proposals the runs logged keep their verdict.
Four change:

- study B's rerun, Sonnet's first mapping (Table 3 by gene symbol, 279 of 303 rows): refused before, accepted now;
- study C, Sonnet's first (Signal by gene symbol, 95 of 102): refused before, accepted now;
- study C, Sonnet's second (`name:gene`): accepted before, refused now, with 0 of 102 rows giving an edge of the graph;
- study D, Sonnet's first (below, 508 of 539): accepted, and refused under the old rule.

Each logged verdict equals the recheck under the rule it was logged under.

**Study C again, as development.** Checked from the log, with no model called, Sonnet's first mapping is accepted and
links 95 graph genes to type 1 diabetes mellitus. Of the 7 rows that give no graph gene, five name a chromosome band
(14q32, 16q23, 17q24, 2p24, 6q27) and two a symbol the graph does not hold (FBXSIK68, THEM131L). The Catalog maps no
gene for the study, so neither Catalog check can be made. It was not joined.

**Study D, unseen.** `add candidates` ranks four Catalog studies, and the first, GCST90841196, is Cheng et al. From the
file's 15 sheets Haiku named "Supp Table 3", the study's loci with their nearest protein-coding gene, and its header row.
Its 539 rows lie within 466 to 1,864, so by the stated check it is a table of the study's loci. Code refused both of
Haiku's mappings (an abstract predicate; then a predicate whose domain is a gene) and accepted Sonnet's first. It reads
the nearest protein-coding gene column by gene symbol, links each gene to refractive error with
`biolink:genetically_associated_with` in a `biolink:Association`, carries the lead variant and its p-value, and keeps the
15 rows whose gene is "NA" as gaps. 508 of 539 rows count (94.2%), and it joins 500 links from 500 graph genes. The other
16 rows that give no graph gene name a symbol the rule does not match, among them LPHN3 and KIRREL, older names of ADGRL3
and KIRREL1: the rule reads current symbols only.

Against the Catalog's positional mapping, 327 genes are shared: 65% of the table's genes [61%, 69%] and 44% of the
Catalog's [41%, 48%], so the two disagree by the stated reading. Part of the difference is one of grain. The Catalog's
932 associations give 1,112 positional genes, 737 of them graph genes, against the table's 500 from 539 loci, so the
table could hold at most 68% of the Catalog's genes. The rest is loci where the table's nearest protein-coding gene is
not among the Catalog's genes: 173 of the table's 500. The mapping was not joined; it joins only on Ian's word.

**What this shows.** On a study no part of the entry had read, the stages ended on the loci table with a mapping that
reads its gene column, without going back to the sheet. The old check would have refused that mapping, as it refused the
same kind of mapping for B and C, and the new one refuses the ungrounded mapping C ended on.

**Limits.** Study D ran once and is one study. Check 8's bar was set after three tables had been seen, and D is the first
table it met unseen. The Catalog check reads positional genes, not what the paper named, and on D the two lists differ in
grain. The symbol rule has no table of previous symbols.
