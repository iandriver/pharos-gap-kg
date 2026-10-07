# What the Pharos snapshots can date

6 Oct 2026. Every figure is printed by a committed command and stored under `docs/data/`.

| command | writes |
|---|---|
| `python -m pharos_gap_kg.sqldump inventory pharos316` (and `pharos319`) | `pharos316_tables.json`, `pharos319_tables.json` |
| `python -m pharos_gap_kg.snapshots pharos316 pharos319` | `snapshots_pharos316_pharos319.json` |
| `python -m pharos_gap_kg.levels pharos316 pharos319` | `levels_pharos316_pharos319.json` |
| `python -m pharos_gap_kg.levels pharos319` | `levels_pharos319_pharos400.json` |
| `python -m pharos_gap_kg.pubyear pharos319` then `python -m pharos_gap_kg.dated pharos319` | `dated_links_pharos319.json` |

Inputs: `pharos316.sql.gz`, `pharos319.sql.gz`, the Pharos 4.0.0 canonical table and its release notes.
`src/pharos_gap_kg/registry.py` pins the SHA-256 of all four.

## 1. The April 2023 and April 2024 databases hold the same links

The two dumps were published a year apart. They have 89 and 92 tables; 70 are identical by digest
(same rows, same ids, same order), and those 70 include every table a link would come from:
`disease`, `cmpd_activity`, `drug_activity`, `ncats_ligands`, `ncats_ligand_activity`, `ncats_ppi`,
`pathway`, `goa`, `phenotype`, `tiga`, `tiga_provenance`, `gwas`. Their `dbinfo` and `dataset` tables
are also identical: data version 6.13.4, newest source load 8 September 2022.

Nineteen tables differ. What moved in that year is literature and expression: `protein2pubmed`
33,283,707 to 43,188,393 rows, `generif` 861,171 to 888,732, `expression` 10,309,294 to 10,563,981.
Three tables exist only in the later dump (`go`, `go_parent`, `mjff_target`).

The development level changed for 339 of 20,412 targets: 223 rose from Tdark to Tbio, 60 fell from
Tbio to Tdark, and 56 fell from Tchem to Tbio. No target gained a ligand grade.

**A split that builds from April 2023 and scores on April 2024 has no link to score.**

## 2. In 3.19 the grade and the ligand table disagree for 56 targets

The ligand activity table names 2,675 targets. 2,619 targets are graded Tchem or Tclin. The other 56
have qualifying activity rows and are graded Tbio. They are the 56 that fell between 3.16 and 3.19,
while the ligand table did not change by a byte. No target holds the grade without rows.

So a label defined by grade would treat these 56 as targets without a ligand. Fourteen of them are
graded Tchem or Tclin in Pharos 4.0.0 and would have counted as gains.

## 3. From 3.19 to Pharos 4.0.0, 1,359 targets gained a ligand

20,191 accessions are in both; 221 only in 3.19; 463 only in 4.0.0, as the release notes say.

| 3.19 \ 4.0.0 | Tdark | Tbio | Tchem | Tclin |
|---|---|---|---|---|
| Tdark | 4,526 | 759 | 46 | 0 |
| Tbio | 118 | 10,797 | 1,317 | 10 |
| Tchem | 3 | 26 | 1,867 | 18 |
| Tclin | 0 | 0 | 0 | 704 |

2,150 targets rose and 147 fell. The release notes report 2,137 rises and no falls. Three of the five
upward counts agree exactly; Tdark to Tbio is 759 here against 743 in the notes, and Tbio to Tchem is
1,317 against 1,320. The difference is unexplained.

**The label.** Of 17,517 targets with no activity row in 3.19, 1,359 (7.8%) are graded Tchem or Tclin
in 4.0.0. The rate differs by family, counted over the 17,573 candidates by grade: epigenetic 24 of
147, kinase 26 of 186, enzyme 359 of 3,246, no family 869 of 11,546, transcription factor 40 of 1,353,
olfactory receptor 0 of 420.

Two limits on reading a rise as new knowledge:

- Pharos 4.0.0 is a rebuilt pipeline. 29 of 2,618 ligand holders lost the grade (1.1%), which shows
  the pipeline moves grades on its own. It does not bound how many rises are the pipeline's.
- The public 4.0.0 table carries counts, not activities, so it cannot say whether a gained ligand was
  published after 3.19's evidence ends. That needs the first qualifying activity for each of the
  1,359 targets from ChEMBL, with its year. Not done.

## 4. Links dated by the first paper they cite

3.19's links cite 29,936 distinct PubMed ids. NCBI returned a year for every one (1971 to 2023).

**Target and ligand.** 473,207 pairs over 2,675 targets and 355,932 compounds. 343,641 have a year,
334,788 have a year on every row, 129,566 have none. The last year is 2021.

| cutoff | train | test | both ends seen | targets with such a pair | targets first seen after |
|---|---|---|---|---|---|
| 2016 | 261,858 | 80,016 | 1,326 | 523 | 316 |
| 2018 | 297,329 | 45,491 | 632 | 271 | 141 |
| 2019 | 313,659 | 29,431 | 389 | 194 | 87 |

After any cutoff nearly every new pair involves a compound that did not exist in the graph before:
44,817 of 45,491 at 2018. Predicting these pairs from graph structure is not possible; it would
need the compound's chemistry.

**Gene and trait** (GWAS associations). 295,327 pairs over 18,005 genes and 2,247 traits. Every pair
is dated, 2005 to 2023.

| cutoff | train | test | both ends seen | genes with such a pair | genes first seen after |
|---|---|---|---|---|---|
| 2016 | 57,411 | 237,916 | 115,925 | 10,255 | 7,357 |
| 2018 | 99,435 | 195,892 | 130,369 | 13,731 | 3,521 |
| 2019 | 133,749 | 161,578 | 116,965 | 14,487 | 2,311 |
| 2020 | 197,727 | 97,600 | 72,208 | 13,951 | 910 |

## What follows for the test

- **Gene and trait links, split by study year,** are the only pair-level prospective test these files
  support: fully dated, and at a 2018 cutoff 130,369 later links join two ends already in the graph.
- **Ligand gain per target, 3.19 to 4.0.0,** is a second, target-level test: 17,517 candidates,
  1,359 positives, features strictly earlier than labels.
- Pair-level target and ligand prediction is dropped, and so is 3.16 against 3.19.

Three things the gene and trait test does not yet control, each a control to build before a model is
scored. GWAS assigns a signal to a gene by position, so neighbours in one region share links for a
reason unrelated to biology. The count of associations per year follows cohort size (4,101 in 2010,
63,978 in 2020), so a later link is not a random later link. And only these links carry dates: the
other tables of the graph (interactions, annotations, disease associations from other sources) are as
of 2022, so used as context they can hold what the held-out link later made known.
