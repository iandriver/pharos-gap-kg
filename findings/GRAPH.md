# One graph, and each source exported as what it asserts

13 Oct 2026. `python -m pharos_gap_kg.kg build pharos319`, then
`python -m pharos_gap_kg.bundle pharos319 --context-only`. The definitions are the two notes of 13 Oct
in `docs/PROTOCOL.md`, written before the build. The record is `docs/data/graph_pharos319.json`.
Nothing here is a pruning result or a prediction result.

## What was changed, and why

Three KGX exports existed, each with its own node file: the dated GWAS links as of 2018, the eleven
disease sources, and the target-level gaps. The dated links after 2018 were in no KGX file, and the
questions read only the dated links and the target gaps, so the disease sources and the withheld rows
could not be asked about. Two of the sources were not exported as what the release says they are.

| source | what the release holds | how it was exported | how it is exported now |
|---|---|---|---|
| DrugCentral Indication | 13,919 rows, every one naming a drug; each is a row of the drug-activity table joined to one of that drug's indications | a direct edge from target to disease | two edges through a drug node: drug to target, drug to indication |
| Expression Atlas | 159,846 rows, every one with a log2 fold change of at least 1.0 in size and a p-value below 0.05 | an edge with neither | an observation with a direction, the fold change, its floor and the p-value |

## What the build counted

**DrugCentral.** Of the 13,919 indication rows, 10,327 run through a target the release marks as the
drug's mechanism of action and 3,592 through another activity of the drug. 7,738 rows name the
indication without an identifier. The drug-activity table holds 1,782 drugs, of which 1,566 have an
indication, in 4,436 rows that are 4,295 pairs of a drug and a target. The indication rows are 5,943
pairs of a drug and an indication over 1,444 indications, of which 979 have no identifier and are
nodes under the local prefix, flagged as not grounded. Every indication row of the release is a path
of two edges in the graph, and the build refuses if one is not.

**Expression Atlas.** 79,696 rows are higher in the condition and 80,150 lower. The smallest fold
change is 1.0 and the largest p-value 0.049996. The source names 142 conditions, which are 107 by
identifier. The floor written on each edge is 0.999: the release holds no row below 1.0 at three
decimals, so a pair that is absent may have a fold change up to that size. It is the release's
inclusion rule and says nothing of an experiment's power.

## The merged graph

650,258 input rows became 620,168 edges and 30,484 gaps over 39,060 nodes, and the export passes 10
of 10 conformance checks with every input row an edge or a gap.

| part | edges | gaps |
|---|---|---|
| dated GWAS links, every year of the release | 295,327 | |
| assertions of a target with a condition, nine sources | 154,363 | 7,248 |
| expression observations | 159,846 | |
| drugs: 4,295 drug to target, 5,943 drug to indication | 10,238 | |
| a retired trait id is the same as its current term | 394 | |
| target-level gaps | | 23,236 |

The parts agree with the earlier records (`tests/test_merged_graph.py`): the dated links are the
bundle's 295,327 under the same edge ids; the nine sources' rows are the disease-source export less
the two remodelled sources, row for row; the 7,248 withheld rows are the ones the lift cut and the
in-graph rule withheld, by source and by rule; the target-level gaps are the 23,236 of their own
export. 141,553 of the 154,363 assertion edges carry the source's own evidence field.

## Joining the dated links to the rest

The first build showed that of the 2,247 GWAS traits only 41 share an identifier with a condition in
the other parts. The GWAS table names traits by EFO ids and the disease sources mostly by MONDO ids.
The ontology has retired 720 of the 2,247 trait ids and names a replacement for each.

- 394 `biolink:same_as` edges tie a retired trait id to its current term, on the ontology's own
  statement. With them 346 traits reach the other parts: 41 under their own id and 305 through their
  current term.
- 326 retired ids are replaced by a term in a vocabulary with no declared kind here (324 in OBA, the
  ontology of measured attributes, and 2 in GSSO) and are not tied. They are measurements.
- 13 traits are tied to a condition of the other parts only by a cross-reference. They are counted
  and not joined.

The questions go one step further than the graph, and say so: a condition is also read under an id
that the sources give the same name. 94 GWAS traits are read with such a condition, and 1,262
conditions of the other parts sit under more than one id with one name. That tie is a match on the
name and is not an edge in the graph. An answer lists every id it read and how each is tied.

## What the questions now answer

`known`, `why` and `connect` read every part, and `find` resolves conditions and drugs
(`docs/ASK.md`, "Beside the dated links"). Three examples from the release:

- ACHE and Alzheimer disease have no GWAS link in the release. Four sources assert the pair, and four
  drugs that act on ACHE as their mechanism target are indicated for it.
- MT-RNR1 and tinnitus: one source asserts the pair, and one CTD row is withheld by the lift cut with
  no paper mentioning both. The answer gives the rule, the count and what would restore the row.
- SH3D21 and malignant mesothelioma: one expression observation and nothing else, reported as an
  observation.

Naming conditions and drugs left every gene symbol but one resolving to its target. `HBD` is also,
letter for letter, the name of a condition, and now needs an id.

## What this does not show

- It prunes nothing. Expression Atlas, DrugCentral Indication and AmyCo remain out of reach of every
  pruning rule tested. Whether a target is the drug's mechanism of action is the release's own mark,
  and no planted-pair test of it has been run: a planted pair has no drug.
- The questions' new answers have not been put to a model or a person. The three usability tests on
  record are of the question layer before this change.
- The bundle the questions read grew from 12 MB to 20 MB.
- A tie by name can join two terms the ontology keeps apart where the release names them alike.
  MONDO:0007088 is read with Alzheimer disease for that reason, and the answer says it was.
- The conditions under UMLS and MESH ids are not mapped to MONDO, and the 326 measurement traits
  replaced by OBA terms are not tied to them.
- The other parts carry no dates. An answer as of an earlier year dates the GWAS links only.
- The graph as of 2018, the disease-source export and the target-gap export are unchanged. So are the
  pruning tests, the fitted models and every scored split.

Added 19 Oct 2026, the remodel for Biolink 4.4.3 (`docs/PROTOCOL.md`, amendments of 19 Oct;
`findings/STANDARDS.md`). The merged graph is now 640,321 edges and 30,484 gaps over 59,071 nodes. The
nine disease sources' 154,363 assertions and 7,248 withheld rows run from the gene (`NCBIGene`) of the
Pharos protein the row names, and name that protein as `pharosgapkg:gene_product`; 20,153 edges of a
new part, `genes`, make each protein the product of its gene (20,011 genes; 259 proteins have no gene
in the release). The 4 assertions on a protein with no gene stay on it, under `biolink:associated_with`.
Indications are `biolink:applied_to_treat`; the GWAS links are `biolink:Association`; expression
observations are `biolink:CorrelatedGeneToDiseaseAssociation`, with their direction as the property
`pharosgapkg:expression_direction`. The ids of the assertions, their gaps and the indications changed;
those of the links, expression observations, drug targets and target-level gaps did not. Compared with
the bundle of commit 054853a, every table the questions read is equal row for row but for those ids and
a new `gene` column, the census of disease gaps is unchanged, and fourteen typed answers are the same
but for edge ids (`docs/data/remodel_check_pharos319.json`).
