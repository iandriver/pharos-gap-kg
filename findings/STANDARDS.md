# Standards

15 Oct 2026. `src/pharos_gap_kg/standards.py`, run in an environment of its own as its docstring
says. The record is `docs/data/standards_pharos319.json`.

The merged graph as exported (`data/graph/pharos319/nodes.tsv` and `edges.tsv`) was read by the
reference KGX validator, kgx 2.7.0, against the Biolink model 4.4.3. The model's schema and prefix map
were read from the installed `biolink-model` package, and every network connection was refused. kgx
checks each node and edge for its required properties, the types of its properties, the prefixes of
its identifiers against Biolink's prefix map, node categories, and edge predicates.

| | |
|---|---|
| nodes read | 39,060 |
| edges read | 620,168 |
| missing required properties, wrong property types, unknown predicates | none |
| nodes with an error | 1,042 (2.7%) |
| subject-object pairs with an error | 1,062 |

**No edge fails on its predicate or a required property, and no node on a required property.** Every
error is an identifier or category Biolink does not hold, from four causes:

- **979 indications that DrugCentral names without an identifier.** The graph gives each an id under
  its own prefix, `pharosgapkg`, which Biolink does not register. The 979 edges that point at them fail
  with them.
- **39 traits and conditions written `Orphanet:`,** where Biolink writes `orphanet:`. 30 are GWAS
  traits whose ids are those of the frozen split, as the GWAS Catalog writes them. Renaming them would
  change the split's digests and every edge id built on them, so they are kept. 60 subject-object pairs
  fail with them.
- **23 conditions named by JensenLab's AmyCo under `AmyCo:`,** a prefix Biolink does not register.
- **The three nodes that kgx-gaps adds for kinds of gap** (`pharosgapkg:GAP:...`). Their category,
  `pharosgapkg:KnowledgeGap`, is not a Biolink category, and kgx does not read their id as a CURIE.
  kgx reports each of these two errors once.

**What this does not show.**

- That the qualifier values are valid. kgx 2.7 checks the prefix of a subject and an object, not the
  value of a qualifier against Biolink's enumerations.
- That edges point at nodes, or anything about the gaps file. The gaps file is kgx-gaps' extension and
  not a KGX file. kgx-gaps' conformance suite checks both, and the graph passes its 10 checks.
- Conformance to a declared Biolink version. The export declares none; 4.4.3 is the version of the
  installed package.

Added 18 Oct 2026, a count. The nodes with a flagged identifier are 1,044: 979, 39 and 23 by prefix,
and the three nodes for kinds of gap. kgx reports the errors of those three once each, so its own
count, in the record, is 1,042.

Added 19 Oct 2026, what kgx does not check. Biolink 4.4.3's own constraints, checked by the ingest
stage's validator on every kind of edge in the export (`python -m pharos_gap_kg.ingest graphcheck
pharos319`, `docs/data/biolink_checks_pharos319.json`, `findings/INGEST.md`): 5,768 of the 620,168
edges are of a kind that passes, and 614,400 fail. The 294,248 GWAS links use
`biolink:GeneToDiseaseOrPhenotypicFeatureAssociation` and the 5,943 drug indications
`biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation`, neither of which is a 4.4.3 class, and 4.4.3
makes `biolink:treats` a mixin. The 154,363 disease-source assertions use
`biolink:gene_associated_with_condition`, whose domain is a gene, with a protein as the subject. The
159,846 expression observations carry `subject_direction_qualifier`, which
`biolink:GeneToDiseaseAssociation` does not carry. kgx 2.7.0 checks that a predicate and a node
category exist, and none of these.

Added 19 Oct 2026, after the remodel (`findings/GRAPH.md`). Checks 2 to 6 pass on every kind of edge:
the merged graph's 640,321 edges of 19 kinds, the disease sources' graph's 317,634 of 5 and the dated
graph of 2018's 99,435 of 4 (`docs/data/biolink_checks_pharos319.json` and its two siblings; the
figures before the remodel are in `docs/data/biolink_checks_pharos319_before_fix.json`). The graph's
own mappings pass the ingest validator on 5 of the 9 scored cards; the other four fail check 8 on
identifiers as the release gives them: DOIDs the release's own DOID table does not hold, AmyCo's prefix,
and GO, Orphanet and PATO traits under one category in the form the validator reads. The reference
validator, run again, reads 59,071 nodes and 640,321 edges and flags what it flagged before, 1,042
nodes by its count, for the same four causes; the gene nodes raise nothing.
