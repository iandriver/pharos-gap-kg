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
