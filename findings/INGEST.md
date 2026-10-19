# Ingest: a model proposes the mapping for a table, and code decides

19 Oct 2026. `python -m pharos_gap_kg.ingest`, under `docs/PROTOCOL.md`, "Ingest", written before any
card was made, and its three amendments of the same day. The report is `docs/data/ingest_pharos319.json`;
the cards, every reply with the validator's verdict on it, the spend ledger and the exported views'
record are under `docs/data/ingest/`.

**The models, on the nine scored cards** (the twelve rules the graph's own code applies to them):

| model | cards accepted (first try) | columns as the hand mapping | predicate exact, or exact or close | cost | seconds a card |
|---|---|---|---|---|---|
| qwen3.6 35B, local | 3 (2) | 2 of 12 | 0, 9 of 12 | $0 | 127 |
| Claude Haiku 5.5 | 9 (4) | 12 of 12 | 9, 11 of 12 | $0.020 | 12 |
| Claude Sonnet 5.5 | 8 (4) | 12 of 12 | 10, 10 of 12 | $0.33 | 12 |
| Claude Opus 5.5 | 8 (2) | 11 of 12 | 8, 10 of 12 | $0.82 | 22 |

Haiku, Sonnet and Opus qualify for new tables (at least 8 cards accepted, columns on at least 11
rules, predicate on at least 9); the local model does not. Of its six refused proposals, four use
classes or predicates Biolink 4.4.3 does not have (such as `biolink:DiseaseToGeneOrGeneProductAssociation`),
one uses identifiers outside the category's prefixes, and one reply holds no JSON object. On 10 of the
12 rules its subject column is not the hand mapping's.

Where the three qualifying models differ from the hand mapping:

- **The subject of a disease source's row.** All three map it to a gene (`NCBIGene`, through the
  Pharos protein's gene), on every disease rule; the hand mapping gives the protein (`UniProtKB`),
  which the predicate's domain does not admit (below).
- **Knowledge level and agent type for DisGeNET, eRAM and Monarch.** All three say
  `knowledge_assertion`, and give a manual or automated agent (except Opus for Monarch); the hand
  mapping says `not_provided` for both, because these sources mix curated and mined rows and the
  release does not say which a row is.
- **One predicate each.** Haiku: `gene_associated_with_condition` for Expression Atlas, where the hand
  mapping says `correlated_with` for an observation of expression. Sonnet: `causes` for UniProt's
  disease rows. Opus: `occurs_together_in_literature_with` for JensenLab's text mining. Sonnet and
  Opus split the DrugCentral activities by action type and map most of them to `affects`.

**The control.** The hand mappings, written in the same form, pass the validator on 1 of the 9 scored
cards, `drug_activity`. On the other eight:

- the nine disease sources' `gene_associated_with_condition` has a gene as its domain, and the graph
  gives a protein as the subject (eight cards);
- Expression Atlas's `subject_direction_qualifier` is not a qualifier `GeneToDiseaseAssociation`
  carries;
- `tiga`'s association class, `GeneToDiseaseOrPhenotypicFeatureAssociation`, is not a class in
  Biolink 4.4.3;
- of the DOID identifiers the eRAM, Monarch and JensenLab rows use where they have no MONDO id, 12 of
  35, 1 of 2 and 7 of 148 are in the release's own DOID table, and 23 of the 75 AmyCo conditions are
  under a prefix Biolink does not register.

**The whole graph.** The same checks 2 to 6, run on every kind of edge the graph exports
(`ingest graphcheck`, `docs/data/biolink_checks_pharos319.json`): 5,768 of the 620,168 edges are of a
kind that passes, and 614,400 fail at least one check.

| edges | what fails |
|---|---|
| 294,248 | GWAS links: `GeneToDiseaseOrPhenotypicFeatureAssociation` is not a Biolink 4.4.3 class |
| 159,846 | expression observations: `GeneToDiseaseAssociation` does not carry `subject_direction_qualifier` |
| 154,363 | disease sources: `gene_associated_with_condition`'s domain is gene, and the subject is a protein |
| 5,943 | drug indications: `treats` is a mixin in 4.4.3, and `ChemicalToDiseaseOrPhenotypicFeatureAssociation` is not a class |

The reference validator, kgx 2.7.0, reported none of these. It checks that a predicate and a node
category exist; it does not check an edge's category, a predicate's domain and range, or qualifiers.

**New tables.** Four Pharos provenance views the graph does not hold. All three qualifying models had
their mappings for `ortholog_disease`, OMIM phenotypes and JAX/MGI phenotypes accepted, and none had
one for `goa` (Gene Ontology annotations) accepted.

| view | rows | edges | gaps | nodes | mapped as |
|---|---|---|---|---|---|
| `ortholog_disease` | 37,852 | 37,847 | 0 | 9,468 | gene, `gene_associated_with_condition`, disease; a prediction by a computational model |
| `phenotype`, OMIM | 14,147 | 3,710 | 10,083 | 17,057 | gene, `gene_associated_with_condition`, the phenotype's MIM; a row with no phenotype is a gap of its gene |
| `phenotype`, JAX/MGI | 58,398 | 57,873 | 0 | 10,189 | gene, `has_phenotype`, MP term |

The mappings exported are Haiku's, the cheapest qualifying model. They are under
`data/graph/pharos319/ingested/`, and nothing reads them until Ian decides that they join. The
reference validator finds no error in them except OMIM's node for the kind of gap, as in the whole
graph (`docs/data/ingest/standards_*_pharos319.json`). Where Sonnet and Opus differ from Haiku: for
orthologs, Sonnet's agent type and its knowledge source's spelling (`infores:monarch-initiative`), and
Opus's predicate (`associated_with`); for OMIM, Sonnet's predicate (`genetically_associated_with`) and
both models' pattern for the phenotype number and condition; for JAX/MGI, both models' agent type and
Opus's knowledge level.

- **OMIM holds one phenotype per row.** 1,281 of the 3,808 rows with a phenotype number name more
  than one, 6,141 in all; each accepted mapping exports one per row, 3,710 gene and phenotype pairs.
  The resolvers give one identifier per cell, and the card cuts each value to 80 characters, so its
  examples of `trait` show at most one phenotype; no model's notes mention rows with several, and the
  validator cannot see it. The view needs a resolver that gives every match before it joins.
- **`goa` cannot pass as Biolink 4.4.3 is written.** Sonnet's last proposal fails only on molecular
  function: `enables` has `PhysicalEntity` as its domain, and neither `Gene` nor `Protein` descends
  from it, while Biolink's `MacromolecularMachineToMolecularActivityAssociation` takes both as its
  subject. `participates_in`, which Haiku and Opus used for processes, has `Occurrent` as its domain.

**Spend.** $1.513 billed, and $0.0035 reserved for a call stopped before its cost was written: $1.517
against the cap of $1.75 (raised from $1.00 on Ian's word when Opus had 7 of its 9 scored cards). It
includes $0.017 on replies set aside (below) and $0.002 on a reply lost when the first run was stopped.

**Set aside before the results.** Two defects of the validator were fixed while the models ran. The
reply parser took an object nested in the proposal for the proposal, so Haiku's first six replies were
read as empty and refused with wrong reasons; they are under `docs/data/ingest/discarded_pharos319/`
and the run was restarted (commit 81dab73). Check 8 judged the objects of rows a proposal's condition
had made gaps; Haiku's `tiga` card, refused for that alone, was set aside and asked again (3ee489b).
`ingest stale` revalidates every recorded attempt with the current code and finds none other.

**What this does not show.**

- That agreement with the hand mapping is correctness. The control shows the hand mappings fail
  Biolink's own constraints more often than the models' proposals do.
- That an accepted mapping is right. OMIM's is accepted and exports 3,710 of the 6,141 phenotype
  numbers its rows name. A person reads each accepted mapping before it joins.
- That the knowledge sources exist: the infores registry is not checked offline.
- How the models do on another release or another source.
