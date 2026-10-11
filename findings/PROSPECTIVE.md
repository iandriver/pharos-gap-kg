# The prospective test: links first reported from 2023

10 Oct 2026. The test, its link rule and what is scored on it were committed to `docs/PROTOCOL.md`
("A prospective test", and three dated additions of 10 Oct) before the steps they govern.

## The split

**The source.** The GWAS Catalog's release of 5 Oct 2026, registered and pinned by digest:
1,198,394 associations from 76,336 studies in 6,905 papers, published 2005 to 2026
(`python -I -m pharos_gap_kg.gwascat extract`). Pharos 3.19 holds GWAS links to 2022.

**The rule declared first failed.** Four rules were declared, built on the Catalog's own gene columns
and on trait ids as written. On the 12,260 studies that Pharos and the Catalog both hold, all published
by 2022, against Pharos's 720,702 (gene, trait, study) triples (`gwascat first-rule`,
`docs/data/gwascat_rule_first_pharos319.json`):

| rule | triples produced | of them Pharos's | of Pharos's, recovered |
|---|---|---|---|
| genes the variant lies in, p at most 5e-8 | 119,752 | 68.1% | 11.3% |
| genes the variant lies in, any p | 135,726 | 60.3% | 11.4% |
| those and the nearest gene each side, p at most 5e-8 | 175,035 | 67.1% | 16.3% |
| those and the nearest gene each side, any p | 199,731 | 59.0% | 16.4% |

The protocol had said that under 90% the test would go ahead and no further rule would be tried. That
sentence was set aside in a dated amendment, with this evidence (`gwascat look`,
`docs/data/gwascat_look_pharos319.json`), all of it from the same studies:

- **The Catalog has replaced half of its trait terms since Pharos's snapshot.** The terms differ for
  6,433 of the 12,260 studies, mostly an EFO term replaced by an OBA or MONDO term. As written, 50.0%
  of the Catalog's (trait, study) pairs are Pharos's. Through the ontology's live term, an obsolete
  term standing for its replacement, 97.1% are, and 98.7% of Pharos's are recovered.
- **Pharos's source took every gene near the variant, not the nearest.** Pharos lists 3 times as
  many genes for a study, at the median, as the Catalog's columns give. Taking every gene within a
  distance of a variant significant at 5e-8:

| distance from the variant | of the rule's (gene, study) pairs, Pharos's | of Pharos's, recovered |
|---|---|---|
| 0 bases | 98.5% | 16.6% |
| 10,000 bases | 98.4% | 24.8% |
| 25,000 bases | 98.4% | 35.2% |
| 50,000 bases | 98.3% | 50.7% |
| 100,000 bases | 98.1% | 79.3% |
| 250,000 bases | 50.9% | 81.7% |

  What the rule produces stays Pharos's up to 100,000 bases and stops being so beyond. Without the
  p-value filter, 87.7% of the pairs at 100,000 bases are Pharos's.

**The amended rule.** A link joins a trait, as the live term of its ontology id, to every target with
a gene within 100,000 bases of a variant associated with it at 5e-8 or smaller. At the level a link is
made, (target, live trait, study), on the shared studies (`gwascat rule`,
`docs/data/gwascat_rule_pharos319.json`): 95.1% of the rule's 636,104 triples are Pharos's, and 92.8% of
Pharos's 651,671 are recovered. For pairs it is 95.3% and 92.0%. Both are over the 90% asked for.
The amendment had expected about 77% recovered, a figure taken gene by gene that left out targets
Pharos ties only to gene ids on alternate assemblies; a dated note corrects it.

**The split** (`gwascat build`, `docs/data/split_gene_trait_2022.json`), frozen by code that reports
counts and runs no scorer:

| | pairs | |
|---|---|---|
| training | 294,954 | Pharos 3.19 by live trait term: 18,005 targets, 2,225 traits. As written they were 295,327 pairs over 2,247 traits |
| prospective test | 194,017 | first reported in 2023 or later, between a training target and a training trait: 1,110 traits, 17,456 targets |
| cold start | 394,383 | 377,156 to a new trait (8,425 new traits), 6,729 to a new target (1,529), 10,498 to both |
| set aside | 228,538 | pairs Pharos lacks that were first reported by 2022; neither positives nor candidates. 139,482 are between a training target and trait |

The test links by first year: 2023 34,067, 2024 58,594, 2025 85,222, 2026 16,134. The Catalog gives 93.3% of Pharos's
training pairs back under the rule.

- **The set-aside pairs are many.** They are pairs from studies published by 2022 that Pharos does not
  hold: the Catalog curated the study after Pharos's snapshot, or Pharos's source left the study out.
  Which, this protocol cannot see, and that is why they are counted neither way.
- **The last year is partial.** 2026 holds what had been curated by 5 Oct.
