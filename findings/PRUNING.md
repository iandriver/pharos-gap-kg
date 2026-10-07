# Planted false pairs: does a corroboration rule remove them

7 Oct 2026. `python -m pharos_gap_kg.planted pharos319` writes `docs/data/planted_pharos319.json`.
The planting, the rule, the measures and the bar were committed to `docs/PROTOCOL.md` before any pair
was planted.

The table is the 276,642 (target, disease) pairs of the eleven undated disease sources; 248,518 are
asserted by one source only. (The protocol gives 276,658 and 248,550, counted before `MIM:` and
`OMIM:` identifiers were written one way.) 16,250 false pairs were planted, 5% of each source's pairs,
by crossing two of the source's own pairs.

The rule keeps a pair that has a second source, a lineage relative for the same target, a GWAS link
in the lineage, or an interaction partner with the same disease, and withholds the rest.

## The rule does not work

| source | pairs | real, one source | planted | planted withheld | real one-source withheld | rule applied |
|---|---|---|---|---|---|---|
| Expression Atlas | 159,748 | 158,479 | 7,988 | 19.1% | 11.7% | no |
| DisGeNET | 82,453 | 59,466 | 4,123 | 55.6% | 27.5% | no |
| CTD | 34,797 | 12,860 | 1,740 | 52.2% | 23.7% | no |
| eRAM | 14,263 | 4,962 | 714 | 54.9% | 11.3% | no |
| Monarch | 9,424 | 1,692 | 472 | 86.0% | 30.0% | yes |
| JensenLab Text Mining | 7,394 | 4,994 | 370 | 53.2% | 8.1% | no |
| UniProt Disease | 5,549 | 401 | 278 | 91.4% | 58.1% | no |
| JensenLab Knowledge MedlinePlus | 4,107 | 972 | 206 | 74.3% | 23.7% | no |
| DrugCentral Indication | 3,523 | 3,145 | 177 | 57.1% | 19.2% | no |
| JensenLab Knowledge UniProtKB-KW | 3,364 | 1,329 | 169 | 54.4% | 15.7% | no |
| JensenLab Knowledge AmyCo | 260 | 218 | 13 | 23.1% | 24.3% | no |
| all sources | 276,642 | 248,518 | 16,250 | 38.9% | 16.4% | |

The bar was 80% of planted pairs withheld at a cost of at most 50% of real one-source pairs. One
source of eleven clears it. Over all sources the rule withholds 38.9% of the false pairs it was meant
to catch.

**Why.** Each criterion alone, over all sources, corroborates this share of planted and of real
one-source pairs:

| criterion | planted | real |
|---|---|---|
| an interaction partner has the same disease | 57.9% | 78.9% |
| the target has a disease in the same lineage | 14.1% | 29.6% |
| a GWAS link in the lineage | 1.7% | 2.2% |

The partner criterion corroborates more than half of pairs made at random. A target has a median of
22 interaction partners and a common disease is asserted for hundreds of targets, so some partner
nearly always shares it. That criterion is not evidence, and it is what lets the false pairs through.
The other two tell real from planted a little and apply to few pairs.

## A score does better than the rule, and not well enough

The neighbours score of a pair, computed on the table with the pair itself left out, ranks a real
one-source pair above a planted one 70.4% of the time over all sources, and 72% to 79% of the time in
eight of the eleven. Set to withhold 80% of the planted pairs, a cut on that score would also withhold
54.8% of the real one-source pairs over all sources (27% for Monarch, 37% for DisGeNET, 66% for
Expression Atlas). This figure was added after the first run, to say what any cut on the score would
cost.

**Expression Atlas pairs are close to indistinguishable from random crossings.** The score separates
them 60.3% of the time, the lowest of the eight sources with at least 1,000 one-source pairs, and the rule
withholds 19% of what was planted there. `findings/CONTEXT.md` reached the same place another way:
Expression Atlas pairs become later GWAS links at less than the base rate.

## A fault in the design

The protocol expected the product of a pair's two degrees to separate real from planted about half
the time, since a planted pair is built from the source's own targets and diseases. It does for the
sources whose pairs mostly have one source (0.50 to 0.59 for DisGeNET, CTD, eRAM, Expression Atlas).
It does not where most pairs have a second source: 0.14 for UniProt Disease and 0.33 for UniProtKB-KW.
Planted pairs were drawn from all of a source's pairs and compared with its one-source pairs, which
in those sources are the rarer diseases. For those sources the comparison is confounded by how common
the disease is, and their rows above should be read with that in mind.

## What was pruned

The rule was applied to Monarch, as the bar said: 508 real pairs are withheld and written out with
the rule's name and what would restore them. 324,374 of 324,882 rows are kept.

## What follows

- Corroboration from inside the graph cannot clean these sources at an acceptable cost. What tells a
  real pair from a planted one has to come from outside it: the source's own evidence fields, or the
  literature.
- The partner criterion should not be used as corroboration for anything.
- A second round would plant among one-source pairs only, so that the comparison is matched, and
  would test a rule that reads evidence the sources carry (a confidence, a z-score, a p-value) on a
  new seed.
