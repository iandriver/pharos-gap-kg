# Kept and withheld rows against ClinGen

26 Oct 2026. `python -m pharos_gap_kg.clingen audit pharos319` writes
`docs/data/clingen_audit_pharos319.json`. The file, the classes, the matching, the measures and how
the result would be read were committed to `docs/PROTOCOL.md` ("Kept and withheld rows against
ClinGen", f3921cd) before any row was matched to a classification.

ClinGen's gene-disease validity file as served on 26 Oct 2026 holds 3,715 curations of 3,066 genes,
which make 3,656 gene-disease pairs: 2,830 supported (a Definitive, Strong or Moderate curation), 535
limited, and 291 against (every curation Disputed, Refuted or No Known Disease Relationship). Of the
161,611 rows of the nine disease sources, 51,862 are in scope, on a gene ClinGen has curated and a
MONDO condition, and 1,503 of those are withheld. A row is judged when ClinGen curated its gene for
its MONDO id.

## The kept rows

ClinGen has judged 6,562 of the 50,359 kept rows in scope. It supports 86.7% [85.8%, 87.5%] of them
and is against 4.6% [4.1%, 5.1%]; the rest, 575, it rates limited.

| source | kept rows in scope | judged | supported | against |
|---|---:|---:|---|---|
| DisGeNET | 16,810 | 1,304 | 82.6% [80.4%, 84.5%] | 7.2% [5.9%, 8.7%] |
| UniProt Disease | 3,643 | 1,276 | 90.4% [88.7%, 91.9%] | 1.4% [0.9%, 2.2%] |
| Monarch | 5,285 | 1,249 | 89.5% [87.7%, 91.1%] | 3.3% [2.4%, 4.4%] |
| CTD | 10,193 | 1,064 | 89.6% [87.6%, 91.3%] | 3.2% [2.3%, 4.4%] |
| JensenLab Knowledge MedlinePlus | 2,361 | 581 | 85.7% [82.6%, 88.3%] | 4.3% [2.9%, 6.3%] |
| eRAM | 6,858 | 572 | 80.6% [77.1%, 83.6%] | 7.9% [5.9%, 10.4%] |
| JensenLab Text Mining | 3,156 | 316 | 91.1% [87.5%, 93.8%] | 6.0% [3.9%, 9.2%] |
| JensenLab Knowledge UniProtKB-KW | 1,998 | 195 | 69.2% [62.4%, 75.3%] | 11.8% [8.0%, 17.1%] |
| JensenLab Knowledge AmyCo | 55 | 5 | 80.0% [37.5%, 96.4%] | 0.0% [0.0%, 43.5%] |
| all nine | 50,359 | 6,562 | 86.7% [85.8%, 87.5%] | 4.6% [4.1%, 5.1%] |

The 299 kept rows ClinGen is against are published claims that its panels classed as disputed,
refuted or with no known disease relationship. CTD, for example, ties BRAF, KRAS and MAP2K1 to
Costello syndrome, MYH6 to hypertrophic cardiomyopathy, and KCNE2 to long QT syndrome. 141 of the 299
come from the five pruned sources. The rule did not test 126 of those, because another source
asserts the same pair and the rule tests only pairs one source asserts. It searched 14 and kept them,
with 1 to 64 papers naming the gene with the disease. The last, SCN9A with epilepsy from JensenLab
Knowledge UniProtKB-KW, is not among the pairs it searched. A second source
did not make these rows right: 126 of the 141 are asserted by two or more sources.

## Kept against withheld

ClinGen has judged 26 of the 1,503 withheld rows in scope, against 6,562 of 50,359 kept. Of the 26 it
supports 15, rates 6 limited and is against 5.

| pruned source | withheld rows in scope | judged: supported / limited / against | supported, kept minus withheld | secondary reading |
|---|---:|---|---|---|
| CTD | 946 | 1: 0 / 1 / 0 | +0.896 [+0.102, +0.913], too few for a verdict | +0.296 [+0.064, +0.584] on 10 |
| Monarch | 353 | 8: 4 / 1 / 3 | +0.395 [+0.110, +0.680], too few for a verdict | +0.378 [+0.218, +0.546] on 32 |
| JensenLab Knowledge UniProtKB-KW | 110 | 6: 1 / 3 / 2 | +0.526 [+0.123, +0.675], too few for a verdict | +0.550 [+0.255, +0.727] on 11 |
| JensenLab Knowledge MedlinePlus | 48 | 0 | none judged | none judged |
| UniProt Disease | 46 | 11: 10 / 1 / 0 | -0.005 [-0.081, +0.282] | -0.028 [-0.075, +0.175] on 18 |
| pooled, weighted by judged withheld rows | 1,503 | 26 | +0.275 [+0.127, +0.425] | +0.290 [+0.186, +0.390] |
| pooled without Monarch | | | +0.222 [+0.082, +0.356] | +0.218 [+0.094, +0.338] |

By the reading stated in advance, the pooled interval lies above 0: the rows the rule withholds are
ones ClinGen supports less often than the rows it keeps, with Monarch or without it. The figure rests
on 26 judged rows. UniProt Disease is the one source with enough for a verdict of its own, and there
no difference is detected: ClinGen supports 10 of the 11 rows the rule withheld from it. The
secondary reading, which also counts a disease within two is_a steps of a curated one, judges 71
withheld rows and points the same way; it decides nothing.

Confirmation, supported over all rows in scope, is 11.3% [11.0%, 11.6%] for kept rows and 1.0% [0.6%,
1.6%] for withheld ones. As the protocol says, ClinGen curates well-published pairs first and the rule
keeps well-published pairs, so this difference is expected whether or not a withheld row is wrong, and
it is not read as a verdict.

## What the rule withheld that ClinGen supports

15 withheld rows are pairs ClinGen supports: UniProt Disease 10, Monarch 4, JensenLab Knowledge
UniProtKB-KW 1. For 13 of them the PubMed search found no paper naming the gene with the disease, so
the lift was 0 and the cut withheld the row. They are Mendelian subtypes under numbered names, such as
CFAP43 with spermatogenic failure 19, KIF5A with amyotrophic lateral sclerosis 25, and OGT with
X-linked intellectual disability 106. The other two are Monarch rows withheld for having no second
source, PRKG1 with familial thoracic aortic aneurysm and dissection and TRAPPC9 with its intellectual
disability syndrome. The search found no paper for 3,672 of the 7,248 withheld rows, 745 of them in
scope.

## Limits

- 109,749 of the 161,611 rows are out of scope, and ClinGen curates mostly Mendelian genes, so these
  figures are agreement with expert curation where it exists, not the precision of the graph.
- 26 judged withheld rows carry the comparison, 11 of them from one source.
- Monarch's associations are assembled from curated sources that may include ClinGen, so its
  agreement is not independent of ClinGen; the pooled figures without Monarch are above.
- ClinGen's curations run to October 2026, after the release, so the rule could not have known them.
- Nothing in the graph changed on this result. Whether to restore the withheld rows ClinGen supports,
  or to mark the kept rows it is against, would be a new protocol entry.
