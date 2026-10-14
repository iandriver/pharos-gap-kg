# Gaps for a disease

13 Oct 2026. `python -I -m pharos_gap_kg.studysize build pharos319`,
`python -m pharos_gap_kg.bundle pharos319 --sizes-only`, `python -m pharos_gap_kg.ask gaps "Parkinson disease"`
and `python -m pharos_gap_kg.diseasegaps census pharos319`. The definitions are in
`docs/PROTOCOL.md`, "Gaps for a disease", written before the question was built. The records are
`docs/data/study_sizes_pharos319.json` and `docs/data/disease_gaps_pharos319.json`. No parameter was
fitted and no test split was scored.

## What the question answers

`gaps` with a trait or a condition lists the genes the graph ties to it and what they lack. A gene is
tied by a dated GWAS link, by two or more disease sources, or by being the mechanism target of a drug
indicated for the condition. Six lists follow, each with a count, its strongest rows, what would close
the gap, and a rate where an earlier section measured one.

| list | closes with | rate quoted |
|---|---|---|
| `no_compound` | a compound at the family's potency cutoff | 8% of such targets held one in Pharos 4.0.0 |
| `dark` | the target's nearest route to Tbio | the share of its band promoted in Pharos 4.0.0 |
| `locus_not_resolved` | evidence for one gene of the locus | not measured |
| `no_genetic_support` | a genome-wide significant association at the gene | by source, 2.1% to 18.2% of named pairs were GWAS-linked later, against 0.6% of all pairs |
| `withheld` | what the row's rule says restores it | not measured |
| `observed_only` | an assertion or a genetic association | not measured |

## How large the studies were

The GWAS Catalog file gives each study's initial sample as a line of text. Of its 76,336 studies,
76,011 have a line with a number. Of the release's 12,335 studies, 12,260 are in the file and 11,952
were read (96.9%): 2,413 with cases and controls, at a median of 2,098 cases and 8,571 controls, and
9,539 with individuals only, at a median of 5,367.

Reading the largest sizes first showed two misreadings, and the reader was corrected before any
answer used a size. A remark in brackets held a year and a PubMed id that were being counted as
people, and "36,989 cases, controls" was being read as controls. What stands in brackets is no longer
read, and the first of the two words after a number decides. Two limits remain. A line that lists
cases of more than one condition is added up. A number that is not a count, such as the year in
"before the 2009 H1N1 pandemic", is read as people; 15 lines of the 76,011 hold something that looks
like a year.

## Every disease trait, counted once

The release holds 513 disease traits with a GWAS link. Asked one at a time:

| | summed over the 513 traits | traits with at least one |
|---|---|---|
| genes tied | 47,610 | |
| of them GWAS-linked | 39,271 | |
| `no_compound` | 36,823 | 507 |
| `dark` | 5,909 | 370 |
| GWAS-linked genes in a locus with another linked gene | 36,452, in 7,042 loci | 461 |
| `no_genetic_support` | 8,339 | 314 |
| `withheld` rows | 1,529 | 122 |
| `observed_only` | 56,220 | 51 |

A gene tied to two traits is counted under each. The median trait has 25 tied genes, 13 of them
GWAS-linked.

- **92.8% of GWAS-linked genes share their locus with another gene that carries the same link.** The
  release's source assigns a signal to every gene within 100,000 bases, so a link names a place more
  often than a gene. Other evidence in the graph singles out one gene in 647 of the 7,042 loci (9.2%).
- **77.3% of tied genes have no drug and no compound at their family's cutoff.**
- **What an absent link could have seen varies widely.** The largest study behind a trait has a size
  for all 513 traits, and for 457 it has cases and controls. For those, the smallest odds ratio a
  variant carried by 20% of chromosomes needed to be found has a median of 1.26 (quartiles 1.12 and
  1.55), and for 124 traits it is above 1.5.

## One disease

For Parkinson disease the graph ties 345 genes: 248 by a GWAS link, 88 by two or more sources and 24
as the mechanism target of an indicated drug. 251 have no drug and no compound at their family's
cutoff, and 22 are still Tdark. 228 of the 248 GWAS-linked genes sit in 39 loci that carry the link
on more than one gene, and other evidence singles out one gene in 6 of them. 97 genes are tied with
no GWAS link; the first two listed, ATP13A2 and PARK7, are each asserted by six sources and rank
2,565 and 8,820 among the GWAS candidates. The largest study behind the trait's links had 20,184
cases and 397,324 controls, so a variant carried by 20% of chromosomes needed an odds ratio of about
1.08 to be found.

## Checks

`tests/test_disease_gaps.py` recounts every list for seven traits, one of them at an earlier year,
and for one condition that is not a GWAS trait, from the bundle's tables without calling the question
layer, and finds loci a second way. The first run disagreed on one list: `observed_only` counted
Pharos entries where every other list counts genes. It now counts genes, and the two agree on all
eight.

## What this does not show

- That a listed gene matters in the disease, or that a gap is worth closing. A gap is evidence the
  graph does not hold, and every answer says so.
- That a gene tied by many sources and lacking a GWAS link is a false lead. The floor is for a
  variant carried by 20% of chromosomes; a gene whose variants are rare is not visible to it.
- The rates. Each was measured by an earlier section on what the graph's own sources later held, and
  the sources other than GWAS carry no dates. `locus_not_resolved`, `withheld` and `observed_only`
  have no measured rate.
- The threshold of two sources is a choice, stated in the protocol and not tuned.
- The count of unresolved loci follows from the 100,000-base rule of the release's source and from
  the 500,000-base neighbour relation used here. Another window gives another count.
- The floor is rough: one variant frequency, an additive model, no covariates, and only studies that
  reported some association for the trait.
- The question has not been put to a model or a person.

## Added 14 Oct 2026: the study each gap calls for

`python -m pharos_gap_kg.diseasegaps census pharos319`. The definitions are in `docs/PROTOCOL.md`,
"The study each gap calls for", written before it was built. The record is
`docs/data/disease_gaps_pharos319.json`. Nothing was fitted and no test split was scored. The lists,
their members and their counts are as above.

Two lists and `why` now end in a stated study.

- `no_genetic_support` gives the cases a study of the trait needs to find a common variant (20% of
  chromosomes, p = 5e-8, 80% power) at odds ratios of 1.05, 1.1 and 1.2, with controls at the ratio
  of the trait's largest study, and marks which of the three that study already had the cases for.
  It also gives the cases for a gene-based test of rare variants (1 in 1,000 chromosomes combined,
  odds ratio 3, p = 2.5e-6), which needs sequence data and which a genotyped GWAS does not run.
- `locus_not_resolved` lists each locus's genes in the order of what the graph holds for them besides
  the link (a tie, one source, an expression observation, nothing), names what each holds, and gives
  the number of genes a study would test.
- `why` with no link gives the trait's floor and the same study.

Over the 457 disease traits whose largest study has cases and controls (median 5.23 controls per case):

| odds ratio | median cases needed | traits whose largest study had them |
|---|---|---|
| 1.05 | 61,921 | 28 (6.1%) |
| 1.1 | 16,227 | 93 (20.4%) |
| 1.2 | 4,435 | 187 (40.9%) |

The rare-variant test needs a median of 15,213 cases.

- **For 270 of the 457 traits (59.1%), the largest study lacked the cases to find a common variant
  with an odds ratio of 1.2.** For those traits an absent link does not yet weigh against an effect
  of that size.
- **Of the 6,395 loci that no tie singles out, one source or an expression observation points to a
  gene in 689 (10.8%).** In the other 5,706, no gene holds anything besides the link, or more than
  one gene holds a tie. Such a locus holds a median of 3 genes, 32,575 in all.

For Parkinson disease the largest study (20,184 cases, 397,324 controls) had the cases to find a
common variant at 1.1 or 1.2. At 1.05 it would take about 55,000 cases at its ratio, and a
rare-variant test about 13,000. Of the 33 loci no tie singles out (165 genes, a median of 4), one
source or an expression observation points to a gene in 3, among them APOE and SLC41A1 (one source
each). In two more, two genes hold ties equally: GAK and TMEM175, and HLA-DRB5 and HLA-DRA.

**Checks.** `tests/test_disease_gaps.py` reads all 2,413 sized studies with cases and controls: at
its own floor's odds ratio each needs its own cases, and the marks agree with the floor. `gaptables`
computes the order of genes and the study sizes with its own code; it agrees with the question on
the 384 locus rows shown for 45 disease traits (five named, forty drawn with a fixed seed) and with
the census over all 513.

**What this does not show.**

- That a study of these sizes would find an association. The gene's variants may be rarer or their
  effects smaller, and the sizes are for one variant at one frequency, with no covariates and no
  staged design.
- That the rare-variant parameters fit any particular gene. They give the scale of such a study.
- That the first gene of a locus is the causal one. The order says where the graph's other evidence
  points. How often the first gene proves to be the one is not measured.
- The proposals have not been put to a model or a person.
