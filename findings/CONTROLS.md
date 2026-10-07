# Lineage, position and locus: how much of a later link is predictable without biology

7 Oct 2026. `python -m pharos_gap_kg.controls pharos319 --cutoff 2018` writes
`docs/data/controls_gene_trait_2018.json`. The three controls were defined in `docs/PROTOCOL.md` and
committed before any was computed. No fitted model has been scored.

Inputs beyond the frozen split: the Experimental Factor Ontology v3.94.0 and the Ensembl 113 gene
annotation, both pinned by digest in the registry. 1,402 of the 1,404 training traits resolve in the
ontology, and 14,470 of the 14,484 training targets are placed on a primary chromosome.

## A fifth of the positives are linked by lineage or by position

Of 130,369 positives:

| | positives | share |
|---|---|---|
| lineage-linked: the target already had a link to the same trait under a broader, narrower or second name | 12,227 | 9.4% |
| position-linked: a gene within 500,000 bases already had a link to the same trait | 15,868 | 12.2% |
| both | 1,749 | 1.3% |
| either | 26,346 | 20.2% |
| neither | 104,023 | 79.8% |

A lineage-linked candidate turns out positive 3.5% of the time and a position-linked one 6.8%,
against 0.64% for a candidate at large.

The positives are 31,885 loci, not 130,369 independent findings: one GWAS signal is assigned to
several neighbouring genes, 4.1 on average.

## What each scorer keeps

Primary measure (mean average precision over traits), with 95% intervals in the stored file.

| scorer | all positives | new links, linked positives removed | new links, every linked candidate removed | each locus counted once |
|---|---|---|---|---|
| B0 chance | 0.0178 | 0.0151 | 0.0162 | 0.0123 |
| B1 target degree | 0.0733 | 0.0496 | 0.0556 | 0.0255 |
| B3 neighbours | 0.0969 | 0.0539 | 0.0662 | 0.0512 |
| B4 lineage only | 0.0459 | 0.0143 | 0.0162 | 0.0339 |
| B5 position only | 0.0468 | 0.0134 | 0.0162 | 0.0310 |

B2 equals B1 in every column. New links are scored over 516 traits.

- **The ontology alone, and the chromosome alone, each rank later links at about two and a half times
  chance** (0.0459 and 0.0468 against 0.0178). Neither knows anything about the trait's biology.
- Both are at chance on new links, as they must be: they score nothing there.
- **Counted by locus, target degree falls to twice chance** (0.0255 against 0.0123), below the ontology
  alone (0.0339) and level with the chromosome alone (0.0310). Most of what degree scored on all
  positives was one region counted many times: Pharos holds the HLA genes as 83 separate entries
  with the same links, and they sit at the top of a ranking by degree for every trait.

## Neighbours against target degree, paired over traits

| positives | difference | 95% interval |
|---|---|---|
| all | +0.0236 | +0.0130 to +0.0359 |
| new links, linked positives removed (a lower bound for a scorer that uses lineage) | +0.0043 | −0.0019 to +0.0102 |
| new links, every linked candidate removed | +0.0106 | +0.0037 to +0.0174 |
| each locus counted once | +0.0257 | +0.0149 to +0.0380 |

**On new links the neighbours score is not shown to beat target degree.** The protocol requires both
readings, and under the first the interval includes zero. More than half of the neighbours gain on
all positives (0.0236) is gone under the second reading (0.0106) and four fifths under the first
(0.0043). Counted by locus the gain is +0.0257, and neighbours scores twice what degree does.

So the figures a fitted model has to beat depend on what is being claimed:

| claim | best parameter-free scorer |
|---|---|
| ranks later links | 0.0969 (neighbours) |
| ranks later links that lineage and position do not explain | 0.0662 with linked candidates removed, 0.0539 with only linked positives removed (neighbours in both) |
| ranks later signals, each locus once | 0.0512 (neighbours) |

## Corrections

**7 Oct, first run.** The first run of this command found 7,925 lineage-linked positives and both
new-link intervals for neighbours against degree above zero (+0.0097, +0.0015 to +0.0177 on the
lower-bound reading). It resolved only 911 of the 1,404 traits. The ontology file names the
replacement of an obsolete term by its full web address, the reader expected a compact id, and so 505
retired EFO disease ids stood for nothing and had no lineage. With the reader fixed the
lineage-linked count is 12,227 and the lower-bound interval includes zero. That run's output is in
the history at commit 87e23b2.

**7 Oct, second run.** 55 training targets had no position: Pharos cross-references them only to gene
ids on alternate haplotype assemblies. Most are HLA allele entries. Each therefore counted as a locus
of its own, and a ranking by degree puts them first for every trait. They are now placed by the
primary gene that carries their symbol (14,470 of 14,484 placed, was 14,415). Position-linked
positives rose from 15,559 to 15,868. The locus-collapsed column moved most: target degree from
0.0492 to 0.0255 and neighbours from 0.0723 to 0.0512, over 31,885 loci where there had been 34,149.
The new-link differences moved from +0.0121 and +0.0060 to +0.0106 and +0.0043, with the same
conclusion. That run's output is in the history at commit 1ec284f.

## What this does not settle

- The window is 500,000 bases and was not varied. A wider window links more positives.
- Lineage follows `is_a` only. Two measurements of one underlying quantity that the ontology files
  under different parents (two cholesterol fractions, say) are not linked here.
- 505 of the 1,404 trait ids in the Pharos table are ids the ontology has since retired. The graph
  still carries them as published.
