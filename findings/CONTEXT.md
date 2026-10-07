# Disease traits, and what the undated disease table already holds

7 Oct 2026. `python -m pharos_gap_kg.context pharos319 --cutoff 2018` writes
`docs/data/context_gene_trait_2018.json`. Both reports were defined in `docs/PROTOCOL.md` and committed
before they were computed.

## Nine tenths of the later links are to measurements, not diseases

| kind of trait | scored traits | positives | chance | target degree | neighbours | neighbours minus degree |
|---|---|---|---|---|---|---|
| disease-like (MONDO, Orphanet or HP term) | 172 | 14,227 | 0.0063 | 0.0671 | 0.0969 | +0.0298 (+0.0096 to +0.0506) |
| measurement or other | 368 | 116,142 | 0.0230 | 0.0762 | 0.0969 | +0.0207 (+0.0069 to +0.0350) |

The largest traits in this data are quantities measured on healthy cohorts (educational attainment,
blood protein levels, height). A test described as "target and disease" is, by count of links, mostly
a test of target and measurement. On disease-like traits alone, how well-studied a target was ranks
later links at more than ten times chance.

## One of the twelve undated sources is the held-out data

Over the 161 scored traits whose live term is a MONDO id: 13,989 positives among 2,314,575
candidates, a base rate of 0.60%. On these traits target degree scores 0.0715 and neighbours 0.1032.

| source | candidates it asserts | of those, positives | precision | share of all positives | training links it also holds | score if used as a predictor |
|---|---|---|---|---|---|---|
| JensenLab Experiment TIGA | 6,302 | 5,569 | 88.4% | 39.8% | 9,642 | 0.6000 |
| DisGeNET | 5,674 | 185 | 3.3% | 1.3% | 1,047 | 0.0120 |
| CTD | 4,503 | 158 | 3.5% | 1.1% | 1,074 | 0.0130 |
| eRAM | 2,092 | 131 | 6.3% | 0.9% | 926 | 0.0134 |
| Expression Atlas | 26,091 | 102 | 0.4% | 0.7% | 464 | 0.0067 |
| JensenLab Knowledge MedlinePlus | 379 | 69 | 18.2% | 0.5% | 420 | 0.0111 |
| JensenLab Text Mining | 1,153 | 55 | 4.8% | 0.4% | 145 | 0.0084 |
| Monarch | 353 | 35 | 9.9% | 0.3% | 213 | 0.0093 |
| DrugCentral Indication | 891 | 17 | 1.9% | 0.1% | 39 | 0.0073 |
| JensenLab Knowledge UniProtKB-KW | 721 | 15 | 2.1% | 0.1% | 40 | 0.0069 |
| UniProt Disease | 91 | 10 | 11.0% | 0.1% | 34 | 0.0070 |
| JensenLab Knowledge AmyCo | 31 | 2 | 6.5% | 0.0% | 5 | 0.0066 |
| any of the twelve | 41,881 | 5,739 | 13.7% | 41.0% | 10,441 | 0.2945 |

- **The source named "JensenLab Experiment TIGA" is the same GWAS aggregation the test links come
  from, loaded into the disease table in 2022.** It asserts 5,569 of the 13,989 positives. A model
  given the disease table as context would score 0.60 on these traits by reading the answers, 5.8
  times what the best parameter-free scorer reaches. This source cannot be context for this test.
- **The other eleven together hold 170 further positives** (5,739 less 5,569). Most are enriched over
  the base rate, some strongly (MedlinePlus 18.2%, UniProt 11.0%, Monarch 9.9%), on few assertions.
  Whether a given row is prior knowledge or a later finding re-entered through curation cannot be
  told from this table.
- **Expression Atlas asserts the most candidates and is below the base rate** (0.4% against 0.60%).
  A gene differentially expressed in a disease is not, by this measure, more likely to gain a genetic
  association with it.

## What this does not settle

- Only an exact match of the live term is counted. A source asserting a broader or narrower disease
  is not.
- The measurement traits have no counterpart in the disease table, so nine tenths of the positives
  are outside this audit.
- Precision here is against later GWAS links only. A source row that is not a later GWAS link may
  still be true.
