# Does the rest of Pharos add to the neighbours score

7 Oct 2026. `python -m pharos_gap_kg.ctxtables pharos319` writes the context tables and
`python -m pharos_gap_kg.enrich pharos319 --cutoff 2018` writes `docs/data/enrich_gene_trait_2018.json`.
The five scorers, how they combine with the neighbours score, and how the result would be read were
committed to `docs/PROTOCOL.md` before any was computed.

Context is what the 3.19 database holds besides the dated GWAS links: 494,130 interaction pairs,
2,902 pathways, 12,634 GO terms, 71,756 compounds active on two or more targets, and 276,642
target-disease pairs from eleven sources. None of it is dated except the compounds.

## Each context table alone

Primary measure and its three variants. Chance is 0.0178, 0.0063, 0.0162 and 0.0123 in the four columns.

| scorer | all traits | disease-like traits | new links | each locus once |
|---|---|---|---|---|
| B1 target degree | 0.0733 | 0.0671 | 0.0556 | 0.0255 |
| B3 neighbours | 0.0969 | 0.0969 | 0.0662 | 0.0512 |
| X1 interaction partners | 0.0290 | 0.0267 | 0.0241 | 0.0215 |
| X2 pathways | 0.0315 | 0.0366 | 0.0222 | 0.0264 |
| X3 GO terms | 0.0294 | 0.0282 | 0.0229 | 0.0252 |
| X4 shared compounds | 0.0185 | 0.0067 | 0.0170 | 0.0142 |
| X5 disease sources | 0.0212 | 0.0167 | 0.0190 | 0.0186 |

- Interactions, pathways and GO terms each carry some signal: on disease-like traits they score four
  to six times chance. Every one is well below the neighbours score, and below target degree except
  when loci are counted once.
- Shared compounds are at chance.
- The eleven disease sources score 0.0167 on disease-like traits. They assert 54,681 candidates, of
  which 464 are later links.

## Combined with the neighbours score

Each context scorer averaged with B3 by percentile rank, and the difference from B3 alone.

| | all traits | 95% interval | each locus once | 95% interval |
|---|---|---|---|---|
| B3 with partners | −0.0301 | −0.0384 to −0.0233 | −0.0061 | −0.0153 to +0.0015 |
| B3 with pathways | −0.0354 | −0.0447 to −0.0276 | −0.0059 | −0.0158 to +0.0021 |
| B3 with GO terms | −0.0419 | −0.0517 to −0.0333 | −0.0027 | −0.0130 to +0.0063 |
| B3 with shared compounds | −0.0290 | −0.0382 to −0.0219 | −0.0121 | −0.0194 to −0.0060 |
| B3 with disease sources | −0.0083 | −0.0146 to −0.0031 | −0.0016 | −0.0060 to +0.0026 |
| B3 with all five | −0.0571 | −0.0685 to −0.0476 | −0.0180 | −0.0284 to −0.0105 |

**Context does not add to the neighbours score.** The stated test was that "B3 with all five" beat B3
on both measures. It is worse on both, and no single table improves B3 on either.

An equal average is a blunt way to combine a strong scorer with a weak one, so a bound was computed as
well: the weight on B3 was chosen on the test links themselves, among 0.5, 0.7, 0.8, 0.9, 0.95, 0.98
and 1. For every one of the five the best weight is 1, which is B3 alone. A weight tuned on the test
is an upper bound and not a result, and here the upper bound is no gain.

## What this shows and what it does not

- For links of this kind, pulling in the rest of Pharos and blending it by rank does not improve on a
  score built from the dated links alone. That holds even though the context is as of 2022 and the
  links being predicted are from 2019 on.
- It does not show that a fitted model cannot use these tables. It shows that the gain, if there is
  one, is not available to any rank blend of these five scorers with B3.
- Counted by locus, four of the five single-table differences have intervals that include zero. On
  that measure the tables do no measurable harm and no measurable good.
- The scorers are one way of reading each table. A pathway could be used differently.
