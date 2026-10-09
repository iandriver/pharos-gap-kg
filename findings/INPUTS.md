# Trait lineage and gene position as model inputs

9 Oct 2026. The models, the forms, the fitting and what counts were committed to `docs/PROTOCOL.md`
("Lineage and position as inputs") before anything was fitted. Every model here was **added after
`mf` and `lightgcn` had been scored on the test links**, and the decision to use these two inputs was
taken after seeing the controls, which are computed on the test links.

**No test link has been read for any model in this file.** Everything below is on the inner split of
the training links: fitted on the 57,411 links first reported in or before 2016, validated on the
22,610 links of 2017 and 2018 over 284 traits.

## The two inputs on the validation links

Two traits are in one lineage if one is an ancestor of the other in the Experimental Factor Ontology
v3.94.0. Two targets are neighbours if genes of theirs lie within 500,000 bases on one chromosome
(Ensembl 113). On the inner split there are 855 lineage pairs among the 850 traits and 76,011
neighbour pairs among the 10,648 targets.

| | candidates | validation positives among them | share of the 22,610 positives |
|---|---|---|---|
| every candidate of the 284 traits | 2,994,372 | 22,610 (0.76%) | |
| lineage-linked (B4 above zero) | 35,854 | 1,270 (3.5%) | 5.6% |
| position-linked (B5 above zero) | 58,546 | 2,367 (4.0%) | 10.5% |
| either | | 3,477 | 15.4% |

A linked candidate is about five times as likely to be a later link as a candidate at random. Alone,
the two scores reach 0.0331 (B4 lineage) and 0.0372 (B5 position) on validation, beside 0.0394 for
target degree and 0.0712 for the neighbours baseline B3.

## The combiners as first declared

A combiner scores a candidate by a weighted sum of parameter-free scores. C0 combines target degree
(B1) and neighbours (B3) and uses neither input; C1 adds position (B5), C2 adds lineage (B4), C3 adds
both. As first declared the weights maximised a likelihood: the probability of a trait's validation
positives under a softmax over its candidates (code at commit `86b57d5`,
`docs/data/combine_gene_trait_2018_first.json`). Each figure is cross-validated: a trait is scored by
weights fitted on the other four fifths of the traits.

The form was chosen for C3 among three: F1 0.0730, F2 0.0774, F3 0.0715. Under F2, each score enters
as its percentile rank among the trait's candidates.

| | scores | cross-validated | against neighbours, 0.0712 |
|---|---|---|---|
| C0 `combine` | B1, B3 | 0.0511 | -0.0201 [-0.0345, -0.0082] |
| C1 `combine_pos` | B1, B3, B5 | 0.0717 | +0.0005 [-0.0139, +0.0133] |
| C2 `combine_lin` | B1, B3, B4 | 0.0556 | -0.0156 [-0.0318, -0.0016] |
| C3 `combine_lin_pos` | B1, B3, B4, B5 | 0.0774 | +0.0062 [-0.0086, +0.0207] |

- **Fitted this way, a combiner that holds the neighbours score is worse than the neighbours score.**
  A weight of zero on degree would make C0 rank exactly as B3 does. The fitted weights put it 0.0201
  below. The likelihood counts a positive the same wherever it ranks, and average precision counts the
  top of the list most.
- **The comparison this was built for would have misled.** C3 is 0.0264 above C0 [+0.0144, +0.0405],
  which reads as a large gain from the two inputs, while C3 is 0.0062 above B3 with an interval that
  spans zero.

The protocol was amended on this evidence, before any model of this section was scored on the test
links. This fit stays on record as a validation result and is not scored on the test links.

## The combiners as amended

The weights now maximise the primary measure itself on the validation positives. The weight on B3 is
fixed at 1 and the others start at 0, so the search starts at B3 alone. It moves one weight at a time
over 23 fixed values (0, and plus or minus each power of two from 1/64 to 16), and nothing in it is
random (code at commit `74fe460`, `docs/data/combine_gene_trait_2018.json`, frozen at `4e5859a`).

The form chosen for C3: F1 0.0757, F2 0.0751, **F3 0.0810**. Under F3 a score enters as log(1 + score),
and again multiplied by log(1 + the number of links the trait has), so a weight can differ between
small and large traits.

| | scores | cross-validated | against neighbours, 0.0712 | against C0 | fitted on all 284 traits |
|---|---|---|---|---|---|
| C0 `combine` | B1, B3 | 0.0732 | +0.0020 [-0.0030, +0.0069] | | 0.0759 |
| C1 `combine_pos` | B1, B3, B5 | 0.0813 | +0.0101 [+0.0036, +0.0173] | +0.0081 [+0.0006, +0.0160] | 0.0853 |
| C2 `combine_lin` | B1, B3, B4 | 0.0757 | +0.0045 [-0.0023, +0.0116] | +0.0025 [-0.0029, +0.0092] | 0.0824 |
| C3 `combine_lin_pos` | B1, B3, B4, B5 | 0.0810 | +0.0098 [+0.0019, +0.0177] | +0.0077 [-0.0005, +0.0165] | 0.0901 |

- **Position adds to the neighbours score on validation.** C1 is 0.0101 above B3 and 0.0081 above C0,
  both intervals above zero. Its weight on log(1 + B5) is 0.5, beside 1 on log(1 + B3).
- **Lineage is not shown to add.** C2 is 0.0025 above C0 with an interval that spans zero, and C3
  (0.0810) is level with C1 (0.0813). The half-width of C2 against C0, 0.0060, is the smallest gain
  from lineage this comparison would have resolved.
- **Fitting weights to degree and neighbours alone does not beat neighbours.** C0 is +0.0020
  [-0.0030, +0.0069], so the gains of C1 and C3 come from position and not from the fitting.
- **Each fit flatters itself, more with more weights.** Scored on the traits it was fitted on, a
  combiner reads 0.0027 (C0, 3 free weights), 0.0039 (C1, 5), 0.0067 (C2, 5) and 0.0092 (C3, 7) above
  its cross-validated figure. For C3 that gap is as large as its margin over B3.
- **The chosen form's figure is a little high too**, because F3 was picked on C3's cross-validated
  score. F1 and F2 are also above B3's 0.0712; no interval was computed for them.

The fitted weights. A term not listed has weight 0. "size" is log(1 + the trait's links).

| | weights |
|---|---|
| C0 | B3 1; B1 1/64; B1 x size -1/64; B3 x size -1/8 |
| C1 | B3 1; B1 1/64; B5 1/2; B3 x size 1/8; B5 x size 1/64 |
| C2 | B3 1; B1 1/32; B4 1/16; B1 x size -1/64; B3 x size -1/8; B4 x size 1/32 |
| C3 | B3 1; B1 1/32; B4 -1/2; B5 1/8; B1 x size -1/64; B3 x size -1/8; B4 x size 1/4 |

The sign of the weight on B3 x size differs between C1 and the other three, and C3 sets B5 x size to
0 where C1 uses it. Weights that move this much between neighbouring models are not individually
meaningful; the cross-validated score is what the fit is judged on.
