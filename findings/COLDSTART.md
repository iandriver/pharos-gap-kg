# Cold start: a trait or a target with no link

10 Oct 2026. `python -m pharos_gap_kg.coldstart pharos319 --cutoff 2018` writes
`docs/data/coldstart_gene_trait_2018.json`. The tasks, the scorers and the rule for reading them were
committed to `docs/PROTOCOL.md` ("Cold start: new traits and new targets") before any scorer was run
on these links. No scorer here has a fitted parameter.

**The links.** 65,523 links first reported after 2018 were left out of the test because their trait or
their target has no training link: 52,601 join a new trait to a training target, 10,008 join a training
trait to a new target, and 2,914 join a new trait to a new target. The neighbours baseline has no
score for any of them. The last group is scored by nothing here.

## New traits

For a trait with no training link, rank the 14,484 training targets. 829 new traits have a later link
to a training target. 543 of them have a lineage relative among the training traits, a trait above or
below them in the ontology, and those hold 37,644 of the 52,601 links. The protocol's count of 843 new
traits included 14 whose only later links are to new targets.

Mean average precision over new traits, with 95% intervals:

| | all new traits | with a relative in training | without | by locus, all |
|---|---|---|---|---|
| N0 chance | 0.0051 [0.0039, 0.0064] | 0.0056 | 0.0042 | 0.0054 |
| N1 target degree | 0.0610 [0.0517, 0.0716] | 0.0653 | 0.0530 | 0.0297 |
| N2 lineage transfer | 0.0427 [0.0353, 0.0510] | 0.0629 | 0.0042 | 0.0445 |
| N3 lineage neighbours | 0.0683 [0.0576, 0.0808] | 0.1021 | 0.0042 | 0.0427 |
| N4 lineage transfer, ties by degree | 0.1019 [0.0886, 0.1162] | 0.1276 | 0.0530 | 0.0671 |
| N5 lineage neighbours, ties by degree | 0.0853 [0.0738, 0.0984] | 0.1023 | 0.0530 | 0.0514 |

- **Lineage transfer works for new traits.** Ranking targets by how many of the new trait's relatives
  they are linked to, ties broken by degree (N4), scores 0.1019 against 0.0610 for degree alone:
  +0.0408 [+0.0300, +0.0519] over all 829 traits. That is the protocol's test and it is passed.
- **Where a relative exists the score nearly doubles**: 0.1276 against 0.0653,
  +0.0623 [+0.0472, +0.0788] over 543 traits. For the 286 traits with no relative the two
  scorers are the same by construction, and both score 0.0530.
- **Counted by locus the gain holds**: 0.0671 against 0.0297, +0.0375 [+0.0277, +0.0476].
- **Degree is what makes it work as a ranking.** Lineage transfer alone (N2) ties every target that
  is linked to no relative, and scores 0.0427, below degree. On the traits with a relative it is level
  with degree (-0.0023 [-0.0148, +0.0102]). The two together are what beats either.
- **Passing the relatives' neighbours scores on does not add.** N5 against N4 is
  -0.0166 [-0.0233, -0.0104]: counting a relative's links is better than also counting what the
  neighbours baseline predicts for the relative. In the earlier period the two are level.
- **The size of it.** 0.1019 for a trait with no link at all is about what the neighbours baseline scores
  on the test links for traits that have links (0.0969). The two are different tasks with different
  candidates, so that is a comparison of scale and not a test.

## New targets

For a training trait, rank the 5,928 targets of the release that have no training link. 322 training
traits have a later link to such a target, 10,008 links over 3,232 new targets. Read the other way, for
each of those new targets the 1,404 training traits are ranked.

| | by trait: ranking the new targets | | by new target: ranking the traits |
|---|---|---|---|
| P0 chance | 0.0063 [0.0053, 0.0075] | Q0 chance | 0.0070 [0.0066, 0.0074] |
| P1 locus degree | 0.0102 [0.0085, 0.0120] | Q1 trait degree | 0.0449 [0.0433, 0.0466] |
| P2 position transfer | 0.0216 [0.0146, 0.0300] | Q2 position transfer | 0.0336 [0.0305, 0.0368] |
| P3 position transfer, ties by locus degree | 0.0250 [0.0172, 0.0349] | Q3 position transfer, ties by trait degree | 0.0611 [0.0576, 0.0648] |

- **Position transfer works for new targets.** The number of a new target's chromosome neighbours
  linked to the trait, ties broken by how many links its neighbours hold in all (P3), scores 0.0250
  against 0.0102 for the tie-breaker alone: +0.0148 [+0.0074, +0.0243] over 322 traits. The protocol's
  test is passed.
- **The score is low in absolute terms**: 0.0250 is four times chance (0.0063). 1,048 of the 10,008
  links (10.5%) are to a target beside one already linked to the trait, and position can
  rank no other above chance.
- **Read by new target it also holds.** Ranking the traits for one new target, position transfer with
  ties by the trait's degree scores 0.0611 against 0.0449 for the trait's degree alone,
  +0.0162 [+0.0130, +0.0197] over 3,232 new targets. Position alone is below the trait's degree
  (-0.0113 [-0.0146, -0.0079]).

## The earlier period

The same scorers with a cutoff of 2016, on the 19,414 links of 2017 and 2018 whose trait or target was
new in 2016. This is inside the training years and decides nothing. It agrees.

| | cutoff 2016 | cutoff 2018 |
|---|---|---|
| new traits scored, and those with a relative | 497, 346 | 829, 543 |
| N1 target degree | 0.0497 | 0.0610 |
| N4 lineage transfer, ties by degree | 0.0663 | 0.1019 |
| N4 against N1 | +0.0166 [+0.0078, +0.0268] | +0.0408 [+0.0300, +0.0519] |
| N5 against N4 | +0.0009 [-0.0037, +0.0053] | -0.0166 [-0.0233, -0.0104] |
| traits with a link to a new target | 210 | 322 |
| P1 locus degree | 0.0049 | 0.0102 |
| P3 position transfer, ties by locus degree | 0.0236 | 0.0250 |
| P3 against P1 | +0.0187 [+0.0072, +0.0336] | +0.0148 [+0.0074, +0.0243] |
| by new target, Q3 against Q1 | +0.0261 [+0.0219, +0.0305] | +0.0162 [+0.0130, +0.0197] |

The gain from lineage is smaller in the earlier period, where a new trait has one relative at the
median and 915 of 11,543 links are to a target linked to a relative; in the later period the median
is two and it is 9,024 of 52,601.

## What this is and is not

- It is a ranked list for a trait or a target that has no link, with a measured precision, on links no
  scorer had been run on. The neighbours baseline cannot produce one.
- It is not a claim about biology. A new trait with a relative is often a narrower or broader form of
  a trait already studied, and a new target beside a linked gene is often the same GWAS signal assigned
  to one more gene.
- The ontology is the 2026 release. A new trait's term and its place in the ontology were fixed after
  2018, in some cases because the trait had been studied. The placement does not use gene links.
- 286 of the 829 new traits have no relative in training, and for them the best scorer here is degree.
  Relations wider than lineage, such as two traits under one parent, were not tried.
- No fitted model was tried on either task.
- The question layer does not serve these rankings yet. A question about a trait with no link is still
  answered as having none.

## Added 12 Oct 2026: the position rankings for targets with no link are served

On Ian's decision the question layer now ranks a target that has no link, by position alone. What is
served and how it is calibrated were written into `docs/PROTOCOL.md` first. The lineage rankings for
new traits are not served: they held on the 2018 split and failed on the prospective one
(`findings/PROSPECTIVE.md`).

**What is served.** For a trait, `predict(trait, unlinked=True)` lists the targets with no link that
have a gene within 500,000 bases linked to the trait, ordered by the number of such genes, ties broken
by the number of links those genes hold in all (P3). For a target with no link, `predict(target)`
lists the traits that a gene within 500,000 bases of it is linked to, ordered by the number of such
genes, ties broken by the number of targets linked to the trait (Q3). `why` gives such a target its
place in the trait's list. Nothing is fitted. The list for a trait is its own list; it is not merged
into the ranking of targets that have links, because the two were never scored together.

**Why these two.** Each beat its own tie-breaker with the interval above zero in three periods:

| | cutoff 2016 | cutoff 2018 | prospective split |
|---|---|---|---|
| for a trait, P3 against P1 | +0.0187 [+0.0072, +0.0336] | +0.0148 [+0.0074, +0.0243] | +0.0250 [+0.0087, +0.0418] |
| for a target, Q3 against Q1 | +0.0261 [+0.0219, +0.0305] | +0.0162 [+0.0130, +0.0197] | +0.0249 [+0.0182, +0.0318] |

**How often a listed pair held** (`python -m pharos_gap_kg.prospect calibrate pharos319`,
`docs/data/calibration_prospect_2022.json`, under `no_link`), on the prospective split's links first
reported in 2023 to 2026, counted over every trait and every target that has a list:

| ranks | for a trait: targets with no link | listed at those ranks | for a target with no link: traits | listed at those ranks |
|---|---|---|---|---|
| 1-5 | 5.0% | 5,440 | 3.3% | 8,102 |
| 6-20 | 2.9% | 7,219 | 2.6% | 15,730 |
| 21-50 | 3.0% | 7,520 | 2.4% | 12,930 |

- **The rates are a few in a hundred.** In the first five for a trait, 5.0% of the listed targets were
  linked to it within four years; in the first five for a target, 3.3% of the listed traits.
- **Being listed is most of the signal.** 2.8% of all 42,982 listed pairs were linked, against 0.13%
  of all pairs of a target with no link and a trait: 22 times. Within a list the order adds
  little below the first five.
- **Position reaches a minority of the later links.** 1,760 of the 2,407 targets with no link are listed for some
  trait, and 1,199 of the 6,729 later links to such targets (17.8%) are to a listed pair. The rest
  cannot be ranked this way at all.
- **A listed target is often the same signal.** A target beside a linked gene is often one more gene
  the same association is assigned to. Every answer says so, and one in the MHC region carries the
  caution the other rankings carry.

**The served lists are the scored lists** (`tests/test_served_no_link.py`). For 60 traits the
prospective split did not merge with another, the number of linked genes beside each of the 2,407
targets with no link is the same from the bundle as in the scored matrix. The tie-breaker is the same
for the 1,190 of 1,763 targets with a linked gene nearby whose neighbouring genes hold no merged
trait, and can only be higher for the others, because the bundle counts a merged trait under each of
its ids. Where every listed target has the scored tie-breaker the order is the scored order.

**Not done.** No usability test of these rankings. A target with no link and no linked gene within
500,000 bases is not ranked. No ranking is offered for a trait that has no link.
