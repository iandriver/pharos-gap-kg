# Asking the graph

Eight questions, answered from one bundle of files with their digests. They can be asked from the
command line, or by any MCP client through `python -m pharos_gap_kg.server`. The bundle is built by
the README's Run steps, or downloaded, checked and unpacked by `python -m pharos_gap_kg.release get
pharos319` (a 16 MB archive attached to each release).

| question | what it answers |
|---|---|
| `find(text)` | which target, trait, condition or drug a loosely given name means |
| `known(node)` | what is linked to a target or a trait, and since when; and what the other sources say of it |
| `why(target, trait)` | the studies behind one link, or why there is none; and who else asserts the pair |
| `connect(a, b)` | what ties two things together when no link does, a drug included |
| `predict(node)` | what may be linked next, and how often that has held |
| `whats_new(node, since)` | what was first reported after a year |
| `gaps(target)` | what a target lacks for the next Pharos level, the nearest way to close it, and how often that closed |
| `gaps(disease)` | the genes tied to a disease, what each lacks, what would supply it, and how often that has been seen |
| `status()` | what the graph is, where it came from, what it leaves out |

Three rules hold for every answer.

- **Names are resolved against the graph.** A name that matches nothing returns "not in the graph".
  A name that matches several things returns the list and asks which. Nothing is guessed.
- **Every answer is as of a year.** A link exists from the year of the first study behind it.
  `--as-of 2018` answers from what had been reported by the end of 2018.
- **No link is not a negative.** It means no study in this source had reported one by that year, and
  every answer that reports an absence says so.

## From the command line

Every example below is the output of the command shown, on the bundle built from Pharos 3.19.

```bash
python -m pharos_gap_kg.ask status
```

> A graph of 295,327 GWAS links between 18,005 of Pharos's 20,412 protein targets and 2,247 traits,
> each dated by the first study behind it (2005 to 2023), from 12,335 studies in 3,141 papers.
> Beside them, as of the release and undated: 154,363 assertions of a target with a condition from 9
> disease sources, 7,248 more rows withheld by a pruning rule and kept as gaps, 1,782 drugs with
> 4,295 target and 5,943 indication edges, and 159,846 expression observations, which are kept apart
> from the assertions.

### Why is A linked to B

```bash
python -m pharos_gap_kg.ask why APOE "Alzheimer disease" --as-of 2010
```

```json
{
 "found": true,
 "edge": "pharosgapkg:e2e3b1a3f2814",
 "answer": "APOE is linked to Alzheimer disease because 4 GWAS studies in 4 paper(s), the first in 2007, reported a genome-wide significant association assigned to this gene. 25 more studies came after 2010. Apart from GWAS, as of the release: 8 sources assert APOE with this condition (CTD, DisGeNET, JensenLab Knowledge AmyCo, JensenLab Knowledge MedlinePlus, JensenLab Knowledge UniProtKB-KW, JensenLab Text Mining, Monarch, eRAM). This was also read under MONDO:0007088, which has the same name in the sources: a match on the name, not an ontology statement.",
 "claim": {
  "subject": "UniProtKB:P02649",
  "predicate": "biolink:genetically_associated_with",
  "object": "MONDO:0004975"
 },
 "studies": [
  {
   "study": "GCST000018",
   "pmid": 17474819,
   "year": 2007
  },
  {
   "study": "GCST000479",
   "pmid": 19734902,
   "year": 2009
  },
  {
   "study": "GCST000682",
   "pmid": 20460622,
   "year": 2010
  },
  {
   "study": "GCST000808",
   "pmid": 20885792,
   "year": 2010
  }
 ],
 "studies_after_as_of": 25,
 "knowledge": {
  "primary_source": "infores:tiga",
  "aggregator": "infores:pharos",
  "level": "statistical_association",
  "agent": "data_analysis_pipeline"
 },
 "admitted_by": "a GWAS association the source lists for this target and trait, with a study dated by the as-of year",
 "caveat": "GWAS assigns a signal to genes by position. A link says an associated variant lies at this gene, not that this gene is the cause.",
 "source": {
  "release": "pharos319",
  "file": "pharos319.sql.gz",
  "published": "2024-04-29",
  "sha256": "86ebfffdd5fd03bcc23ef43a31251d40572ccd7feec62dadd3215f276f8b0953",
  "tables": [
   "tiga",
   "tiga_provenance",
   "protein"
  ],
  "as_of": 2010,
  "data_ends": 2023
 },
 "asserted_elsewhere": true,
 "other_sources": {
  "read_under": [
   {
    "id": "MONDO:0004975",
    "tie": "as published"
   },
   {
    "id": "MONDO:0007088",
    "tie": "the sources give it the same name"
   }
  ],
  "asserted_by": [
   {
    "source": "CTD",
    "edge": "pharosgapkg:e65eda039a02f",
    "source_evidence": "therapeutic=0"
   },
   {
    "source": "DisGeNET",
    "edge": "pharosgapkg:e4bd3c8d4920c",
    "source_evidence": "pubmed_ids=5; score=0.7; upstream_sources=1"
   },
   {
    "source": "JensenLab Knowledge AmyCo",
    "edge": "pharosgapkg:e368408d871fe"
   },
   {
    "source": "JensenLab Knowledge MedlinePlus",
    "edge": "pharosgapkg:e3951382a1d4b"
   },
   {
    "source": "JensenLab Knowledge UniProtKB-KW",
    "edge": "pharosgapkg:eee4a17c27c61"
   },
   {
    "source": "JensenLab Text Mining",
    "edge": "pharosgapkg:e2d62b048b8de",
    "source_evidence": "zscore=8.483"
   },
   {
    "source": "Monarch",
    "edge": "pharosgapkg:e2250f2aecb8e",
    "source_evidence": "s2o=48.68"
   },
   {
    "source": "eRAM",
    "edge": "pharosgapkg:ee1d51be716a7",
    "source_evidence": "upstream_sources=3"
   }
  ],
  "withheld": [],
  "through_a_drug": [],
  "drug_paths": 0,
  "expression_observations": [],
  "dated": false,
  "answer": "Apart from GWAS, as of the release: 8 sources assert APOE with this condition (CTD, DisGeNET, JensenLab Knowledge AmyCo, JensenLab Knowledge MedlinePlus, JensenLab Knowledge UniProtKB-KW, JensenLab Text Mining, Monarch, eRAM). This was also read under MONDO:0007088, which has the same name in the sources: a match on the name, not an ontology statement."
 }
}
```

The edge id is a digest of the claim, so it is the same in every build and in the KGX export. The
digest under `source` is that of the published file the studies were read from.

### Why is there no link

```bash
python -m pharos_gap_kg.ask why TREM2 "Alzheimer disease" --as-of 2012
```

> As of 2012 there was no link between TREM2 and Alzheimer disease. One was first reported in
> 2013. Nothing nearby in the graph points either way. The largest study behind this trait's links
> (GCST001025, 2011): it had 6,688 cases and 13,685 controls, so a variant carried by 20% of
> chromosomes needed an odds ratio of about 1.18 to be found. An absent link is no effect of about
> that size found at a common variant; smaller effects and rarer variants were not visible. This
> is a rough floor, and a study that reported nothing for the trait is not in the source. The
> study an absent link calls for: for a common variant, about 77,000 cases for an odds ratio of
> 1.05, 20,000 for 1.1 and 5,500 for 1.2, with the largest study's 2.0 controls per case; the
> largest study had enough for 1.2. For rare variants (1 in 1,000 chromosomes combined, odds ratio
> 3), a gene-based test on sequence data would need about 19,000 cases. Apart from GWAS, as of the
> release: 5 sources assert TREM2 with this condition (CTD, DisGeNET, JensenLab Text Mining,
> Monarch, eRAM). This was also read under MONDO:0007088, which has the same name in the sources:
> a match on the name, not an ontology statement.

```bash
python -m pharos_gap_kg.ask why HLA-B "Parkinson disease"
```

> There is no link between HLA-B and Parkinson disease in this graph (data to 2023). No link here
> means no study in this source had reported one by that year. It is not evidence that there is no
> association. What nearby evidence exists: by this graph's ranking (shared-trait neighbours,
> traits already linked, and linked genes nearby) HLA-B ranks 16 of 20,152 candidates for this
> trait. The largest study behind this trait's links (GCST010991, 2021): it had 20,184 cases and
> 397,324 controls, so a variant carried by 20% of chromosomes needed an odds ratio of about 1.08
> to be found. An absent link is no effect of about that size found at a common variant; smaller
> effects and rarer variants were not visible. This is a rough floor, and a study that reported
> nothing for the trait is not in the source. The study an absent link calls for: for a common
> variant, about 55,000 cases for an odds ratio of 1.05, 14,000 for 1.1 and 3,900 for 1.2, with
> the largest study's 19.7 controls per case; the largest study had enough for 1.1 and 1.2. For
> rare variants (1 in 1,000 chromosomes combined, odds ratio 3), a gene-based test on sequence
> data would need about 13,000 cases. HLA-B is in the MHC region of chromosome 6, where one
> association signal spans many genes, so a high rank here says little about this gene in
> particular. Apart from GWAS, nothing in this graph ties HLA-B to this condition: no disease
> source asserts the pair, no drug path runs between them, and no expression observation is on
> record. That is not evidence that there is no association.

### How are two things related

```bash
python -m pharos_gap_kg.ask connect "Parkinson disease" "Crohn's disease"
```

> As of 2023, Parkinson disease and Crohn's disease share 45 target(s) (of 260 and 650).

### What might come next, and how often has that held

```bash
python -m pharos_gap_kg.ask predict "type 2 diabetes mellitus" --top 5
```

> The 5 highest-ranked candidate targets for type 2 diabetes mellitus as of 2023, out of 17,891. Genes
> within 500,000 bases of each other are shown as one lead. In the MHC region of chromosome 6, where
> neighbouring candidates are best read as one lead: 2 of them. A target with no link is not in this
> list; 431 of them lie beside a gene with this link, and unlinked=true ranks those.

| rank | lead | shared-trait neighbours | traits already linked | genes nearby with this link |
|---|---|---|---|---|
| 1 | GPANK1, with C6orf47, CSNK2B, LY6G5B and others (MHC) | 19.214 | 104 | 90 |
| 16 | HIST1H1T, with HFE (MHC) | 19.598 | 103 | 11 |
| 17 | MYBPC3, with SPI1, MADD | 19.181 | 97 | 12 |
| 20 | MARK4 | 15.893 | 73 | 20 |
| 23 | CTDNEP1, with ELP5, GABARAP, PHF23, CLDN7 | 14.547 | 74 | 23 |

Beside each candidate the answer states how often its rank has held:

> In the test on links first reported from 2023 to 2026, 11% of candidates at ranks 1-10 that an
> earlier link explains by trait wording or by position were reported by 2026, counted over every
> trait ranked (base rate 0.5%). A candidate not reported is not known to be absent.

The share is counted over every trait the graph ranks, whether or not the trait had any later link.
Until 12 Oct 2026 the answers quoted a share counted only over the traits that did, which was about
twice as high (21% here) and did not say so.

**What the ranking is.** A candidate's score has three parts, each computed from the links as of the
year asked and each entering as log(1 + part): shared-trait neighbours, weighted 1 (a target scores
high for a trait when the traits it already has share many targets with that trait); the number of
traits the target is already linked to, weighted 1/32; and the number of genes within 500,000 bases of
it that already have the link, weighted 1/4. Every candidate shows its three parts. The two weights
were fitted on links to 2022.

**How it has held.** On 194,017 links first reported from 2023 to 2026, with every link in the release
as training, this ranking scored 0.101 (mean average precision over 1,110 traits) against 0.095 for
shared-trait neighbours alone, +0.0066 with a 95% interval of +0.0031 to +0.0102
(`findings/PROSPECTIVE.md`). The gain is on candidates that an earlier link already explains by gene
position or trait wording, which is why the five above are all genes beside a gene that has the link.
On candidates no earlier link explains it was level with neighbours alone. The `method` field of every
answer says this. Until 11 Oct 2026 the ranking was shared-trait neighbours alone.

A target with no link as of the year is not in this list, because the ranking was tested on targets
that have one. Such targets are ranked separately, by position alone, as shown under "Targets with no
link" below. Ranking the traits for one target uses the same score; that reading was checked on the
same later links after the test (0.063 against 0.062 for neighbours alone) and was not declared before
it.

```bash
python -m pharos_gap_kg.ask predict "type 2 diabetes mellitus" --as-of 2018 --top 8 --unexplained-only
```

> The 8 highest-ranked candidate targets for type 2 diabetes mellitus as of 2018, out of 19,260.
> Genes within 500,000 bases of each other are shown as one lead. In the MHC region of chromosome 6,
> where neighbouring candidates are best read as one lead: 1 of them. Reported between 2019 and 2023:
> 3 of these 8. A target with no link is not in this list; 588 of them lie beside a gene with this
> link, and unlinked=true ranks those.

| rank | lead | reported later |
|---|---|---|
| 41 | ZPR1, with APOA5, BUD13, SIK3, APOA1, APOC3, APOA4 | |
| 49 | AIM2, with CADM3, FCER1A, ACKR1 | |
| 54 | HIST1H4A, with HIST1H3A, HIST1H2BC and others | 2020 |
| 60 | TMEM258, with MYRF, FADS1, FADS2, FEN1, FADS3, DAGLA | 2020 |
| 65 | HLA-A | 2020 |
| 101 | SLC22A5, with IRF1 | |
| 116 | MYB, with HBS1L | |
| 132 | CTCFL | |

`--unexplained-only` drops candidates that an earlier link already explains by trait wording or by
gene position. What remains starts at rank 41 here, because the ranking puts genes beside a linked gene
first. Asking as of a past year shows which candidates were in fact reported afterwards. It is an
illustration and not a test: the two weights were fitted on links later than that year, and the
answer's `method` field says so. The test is the one on links first reported from 2023.

### Targets with no link

About one target in eight has no GWAS link in the release (2,407 of 20,412), and the ranking above
does not cover them. They are ranked by position alone, with nothing fitted.

```bash
python -m pharos_gap_kg.ask predict "Parkinson disease" --unlinked --top 5
```

> The 5 highest-ranked targets with no link as of 2023 that lie beside a gene linked to Parkinson
> disease, out of 38. Targets within 500,000 bases of each other are shown as one lead.

| rank | target | level | genes nearby with this link | nearest of them |
|---|---|---|---|---|
| 1 | RUSC1-AS1 | Tdark | 27 | FDPS, RUSC1 (overlapping), ASH1L |
| 2 | PPP1R2P1 (MHC) | Tdark | 20 | HLA-DQB2, HLA-DQA2, HLA-DQB1 |
| 3 | FBXL19-AS1, with PYCARD-AS1 | Tdark | 16 | FBXL19 (overlapping), ORAI3, SETD1A |
| 5 | FAM223B, with FAM223A | Tdark | 13 | SLC10A3, UBL4A, LAGE3 |
| 7 | TMEM271 | Tdark | 9 | PCGF3, CPLX1, GAK |

A target is listed when a gene within 500,000 bases of it has the link, and the order is by how many
such genes there are, ties broken by how many links those genes hold in all. Asked about a target that
has no link, `predict` gives the traits its neighbouring genes are linked to:

```bash
python -m pharos_gap_kg.ask predict TMEM271 --top 4
```

> TMEM271 has no link as of 2023. The 4 highest-ranked traits that a gene within 500,000 bases of it is
> linked to, out of 65. Traits the ontology files under one another are shown as one lead.

| rank | trait | genes nearby with this link | nearest of them |
|---|---|---|---|
| 1 | appendicular lean mass | 13 | PDE6B, ATP5ME, MYL5 |
| 2 | heel bone mineral density | 13 | PDE6B, ATP5ME, MYL5 |
| 3 | BMI-adjusted waist-hip ratio | 13 | PDE6B, ATP5ME, MYL5 |
| 4 | high density lipoprotein cholesterol measurement | 9 | PCGF3, CPLX1, GAK |

Beside each row the answer states how often that has held:

> In the test on links first reported from 2023 to 2026, 3.3% of traits listed at ranks 1-5 for a
> target with no link were linked to it by 2026. Of all listed pairs 2.8% were, and of all pairs of a
> target with no link and a trait 0.13%. A pair not reported is not known to be absent.

**How to read it.** The rates are low. Being listed at all is most of the signal: a listed pair was
linked within four years 2.8% of the time, about twenty times the rate for any pair, and the top of a
list is a little better than the rest (5.0% in the first five for a trait, 3.3% for a target). Position
reaches only a pair where a gene nearby already has the link, which was 1,199 of the 6,729 later links
to such targets. And a target beside a linked gene is often the same association signal assigned to one
more gene, so a link reported later would say little about that gene in particular. Every answer says
so. `why` gives such a target its place in the list too:

```bash
python -m pharos_gap_kg.ask why TMEM271 "Parkinson disease"
```

> There is no link between TMEM271 and Parkinson disease in this graph (data to 2023). No link
> here means no study in this source had reported one by that year. It is not evidence that there
> is no association. What nearby evidence exists: PCGF3, 129,453 bases away, was linked to this
> trait; TMEM271 has no link of its own, and by position it ranks 7 of the 38 targets with no link
> that lie beside a gene with this one. The largest study behind this trait's links (GCST010991,
> 2021): it had 20,184 cases and 397,324 controls, so a variant carried by 20% of chromosomes
> needed an odds ratio of about 1.08 to be found. An absent link is no effect of about that size
> found at a common variant; smaller effects and rarer variants were not visible. This is a rough
> floor, and a study that reported nothing for the trait is not in the source. The study an absent
> link calls for: for a common variant, about 55,000 cases for an odds ratio of 1.05, 14,000 for
> 1.1 and 3,900 for 1.2, with the largest study's 19.7 controls per case; the largest study had
> enough for 1.1 and 1.2. For rare variants (1 in 1,000 chromosomes combined, odds ratio 3), a
> gene-based test on sequence data would need about 13,000 cases. Apart from GWAS, nothing in this
> graph ties TMEM271 to this condition: no disease source asserts the pair, no drug path runs
> between them, and no expression observation is on record. That is not evidence that there is no
> association.

The ranking beat its own tie-breaker in three periods (`findings/COLDSTART.md`,
`findings/PROSPECTIVE.md`). Because nothing is fitted, asking as of an earlier year is a backtest in
the plain sense. A target with no link and no linked gene within 500,000 bases is not ranked.

### Beside the dated links: other sources, drugs and expression

The dated GWAS links are one part of the graph. `known`, `why` and `connect` also read the rest, and
each of their answers ends with it: which of nine disease sources assert a target with a condition,
which rows a pruning rule withheld and on what count, the drugs that act on a target and what they
treat, and where a gene's expression differs in a condition. Four things hold for that part.

- **It carries no dates.** It is as of the release, whatever year was asked, and the answer says so.
- **An expression observation is not an assertion.** It says a gene's expression was higher or lower
  in a condition, with the fold change and the p-value. It is listed apart from what sources assert.
- **A drug ties a target to a condition through two edges**, drug to target and drug to indication.
  The answer says whether the target is the one the release marks as the drug's mechanism of action.
- **A condition is read under every id that names it**, and the answer lists them with how each is
  tied: its own id, the current term of a retired id (the ontology's statement), or an id the sources
  give the same name (a match on the name alone, which the answer points out).

```bash
python -m pharos_gap_kg.ask why ACHE "Alzheimer disease"
```

> There is no link between ACHE and Alzheimer disease in this graph (data to 2023). No link here
> means no study in this source had reported one by that year. It is not evidence that there is no
> association. What nearby evidence exists: NYAP1, 395,192 bases away, was linked to this trait; by
> this graph's ranking (shared-trait neighbours, traits already linked, and linked genes nearby)
> ACHE ranks 2,147 of 19,794 candidates for this trait. The largest study behind this trait's links
> (GCST007320, 2019): it had 24,087 cases and 383,378 controls, so a variant carried by 20% of
> chromosomes needed an odds ratio of about 1.08 to be found. It also had 47,793 proxy cases (people
> counted by an affected relative); counted as cases, about 1.05. A proxy carries less information
> than a case, so the smallest effect it could find lies between 1.05 and 1.08. An absent link is no
> effect of about that size found at a common variant; smaller effects and rarer variants were not
> visible. This is a rough floor, and a study that reported nothing for the trait is not in the
> source. The study an absent link calls for: for a common variant, about 55,000 cases for an odds
> ratio of 1.05, 14,000 for 1.1 and 4,000 for 1.2, with the largest study's 15.9 controls per case;
> the largest study had enough for 1.1 and 1.2, and with its 47,793 proxy cases counted as cases,
> enough for all three. For rare variants (1 in 1,000 chromosomes combined, odds ratio 3), a
> gene-based test on sequence data would need about 14,000 cases. Apart from GWAS, as of the
> release: 4 sources assert ACHE with this condition (CTD, DisGeNET, JensenLab Text Mining, eRAM); 4
> drugs that act on ACHE are indicated for it (donepezil, galantamine, rivastigmine, tacrine; for 4
> of them ACHE is the mechanism target). This was also read under MONDO:0007088, which has the same
> name in the sources: a match on the name, not an ontology statement.

Each path names both of its edges:

```json
{
 "drug": "donepezil",
 "drug_id": "DrugCentral:946",
 "action_on_target": "INHIBITOR",
 "activity": "Ki=7.8",
 "target_is_the_drugs_mechanism": true,
 "edges": [
  "pharosgapkg:e500c7efb5cdf",
  "pharosgapkg:e14465e7870f1"
 ]
}
```

A row a pruning rule withheld is reported with the rule, what was counted and what would restore it:

```bash
python -m pharos_gap_kg.ask why MT-RNR1 Tinnitus
```

> There is no link between MT-RNR1 and Tinnitus in this graph (data to 2023). No link here means
> no study in this source had reported one by that year. It is not evidence that there is no
> association. Nothing nearby in the graph points either way. The largest study behind this
> trait's links (GCST011785, 2020): it had 172,995 individuals, so a variant carried by 20% of
> chromosomes needed an effect of about 0.027 standard deviations per allele to be found. An
> absent link is no effect of about that size found at a common variant; smaller effects and rarer
> variants were not visible. This is a rough floor, and a study that reported nothing for the
> trait is not in the source. The study an absent link calls for: for a common variant, about
> 310,000 individuals for an effect of 0.02 standard deviations per allele, 50,000 for 0.05 and
> 12,000 for 0.1; the largest study had enough for 0.05 and 0.1. Apart from GWAS, as of the
> release: 1 source asserts MT-RNR1 with this condition (DisGeNET); 1 row from CTD is withheld by
> a pruning rule (PubMed mentions this gene with this disease no more often than the cut for its
> source; 0 papers mention both, a lift of 0 against a cut of 2.63). This was also read under
> MESH:D014012, UMLS:C0040264, which have the same name in the sources: a match on the name, not
> an ontology statement.

```json
{
 "source": "CTD",
 "gap": "pharosgapkg:g64ff2fd4c55c",
 "rule": "literature_lift",
 "why": "PubMed mentions this gene with this disease no more often than the cut for its source",
 "restored_by": "a lift above the cut for its source, or a second source",
 "papers_mentioning_both": 0,
 "lift": 0.0,
 "cut_for_this_source": 2.626
}
```

An expression observation, for a condition that is not a GWAS trait:

```bash
python -m pharos_gap_kg.ask why SH3D21 "malignant mesothelioma"
```

> malignant mesothelioma is not a GWAS trait in this graph, so there is no dated link to ask about.
> Apart from GWAS, as of the release: its expression was higher in this condition in 1 observation
> and lower in 0. An expression observation says the gene's expression differed in the condition. It
> does not say the gene is involved in it.

```json
{
 "condition_as_the_source_names_it": "malignant mesothelioma",
 "direction": "higher",
 "log2_fold_change": 1.6,
 "p": 2.12e-07,
 "edge": "pharosgapkg:eb3aeed147972"
}
```

A drug, and a drug with a target:

```bash
python -m pharos_gap_kg.ask known donepezil
python -m pharos_gap_kg.ask connect donepezil ACHE
```

> As of the release, donepezil acts on 1 target (ACHE as its mechanism of action) and is indicated
> for 2 conditions.

> donepezil acts on ACHE as its mechanism of action. Of the 2 condition(s) donepezil is indicated
> for, a disease source also asserts ACHE with 1. These sources carry no dates: this is as of the
> release, whatever year was asked.

An indication the release names without an identifier is in the graph under a local id, and says so:

```bash
python -m pharos_gap_kg.ask known "Anesthesia for cesarean section"
```

> Anesthesia for cesarean section is not a GWAS trait in this graph, so no dated link is held for
> it. Apart from GWAS, as of the release: 3 drugs are indicated for it.

`predict` and `whats_new` read the dated links only, and say so when asked about a drug or a
condition that has none.

### What does a target lack

Pharos grades each target by explicit evidence, so what a target lacks for the next level is on
record. `gaps` reads it, names the nearest way to close it, and says how often targets that close were
promoted between two releases (`findings/TGAPS.md`).

```bash
python -m pharos_gap_kg.ask gaps LSM14B
```

> LSM14B is Tdark. To leave it needs 2 criteria to clear; the nearest route is: PubMed score up by
> 0.085 (now 4.915, clears at 5); and 1 more antibody (now 50, clears at 51). Of 1,021 targets in this
> band in the Pharos 3.19 database, 11% (10% to 13%) were no longer Tdark in Pharos 4.0.0 (May 2026).
> In the mirror band, 0% of 1,589 targets fell the other way, which is what the pipeline alone moves
> at this distance.

Each gap in the answer names the graph edge that states it (`pharosgapkg:g46b1570d5751` for this one),
in `data/graph/pharos319_target_gaps/gaps.tsv`. Every answer carries the same caution: a gap is what
the grade lacks, and closing it is not the same as learning the protein's function. A target with biological evidence and no compound has a different
gap, without a distance:

```bash
python -m pharos_gap_kg.ask gaps TREM2
```

> TREM2 is Tbio: it has biological evidence and no qualifying compound. It lacks a compound at 1 uM or
> better, the potency cutoff for its family. Of 17,517 targets with no qualifying compound in the
> Pharos 3.19 database, 1,359 (8%) held one in Pharos 4.0.0 (May 2026).

With no target it lists the Tdark targets nearest to leaving, optionally within one family:

```bash
python -m pharos_gap_kg.ask gaps --family Kinase --top 3
```

> The 3 Tdark targets in the family Kinase nearest to leaving Tdark, of 5,499 with a gap. Nearest
> first.

SBK2 (6 more antibodies), NRK (12 more antibodies) and RPS6KC1 (1 more GeneRIF).

### What is missing for a disease

`gaps` with a trait or a condition lists the genes the graph ties to it and what they lack. A gene is
tied by a dated GWAS link, by assertions from two or more disease sources, or by being the mechanism
target of a drug indicated for the condition. One source alone, an expression observation and a
withheld row are counted and are not ties. Six lists follow, each with a count, its strongest rows,
what would close the gap and how often that has been seen:

| list | which genes | closes with |
|---|---|---|
| `no_compound` | tied, with no drug acting on it and no compound at its family's potency cutoff | a compound at that cutoff |
| `dark` | tied and still Tdark | its nearest route to the next level |
| `locus_not_resolved` | GWAS-linked genes within 500,000 bases of one another: one signal on several genes | evidence for one gene of the locus; the genes are listed in the order to test them |
| `no_genetic_support` | tied by sources or by a drug, with no GWAS link | a genome-wide significant association at the gene; the study size that would find one is given |
| `withheld` | rows a pruning rule withheld for the condition | what the row's rule says restores it |
| `observed_only` | an expression observation and nothing else | an assertion or a genetic association |

```bash
python -m pharos_gap_kg.ask gaps "Parkinson disease" --top 2
```

> The graph ties 345 gene(s) to Parkinson disease, some in more than one way: 248 by a GWAS link,
> 88 by two or more disease sources, 24 as the mechanism target of an indicated drug. Of them, 251
> have no drug and no compound at their family's cutoff; 22 are still Tdark; 228 GWAS-linked genes
> sit in 39 loci that carry the link on more than one gene, and other evidence singles out one
> gene in 6 of them (in 3 more, one source or an expression observation points to a gene); 97 are
> tied with no GWAS link. 0 row(s) are withheld by a pruning rule, and 79 gene(s) have only an
> expression observation. The largest study behind this trait's links (GCST010991, 2021): it had
> 20,184 cases and 397,324 controls, so a variant carried by 20% of chromosomes needed an odds
> ratio of about 1.08 to be found. An absent link is no effect of about that size found at a
> common variant; smaller effects and rarer variants were not visible. This is a rough floor, and
> a study that reported nothing for the trait is not in the source. The study an absent link calls
> for: for a common variant, about 55,000 cases for an odds ratio of 1.05, 14,000 for 1.1 and
> 3,900 for 1.2, with the largest study's 19.7 controls per case; the largest study had enough for
> 1.1 and 1.2. For rare variants (1 in 1,000 chromosomes combined, odds ratio 3), a gene-based
> test on sequence data would need about 13,000 cases. These sources carry no dates: all but the
> GWAS links are as of the release, whatever year was asked.

The answer leads with that summary and the count of every list, then shows each list's two strongest
rows. One list is shown in full by naming it: `--kind no_genetic_support --top 10`.

A row names what ties the gene and the gap edge it rests on:

```json
{
 "id": "UniProtKB:Q9BSA9",
 "symbol": "TMEM175",
 "pharos_level": "Tbio",
 "tied_by": {
  "gwas_link": {
   "studies": 10,
   "first_reported": 2010,
   "edge": "pharosgapkg:edcf419e0970e"
  },
  "sources": [
   "JensenLab Knowledge UniProtKB-KW",
   "UniProt Disease"
  ]
 },
 "lacks": "a compound at 1 uM or better, the potency cutoff for its family",
 "graph_edge": "pharosgapkg:g306d5ab715f8"
}
```

A locus lists its genes in the order to test them, what each holds besides the link, any gene that other
evidence singles out, and how many genes a study would test:

```json
{
 "genes": [
  "BCKDK",
  "PRSS53",
  "STX1B",
  "STX4",
  "VKORC1",
  "ZNF646",
  "ZNF668",
  "FUS"
 ],
 "genes_with_the_link": 16,
 "what_each_holds_besides_the_link": {},
 "other_evidence_for": {},
 "singled_out_by_other_evidence": null,
 "in_the_MHC": false,
 "proposal": "16 genes to test. Nothing else in the graph points to one of them.",
 "genes_more": 8
}
```

An absent GWAS link is given with what the largest study could have seen, the study it calls for (the
cases to find a common variant at odds ratios of 1.05, 1.1 and 1.2, and a gene-based test of rare
variants), and the rate at which pairs a source named were GWAS-linked later:

```json
{
 "kind": "no_genetic_support",
 "targets": 97,
 "closes_with": "a genome-wide significant association at the gene, in a study of this trait. The size of study that would find one is under what_an_absent_link_could_have_seen, the_study_that_would_close_it: the same for every gene of the trait",
 "how_often_this_closed": "In the test on links after 2018, over 161 disease traits, a pair one of these sources named with no GWAS link was GWAS-linked later between 2.1% of the time (JensenLab Knowledge UniProtKB-KW) and 18.2% (JensenLab Knowledge MedlinePlus), by source, against 0.6% of all pairs. The sources carry no dates, so a source may have named a pair because of the very studies that linked it later.",
 "what_an_absent_link_could_have_seen": {
  "studies_behind_its_links": 29,
  "with_a_size": 29,
  "largest": {
   "study": "GCST010991",
   "year": 2021,
   "design": "cases and controls",
   "cases": 20184,
   "controls": 397324,
   "smallest_odds_ratio_it_could_find": 1.084
  },
  "the_study_that_would_close_it": {
   "design": "cases and controls",
   "controls_per_case": 19.69,
   "largest_study_cases": 20184,
   "common_variant": [
    {
     "odds_ratio": 1.05,
     "cases": 54627,
     "controls": 1075338,
     "the_largest_study_had_them": false
    },
    {
     "odds_ratio": 1.1,
     "cases": 14316,
     "controls": 281812,
     "the_largest_study_had_them": true
    },
    {
     "odds_ratio": 1.2,
     "cases": 3912,
     "controls": 77009,
     "the_largest_study_had_them": true
    }
   ],
   "rare_variants": {
    "combined_frequency": 0.001,
    "odds_ratio": 3.0,
    "p": 2.5e-06,
    "cases": 13421,
    "controls": 264194,
    "data": "sequence, which a genotyped GWAS is not"
   }
  }
 },
 "top": [
  {
   "id": "UniProtKB:Q9NQ11",
   "symbol": "ATP13A2",
   "pharos_level": "Tbio",
   "tied_by": {
    "sources": [
     "CTD",
     "DisGeNET",
     "JensenLab Knowledge MedlinePlus",
     "JensenLab Knowledge UniProtKB-KW",
     "JensenLab Text Mining",
     "Monarch"
    ]
   },
   "rank_among_gwas_candidates": 2565
  }
 ]
}
```

For a condition that is not a GWAS trait, the lists that need a dated link are left out and the
answer says why:

```bash
python -m pharos_gap_kg.ask gaps "malignant mesothelioma"
```

> The graph ties 105 gene(s) to malignant mesothelioma, some in more than one way: 105 by two or
> more disease sources, 0 as the mechanism target of an indicated drug. Of them, 65 have no drug and
> no compound at their family's cutoff; 1 is still Tdark. 0 row(s) are withheld by a pruning rule,
> and 3,127 gene(s) have only an expression observation. malignant mesothelioma is not a GWAS trait
> in this graph, so no dated link is present or absent for it: its count of loci and its count of
> genes with no GWAS link are null, not zero, and those two lists are not given. These sources carry
> no dates: all but the GWAS links are as of the release, whatever year was asked.

Every answer ends with how to read it:

> A gap is evidence this graph does not hold. It is not a finding that the target matters in the
> disease, and closing it may show that it does not. Two or more sources, a GWAS link or a drug's
> mechanism tie a target to the condition; one source alone does not.

### A study added after the release

The graph's GWAS links stop where Pharos 3.19 stops. One study published after it is added by hand:
the Alzheimer disease GWAS of Uffelmann et al. (Nature, 28 Oct 2026). Its Supplementary Data lists
the genes the study predicts act at each of its 127 loci; 181 of its 191 gene symbols name a graph
gene, and each becomes a link to Alzheimer disease with its locus and the year 2026. The file is
licensed CC BY-NC-ND 4.0, so the table is built from your own copy and is in neither the repository
nor the release archive:

```bash
python -m pharos_gap_kg.adstudy fetch uffelmann2026_supp
python -I -m pharos_gap_kg.later build pharos319
```

The table lands in `data/later/pharos319/`, beside the bundle, and the command line and the server
read it from there. The release's links, rankings, lists and counts do not change. `known` and
`whats_new` for Alzheimer disease, or for a gene the study names, add the study's links apart from the
release's; `why` says whether the study names the gene; and `gaps` for Alzheimer disease marks the
genes with no GWAS link, and the loci, where the study names a gene:

```bash
python -m pharos_gap_kg.ask why ABI3 "Alzheimer disease"
```

> Added after the release: Uffelmann et al., Nature (2026) reports a link the release does not hold:
> it names ABI3 as a predicted effector gene, at locus 109 (rs12185242, reported before that study;
> its FLAMES, Exome and Literature columns also name it). The absence above stands as of the release.

That sentence ends the answer; the rest is the answer without the study. What was built is recorded
in `docs/data/later_pharos319.json`, and the rule is "A study added after the release" in
`docs/PROTOCOL.md`.

## From an MCP client

```bash
uv pip install -e ".[server]" -c constraints.txt     # from a checkout; the package is not on PyPI
python -m pharos_gap_kg.server                # stdio; reads data/bundle
python -m pharos_gap_kg.server --bundle PATH  # or name the bundle
```

The server holds no model and makes no network call. It checks each bundle file against the
manifest's digest when it starts and refuses a bundle that fails. A client's model chooses the
question and words the reply; every fact in the reply comes from the bundle.

For Claude Code:

```bash
claude mcp add pharos-gap-kg -- /path/to/pharos-gap-kg/.venv/bin/python -m pharos_gap_kg.server
```

For a client configured by JSON:

```json
{"mcpServers": {"pharos-gap-kg": {"command": "/path/to/pharos-gap-kg/.venv/bin/python",
                                    "args": ["-m", "pharos_gap_kg.server"]}}}
```

### Claude Desktop on macOS

Quit the app (Cmd-Q), then from the checkout with its `.venv` active:

```bash
python -m pharos_gap_kg.desktop install
```

Open the app and ask in a new chat "Using pharos-gap-kg, what are the gaps for Parkinson disease?". The
reply shows a call to the `gaps` tool before the answer; a reply without one means the server is not
connected, and `~/Library/Logs/Claude/mcp-server-pharos-gap-kg.log` says why. `desktop check` says
whether the app has started the server, and `desktop remove` takes the entry out of the config.

`install` does three things, each for a failure the builder's own setup met:

- Claude Desktop could not start the server from a checkout under `~/Downloads`. macOS protects
  Downloads, Documents and Desktop, and the app starts a server as a separate program that the app's
  own permissions do not cover. So `install` puts a copy of the package, at the versions in
  `constraints.txt`, in an environment of its own, and the release's bundle beside it, in
  `~/.pharos-gap-kg` (`--home` names another place). The bundle comes from the checkout's archive if it
  holds the committed one, and is downloaded otherwise; every digest is checked either way.
- The app keeps its config in memory and writes it back when its settings change, which drops an entry
  added while it runs. `install` writes the entry only while the app is closed, after a dated backup,
  and keeps every other key. Run with the app open, it installs the copy and the bundle and says to
  quit the app and run it again.
- Before it finishes, it starts the installed server over MCP in a sandbox that cannot read Downloads,
  Documents or Desktop, and asks it the Parkinson disease question. A server run from a checkout under
  `~/Downloads` does not start in that sandbox.

The entry it writes, for anyone adding it by hand (with the app quit):

```json
{"mcpServers": {"pharos-gap-kg": {"command": "/Users/<you>/.pharos-gap-kg/venv/bin/python",
  "args": ["-m", "pharos_gap_kg.server", "--bundle", "/Users/<you>/.pharos-gap-kg/data/bundle/pharos319"],
  "env": {"PHAROS_GAP_KG_DATA": "/Users/<you>/.pharos-gap-kg/data"}}}}
```

### Over the web

`--http` serves the same questions as MCP's streamable HTTP at `/mcp`, for a client that connects by
URL, such as a custom connector in Claude:

```bash
python -m pharos_gap_kg.server --http --port 8000                 # http://127.0.0.1:8000/mcp
docker build -t pharos-gap-kg . && docker run --rm -p 8080:8080 pharos-gap-kg
```

It answers each request on its own and keeps no session, so it can be restarted or run as several
copies. It serves plain HTTP and is meant to sit behind a host that terminates TLS. In front of the
questions it:

- limits each client to 60 requests a minute and all clients together to 1,200, answering 429 with
  Retry-After beyond that. With `--behind-proxy` a client is the address the proxy names, which a
  client can set itself, so the per-client limit holds only for clients that do not, and the total is
  the cap that always holds;
- answers at most 4 requests at a time, and refuses a request body over 64 KB;
- answers only requests addressed to a name given with `--public-host` (421 otherwise), or to any name
  with `--public-host '*'`. Bound to a public address, it will not start without one;
- answers `/health` with the bundle's name.

The image (`Dockerfile`) unpacks the release's bundle, checking every file against the digests
committed in `docs/data/bundle_archive_pharos319.json`, and runs as an unprivileged user.

A copy runs at `https://pm4ie6tprh.us-east-1.awsapprunner.com/mcp`: this image on AWS App Runner, one
instance of 1 vCPU and 2 GB, health-checked at `/health`. In Claude, add it as a custom connector
(Settings, Connectors). `scripts/deploy_apprunner.sh` builds the image from the current commit and
creates or updates that service.

## The subgraph behind an answer

`python -m pharos_gap_kg.subgraph` takes the same question and arguments as `ask` and gives the nodes
and edges the answer rests on as KGX JSON, to load into another tool or check with a KGX validator.
Over MCP the same is the `subgraph` tool.

```bash
python -m pharos_gap_kg.subgraph why HMGCR hypercholesterolemia --out hmgcr.json
python -m pharos_gap_kg.subgraph gaps --disease "Parkinson disease" --out parkinson_gaps.json
```

It holds every edge and gap the answer cites, every edge and gap between two nodes the answer names,
and the nodes at their ends. A disease source's assertion runs from the gene of the protein it names,
so the subgraph also holds the gene of each target it names and of each assertion it holds, joined to
the protein by `biolink:has_gene_product`. A GWAS link enters only if it was first reported by the year the answer
is as of; the other sources' edges are undated and enter as of the release, as the answer reads them.
An edge row is the release's: a link's publications run to the release whatever year was asked. Gaps
(rows a pruning rule withheld, and what a target lacks for its next level) are under `gaps`, the
kgx-gaps extension, which a KGX reader ignores. `graph` records the question, the release and its
digest, what was included and the counts.

It is built from the bundle, so it needs no graph files, and its rows are the export's.
`python -m pharos_gap_kg.subgraph --check pharos319` compares every edge and gap the bundle can give
with the merged graph's files: 639,927 edges (all but the ontology's 394 `same_as` edges, which no
answer cites) and 30,484 gaps, with no difference in id, subject, predicate, object, category,
knowledge source, knowledge level or agent type (`docs/data/subgraph_check_pharos319.json`). The
reference KGX validator, run on the subgraphs of ten answers, reports only what it reports for the
whole graph: identifiers under the project's own prefix, and the node that names a kind of gap
(`docs/data/subgraph_standards_pharos319.json`).

## What the bundle is

`python -m pharos_gap_kg.bundle pharos319` writes `data/bundle/pharos319/`, about 20 MB in all: nodes,
links, the study behind each link, trait lineage, neighbouring genes, the scorer the questions rank by
with its test on later links and its calibration by rank, the target-level gaps where they have been
computed, the merged graph's other parts where it has been built (`python -m pharos_gap_kg.kg build
pharos319`, then `bundle pharos319 --context-only`), the size of each study behind the links where
they have been read (`python -I -m pharos_gap_kg.studysize build pharos319`, then `bundle pharos319
--sizes-only`), and a manifest with the digest of every file and of every input. Each row of the other parts carries the id of its edge or gap in the merged graph. The links are every link in the release;
"as of 2018" is a filter the questions apply. Where a dated KGX graph has been built, its edges must
be in the bundle under the same ids or the build refuses.
