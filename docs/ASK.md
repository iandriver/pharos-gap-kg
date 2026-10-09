# Asking the graph

Eight questions, answered from one bundle of files with their digests. They can be asked from the
command line, or by any MCP client through `python -m pharos_gap_kg.server`.

| question | what it answers |
|---|---|
| `find(text)` | which target or trait a loosely given name means |
| `known(node)` | what is linked to a target or a trait, and since when |
| `why(target, trait)` | the studies behind one link, or why there is none |
| `connect(a, b)` | what ties two things together when no link does |
| `predict(node)` | what may be linked next, and how often that has held |
| `whats_new(node, since)` | what was first reported after a year |
| `gaps(target)` | what a target lacks for the next Pharos level, the nearest way to close it, and how often that closed |
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

### Why is A linked to B

```bash
python -m pharos_gap_kg.ask why APOE "Alzheimer disease" --as-of 2010
```

```json
{
 "found": true,
 "edge": "pharosgapkg:e2e3b1a3f2814",
 "answer": "APOE is linked to Alzheimer disease because 4 GWAS studies in 4 paper(s), the first in 2007, reported a genome-wide significant association assigned to this gene. 25 more studies came after 2010.",
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
  "file": "pharos319.sql.gz",
  "published": "2024-04-29",
  "sha256": "86ebfffdd5fd03bcc23ef43a31251d40572ccd7feec62dadd3215f276f8b0953",
  "as_of": 2010
 }
}
```

The edge id is a digest of the claim, so it is the same in every build and in the KGX export. The
digest under `source` is that of the published file the studies were read from.

### Why is there no link

```bash
python -m pharos_gap_kg.ask why TREM2 "Alzheimer disease" --as-of 2012
```

> As of 2012 there was no link between TREM2 and Alzheimer disease. One was first reported in 2013.
> Nothing nearby in the graph points either way.

```bash
python -m pharos_gap_kg.ask why HLA-B "Parkinson disease"
```

> There is no link between HLA-B and Parkinson disease in this graph (data to 2023). No link here
> means no study in this source had reported one by that year. It is not evidence that there is no
> association. What nearby evidence exists: by shared-trait neighbours HLA-B ranks 1 of 20,152
> candidates for this trait. HLA-B is in the MHC region of chromosome 6, where one association signal
> spans many genes, so a high rank here says little about this gene in particular.

### How are two things related

```bash
python -m pharos_gap_kg.ask connect "Parkinson disease" "Crohn's disease"
```

> As of 2023, Parkinson disease and Crohn's disease share 45 target(s) (of 260 and 650).

### What might come next, and how often has that held

```bash
python -m pharos_gap_kg.ask predict "type 2 diabetes mellitus" --as-of 2018 --top 8 --unexplained-only
```

> The 8 highest-ranked candidate targets for type 2 diabetes mellitus as of 2018, out of 19,260.
> Genes within 500,000 bases of each other are shown as one lead. In the MHC region of chromosome 6,
> where neighbouring candidates are best read as one lead: 1 of them. Reported between 2019 and 2023:
> 3 of these 8.

| rank | lead | reported later |
|---|---|---|
| 2 | ZPR1, with APOA5, BUD13, SIK3, APOA1, APOC3 | |
| 5 | AIM2, with CADM3, FCER1A, ACKR1 | |
| 6 | HIST1H4A, with HIST1H3A, HIST1H2BC and others | 2020 |
| 11 | TMEM258, with MYRF, FADS1, FADS2, FEN1, FADS3 | 2020 |
| 13 | HLA-A | 2020 |
| 50 | SLC22A5, with IRF1 | |
| 60 | MYB, with HBS1L | |
| 69 | CTCFL | |

Beside each candidate the answer states how often its rank has held:

> In the 2018 backtest, 8% of candidates at ranks 1-10 for a disease, that no earlier link explains by
> wording or position were reported by 2023 (base rate 0.6%). A candidate not reported is not known to
> be absent.

Asking as of a past year is a backtest anyone can run: the answer says which of its candidates were
in fact reported afterwards. `--unexplained-only` drops candidates that an earlier link already
explains by trait wording or by gene position. The ranking is the parameter-free neighbours score
measured in `findings/BASELINES.md` and `findings/CONTROLS.md`; the `method` field of every answer
quotes those results, including that on links wording and position do not explain it was not shown to
beat ranking by how well-studied a target is.

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

Every answer carries the same caution: a gap is what the grade lacks, and closing it is not the same
as learning the protein's function. A target with biological evidence and no compound has a different
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

## From an MCP client

```bash
pip install "pharos-gap-kg[server]"
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

## What the bundle is

`python -m pharos_gap_kg.bundle pharos319` writes `data/bundle/pharos319/`, about 11 MB: nodes,
links, the study behind each link, trait lineage, neighbouring genes, the backtest summary, the
target-level gaps where they have been computed, and a manifest with the digest of every file and of
every input. The links are every link in the release;
"as of 2018" is a filter the questions apply. Where a dated KGX graph has been built, its edges must
be in the bundle under the same ids or the build refuses.
