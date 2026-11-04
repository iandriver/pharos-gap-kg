# Demo video script

The organisers asked for a hosted demo for 17 November, with a video as the fallback, and want it
to show usability. Ian narrates. The script runs under five minutes: about 600 words of narration
plus the pauses for screen actions. Cut the wait for Claude's reply in scene 6 when editing.

Record it against the hosted server once it has been redeployed from 2157c4e or later. Each figure
in the narration is one the server shows on screen. If an answer on screen differs, say what the
screen says.

## Before recording

- Open three browser tabs, light mode, about 1280 pixels wide, zoomed to 125%:
  1. `https://pm4ie6tprh.us-east-1.awsapprunner.com/`
  2. `https://pm4ie6tprh.us-east-1.awsapprunner.com/add?run=D`
  3. A new chat in Claude with the pharos-gap-kg connector switched on.
- Check that the redeploy is live: the old deploy has no web pages, so if tab 1 shows the page,
  it is the new one.
- Do one practice pass of scene 6 in Claude, then start a fresh chat for the recording.

## 1. The page (0:00 to 0:30)

On screen: the Ask page, from the top.

> This is the Pharos gap graph. I built it from NIH Pharos 3.19, and it holds 295,327 GWAS links
> between genes and traits, each dated by the first study that reported it. It's for a scientist
> who wants to know two things about a disease: what is known, and what is missing. Code writes
> every answer on this page from the graph, with no language model.

## 2. What is missing for Parkinson disease (0:30 to 1:35)

On screen: click the "Parkinson disease" chip. The exact match opens with the "What is missing"
tab.

> I'll look up Parkinson disease. The first tab is what is missing. The graph ties 345 genes to
> Parkinson disease, and the first sentence counts what they lack: 251 have no drug and no potent
> compound, and 22 are still Tdark, the least-studied level in Pharos.

On screen: keep reading down to "The largest study behind this trait's links".

> The answer also bounds what a missing link could hide. The largest Parkinson study here had about
> 20,000 cases, so it could find a common variant with an odds ratio of about 1.08. Then it gives
> the size of study that would find a smaller effect.

On screen: scroll to the first list. Point to TMEM175.

> Each kind of gap has its own list. TMEM175 has a GWAS link from ten studies since 2010, and still
> no compound at the potency cutoff for its protein family. Below the list is how often this kind
> of gap has closed: 8% of such targets had gained a compound by Pharos 4.0.

## 3. A gene with no GWAS link (1:35 to 2:10)

On screen: scroll to the list "tied by two or more sources or by a drug's mechanism, with no GWAS
link". Point to PARK7. Then type PARK7 and Parkinson disease into the box "Is a gene linked to a
disease" and click "Ask why".

> This list runs the other way. Mutations in PARK7 cause a rare inherited form of Parkinson
> disease, and six sources name it, but no GWAS study has linked it. If I ask why, the answer says
> there is no GWAS link, and that this is not evidence of no association. It ranks PARK7 at 8,820
> out of 20,152 candidates for the trait.

## 4. Candidate genes (2:10 to 2:35)

On screen: scroll back up to the Parkinson disease card and click "Candidate genes".

> Candidate genes ranks every target for the trait from what the graph already holds. The ranking
> was tested on 194,017 GWAS Catalog links first reported from 2023 to 2026, and scored 0.101
> against 0.095 for the baseline it had to beat. The page calls these leads for a person to check.

## 5. Adding a table from a paper (2:35 to 3:50)

On screen: tab 2, the "Adding a table" page, with "Refractive error" selected.

> The second page shows how a scientist adds a table from a paper the graph doesn't hold. This
> run is a refractive error study from Nature Genetics, published in April 2026, that nobody on the
> project had read. The person gives only the file and a few lines about the paper.

On screen: open "The file's 15 sheets", then scroll through "1. The sheet" and "2. The mapping".

> A model picked the sheet, Supplementary Table 3 of fifteen, and proposed how each row becomes an
> edge. Code refused its first two proposals and said why: the first used a Biolink predicate that
> is deprecated, and the second put a protein where Biolink allows only a gene. A second model
> started over, and its proposal passed.

On screen: scroll to "3. The person" and "The check fixed before the run".

> The person sees what would join: 500 links from 500 genes, and the 31 rows that gave no gene.
> Before the run I fixed a check against the genes the GWAS Catalog maps to the same study. The
> table names one gene per locus, so it could hold at most 68% of the Catalog's genes, and it holds
> 44%. The table stays out of the graph until the person says to join it.

## 6. Asking from Claude (3:50 to 4:40)

On screen: tab 3, Claude. Type: "Using pharos-gap-kg, was TREM2 linked to Alzheimer disease as of
2012?"

> The same questions are served over MCP, so Claude can ask them too. I'll ask whether TREM2 was
> linked to Alzheimer disease as of 2012. Claude calls the why tool with that year, and the graph
> answers that there was no link in 2012 and that the first was reported in 2013. The largest
> Alzheimer study by then had about 6,700 cases, too few to find a common variant with an odds
> ratio under about 1.18.

On screen: expand the tool call so the JSON answer shows.

> Claude writes the sentences, and the facts come from this answer, which anyone can open and
> check.

## 7. Close (4:40 to 4:50)

On screen: the Ask page, then the GitHub repository.

> The code is in the repository, with the protocol for each test, written before the test ran.
> The page is at the address on screen.
