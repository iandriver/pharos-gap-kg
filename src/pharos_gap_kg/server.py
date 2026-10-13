"""The eight questions as an MCP server, so any MCP client can put them to the graph.

The server holds no model and makes no network call. It loads one bundle, checks each file against
the manifest's digest, and answers from it. A client's language model chooses which question to ask
and words the reply; every fact in the reply comes from here.

The tool descriptions say which question answers which kind of request. That one paragraph is what
lets a model pick the second call of a two-step question, so it is written for the model, not for a
person reading the code.

Run:  python -m pharos_gap_kg.server                 # stdio, the bundle under data/bundle
      python -m pharos_gap_kg.server --bundle PATH
      PHAROS_GAP_KG_BUNDLE=PATH python -m pharos_gap_kg.server

Needs the optional dependency:  pip install "pharos-gap-kg[server]"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pharos_gap_kg import ask

INSTRUCTIONS = """A graph of dated GWAS links between human protein targets (from NIH Pharos) and traits, joined to what nine disease sources assert, to drugs with their targets and indications, and to expression observations. Eight tools; pick by what is asked.

- Someone names a gene, protein or disease loosely -> find. If it reports an exact match, use that id. If it reports no exact match and several names, ask which is meant; do not pick one.
- "What is known about X" -> known. X may be a gene, a trait, a disease the other sources name, or a drug.
- "Why is X linked to Y", "what is the evidence", "is X linked to Y", "which sources say X is involved in Y", "is there a drug for Y that acts on X" -> why. It also answers when there is NO link, and says what nearby evidence exists. Its answer ends with what the other sources say of the pair: who asserts it, any row a pruning rule withheld and why, drugs that act on X and treat Y, and expression observations.
- "Where does X rank for Y", "how strong a candidate is X for Y", "how often has a rank like X's held" -> why, not predict. With no link, why gives X's own rank among Y's candidates and the rate that goes with that rank. predict shows only the top of a list, with neighbouring genes folded into one row.
- "How are X and Y related" (two genes, two diseases, a gene and a disease with no direct link, or a drug with a gene or a disease) -> connect.
- "What might be linked to X next", "candidate targets for Y" -> predict. Quote its how_often_this_held text with any candidate you pass on.
- "Which unstudied genes might be linked to Y", "targets with no associations yet for Y" -> predict with unlinked=true. For a gene that has no link at all, predict on the gene gives the traits its neighbouring genes are linked to. Both are by position alone, and the rates they quote are a few in a hundred; say so.
- "What changed since YEAR" -> whats_new.
- "What is missing for X", "what would it take to characterise X", "which understudied targets are closest" -> gaps. Quote its how_often_this_closed text with any gap you pass on, and its reading.
- "What is this, where is it from, can I trust it" -> status.

Rules for using the answers:
- Every tool takes names or ids. If a tool returns found=false with candidates, show them and ask which was meant; do not choose.
- If a tool says something is not in the graph, say so. Do not fill the gap from your own knowledge.
- "No link" means no study in this source reported one by that year. It is never evidence of no association; say this when reporting an absence.
- A link places an associated variant at a gene. It does not say the gene causes the trait.
- as_of=YEAR answers from what had been reported by the end of that year. Leave it out for everything in the release.
- Only the GWAS links are dated. What the other sources, drugs and expression observations say is as of the release whatever year is asked; say so when a year matters.
- An expression observation says a gene's expression was higher or lower in a condition. Report it as an observation, never as "X is involved in Y".
- A drug that acts on a gene and treats a disease ties the two only as strongly as the gene's part in the drug's action: pass on whether the gene is the drug's mechanism target.
- When an answer says it was also read under another id with the same name, pass that on."""


def build(bundle_dir: Path):
    try:
        from mcp.server.mcpserver import MCPServer
    except ModuleNotFoundError as e:
        raise SystemExit('the server needs the MCP library:  pip install "pharos-gap-kg[server]"') from e
    G = ask.Graph(bundle_dir)
    app = MCPServer("pharos-gap-kg", instructions=INSTRUCTIONS)

    @app.tool(description="Resolve a loosely given name to targets, traits, conditions or drugs in the graph. kind is 'any', 'target', 'trait' "
                          "(GWAS traits and the conditions the other sources name), 'condition' or 'drug'. "
                          "Says whether any name matches exactly (exact_match) and whether to ask which is meant (ask_which_is_meant), "
                          "and returns each match's id, name and how many links or rows it has. Use it first when a name could mean several things.")
    def find(text: str, kind: str = "any", limit: int = 8) -> dict:
        return G.find(text, kind if kind in ("any", "target", "trait", "condition", "drug") else "any", max(1, min(int(limit), 25)))

    @app.tool(description="What is linked to one target or one trait, as of a year: how many links, split by disease and measurement, "
                          "and the best-supported ones with the year each was first reported. A target with no link is reported as such. "
                          "Under other_sources, undated: the conditions disease sources assert it with, rows a pruning rule withheld, drugs that act on it "
                          "(and whether it is their mechanism target), and expression observations, which are listed apart. For a drug: its targets and "
                          "indications. For a condition that is not a GWAS trait: the targets sources assert with it and the drugs indicated for it.")
    def known(node: str, as_of: int | None = None, limit: int = 10) -> dict:
        return G.known(node, as_of, max(1, min(int(limit), 50)))

    @app.tool(description="The evidence behind one target-trait link: the GWAS studies, their papers and years, and the published file "
                          "they were read from with its digest. If there is no link as of the year, says so, says whether one was reported "
                          "later, and lists nearby evidence (the same trait under another wording, a neighbouring gene with the link). With no link "
                          "it also gives the target's own rank among the trait's candidates and how often candidates at that rank were later "
                          "reported: use this, not predict, to place one named target. Under other_sources, undated: which disease sources assert the "
                          "pair (asserted_by), rows a pruning rule withheld with the rule, the count behind it and what would restore the row (withheld), "
                          "drugs that act on the target and treat the condition with the two edges of each path (through_a_drug), and expression "
                          "observations with direction, fold change and p-value. The second argument may be a disease the other sources name.")
    def why(target: str, trait: str, as_of: int | None = None) -> dict:
        return G.why(target, trait, as_of)

    @app.tool(description="What ties two things together, in either order: a target and a trait, two targets, two traits, or a drug with a target "
                          "or a condition (does it act on the target, is it indicated, which of its targets are tied to the condition). For a target and "
                          "a trait it reports any direct link and the indirect routes (trait wording, neighbouring genes, shared traits). "
                          "For two targets, the traits they share and whether their genes are neighbours. For two traits, shared targets "
                          "and the ontology's relation between them.")
    def connect(a: str, b: str, as_of: int | None = None) -> dict:
        return G.connect(a, b, as_of)

    @app.tool(description="The top of the ranked candidates for what may be linked next: targets for a trait, or traits for a target. It lists "
                          "only the highest, with neighbouring genes folded into one row; for the rank of one named target use why. Each candidate says "
                          "the three parts of its score, whether an earlier link already explains it, how often candidates at that rank were "
                          "reported in the test on later links, and (when as_of is in the past) whether it was in fact reported later. The "
                          "ranking favours genes beside a gene that already has the link; unexplained_only=true drops candidates that trait "
                          "wording or gene position already explains. A target with no link at all is not in that list: for a trait, "
                          "unlinked=true ranks the targets with no link that lie beside a gene that has this one, and for a target with no link "
                          "predict lists the traits its neighbouring genes are linked to, both by position alone. These are leads, not findings.")
    def predict(node: str, as_of: int | None = None, top: int = 10, unexplained_only: bool = False, unlinked: bool = False) -> dict:
        return G.predict(node, as_of, top, bool(unexplained_only), bool(unlinked))

    @app.tool(description="Links to one target or trait that were first reported after a year, counted by year, best-supported first.")
    def whats_new(node: str, since: int, limit: int = 10) -> dict:
        return G.whats_new(node, since, max(1, min(int(limit), 50)))

    @app.tool(description="What a target lacks for the next Pharos development level, read from the release's own evidence: for a Tdark target, "
                          "how many criteria must clear, the nearest route, and how often targets that close were promoted in the back-test; for a "
                          "target with no qualifying compound, the potency cutoff it lacks. With no target, lists the Tdark targets nearest to "
                          "leaving Tdark (optionally within one family, such as Kinase, GPCR or IC).")
    def gaps(target: str | None = None, family: str | None = None, top: int = 10) -> dict:
        return G.gaps(target, family, top)

    @app.tool(description="What this graph is: its size and years, the published file it was built from and that file's digest, the checks "
                          "it passed, what it leaves out and why, what its candidates are ranked by, and how that ranking held on later links.")
    def status() -> dict:
        return G.status()

    return app


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundle", type=Path)
    a = ap.parse_args(argv)
    build(a.bundle or ask.default_bundle()).run(transport="stdio")
    return 0


if __name__ == "__main__":
    sys.exit(main())
