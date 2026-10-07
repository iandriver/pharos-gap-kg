"""The seven questions as an MCP server, so any MCP client can put them to the graph.

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

INSTRUCTIONS = """A graph of dated GWAS links between human protein targets (from NIH Pharos) and traits. Seven tools; pick by what is asked.

- Someone names a gene, protein or disease loosely -> find. If it reports an exact match, use that id. If it reports no exact match and several names, ask which is meant; do not pick one.
- "What is known about X" -> known.
- "Why is X linked to Y", "what is the evidence", "is X linked to Y" -> why. It also answers when there is NO link, and says what nearby evidence exists.
- "How are X and Y related" (two genes, two diseases, or a gene and a disease with no direct link) -> connect.
- "What might be linked to X next", "candidate targets for Y" -> predict. Quote its how_often_this_held text with any candidate you pass on.
- "What changed since YEAR" -> whats_new.
- "What is this, where is it from, can I trust it" -> status.

Rules for using the answers:
- Every tool takes names or ids. If a tool returns found=false with candidates, show them and ask which was meant; do not choose.
- If a tool says something is not in the graph, say so. Do not fill the gap from your own knowledge.
- "No link" means no study in this source reported one by that year. It is never evidence of no association; say this when reporting an absence.
- A link places an associated variant at a gene. It does not say the gene causes the trait.
- as_of=YEAR answers from what had been reported by the end of that year. Leave it out for everything in the release."""


def build(bundle_dir: Path):
    try:
        from mcp.server.mcpserver import MCPServer
    except ModuleNotFoundError as e:
        raise SystemExit('the server needs the MCP library:  pip install "pharos-gap-kg[server]"') from e
    G = ask.Graph(bundle_dir)
    app = MCPServer("pharos-gap-kg", instructions=INSTRUCTIONS)

    @app.tool(description="Resolve a loosely given name to targets or traits in the graph. kind is 'any', 'target' or 'trait'. "
                          "Says whether any name matches exactly (exact_match) and whether to ask which is meant (ask_which_is_meant), "
                          "and returns each match's id, name and how many links it has. Use it first when a name could mean several things.")
    def find(text: str, kind: str = "any", limit: int = 8) -> dict:
        return G.find(text, kind if kind in ("any", "target", "trait") else "any", max(1, min(int(limit), 25)))

    @app.tool(description="What is linked to one target or one trait, as of a year: how many links, split by disease and measurement, "
                          "and the best-supported ones with the year each was first reported. A target with no link is reported as such.")
    def known(node: str, as_of: int | None = None, limit: int = 10) -> dict:
        return G.known(node, as_of, max(1, min(int(limit), 50)))

    @app.tool(description="The evidence behind one target-trait link: the GWAS studies, their papers and years, and the published file "
                          "they were read from with its digest. If there is no link as of the year, says so, says whether one was reported "
                          "later, and lists nearby evidence (the same trait under another wording, a neighbouring gene with the link).")
    def why(target: str, trait: str, as_of: int | None = None) -> dict:
        return G.why(target, trait, as_of)

    @app.tool(description="What ties two things together, in either order: a target and a trait, two targets, or two traits. For a target and "
                          "a trait it reports any direct link and the indirect routes (trait wording, neighbouring genes, shared traits). "
                          "For two targets, the traits they share and whether their genes are neighbours. For two traits, shared targets "
                          "and the ontology's relation between them.")
    def connect(a: str, b: str, as_of: int | None = None) -> dict:
        return G.connect(a, b, as_of)

    @app.tool(description="Ranked candidates for what may be linked next: targets for a trait, or traits for a target. Each candidate says "
                          "why it ranks, whether an earlier link already explains it, how often candidates at that rank were later reported "
                          "in the backtest, and (when as_of is in the past) whether it was in fact reported later. These are leads, not "
                          "findings. unexplained_only=true drops candidates that trait wording or gene position already explains.")
    def predict(node: str, as_of: int | None = None, top: int = 10, unexplained_only: bool = False) -> dict:
        return G.predict(node, as_of, top, bool(unexplained_only))

    @app.tool(description="Links to one target or trait that were first reported after a year, counted by year, best-supported first.")
    def whats_new(node: str, since: int, limit: int = 10) -> dict:
        return G.whats_new(node, since, max(1, min(int(limit), 50)))

    @app.tool(description="What this graph is: its size and years, the published file it was built from and that file's digest, the checks "
                          "it passed, what it leaves out and why, and how its ranking performed in the backtest.")
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
