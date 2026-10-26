"""The eight questions as an MCP server, so any MCP client can put them to the graph, and the subgraph behind an answer as KGX.

The server holds no model and makes no network call. It loads one bundle, checks each file against
the manifest's digest, and answers from it. A client's language model chooses which question to ask
and words the reply; every fact in the reply comes from here.

The tool descriptions say which question answers which kind of request. That one paragraph is what
lets a model pick the second call of a two-step question, so it is written for the model, not for a
person reading the code.

Over the web (`--http`) the same server answers MCP's streamable HTTP at /mcp, statelessly, so any
request can go to any copy and a restart loses nothing. It is meant to sit behind a host that
terminates TLS. In front of it, `Guard` limits requests per client and in all, caps how many run at
once, and answers /health; the request body is capped at 64 KB, and the Host header must be one the
server was told it answers to.

Run:  python -m pharos_gap_kg.server                 # stdio, the bundle under data/bundle
      python -m pharos_gap_kg.server --bundle PATH
      PHAROS_GAP_KG_BUNDLE=PATH python -m pharos_gap_kg.server
      python -m pharos_gap_kg.server --http --port 8000                      # http://127.0.0.1:8000/mcp
      python -m pharos_gap_kg.server --http --host 0.0.0.0 --public-host kg.example.org --behind-proxy

Needs the optional dependency, from a checkout:  uv pip install -e ".[server]" -c constraints.txt
"""
from __future__ import annotations

import argparse
import asyncio
import json
import sys
import time
from collections import defaultdict, deque
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
- "What is missing for X", "what would it take to characterise X", "which understudied targets are closest" -> gaps with target=X. Quote its how_often_this_closed text with any gap you pass on, and its reading.
- "What is missing for disease Y", "where are the gaps in Y", "which targets for Y have no drug", "which Y loci are unresolved", "what would fill the gaps in Y" -> gaps with disease=Y. It lists the genes tied to Y and six kinds of gap, each with what would close it and how often that has been seen. Pass on its reading: a gap is evidence the graph does not hold, not a finding that the gene matters. When it gives the size of the largest study, pass that on with any "no GWAS link". "What study would close it", "how large a study is needed", "which gene of a locus to test first" -> gaps with disease=Y: it gives the cases a study needs to find a common variant at odds ratios of 1.05, 1.1 and 1.2 and for a rare-variant gene test, and lists each locus's genes in the order to test them; pass on that the order is where the graph's evidence points, not which gene is causal.
- "What is this, where is it from, can I trust it" -> status.
- "Give me the graph behind that answer", "export it as KGX", "which edges does that rest on" -> subgraph, with the same question and arguments.

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
        raise SystemExit('the server needs the MCP library; from a checkout:  uv pip install -e ".[server]" -c constraints.txt') from e
    G = ask.Graph(bundle_dir)
    app = MCPServer("pharos-gap-kg", instructions=INSTRUCTIONS)

    @app.tool(description="Resolve a loosely given name to targets, traits, conditions or drugs in the graph. kind is 'any', 'target', 'drug', or "
                          "'trait' / 'condition', which mean the same: any disease, phenotype or measurement, whether it is a GWAS trait or is named only "
                          "by the other sources. "
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
                          "reported: use this, not predict, to place one named target. With no link it also gives what the trait's largest study could "
                          "have seen and the study an absent link calls for (what_an_absent_link_could_have_seen). Under other_sources, undated: which disease sources assert the "
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

    @app.tool(description="What is missing, for a target or for a disease. target=X: what X lacks for the next Pharos development level, read from the "
                          "release's own evidence: for a Tdark target, how many criteria must clear, the nearest route, and how often targets that close "
                          "were promoted in the back-test; for a target with no qualifying compound, the potency cutoff it lacks. disease=Y (a trait or a "
                          "condition): the genes the graph ties to Y (by a GWAS link, by two or more disease sources, or as the mechanism target of an "
                          "indicated drug) and what they lack, in six lists: no_compound, dark, locus_not_resolved (one GWAS signal on several genes, "
                          "listed in the order to test them, with what each holds besides the link), "
                          "no_genetic_support (tied with no GWAS link, with the size of the largest study, the smallest effect it could have found and the "
                          "cases a study would need, the_study_that_would_close_it), "
                          "withheld (rows a pruning rule withheld) and observed_only. The answer leads with a summary and the count of every list "
                          "(counts), then each list with closes_with, how_often_this_closed and its two strongest rows. To see one list in full, call "
                          "again with kind set to its name (and top up to 25). With neither target nor disease, lists the Tdark targets nearest to "
                          "leaving Tdark (optionally within one family, such as Kinase, GPCR or IC).")
    def gaps(target: str | None = None, family: str | None = None, top: int | None = None, disease: str | None = None, as_of: int | None = None, kind: str | None = None) -> dict:
        return G.gaps(target, family, top, disease, as_of, kind)

    @app.tool(description="What this graph is: its size and years, the published file it was built from and that file's digest, the checks "
                          "it passed, what it leaves out and why, what its candidates are ranked by, and how that ranking held on later links.")
    def status() -> dict:
        return G.status()

    @app.tool(description="The nodes and edges an answer rests on, as KGX JSON, to load into another tool or check with a KGX validator. "
                          "question is known, why, connect, predict, whats_new or gaps; a and b are that tool's first two arguments (the node; the "
                          "target and the trait; the two things; the node and the year for whats_new; the target for gaps), disease is gaps' disease. "
                          "It holds every edge and gap the answer cites, every edge and gap between two nodes the answer names (a GWAS link only if "
                          "first reported by as_of), and their end nodes, with the graph's own ids, predicates and knowledge sources. Gaps are under "
                          "'gaps'. limit caps the edges and gaps returned, cited ones first; graph.counts says how many there are in full.")
    def subgraph(question: str, a: str | None = None, b: str | None = None, as_of: int | None = None, disease: str | None = None,
                 limit: int = 200) -> dict:
        from pharos_gap_kg import subgraph as sg
        if question not in ("known", "why", "connect", "predict", "whats_new", "gaps"):
            return {"error": "question must be one of known, why, connect, predict, whats_new, gaps"}
        args = [x for x in (a, b) if x is not None]
        if question == "whats_new" and len(args) == 2:
            args[1] = int(args[1])
        answer = sg.answer_for(G, question, args, as_of, disease=disease)
        return sg.of(answer, bundle_dir, limit=max(1, min(int(limit), 2000)),
                     question={"question": question, "args": args, "as_of": as_of, "disease": disease})

    return app


LOCAL = ("127.0.0.1", "localhost", "::1")
MAX_BODY = 64 * 1024        # a question's arguments are a few hundred bytes


class Guard:
    """In front of the MCP app: /health, requests per minute per client and in all, and how many run at once.

    A client is the address the request came from (with --behind-proxy, the one the proxy names). A
    refused request gets 429 and Retry-After, and costs the graph nothing."""

    def __init__(self, app, per_client: int, in_all: int, at_once: int, release: str = "", clock=time.monotonic):
        self.app, self.per_client, self.in_all, self.release, self.clock = app, per_client, in_all, release, clock
        self.at_once = asyncio.Semaphore(at_once)
        self.seen: dict[str, deque] = defaultdict(deque)
        self.all: deque = deque()

    def wait(self, who: str) -> float | None:
        """Seconds until this client may ask again, or None and the request is counted."""
        now = self.clock()
        for q in (self.all, self.seen[who]):
            while q and now - q[0] >= 60:
                q.popleft()
        if len(self.seen) > 10_000:
            for k in [k for k, q in self.seen.items() if not q]:
                del self.seen[k]
        for q, cap in ((self.all, self.in_all), (self.seen[who], self.per_client)):
            if len(q) >= cap:
                return max(1.0, 60 - (now - q[0]))
        self.all.append(now)
        self.seen[who].append(now)
        return None

    async def _send(self, send, status: int, body: dict, headers: list | None = None) -> None:
        await send({"type": "http.response.start", "status": status,
                    "headers": [(b"content-type", b"application/json")] + (headers or [])})
        await send({"type": "http.response.body", "body": json.dumps(body).encode()})

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            return await self.app(scope, receive, send)
        if scope["path"] == "/health":
            return await self._send(send, 200, {"status": "ok", "bundle": self.release})
        wait = self.wait((scope.get("client") or ("unknown",))[0])
        if wait is not None:
            return await self._send(send, 429, {"error": "too many requests; try again shortly"}, [(b"retry-after", str(int(wait + 0.999)).encode())])
        async with self.at_once:
            await self.app(scope, receive, send)


def http_app(bundle_dir: Path, host: str, public_hosts: list[str], per_client: int, in_all: int, at_once: int):
    from mcp.server.transport_security import TransportSecuritySettings
    if "*" in public_hosts:
        security = TransportSecuritySettings(enable_dns_rebinding_protection=False)
    elif public_hosts:
        security = TransportSecuritySettings(allowed_hosts=[h for p in public_hosts for h in (p, f"{p}:*")],
                                             allowed_origins=[o for p in public_hosts for o in (f"https://{p}", f"https://{p}:*")])
    else:
        security = None     # the library then allows only the local names, when bound to one
    app = build(bundle_dir).streamable_http_app(stateless_http=True, json_response=True, max_request_body_size=MAX_BODY,
                                                transport_security=security, host=host)
    return Guard(app, per_client, in_all, at_once, release=Path(bundle_dir).name)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundle", type=Path)
    ap.add_argument("--http", action="store_true", help="serve MCP's streamable HTTP at /mcp instead of stdio")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--public-host", action="append", default=[], help="a Host the server answers to; '*' for any")
    ap.add_argument("--behind-proxy", action="store_true", help="take the client's address from the proxy's X-Forwarded-For")
    ap.add_argument("--per-minute", type=int, default=60, help="requests a minute from one client")
    ap.add_argument("--in-all-per-minute", type=int, default=1200)
    ap.add_argument("--at-once", type=int, default=4, help="requests answered at the same time")
    a = ap.parse_args(argv)
    bundle_dir = a.bundle or ask.default_bundle()
    if not a.http:
        build(bundle_dir).run(transport="stdio")
        return 0
    if a.host not in LOCAL and not a.public_host:
        raise SystemExit("REFUSED: bound to a public address, the server needs --public-host for each name it is reached by, or --public-host '*'")
    import uvicorn
    uvicorn.run(http_app(bundle_dir, a.host, a.public_host, a.per_minute, a.in_all_per_minute, a.at_once), host=a.host, port=a.port,
                proxy_headers=a.behind_proxy, forwarded_allow_ips="*" if a.behind_proxy else None, log_level="info")
    return 0


if __name__ == "__main__":
    sys.exit(main())
