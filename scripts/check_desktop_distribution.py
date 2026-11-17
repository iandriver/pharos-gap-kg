"""Use the Desktop copy's Python to check a wheel installation with the pharos319 bundle.

Set PHAROS_GAP_KG_DATA to the installed data directory. The check uses a temporary
CSV and the local MCP server. It reads the bundle without joining a study.
"""
import asyncio
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from mcp import ClientSession, StdioServerParameters, stdio_client
from pharos_gap_kg import paths


async def check():
    if paths.ROOT is not None:
        raise RuntimeError("Use the installed wheel's Python outside the source checkout")
    with TemporaryDirectory() as folder:
        csv = Path(folder) / "genes.csv"
        csv.write_text("gene\nTMEM175\nPARK7\n")
        params = StdioServerParameters(command=sys.executable,
                                       args=["-I", "-m", "pharos_gap_kg.server"],
                                       cwd=folder, env={"PHAROS_GAP_KG_DATA": str(paths.DATA)})
        about = {"citation": "Installation check", "journal": "Fixture", "published": "2026-11-17",
                 "doi": "10.0000/fixture", "licence": "fixture", "what": "Synthetic table",
                 "disease": "Parkinson disease", "add": "Check the gene column"}
        table = {"file": str(csv), "about": about, "sheet": csv.name, "header_row": 1}
        async with stdio_client(params) as streams, ClientSession(*streams[:2]) as session:
            await session.initialize()
            tools = sorted(tool.name for tool in (await session.list_tools()).tools)

            async def call(name, arguments):
                reply = await session.call_tool(name, arguments)
                if getattr(reply, "is_error", getattr(reply, "isError", False)):
                    raise RuntimeError(f"{name} failed: {reply.content}")
                return json.loads(reply.content[0].text)

            gaps = await call("gaps", {"disease": about["disease"]})
            card = await call("add_card", table)
            mapping = {"map": True, "subject": {"columns": ["gene"], "resolver": "gene_symbol",
                                                "category": "biolink:Protein"},
                       "object": {"constant": gaps["condition"]["id"], "category": "biolink:Disease"},
                       "predicate": "biolink:genetically_associated_with", "association_category": "biolink:Association",
                       "knowledge_level": "statistical_association", "agent_type": "data_analysis_pipeline",
                       "qualifiers": {}, "edge_properties": [], "keep_if": None, "split_by": None, "by_value": {}}
            result = await call("add_check", {**table, "stage": "mapping", "proposal": mapping})
            if not ({"add_card", "add_check", "add_join"} <= set(tools) and gaps["found"]
                    and card["stage"] == "mapping" and result["accepted"]
                    and result["would_join"]["graph_genes"] == 2):
                raise RuntimeError(f"Installed server failed its dry run: {result}")
            return {"tools": tools, "disease_found": gaps["found"], "preview_genes": 2,
                    "checked_against": result["checked_against"],
                    "schema_sha256": hashlib.sha256(paths.BIOLINK_SCHEMA.read_bytes()).hexdigest()}


if __name__ == "__main__":
    print(json.dumps(asyncio.run(check()), indent=2))
