"""Preview the published Higginbotham CSF panel through the installed local MCP tools.

Pass aaz9360_Table_S2.xlsx from the authors' supplement as the only argument.
This development example uses the recorded Codex proposal and calls no hosted model.
"""
import asyncio
import hashlib
import json
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
import pandas as pd

from pharos_gap_kg import add, paths


ABOUT = {"citation": "Higginbotham et al. Integrated proteomics reveals brain-based cerebrospinal fluid biomarkers in asymptomatic and symptomatic Alzheimer’s disease",
         "journal": "Science Advances", "published": "2020-10-21", "doi": "10.1126/sciadv.aaz9360",
         "licence": "CC BY-NC 4.0", "what": "Brain and cerebrospinal fluid proteomics in Alzheimer disease.",
         "disease": "Alzheimer disease", "add": "Identify the authors’ brain-linked CSF marker panel. Keep the disease versus control log2 difference."}
CHOICE = {"sheet": "Table S2G", "header_row": 3}
MAPPING = {"map": True, "subject": {"columns": ["Gene Symbol"], "resolver": "gene_symbol", "category": "biolink:Protein"},
           "object": {"constant": "MONDO:0004975", "category": "biolink:Disease"},
           "predicate": "biolink:biomarker_for", "association_category": "biolink:Association",
           "knowledge_level": "statistical_association", "agent_type": "data_analysis_pipeline", "qualifiers": {},
           "edge_properties": ["Brain Module", "CSF Log2 Difference"], "keep_if": None, "split_by": None, "by_value": {}}
SOURCE_SHA = "6bcfdfa517c61387103637e34f18df1aefc997c9d0f2e33832d5a1df8d087f0c"


async def preview(source, directory):
    bundle = directory / "bundle/pharos319"
    shutil.copytree(paths.bundle("pharos319"), bundle)
    params = StdioServerParameters(command=sys.executable, args=["-m", "pharos_gap_kg.server", "--bundle", str(bundle)],
                                   env={"PHAROS_GAP_KG_DATA": str(directory)})
    args = {"file": str(source), "about": ABOUT}
    async with stdio_client(params) as streams, ClientSession(*streams[:2]) as client:
        await client.initialize()

        async def call(tool, more):
            result = await client.call_tool(tool, {**args, **more})
            if result.is_error:
                raise ValueError(result.content)
            return json.loads(result.content[0].text)

        sheets = await call("add_card", {})
        choice = await call("add_check", {"stage": "sheet", "proposal": CHOICE})
        if not choice["accepted"]:
            raise ValueError(choice)
        card = await call("add_card", CHOICE)
        result = await call("add_check", {**CHOICE, "stage": "mapping", "proposal": MAPPING})
    if not result["accepted"]:
        raise ValueError(result["problems"])

    # Count actual measurements independently of the validator. Section labels, repeated
    # headers and the AD-CT direction labels have no numeric CSF difference.
    table = add.read_sheets(source)[CHOICE["sheet"]]
    frame, problems = add.frame(table, CHOICE["header_row"])
    if problems:
        raise ValueError(problems)
    measured = pd.to_numeric(frame["CSF Log2 Difference"], errors="coerce").notna()
    proteins = pd.read_parquet(bundle / "nodes.parquet")
    proteins = proteins.loc[proteins.kind.eq("target"), ["id", "name"]]
    matches = frame.loc[measured].merge(proteins, left_on="Gene Symbol", right_on="name")
    if set(frame.loc[~measured, "Gene Symbol"]) & set(proteins.name):
        raise ValueError("A heading resolves to a graph target")
    report = result["report"]
    distinct = matches.id.nunique()
    assert result["would_join"]["graph_genes"] == distinct
    assert report["source_rows"] == report["accounted_source_rows"] == len(frame)
    assert report["outcomes"]["graph_pair"] == matches["Gene Symbol"].nunique()
    return {"about": ABOUT, "workbook": source.name, "sha256": SOURCE_SHA,
            "download": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7577712/supplementaryFiles",
            "sheets": len(sheets["card"]["sheets"]), **CHOICE, "sheet_check": choice, "profile": card["profile"],
            "mapping": MAPPING, "references": result["checked_against"], "accepted": result["accepted"],
            "rows": report["source_rows"], "accounted": report["accounted_source_rows"],
            "outcomes": report["outcomes"], "issues": report["issue_counts"], "preview": result["would_join"],
            "independent_check": {"numeric_measurement_rows": int(measured.sum()),
                                  "distinct_measured_symbols": frame.loc[measured, "Gene Symbol"].nunique(),
                                  "structural_rows": int((~measured).sum()), "mapped_measured_symbols": matches["Gene Symbol"].nunique(),
                                  "distinct_targets": distinct},
            "joined": False, "scope": "One development example. Codex selected the sheet and proposed the mapping. Read-only preview on the compact release. The authors propose candidate biomarkers."}


if __name__ == "__main__":
    source = Path(sys.argv[1]).resolve()
    if hashlib.sha256(source.read_bytes()).hexdigest() != SOURCE_SHA:
        raise SystemExit("The supplied workbook differs from the published file used for this test.")
    with TemporaryDirectory() as tmp:
        print(json.dumps(asyncio.run(preview(source, Path(tmp))), indent=2))
