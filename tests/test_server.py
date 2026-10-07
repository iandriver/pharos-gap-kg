"""The server, driven through the MCP protocol over stdio, on the six-node graph of the self-test."""
import asyncio
import json
import sys

import pytest

mcp = pytest.importorskip("mcp")
from pharos_gap_kg import ask_selftest  # noqa: E402

NAMES = ["connect", "find", "known", "predict", "status", "whats_new", "why"]


async def _session(bundle, work):
    from mcp import ClientSession, StdioServerParameters, stdio_client
    params = StdioServerParameters(command=sys.executable, args=["-m", "pharos_gap_kg.server", "--bundle", str(bundle)])
    async with stdio_client(params) as streams:
        async with ClientSession(streams[0], streams[1]) as s:
            init = await s.initialize()
            return await work(s, init)


def _payload(result):
    structured = getattr(result, "structured_content", None) or getattr(result, "structuredContent", None)
    return structured if structured else json.loads(result.content[0].text)


def test_the_server_offers_the_seven_questions_and_answers_them(tmp_path):
    ask_selftest.tiny(tmp_path)

    async def work(s, init):
        tools = (await s.list_tools()).tools
        status = _payload(await s.call_tool("status", {}))
        why = _payload(await s.call_tool("why", {"target": "AAA", "trait": "MONDO:1", "as_of": 2015}))
        unknown = _payload(await s.call_tool("why", {"target": "QQQ", "trait": "MONDO:1"}))
        ambiguous = _payload(await s.call_tool("known", {"node": "disease"}))
        predict = _payload(await s.call_tool("predict", {"node": "arrhythmia", "as_of": 2018, "top": 5}))
        return tools, init, status, why, unknown, ambiguous, predict

    tools, init, status, why, unknown, ambiguous, predict = asyncio.run(_session(tmp_path, work))
    assert sorted(t.name for t in tools) == NAMES
    assert all(len(t.description or "") > 80 for t in tools), "each tool says what kind of request it answers"
    assert "do not choose" in (init.instructions or "") and "not in the graph" in (init.instructions or "")
    assert "11 GWAS links" in status["answer"]
    assert why["found"] and [x["study"] for x in why["studies"]] == ["S1", "S2"] and why["source"]["file"] == "tiny.sql.gz"
    assert unknown["found"] is False and unknown["problem"] == "not in the graph"
    assert ambiguous["found"] is False and ambiguous["problem"] == "ambiguous" and len(ambiguous["candidates"]) >= 3
    assert [c["symbol"] for c in predict["candidates"]] == ["AAA", "BBB"]
    assert all("not known to be absent" in c["how_often_this_held"] for c in predict["candidates"])


def test_the_server_refuses_a_bundle_that_was_altered(tmp_path):
    ask_selftest.tiny(tmp_path)
    (tmp_path / "links.parquet").write_bytes((tmp_path / "links.parquet").read_bytes() + b" ")
    import subprocess
    r = subprocess.run([sys.executable, "-m", "pharos_gap_kg.server", "--bundle", str(tmp_path)], capture_output=True, text=True, timeout=60)
    assert r.returncode != 0 and "REFUSED" in (r.stderr + r.stdout)
