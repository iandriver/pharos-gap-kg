"""The Claude Desktop installer: its config handling, and the question it puts to a server over MCP, on the six-node graph of the self-test."""
import asyncio
import os
import sys

import pytest

from pharos_gap_kg import ask_selftest, desktop


def test_the_selftest_passes():
    assert desktop._selftest() == 0


def test_the_check_asks_a_running_server_the_gaps_for_a_disease(tmp_path):
    """The six-node graph has no disease sources, so the server's answer is that it cannot give a condition's gaps, and the check would fail."""
    pytest.importorskip("mcp")
    ask_selftest.tiny(tmp_path)
    r = asyncio.run(desktop._ask(sys.executable, str(tmp_path), dict(os.environ), "heart disease"))
    assert "gaps" in r["tools"] and "subgraph" in r["tools"] and len(r["tools"]) == 9
    assert r["found"] is False and "merged graph's other parts" in r["answer"]
