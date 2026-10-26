"""The server over HTTP: the guard in front of it, and the MCP questions answered through it, on the six-node graph of the self-test."""
import asyncio
import json
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request

import pytest

mcp = pytest.importorskip("mcp")
from pharos_gap_kg import ask_selftest  # noqa: E402
from pharos_gap_kg.server import Guard  # noqa: E402


class Clock:
    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t


def test_the_guard_counts_per_client_and_in_all_over_a_minute():
    c = Clock()
    g = Guard(None, per_client=3, in_all=5, at_once=1, clock=c)
    assert [g.wait("a") for _ in range(3)] == [None, None, None]
    assert g.wait("a") == 60.0, "a fourth request inside the minute waits for the first to age out"
    c.t += 20
    assert g.wait("a") == 40.0 and g.wait("b") is None and g.wait("b") is None
    assert g.wait("c") == 40.0, "five in all inside the minute: the next client waits too"
    c.t += 41
    assert g.wait("a") is None and g.wait("c") is None, "a minute on, the first requests no longer count"


def _port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture
def serve(tmp_path):
    ask_selftest.tiny(tmp_path)
    procs = []

    def start(*extra):
        port = _port()
        p = subprocess.Popen([sys.executable, "-m", "pharos_gap_kg.server", "--http", "--port", str(port), "--bundle", str(tmp_path), *extra],
                             stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        procs.append(p)
        url = f"http://127.0.0.1:{port}"
        for _ in range(150):
            try:
                with urllib.request.urlopen(url + "/health", timeout=1) as r:
                    return url, json.loads(r.read())
            except OSError:
                if p.poll() is not None:
                    raise RuntimeError(p.stderr.read())
                time.sleep(0.2)
        raise RuntimeError("the server did not answer /health")
    yield start
    for p in procs:
        p.terminate()
        p.wait(10)


def _post(url: str, body: bytes, host: str | None = None) -> int:
    h = {"content-type": "application/json", "accept": "application/json, text/event-stream"} | ({"host": host} if host else {})
    try:
        with urllib.request.urlopen(urllib.request.Request(url + "/mcp", data=body, headers=h), timeout=10) as r:
            return r.status
    except urllib.error.HTTPError as e:
        return e.code


INIT = json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                   "params": {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "t", "version": "0"}}}).encode()


def test_the_questions_over_http(serve):
    url, health = serve()
    assert health["status"] == "ok"

    async def go():
        from mcp import ClientSession
        from mcp.client.streamable_http import streamable_http_client
        async with streamable_http_client(url + "/mcp") as streams:
            async with ClientSession(streams[0], streams[1]) as s:
                await s.initialize()
                tools = sorted(t.name for t in (await s.list_tools()).tools)
                r = await s.call_tool("why", {"target": "AAA", "trait": "MONDO:1", "as_of": 2015})
                body = getattr(r, "structured_content", None) or getattr(r, "structuredContent", None) or json.loads(r.content[0].text)
                return tools, body
    tools, why = asyncio.run(go())
    assert "why" in tools and "gaps" in tools and len(tools) == 9
    assert why["found"] and [x["study"] for x in why["studies"]] == ["S1", "S2"]


def test_requests_over_the_limit_wait(serve):
    url, _ = serve("--per-minute", "3")
    codes = [_post(url, INIT) for _ in range(5)]
    assert codes == [200, 200, 200, 429, 429]


def test_a_large_body_and_a_host_it_was_not_told_of_are_refused(serve):
    url, _ = serve("--public-host", "kg.example.org")
    assert _post(url, INIT) == 421, "Host 127.0.0.1 is not the public name"
    assert _post(url, INIT, host="kg.example.org") == 200
    big = json.dumps({"jsonrpc": "2.0", "id": 2, "method": "ping", "params": {"pad": "x" * 70_000}}).encode()
    assert _post(url, big, host="kg.example.org") == 413


def test_a_public_address_needs_its_names():
    p = subprocess.run([sys.executable, "-m", "pharos_gap_kg.server", "--http", "--host", "0.0.0.0", "--bundle", "nowhere"], capture_output=True, text=True, timeout=60)
    assert p.returncode != 0 and "--public-host" in p.stderr
