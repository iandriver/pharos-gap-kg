"""Install the question server for Claude Desktop on macOS in one command, and check that it answers.

Claude Desktop starts an MCP server as a program of its own, which macOS does not let read
~/Downloads, ~/Documents or ~/Desktop, so a server run from a checkout there does not start. `install`
puts a copy of the package in its own environment and the release's bundle beside it, outside those
folders (~/.pharos-gap-kg unless --home names another place), and adds the server to the app's config.
The app keeps its config in memory and writes it back when its settings change, so an entry added while
it runs is lost: the config is written only while the app is closed. Run `install` again after quitting
the app; the copy and the bundle already in place are kept, and the package is reinstalled from the
checkout.

`check` starts the installed server over MCP with those three folders unreadable, as the app would,
asks it the Parkinson disease question, and says whether the config holds the entry and whether the app
has started the server.

Run, from a checkout with its environment active:
      python -m pharos_gap_kg.desktop install     # the copy, the bundle, then the config entry (quit Claude first)
      python -m pharos_gap_kg.desktop check
      python -m pharos_gap_kg.desktop remove      # the config entry only; the copy stays until its folder is deleted
      python -m pharos_gap_kg.desktop --selftest
Options: --home DIR (default ~/.pharos-gap-kg), --release NAME (default pharos319).
"""
from __future__ import annotations

import argparse
import asyncio
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from pharos_gap_kg import paths, release

NAME = "pharos-gap-kg"
RELEASE = "pharos319"
HOME, CONFIG, LOGS, PROTECTED = paths.DESKTOP_HOME, paths.DESKTOP_CONFIG, paths.DESKTOP_LOGS, paths.DESKTOP_UNREADABLE
APP = "/Claude.app/Contents/MacOS/Claude"
QUESTION = "Parkinson disease"


def python_in(home: Path) -> Path:
    return home / "venv" / "bin" / "python"


def bundle_in(home: Path, rel: str) -> Path:
    return home / "data" / "bundle" / rel


def entry(home: Path, rel: str) -> dict:
    return {"command": str(python_in(home)), "args": ["-m", "pharos_gap_kg.server", "--bundle", str(bundle_in(home, rel))],
            "env": {"PHAROS_GAP_KG_DATA": str(home / "data")}}


def app_running() -> bool:
    out = subprocess.run(["ps", "-axo", "comm="], capture_output=True, text=True).stdout
    return any(line.endswith(APP) for line in out.splitlines())


def _run(cmd: list, **kw) -> None:
    print("  $ " + " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def install_copy(home: Path) -> None:
    """An environment of its own at home/venv with the package and the server's library, at the checkout's pinned versions."""
    if paths.ROOT is None:
        raise SystemExit("REFUSED: run this from a checkout of pharos-gap-kg; the copy is installed from it")
    py, uv = python_in(home), shutil.which("uv")
    if not py.is_file():
        home.mkdir(parents=True, exist_ok=True)
        _run([uv, "venv", home / "venv", "--python", sys.executable] if uv else [sys.executable, "-m", "venv", home / "venv"])
    spec, pins = f"{paths.ROOT}[server]", paths.ROOT / "constraints.txt"
    if uv:
        _run([uv, "pip", "install", "--python", py, spec, "-c", pins, "--reinstall-package", NAME])
    else:
        _run([py, "-m", "pip", "install", "--quiet", spec, "-c", pins])
        _run([py, "-m", "pip", "install", "--quiet", "--force-reinstall", "--no-deps", paths.ROOT])


def install_bundle(home: Path, rel: str) -> str:
    """The release's bundle at home/data/bundle/<release>, from the checkout's archive if it has the committed one, else downloaded; every digest checked."""
    if (bundle_in(home, rel) / "manifest.json").is_file():
        return "already in place"
    env = {k: v for k, v in os.environ.items() if k not in ("PHAROS_GAP_KG_BUNDLE",)} | {"PHAROS_GAP_KG_DATA": str(home / "data")}
    local = paths.RELEASE / release.name(rel)
    rec = json.loads(release.record_path(rel).read_text()) if release.record_path(rel).is_file() else {}
    if local.is_file() and rec and release.sha(local.read_bytes()) == rec["sha256"]:
        _run([sys.executable, "-m", "pharos_gap_kg.release", "unpack", rel, local], env=env)
        return f"unpacked from {local}"
    _run([sys.executable, "-m", "pharos_gap_kg.release", "get", rel], env=env)
    return "downloaded"


def write_entry(config: Path, value: dict | None, running: bool) -> str:
    """Put value under mcpServers[NAME] (None removes it), keeping every other key; a dated backup first. Refuses while the app runs."""
    current = json.loads(config.read_text()) if config.is_file() else {}
    if current.get("mcpServers", {}).get(NAME) == value:
        return "already in place" if value else "not there"
    if running:
        raise SystemExit("REFUSED: Claude is open, and it writes its own copy of the config back over this file. "
                         "Quit it (Cmd-Q), run the same command again, then open it.")
    if not config.parent.is_dir():
        raise SystemExit(f"REFUSED: {config.parent} does not exist. Install Claude Desktop, open it once, quit it, and run this again.")
    if config.is_file():
        shutil.copy2(config, config.with_name(f"{config.stem}.backup-{time.strftime('%Y%m%d-%H%M%S')}.json"))
    servers = current.setdefault("mcpServers", {})
    if value is None:
        servers.pop(NAME, None)
    else:
        servers[NAME] = value
    tmp = config.with_name(config.name + ".tmp")
    tmp.write_text(json.dumps(current, indent=2) + "\n")
    os.replace(tmp, config)
    return "added" if value else "removed"


def started(logs: Path) -> dict:
    """When the app last started the server, by its own log, and the last problem it logged after that."""
    lines = [l for l in (logs / "mcp.log").read_text(errors="replace").splitlines() if f"[{NAME}]" in l] if (logs / "mcp.log").is_file() else []
    ok = [l for l in lines if "Server started and connected successfully" in l]
    bad = [l for l in lines if re.search(r"\[(error|warn)\]", l)]
    out = {"started": ok[-1][:24] if ok else None}
    if bad and (not ok or bad[-1] > ok[-1]):
        out["problem"] = bad[-1][:300]
    return out


def sandbox_profile(folders: list[Path]) -> str:
    return "(version 1)(allow default)(deny file-read* " + " ".join(f'(subpath "{f}")' for f in folders) + ")"


async def _ask(py: str, bundle: str, env: dict, disease: str) -> dict:
    from mcp import ClientSession, StdioServerParameters, stdio_client
    params = StdioServerParameters(command=py, args=["-m", "pharos_gap_kg.server", "--bundle", bundle], env=env, cwd=str(Path(py).parents[2]))
    async with stdio_client(params) as streams:
        async with ClientSession(streams[0], streams[1]) as s:
            await s.initialize()
            tools = sorted(t.name for t in (await s.list_tools()).tools)
            r = await s.call_tool("gaps", {"disease": disease})
            body = getattr(r, "structured_content", None) or getattr(r, "structuredContent", None) or json.loads(r.content[0].text)
            return {"tools": tools, "answer": body.get("answer", ""), "found": body.get("found")}


def ask_server(py: Path, bundle: Path, data: Path, disease: str = QUESTION) -> dict:
    """Start the server as the app would and ask it one question; run by the copy's own Python, inside a sandbox that cannot read the protected folders."""
    env = {"PHAROS_GAP_KG_DATA": str(data), "PATH": "/usr/bin:/bin", "HOME": str(Path.home())}
    cmd = [str(py), "-m", "pharos_gap_kg.desktop", "_ask", str(bundle), str(data), disease]
    sandbox = shutil.which("sandbox-exec") if sys.platform == "darwin" else None
    if sandbox:
        cmd = [sandbox, "-p", sandbox_profile(PROTECTED)] + cmd
    p = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(py.parents[2]), timeout=300)
    if p.returncode != 0:
        return {"error": (p.stderr or p.stdout).strip()[-600:], "sandboxed": bool(sandbox)}
    return json.loads(p.stdout.strip().splitlines()[-1]) | {"sandboxed": bool(sandbox)}


def check(home: Path, rel: str, config: Path = CONFIG, logs: Path = LOGS) -> int:
    py, bundle, good = python_in(home), bundle_in(home, rel), True
    print(f"  the copy: {py if py.is_file() else 'missing, run install'}")
    print(f"  the bundle: {bundle if (bundle / 'manifest.json').is_file() else 'missing, run install'}")
    if py.is_file() and (bundle / "manifest.json").is_file():
        r = ask_server(py, bundle, home / "data")
        where = "with Downloads, Documents and Desktop unreadable" if r.get("sandboxed") else "unsandboxed"
        if "error" in r or not r.get("found"):
            good = False
            print(f"  the server, {where}: FAILED\n    {r.get('error') or r}")
        else:
            print(f"  the server, {where}: {len(r['tools'])} tools; asked the gaps for {QUESTION}: {r['answer'].split('; ')[0][:160]}")
    else:
        good = False
    there = (json.loads(config.read_text()) if config.is_file() else {}).get("mcpServers", {}).get(NAME)
    state = "in place" if there == entry(home, rel) else "missing" if there is None else "differs from this install"
    good &= state == "in place"
    print(f"  Claude Desktop's config: {state}")
    s = started(logs)
    print(f"  started by the app: {s['started'] + ' UTC' if s['started'] else 'not yet'}" + (f"\n  the app's last problem: {s['problem']}" if "problem" in s else ""))
    if good and not s["started"]:
        print("  Open Claude Desktop and ask in a new chat: \"Using pharos-gap-kg, what are the gaps for Parkinson disease?\" "
              "The reply should show a call to the gaps tool before the answer.")
    return 0 if good else 1


def install(home: Path, rel: str) -> int:
    if sys.platform != "darwin":
        raise SystemExit("REFUSED: this installer is for Claude Desktop on macOS; docs/ASK.md gives the config entry for other clients")
    print(f"1. The copy, in {home / 'venv'}")
    install_copy(home)
    print(f"2. The bundle: {install_bundle(home, rel)}")
    print(f"3. Claude Desktop's config, {CONFIG}")
    try:
        print(f"   {write_entry(CONFIG, entry(home, rel), app_running())}")
    except SystemExit as e:
        print(f"   {e}")
    print("4. The check")
    return check(home, rel)


def _selftest() -> int:
    import tempfile
    ok = True

    def ck(label: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {label}")

    def refused(fn) -> bool:
        try:
            fn()
            return False
        except SystemExit as e:
            return str(e).startswith("REFUSED")

    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        cfg, home = t / "Claude" / "claude_desktop_config.json", t / "home"
        e = entry(home, "r")
        ck("without the app's folder nothing is written", refused(lambda: write_entry(cfg, e, False)) and not cfg.parent.exists())
        cfg.parent.mkdir()
        ck("a config the app has not written yet is created", write_entry(cfg, e, False) == "added" and json.loads(cfg.read_text()) == {"mcpServers": {NAME: e}})
        cfg.write_text(json.dumps({"preferences": {"x": 1}, "mcpServers": {"other": {"command": "o"}}}))
        before = cfg.read_bytes()
        ck("while the app runs nothing is written and no backup is made", refused(lambda: write_entry(cfg, e, True))
           and cfg.read_bytes() == before and len(list(cfg.parent.iterdir())) == 1)
        ck("the entry is added beside every other key, after a backup", write_entry(cfg, e, False) == "added"
           and json.loads(cfg.read_text()) == {"preferences": {"x": 1}, "mcpServers": {"other": {"command": "o"}, NAME: e}}
           and [p.read_bytes() for p in cfg.parent.glob("*.backup-*.json")] == [before])
        ck("an entry already in place is left, even while the app runs", write_entry(cfg, e, True) == "already in place")
        ck("the entry runs the copy's Python on the copy's bundle", e["command"] == str(home / "venv/bin/python") and e["args"][-1] == str(home / "data/bundle/r")
           and e["env"] == {"PHAROS_GAP_KG_DATA": str(home / "data")})
        time.sleep(1.1)
        ck("remove takes out the entry and keeps the rest", write_entry(cfg, None, False) == "removed"
           and json.loads(cfg.read_text()) == {"preferences": {"x": 1}, "mcpServers": {"other": {"command": "o"}}} and len(list(cfg.parent.glob("*.backup-*.json"))) == 2)
        logs = t / "logs"
        logs.mkdir()
        ck("no log: never started", started(logs) == {"started": None})
        (logs / "mcp.log").write_text(f"2026-10-21T14:25:28.368Z [info] [Other] Server started and connected successfully\n"
                                      f"2026-10-21T14:26:00.000Z [info] [{NAME}] Server started and connected successfully\n"
                                      f"2026-10-21T14:27:00.000Z [error] [{NAME}] Server disconnected\n")
        s = started(logs)
        ck("the app's log gives the last start, and a problem logged after it", s["started"] == "2026-10-21T14:26:00.000Z" and "disconnected" in s["problem"])
        p = sandbox_profile([t / "Downloads", t / "Desktop"])
        ck("the sandbox denies reading each protected folder", p.count("(subpath ") == 2 and f'"{t / 'Downloads'}"' in p and "deny file-read*" in p)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    sys.stdout.reconfigure(line_buffering=True)
    if "--selftest" in argv:
        return _selftest()
    if argv[:1] == ["_ask"]:
        print(json.dumps(asyncio.run(_ask(sys.executable, argv[1], dict(os.environ, PHAROS_GAP_KG_DATA=argv[2]), argv[3]))))
        return 0
    ap = argparse.ArgumentParser(prog="python -m pharos_gap_kg.desktop", description=__doc__.split("\n\n")[0])
    ap.add_argument("command", choices=["install", "check", "remove"])
    ap.add_argument("--home", type=Path, default=HOME)
    ap.add_argument("--release", default=RELEASE)
    a = ap.parse_args(argv)
    home = a.home.expanduser().resolve()
    if any(home == f or f in home.parents for f in PROTECTED):
        raise SystemExit(f"REFUSED: {home} is inside a folder Claude Desktop's servers cannot read; choose another --home")
    if a.command == "install":
        return install(home, a.release)
    if a.command == "check":
        return check(home, a.release)
    print(f"  {write_entry(CONFIG, None, app_running())}; the copy and the bundle stay in {home} until that folder is deleted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
