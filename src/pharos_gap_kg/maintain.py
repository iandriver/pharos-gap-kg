"""Keeping the graph current: which outputs a changed input makes stale, and re-deriving only those.

Each stage of the README's build block reads some files and writes others. `trace` runs the build
once with every file a stage opens recorded: the files Python opens, and the parquet files pyarrow
reads and writes itself. That gives the map of the build (docs/data/build_graph_<release>.json) and
the digest of every file as each stage saw it (data/stamps/<release>.json). A stage is stale when a
file it read has changed since, or when a stage before it that writes one of its inputs is stale.
`rebuild` runs the stale stages in build order and skips a stage whose inputs, once the stages before
it have run, are unchanged; it stamps as it goes. The definitions are in docs/PROTOCOL.md, "Keeping
the graph current".

Run, in a built checkout:
      python -m pharos_gap_kg.maintain trace pharos319      # the build under the tracer: the map and the stamps
      python -m pharos_gap_kg.maintain stale pharos319      # what is stale now, and why
      python -m pharos_gap_kg.maintain whatif pharos319 data:eval/pruning2/pharos319/applied.parquet
      python -m pharos_gap_kg.maintain rebuild pharos319    # only the stale stages, in order
      python -m pharos_gap_kg.maintain digests pharos319 OUT.json   # every mapped file's digest, to compare two checkouts
The demonstration:
      python -m pharos_gap_kg.maintain change pharos319     # PubMed asked again: 20 pairs lose one paper
      python -m pharos_gap_kg.maintain demo pharos319       # change, stale, rebuild, and the record of it
      python -m pharos_gap_kg.maintain --selftest
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from pharos_gap_kg import paths


def stages(release: str) -> list[dict]:
    """The README's build block, in order. `isolated` stages read a downloaded archive and run with -I."""
    s = lambda name, *args, isolated=False: {"name": name, "args": list(args), "isolated": isolated}
    return [s("split", "split", release, "--cutoff", "2018"), s("graph", "graph", "build", release, "--cutoff", "2018"),
            s("ctxtables", "ctxtables", release), s("bundle_links", "bundle", release), s("planted", "planted", release),
            s("evidence", "evidence", release), s("prune2_freeze", "prune2", "freeze", release), s("prune2_ingraph", "prune2", "ingraph", release),
            s("applycut", "applycut", "write", release), s("tgaps", "tgaps", release), s("tgapgraph", "tgapgraph", release),
            s("ctxgraph", "ctxgraph", release), s("kg", "kg", "build", release),
            s("gwascat_extract", "gwascat", "extract", isolated=True), s("gwascat_rule", "gwascat", "rule", release), s("gwascat_build", "gwascat", "build", release),
            s("studysize", "studysize", "build", release, isolated=True), s("bundle", "bundle", release),
            s("census", "diseasegaps", "census", release), s("archive", "release", "archive", release)]


# ---- files as tokens, so a map made in one checkout reads in another -----------------------------------------

def _roots() -> list[tuple[str, Path]]:
    out = [("raw", paths.RAW.resolve()), ("data", paths.DATA.resolve())]
    if paths.ROOT is not None:
        out += [("repo", paths.ROOT.resolve())]
    return out


def token(p) -> str | None:
    """'raw:…', 'data:…' or 'repo:…' for a file this build owns or reads; None for anything else."""
    try:
        q = Path(os.fsdecode(p)).resolve()
    except (TypeError, ValueError, OSError):
        return None
    for name, root in _roots():
        try:
            rel = q.relative_to(root)
        except ValueError:
            continue
        r = rel.as_posix()
        if name == "repo" and not (r.startswith("docs/data/") or r.startswith("recorded/")):
            return None
        if name == "data" and (r.startswith("stamps/") or r.startswith("raw/") or r.startswith("eval/literature/")):
            return None
        return f"{name}:{r}"
    return None


def path_of(tok: str) -> Path:
    name, rel = tok.split(":", 1)
    return dict(_roots())[name] / rel


def sha(tok: str) -> str | None:
    p = path_of(tok)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with p.open("rb") as f:
        while b := f.read(1 << 20):
            h.update(b)
    return h.hexdigest()


# ---- the tracer -------------------------------------------------------------------------------------------

def _run(log: str, module: str, args: list[str]) -> int:
    """Run one stage in this process with every file it opens recorded to `log`."""
    import runpy
    reads, writes = set(), set()

    def note(p, w: bool) -> None:
        t = token(p) if isinstance(p, (str, bytes, os.PathLike)) else None
        if t:
            (writes if w else reads).add(t)

    def hook(event, a):
        if event == "open" and a:
            mode, flags = (a[1] if len(a) > 1 else None), (a[2] if len(a) > 2 else 0)
            w = any(c in mode for c in "wax+") if isinstance(mode, str) else bool((flags or 0) & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND | os.O_TRUNC))
            note(a[0], w)
        elif event in ("os.rename", "os.replace") and len(a) > 1:
            note(a[1], True)

    sys.addaudithook(hook)
    import pyarrow.parquet as pq

    def wrap(f, w, pos):
        def g(*a, **k):
            note(a[pos] if len(a) > pos else k.get("where" if w else "source"), w)
            return f(*a, **k)
        return g

    pq.read_table, pq.read_metadata, pq.read_schema = (wrap(pq.read_table, False, 0), wrap(pq.read_metadata, False, 0), wrap(pq.read_schema, False, 0))
    pq.write_table = wrap(pq.write_table, True, 1)
    _pf_init, _pw_init = pq.ParquetFile.__init__, pq.ParquetWriter.__init__
    pq.ParquetFile.__init__ = lambda self, source, *a, **k: (note(source, False), _pf_init(self, source, *a, **k))[1]
    pq.ParquetWriter.__init__ = lambda self, where, *a, **k: (note(where, True), _pw_init(self, where, *a, **k))[1]
    sys.argv = [module] + args
    code = 0
    try:
        runpy.run_module(f"pharos_gap_kg.{module}", run_name="__main__", alter_sys=True)
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
    finally:
        Path(log).write_text(json.dumps({"reads": sorted(reads - writes), "writes": sorted(writes), "read_and_written": sorted(reads & writes), "exit": code}))
    return code


def run_stage(st: dict, log: Path) -> dict:
    cmd = [sys.executable] + (["-I"] if st["isolated"] else []) + ["-m", "pharos_gap_kg.maintain", "_run", str(log), st["args"][0], *st["args"][1:]]
    t0 = time.time()
    p = subprocess.run(cmd, cwd=paths.ROOT, capture_output=True, text=True)
    out = json.loads(log.read_text()) if log.is_file() else {"exit": p.returncode, "reads": [], "writes": []}
    if out["exit"] != 0:
        raise SystemExit(f"REFUSED: stage {st['name']} failed ({out['exit']}):\n{(p.stdout + p.stderr)[-2000:]}")
    return {**out, "seconds": round(time.time() - t0, 1)}


# ---- the map, the stamps, and what is stale ----------------------------------------------------------------

def map_path(release: str) -> Path:
    return paths.DOCS_DATA / f"build_graph_{release}.json"


def stamps_path(release: str) -> Path:
    return paths.DATA / "stamps" / f"{release}.json"


def producers(gmap: list[dict]) -> dict[str, dict[str, str]]:
    """stage -> {input file: the stage before it that last wrote that file}; external inputs have none."""
    last, out = {}, {}
    for st in gmap:
        out[st["name"]] = {f: last[f] for f in st["reads"] if f in last}
        for f in st["writes"]:
            last[f] = st["name"]
    return out


def stale(gmap: list[dict], stamps: dict, digest=sha) -> dict[str, str]:
    """stage -> why it is stale, for every stale stage, in build order."""
    prod, out = producers(gmap), {}
    for st in gmap:
        n, seen = st["name"], stamps.get(st["name"], {}).get("inputs", {})
        up = sorted({p for f, p in prod[n].items() if p in out})
        changed = [f for f in st["reads"] if f not in prod[n] and digest(f) != seen.get(f)]
        missing = [f for f in st["writes"] if digest(f) is None]
        if up:
            out[n] = "reads what " + ", ".join(up) + " writes"
        elif changed:
            out[n] = "changed: " + ", ".join(changed[:3]) + (f" and {len(changed) - 3} more" if len(changed) > 3 else "")
        elif missing:
            out[n] = "missing: " + ", ".join(missing[:3])
    return out


def whatif(gmap: list[dict], tok: str) -> list[str]:
    """The stages a change to one file would make stale."""
    return list(stale(gmap, {st["name"]: {"inputs": {f: "same" for f in st["reads"]}} for st in gmap}, digest=lambda f: "changed" if f == tok else "same"))


def trace(release: str) -> dict:
    """Run the build under the tracer: write the map, and the stamps of what each stage saw and wrote."""
    logs = paths.DATA / "stamps" / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    gmap, stamps = [], {}
    for st in stages(release):
        r = run_stage(st, logs / f"{st['name']}.json")
        gmap.append({**st, "reads": r["reads"], "writes": r["writes"]})
        stamps[st["name"]] = {"inputs": {f: sha(f) for f in r["reads"]}, "outputs": {f: sha(f) for f in r["writes"]}, "seconds": r["seconds"]}
        print(f"  {st['name']:16s} reads {len(r['reads']):3d}  writes {len(r['writes']):3d}  {r['seconds']:7.1f} s", flush=True)
    map_path(release).write_text(json.dumps({"release": release, "stages": gmap}, indent=1) + "\n")
    stamps_path(release).write_text(json.dumps(stamps, indent=1) + "\n")
    return {"stages": len(gmap)}


def rebuild(release: str) -> list[dict]:
    """Run, in build order, each stage whose inputs differ from what it last saw; stamp as each one runs."""
    gmap = json.loads(map_path(release).read_text())["stages"]
    stamps = json.loads(stamps_path(release).read_text())
    logs, done = paths.DATA / "stamps" / "logs", []
    last = {f: st["name"] for st in gmap for f in st["writes"]}           # the stage whose version of a file the build keeps
    known = {f for st in gmap for f in st["reads"] + st["writes"]}
    for i, st in enumerate(gmap):
        seen = stamps[st["name"]]["inputs"]
        now = {f: sha(f) for f in st["reads"]}
        kept = {f: h for f, h in stamps[st["name"]]["outputs"].items() if last.get(f) == st["name"]}
        if now == seen and all(sha(f) == h for f, h in kept.items()):
            continue
        prev = stamps[st["name"]]["outputs"]
        r = run_stage(st, logs / f"{st['name']}.json")
        later = {f for s2 in gmap[i + 1:] for f in s2["writes"]}
        if set(r["writes"]) != set(st["writes"]) or not set(r["reads"]) <= set(st["reads"]) | later or not set(r["reads"]) <= known:
            raise SystemExit(f"REFUSED: {st['name']} read or wrote files the map does not give it; trace the build again")
        stamps[st["name"]] = {"inputs": now, "outputs": {f: sha(f) for f in r["writes"]}, "seconds": r["seconds"]}
        stamps_path(release).write_text(json.dumps(stamps, indent=1) + "\n")
        done.append({"stage": st["name"], "seconds": r["seconds"], "outputs_changed": sorted(f for f, h in stamps[st["name"]]["outputs"].items() if h != prev.get(f))})
        print(f"  ran {st['name']:16s} {r['seconds']:7.1f} s", flush=True)
    return done


def digests(release: str) -> dict[str, str | None]:
    gmap = json.loads(map_path(release).read_text())["stages"]
    return {f: sha(f) for f in sorted({f for st in gmap for f in st["reads"] + st["writes"]})}


# ---- the demonstration's change ----------------------------------------------------------------------------

CHANGED_PAIRS = 20


def change(release: str) -> dict:
    """PubMed asked again: the CTD one-source pairs outside the sample with the smallest lift above CTD's
    cut each lose one paper that mentions both; lift and verdict are recomputed as `applycut search` would."""
    import pandas as pd

    from pharos_gap_kg import applycut, litcheck, prune2
    W = prune2.work(release)
    a = pd.read_parquet(W / "applied.parquet")
    total = json.loads((W / "applied.json").read_text())["pubmed_records"]
    sampled = set(map(tuple, pd.read_parquet(W / "search.parquet")[prune2.KEY].to_numpy()))
    keyed = list(map(tuple, a[prune2.KEY].to_numpy()))
    ok = (a.source == "CTD") & a.unsearchable.isna() & ~a.lit_withheld.astype(bool) & a.both.fillna(0).gt(0) & ~pd.Series([k in sampled for k in keyed], index=a.index)
    pick = a[ok].sort_values(["lift"] + prune2.KEY).head(CHANGED_PAIRS).index
    a.loc[pick, "both"] = a.loc[pick, "both"] - 1
    a.loc[pick, "ids"] = [("|".join(i.split("|")[:-1])) if isinstance(i, str) else i for i in a.loc[pick, "ids"]]
    a.loc[pick, "lift"] = [litcheck.lift(int(b), int(g), int(d), total) for b, g, d in zip(a.loc[pick, "both"], a.loc[pick, "gene_papers"], a.loc[pick, "disease_papers"])]
    a["lit_withheld"] = applycut.decide(a.lift, a.unsearchable, a.cut)
    tmp = W / "applied.parquet.tmp"
    a.to_parquet(tmp, index=False)
    os.replace(tmp, W / "applied.parquet")
    return {"pairs_changed": int(len(pick)), "now_withheld_by_the_cut": int(a.loc[pick, "lit_withheld"].sum())}


def demo(release: str) -> dict:
    gmap = json.loads(map_path(release).read_text())["stages"]
    before = digests(release)
    what = change(release)
    st_ = stale(gmap, json.loads(stamps_path(release).read_text()))
    t0 = time.time()
    ran = rebuild(release)
    after = digests(release)
    outs = {f for r in ran for st in gmap if st["name"] == r["stage"] for f in st["writes"]}
    changed = sorted(f for f in after if after[f] != before.get(f))
    out = {"release": release, "change": {**what, "file": "data:eval/pruning2/" + release + "/applied.parquet"},
           "stages": len(gmap), "stale": st_, "ran": [r["stage"] for r in ran], "seconds_by_stage": {r["stage"]: r["seconds"] for r in ran},
           "seconds": round(time.time() - t0, 1), "files_mapped": len(after), "files_changed": changed,
           "every_changed_file_is_an_output_of_a_stage_that_ran_or_the_change": all(f in outs or f == "data:eval/pruning2/" + release + "/applied.parquet" for f in changed),
           "digests_after": after}
    (paths.DOCS_DATA / f"maintenance_demo_{release}.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def _selftest() -> int:
    ok = True

    def ck(label: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {label}")

    g = [{"name": "a", "reads": ["raw:x"], "writes": ["data:A"]}, {"name": "b", "reads": ["data:A", "raw:y"], "writes": ["data:B"]},
         {"name": "c", "reads": ["raw:z"], "writes": ["data:C"]}, {"name": "d", "reads": ["data:B", "data:C"], "writes": ["data:D"]},
         {"name": "e", "reads": ["data:A"], "writes": ["data:A2"]}]
    ck("a stage's inputs come from the last stage before it that wrote them", producers(g)["d"] == {"data:B": "b", "data:C": "c"} and producers(g)["a"] == {})
    st = {s["name"]: {"inputs": {f: "h" for f in s["reads"]}} for s in g}
    fresh = lambda f: "h"
    ck("nothing changed, nothing stale", stale(g, st, digest=fresh) == {})
    s1 = stale(g, st, digest=lambda f: "new" if f == "raw:y" else "h")
    ck("a changed input makes its stage stale, and every stage downstream; the others stay fresh", list(s1) == ["b", "d"] and "raw:y" in s1["b"] and "b" in s1["d"])
    ck("what a change to one file would make stale", whatif(g, "raw:x") == ["a", "b", "d", "e"] and whatif(g, "raw:z") == ["c", "d"])
    ck("an output that is missing makes its stage stale, and what reads it", list(stale(g, st, digest=lambda f: None if f == "data:C" else "h")) == ["c", "d"])
    import tempfile
    with tempfile.TemporaryDirectory() as t:
        d = Path(t) / "data"
        (d / "raw").mkdir(parents=True)
        import pandas as pd
        pd.DataFrame({"x": [1, 2]}).to_parquet(d / "raw" / "in.parquet")
        (d / "raw" / "in.txt").write_text("hello")
        log = Path(t) / "log.json"
        env = {**os.environ, "PHAROS_GAP_KG_DATA": str(d), "PHAROS_GAP_KG_RAW": str(d / "raw")}
        p = subprocess.run([sys.executable, "-m", "pharos_gap_kg.maintain", "_run", str(log), "maintain", "_touch"], cwd=paths.ROOT, env=env, capture_output=True, text=True)
        r = json.loads(log.read_text()) if log.is_file() else {}
        ck("the tracer sees a parquet file pyarrow reads, a text file Python reads, and the files written",
           r.get("reads") == ["raw:in.parquet", "raw:in.txt"] and r.get("writes") == ["data:out/a.parquet", "data:out/b.json"] and r.get("exit") == 0)
    return 0 if ok else 1


def _touch() -> None:
    """For the self-test: read one parquet and one text file from RAW, write one of each to DATA."""
    import pandas as pd
    x = pd.read_parquet(paths.RAW / "in.parquet")
    t = (paths.RAW / "in.txt").read_text()
    (paths.DATA / "out").mkdir(exist_ok=True)
    x.to_parquet(paths.DATA / "out" / "a.parquet")
    (paths.DATA / "out" / "b.json").write_text(json.dumps({"t": t}))


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if argv and argv[0] == "_run":
        return _run(argv[1], argv[2], argv[3:])
    if argv == ["_touch"]:
        _touch()
        return 0
    if len(argv) < 2:
        raise SystemExit(__doc__)
    cmd, release = argv[0], argv[1]
    if cmd == "trace":
        print(trace(release))
    elif cmd == "stale":
        s = stale(json.loads(map_path(release).read_text())["stages"], json.loads(stamps_path(release).read_text()))
        print("  nothing is stale" if not s else "\n".join(f"  {k:16s} {v}" for k, v in s.items()))
    elif cmd == "whatif":
        print("  " + ", ".join(whatif(json.loads(map_path(release).read_text())["stages"], argv[2])))
    elif cmd == "rebuild":
        print(f"  {len(rebuild(release))} stage(s) ran")
    elif cmd == "digests":
        Path(argv[2]).write_text(json.dumps(digests(release), indent=1) + "\n")
    elif cmd == "change":
        print(change(release))
    elif cmd == "demo":
        o = demo(release)
        print(f"  stale {list(o['stale'])}; ran {o['ran']} in {o['seconds']} s; {len(o['files_changed'])} of {o['files_mapped']} files changed")
    else:
        raise SystemExit(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
