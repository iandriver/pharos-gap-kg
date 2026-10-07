"""Every location this package reads or writes. No other module names a path.

A checkout keeps its data under `<repo>/data` (git-ignored). The raw snapshots are several GB each, so
`data/raw` is usually a symlink to a larger volume, or `PHAROS_GAP_KG_RAW` names one. An installed copy
has no repo root to stand on, so it requires `PHAROS_GAP_KG_DATA`.

Run:  python -m pharos_gap_kg.paths            # print the resolved locations
      python -m pharos_gap_kg.paths --selftest
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()


def _root() -> Path | None:
    """The repo root when this file sits in a checkout (`<root>/src/pharos_gap_kg/paths.py`), else None."""
    cand = _HERE.parents[2]
    return cand if (cand / "pyproject.toml").is_file() and (cand / "src").is_dir() else None


ROOT = _root()


def _data() -> Path:
    env = os.environ.get("PHAROS_GAP_KG_DATA")
    if env:
        return Path(env)
    if ROOT is None:
        raise RuntimeError("pharos-gap-kg is installed outside a checkout: set PHAROS_GAP_KG_DATA to a data directory")
    return ROOT / "data"


DATA = _data()
RAW = Path(os.environ.get("PHAROS_GAP_KG_RAW", DATA / "raw"))   # snapshots as published, never edited
TABLES = DATA / "tables"                                      # tables/<release>/<table>.parquet
GRAPH = DATA / "graph"                                        # graph/<release>/{nodes,edges,gaps}.tsv
EVAL = DATA / "eval"                                          # frozen splits and scores
BUNDLE = Path(os.environ.get("PHAROS_GAP_KG_BUNDLE", DATA / "bundle"))   # bundle/<release>: what the front end reads
DOCS_DATA = (ROOT / "docs" / "data") if ROOT else DATA / "docs"  # small committed summaries


def tables(release: str) -> Path:
    return TABLES / release


def bundle(release: str) -> Path:
    return BUNDLE / release


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("data directory resolves", isinstance(DATA, Path))
    ck("raw, tables, graph, eval and bundle are five different places", len({RAW, TABLES, GRAPH, EVAL, BUNDLE}) == 5)
    ck("tables('x') sits under TABLES", tables("x").parent == TABLES)
    src = _HERE.parent
    offenders = []
    for p in sorted(src.glob("*.py")):
        if p.name == "paths.py":
            continue
        for i, line in enumerate(p.read_text().splitlines(), 1):
            code = line.split("#", 1)[0]
            if "/Users/" in code or "/Volumes/" in code or "/home/" in code:
                offenders.append(f"{p.name}:{i}")
    ck(f"no module hardcodes a machine path {offenders or ''}", not offenders)
    return 0 if ok else 1


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        raise SystemExit(_selftest())
    for k in ("ROOT", "DATA", "RAW", "TABLES", "GRAPH", "EVAL", "BUNDLE"):
        print(f"{k:7s} {globals()[k]}")
