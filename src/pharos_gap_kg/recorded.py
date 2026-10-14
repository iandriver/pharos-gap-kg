"""Answers from live services, kept as they came, so that a rebuild reproduces the graph without asking again.

Three stages ask a service whose answers move with time. `pubyear` asks NCBI for the year of each cited
paper. `prune2 search` and `applycut search` ask PubMed how many papers mention a gene, a disease and
both, and the lift they compute divides by PubMed's size on the day. Run again later, they withhold a
different set of rows. Their outputs are committed under `recorded/<release>/` with a digest of each
file, and `restore` copies them to where the stages write them. After that `pubyear` finds nothing to
ask, and the two searches need not be run: `applycut write` reads what was recorded.

A file is copied, never rewritten, so its bytes and digest are those recorded. `restore` refuses to
replace a file already on disk that differs from the recorded one: that is an answer asked again.

Run:  python -m pharos_gap_kg.recorded restore pharos319
      python -m pharos_gap_kg.recorded save pharos319       # what the stages wrote, recorded for others
      python -m pharos_gap_kg.recorded --selftest
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

from pharos_gap_kg import paths


def places(release: str) -> dict[str, Path]:
    """recorded file name -> where the stage that asks for it writes it."""
    work = paths.EVAL / "pruning2" / release
    return {"pubmed_year.parquet": paths.TABLES / "pubmed_year.parquet",
            "search.parquet": work / "search.parquet", "search.json": work / "search.json",
            "applied.parquet": work / "applied.parquet", "applied.json": work / "applied.json"}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(release: str) -> dict:
    d = paths.RECORDED / release
    d.mkdir(parents=True, exist_ok=True)
    man = {}
    for name, src in places(release).items():
        if not src.is_file():
            raise SystemExit(f"nothing to record: {src} is not on disk")
        shutil.copyfile(src, d / name)
        man[name] = sha(d / name)
    (d / "manifest.json").write_text(json.dumps({"release": release, "files": man}, indent=1) + "\n")
    return man


def restore(release: str) -> dict:
    d = paths.RECORDED / release
    if not (d / "manifest.json").is_file():
        raise SystemExit(f"no answers are recorded for {release} under {d}")
    man = json.loads((d / "manifest.json").read_text())["files"]
    out = {}
    for name, dest in places(release).items():
        if name not in man or sha(d / name) != man[name]:
            raise SystemExit(f"REFUSED: {d / name} is not the file its manifest records")
        if dest.is_file():
            if sha(dest) != man[name]:
                raise SystemExit(f"REFUSED: {dest} is on disk and differs from the recorded answer; it was asked again. "
                                 f"Move it aside to restore the recorded one.")
            out[name] = "already in place"
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(d / name, dest)
        out[name] = "restored"
    return out


def _selftest() -> int:
    import tempfile
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    old = (paths.EVAL, paths.TABLES, paths.RECORDED)
    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        paths.EVAL, paths.TABLES, paths.RECORDED = t / "eval", t / "tables", t / "recorded"
        try:
            for i, p in enumerate(places("r").values()):
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(bytes([i]) * 10)
            man = save("r")
            ck("every answer is recorded with the digest of its bytes", len(man) == 5 and all(sha(paths.RECORDED / "r" / n) == h for n, h in man.items()))
            ck("what is already in place is left alone", set(restore("r").values()) == {"already in place"})
            for p in places("r").values():
                p.unlink()
            ck("a clean data directory gets every answer back, byte for byte", set(restore("r").values()) == {"restored"}
               and all(sha(p) == man[n] for n, p in places("r").items()))
            places("r")["applied.parquet"].write_bytes(b"asked again")
            try:
                restore("r")
                ck("an answer asked again is not overwritten", False)
            except SystemExit:
                ck("an answer asked again is not overwritten", places("r")["applied.parquet"].read_bytes() == b"asked again")
            places("r")["applied.parquet"].unlink()
            (paths.RECORDED / "r" / "search.json").write_bytes(b"edited")
            try:
                restore("r")
                ck("a recorded file that does not match its manifest is refused", False)
            except SystemExit:
                ck("a recorded file that does not match its manifest is refused", True)
        finally:
            paths.EVAL, paths.TABLES, paths.RECORDED = old
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 2 or argv[0] not in ("restore", "save"):
        raise SystemExit(__doc__)
    out = save(argv[1]) if argv[0] == "save" else restore(argv[1])
    for k, v in out.items():
        print(f"  {k:22s} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
