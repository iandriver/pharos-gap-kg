"""The bundle as one archive, to attach to a published release, and to get back without building anything.

The questions read a bundle of about 20 MB. Building it takes 2.5 GB of published files and a few
minutes; getting it takes one command. `archive` writes the built bundle as a gzipped tar with fixed
names, times and owners, and commits its digest and each member's digest to
`docs/data/bundle_archive_<release>.json`. `get` downloads the archive attached to the repository's
latest release, and `unpack` installs one already on disk. Both refuse an archive whose digest is not
the committed one, any member that is not a plain file named in the record, and a member whose digest
differs, and they write nothing until every member has been checked. The questions then check each file
against the bundle's own manifest, as they do for a bundle built here.

Run:  python -m pharos_gap_kg.release archive pharos319        # the author's side: data/release/<file>
      python -m pharos_gap_kg.release get pharos319            # download, check, unpack into data/bundle
      python -m pharos_gap_kg.release get pharos319 --from URL # the same archive from another place
      python -m pharos_gap_kg.release unpack pharos319 FILE    # the same, for an archive already downloaded
      python -m pharos_gap_kg.release --selftest
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import sys
import tarfile
import urllib.request
from pathlib import Path

from pharos_gap_kg import paths

REPO = "iandriver/pharos-gap-kg"


def name(release: str) -> str:
    return f"pharos-gap-kg-bundle-{release}.tar.gz"


def url(release: str) -> str:
    return f"https://github.com/{REPO}/releases/latest/download/{name(release)}"


def record_path(release: str) -> Path:
    return paths.DOCS_DATA / f"bundle_archive_{release}.json"


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def pack(d: Path, release: str) -> bytes:
    """The bundle directory as a gzipped tar whose bytes depend only on the files' names and contents."""
    raw = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
        with tarfile.open(fileobj=gz, mode="w", format=tarfile.USTAR_FORMAT) as tar:
            for f in sorted(p for p in d.iterdir() if p.is_file()):
                data = f.read_bytes()
                ti = tarfile.TarInfo(f"{release}/{f.name}")
                ti.size, ti.mtime, ti.mode, ti.uid, ti.gid, ti.uname, ti.gname = len(data), 0, 0o644, 0, 0, "", ""
                tar.addfile(ti, io.BytesIO(data))
    return raw.getvalue()


def archive(release: str, check: bool = True) -> dict:
    d = paths.bundle(release)
    if not (d / "manifest.json").is_file():
        raise SystemExit(f"no bundle at {d}; build it first:  python -m pharos_gap_kg.bundle {release}")
    if check:
        from pharos_gap_kg import ask
        ask.Graph(d)                                   # refuses a bundle whose files do not match its manifest
    blob = pack(d, release)
    paths.RELEASE.mkdir(parents=True, exist_ok=True)
    (paths.RELEASE / name(release)).write_bytes(blob)
    rec = {"release": release, "file": name(release), "url": url(release), "bytes": len(blob), "sha256": sha(blob),
           "members": {f.name: sha(f.read_bytes()) for f in sorted(d.iterdir()) if f.is_file()}}
    record_path(release).write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def unpack(release: str, file: Path, check: bool = True) -> dict:
    """Check the archive and every member against the committed record, then write the bundle."""
    if not record_path(release).is_file():
        raise SystemExit(f"REFUSED: no committed record {record_path(release).name}, so there is nothing to check the archive against")
    rec = json.loads(record_path(release).read_text())
    blob = Path(file).read_bytes()
    if sha(blob) != rec["sha256"]:
        raise SystemExit(f"REFUSED: {Path(file).name} is not the archive whose digest is committed in {record_path(release).name}")
    got = {}
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tar:
        for m in tar.getmembers():
            base = m.name[len(release) + 1:] if m.name.startswith(release + "/") else None
            if not m.isreg() or base not in rec["members"] or base in got:
                raise SystemExit(f"REFUSED: the archive holds {m.name!r}, which is not a file the record names")
            data = tar.extractfile(m).read()
            if sha(data) != rec["members"][base]:
                raise SystemExit(f"REFUSED: {m.name} is not the file the record names")
            got[base] = data
    if set(got) != set(rec["members"]):
        raise SystemExit(f"REFUSED: the archive lacks {sorted(set(rec['members']) - set(got))}")
    d = paths.bundle(release)
    clash = [n for n in got if (d / n).is_file() and sha((d / n).read_bytes()) != rec["members"][n]]
    if clash:
        raise SystemExit(f"REFUSED: {d} already holds other versions of {clash}; move that bundle aside first")
    d.mkdir(parents=True, exist_ok=True)
    out = {}
    for n, data in sorted(got.items()):
        if (d / n).is_file():
            out[n] = "already in place"
            continue
        (d / n).write_bytes(data)
        out[n] = "written"
    if check:
        from pharos_gap_kg import ask
        ask.Graph(d)
    return out


def get(release: str, src: str | None = None, check: bool = True) -> dict:
    rec = json.loads(record_path(release).read_text()) if record_path(release).is_file() else None
    if rec is None:
        raise SystemExit(f"REFUSED: no committed record {record_path(release).name}")
    f = paths.RELEASE / name(release)
    if not (f.is_file() and sha(f.read_bytes()) == rec["sha256"]):
        paths.RELEASE.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(src or rec["url"], headers={"User-Agent": "pharos-gap-kg"})
        with urllib.request.urlopen(req, timeout=300) as r:
            data = r.read(rec["bytes"] + 1)
        if len(data) != rec["bytes"]:
            raise SystemExit(f"REFUSED: the download is {len(data):,} bytes and the record says {rec['bytes']:,}")
        f.write_bytes(data)
    return unpack(release, f, check)


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

    old = (paths.BUNDLE, paths.RELEASE, paths.DOCS_DATA)
    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        paths.BUNDLE, paths.RELEASE, paths.DOCS_DATA = t / "bundle", t / "release", t / "docs"
        try:
            paths.DOCS_DATA.mkdir()
            b = paths.bundle("r")
            b.mkdir(parents=True)
            for i, n in enumerate(("manifest.json", "links.parquet", "nodes.parquet")):
                (b / n).write_bytes(bytes([i + 1]) * (100 + i))
            rec = archive("r", check=False)
            first = (paths.RELEASE / name("r")).read_bytes()
            ck("an archive's bytes depend only on the files' names and contents", pack(b, "r") == first and rec["sha256"] == sha(first) and len(rec["members"]) == 3)
            keep = {n: (b / n).read_bytes() for n in rec["members"]}
            for n in keep:
                (b / n).unlink()
            ck("unpacking writes every member back, byte for byte", set(unpack("r", paths.RELEASE / name("r"), check=False).values()) == {"written"}
               and all((b / n).read_bytes() == v for n, v in keep.items()))
            ck("unpacking again leaves what is in place", set(unpack("r", paths.RELEASE / name("r"), check=False).values()) == {"already in place"})
            bad = bytearray(first)
            bad[-20] ^= 1
            (t / "bad.tar.gz").write_bytes(bytes(bad))
            ck("an archive whose digest is not the committed one is refused", refused(lambda: unpack("r", t / "bad.tar.gz", check=False)))
            (b / "links.parquet").write_bytes(b"another version")
            ck("a bundle on disk that differs is not overwritten", refused(lambda: unpack("r", paths.RELEASE / name("r"), check=False))
               and (b / "links.parquet").read_bytes() == b"another version")
            (b / "links.parquet").write_bytes(keep["links.parquet"])
            raw = io.BytesIO()
            with tarfile.open(fileobj=raw, mode="w:gz") as tar:
                ti = tarfile.TarInfo("r/../../escape")
                ti.size = 3
                tar.addfile(ti, io.BytesIO(b"bad"))
            (t / "evil.tar.gz").write_bytes(raw.getvalue())
            json.dump({**rec, "sha256": sha(raw.getvalue())}, open(record_path("r"), "w"))
            ck("a member that is not a file the record names is refused, and nothing escapes", refused(lambda: unpack("r", t / "evil.tar.gz", check=False))
               and not (t / "escape").exists())
            json.dump(rec, open(record_path("r"), "w"))
            (paths.RELEASE / name("r")).unlink()
            for n in keep:
                (b / n).unlink()
            (t / "served").mkdir()
            (t / "served" / name("r")).write_bytes(first)
            ck("get downloads the archive, checks it and unpacks it", set(get("r", src=(t / "served" / name("r")).as_uri(), check=False).values()) == {"written"}
               and (paths.RELEASE / name("r")).read_bytes() == first)
            ck("an archive already downloaded and correct is not asked for again",
               set(get("r", src=(t / "nowhere" / name("r")).as_uri(), check=False).values()) == {"already in place"})
        finally:
            paths.BUNDLE, paths.RELEASE, paths.DOCS_DATA = old
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) >= 2 and argv[0] == "archive":
        r = archive(argv[1])
        print(f"  {paths.RELEASE / r['file']}: {r['bytes']:,} bytes, sha256 {r['sha256']}; recorded in {record_path(argv[1]).name}")
        print(f"  attach it to the release; `get` downloads it from {r['url']}")
        return 0
    if len(argv) >= 2 and argv[0] in ("get", "unpack"):
        src = argv[argv.index("--from") + 1] if "--from" in argv else None
        out = get(argv[1], src) if argv[0] == "get" else unpack(argv[1], Path(argv[2]))
        print(f"  {sum(v == 'written' for v in out.values())} file(s) written, {sum(v != 'written' for v in out.values())} already in place -> {paths.bundle(argv[1])}")
        print(f"  ask:  python -m pharos_gap_kg.ask gaps \"Parkinson disease\"")
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
