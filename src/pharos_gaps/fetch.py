"""Fetch a registered source: ranged, resumable, and refused if it is not the file the registry describes.

The publisher serves about 0.9 MB/s per connection, so a 2.2 GB snapshot is fetched as several byte
ranges at once and joined. Each range resumes from what is already on disk. The joined file's SHA-256,
size, and the server's Last-Modified are written beside it (`<file>.lineage.json`); everything built
from the file cites that digest.

Run:  python -m pharos_gaps.fetch pharos319
      python -m pharos_gaps.fetch pharos400_canonicals --parts 1
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import threading
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from pharos_gaps.registry import SOURCES, Source

CHUNK = 1 << 20
UA = {"User-Agent": "pharos-gaps (research; github.com/iandriver)"}


def head(url: str) -> tuple[int, str]:
    req = urllib.request.Request(url, method="HEAD", headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"]), r.headers.get("Last-Modified", "")


def _range(url: str, start: int, end: int, part: Path, errors: list) -> None:
    """Bytes start..end inclusive into `part`, resuming from the part's current size."""
    want = end - start + 1
    for attempt in range(12):
        have = part.stat().st_size if part.exists() else 0
        if have >= want:
            return
        try:
            req = urllib.request.Request(url, headers={**UA, "Range": f"bytes={start + have}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r, part.open("ab") as f:
                if r.status != 206:
                    raise OSError(f"server answered {r.status} to a range request")
                while True:
                    buf = r.read(CHUNK)
                    if not buf:
                        break
                    f.write(buf)
        except Exception as e:  # network errors resume from the part's size
            time.sleep(min(60, 2 ** attempt))
            last = e
    else:
        errors.append(f"{part.name}: gave up after 12 attempts ({last})")


def fetch(s: Source, parts: int = 4) -> Path:
    dest = s.path
    side = dest.with_name(dest.name + ".lineage.json")
    if dest.is_file() and dest.stat().st_size == s.bytes and side.is_file():
        print(f"{s.filename}: already fetched ({json.loads(side.read_text())['sha256'][:12]})")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    length, modified = head(s.url)
    if length != s.bytes:
        raise SystemExit(f"REFUSED {s.filename}: the server offers {length:,} bytes, the registry describes "
                         f"{s.bytes:,}. The published file changed; update the registry deliberately.")
    parts = max(1, min(parts, length // (8 * CHUNK) or 1))
    step = -(-length // parts)
    spans = [(i * step, min(length, (i + 1) * step) - 1) for i in range(parts)]
    files = [dest.with_name(f"{dest.name}.part{i}") for i in range(parts)]
    errors: list[str] = []
    threads = [threading.Thread(target=_range, args=(s.url, a, b, f, errors), daemon=True)
               for (a, b), f in zip(spans, files)]
    t0 = time.time()
    for t in threads:
        t.start()
    while any(t.is_alive() for t in threads):
        time.sleep(30)
        got = sum(f.stat().st_size for f in files if f.exists())
        print(f"  {s.filename}: {got / 1e9:5.2f} / {length / 1e9:.2f} GB  ({time.time() - t0:5.0f} s)", flush=True)
    if errors:
        raise SystemExit("INCOMPLETE, re-run to resume:\n  " + "\n  ".join(errors))
    h = hashlib.sha256()
    tmp = dest.with_name(dest.name + ".joining")
    with tmp.open("wb") as out:
        for (a, b), f in zip(spans, files):
            if f.stat().st_size != b - a + 1:
                raise SystemExit(f"INCOMPLETE {f.name}: {f.stat().st_size:,} of {b - a + 1:,} bytes; re-run to resume")
            with f.open("rb") as src:
                while buf := src.read(CHUNK):
                    h.update(buf)
                    out.write(buf)
    digest = h.hexdigest()
    if s.sha256 and digest != s.sha256:
        tmp.rename(dest.with_name(dest.name + ".rejected"))
        raise SystemExit(f"REFUSED {s.filename}: sha256 {digest[:16]} is not the pinned {s.sha256[:16]}")
    tmp.rename(dest)
    for f in files:
        f.unlink()
    side.write_text(json.dumps({
        "key": s.key, "url": s.url, "bytes": length, "sha256": digest, "server_last_modified": modified,
        "published": s.published, "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }, indent=1) + "\n")
    print(f"{s.filename}: {length:,} bytes, sha256 {digest}")
    return dest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("keys", nargs="+", choices=sorted(SOURCES))
    ap.add_argument("--parts", type=int, default=4)
    a = ap.parse_args(argv)
    for k in a.keys:
        fetch(SOURCES[k], a.parts)
    return 0


if __name__ == "__main__":
    sys.exit(main())
