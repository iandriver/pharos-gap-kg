"""External data this package reads: where it is published, how big it is, and how to fetch it.

A missing file yields an instruction, not a traceback: `require(key)` names the fetch command.

Sizes and dates are the ones the publisher's directory listing showed on 2026-10-06. `fetch` refuses a
download whose size differs, and a digest that differs from a pinned one: every "this link appeared
later" label is a difference between two of these files, so a snapshot replaced in place would move
the labels without any code changing.

Run:  python -m pharos_gaps.registry          # the table, with what is on disk
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from pharos_gaps import paths

BASE = "https://opendata.ncats.nih.gov/public/pharos"


@dataclass(frozen=True)
class Source:
    key: str
    url: str
    bytes: int
    published: str          # date in the publisher's listing, ISO
    what: str
    sha256: str | None = None   # pinned after the first verified fetch

    @property
    def filename(self) -> str:
        return self.url.rsplit("/", 1)[1]

    @property
    def path(self) -> Path:
        return paths.RAW / self.filename

    @property
    def fetch_command(self) -> str:
        return f"python -m pharos_gaps.fetch {self.key}"


_S = [
    Source("tcrd6124", f"{BASE}/archive/pharos_tcrd6124m.sql.gz", 8_251_845_338, "2022-07-07",
           "TCRD 6.12.4, the full database behind Pharos 3.12 (MySQL dump)"),
    Source("pharos316", f"{BASE}/archive/pharos316.sql.gz", 2_195_610_389, "2023-04-12",
           "Pharos 3.16 database (MySQL dump): the snapshot the graph is BUILT from"),
    Source("pharos317", f"{BASE}/archive/pharos317.sql.gz", 2_198_542_435, "2023-08-11",
           "Pharos 3.17 database (MySQL dump)"),
    Source("pharos318", f"{BASE}/archive/pharos318.sql.gz", 2_206_028_018, "2023-10-24",
           "Pharos 3.18 database (MySQL dump)"),
    Source("pharos319", f"{BASE}/archive/pharos319.sql.gz", 2_257_847_063, "2024-04-29",
           "Pharos 3.19 database (MySQL dump): the snapshot links are SCORED against"),
    Source("pharos400_canonicals", f"{BASE}/pharos400_tdls_canonicals.csv", 44_614_999, "2026-05-28",
           "Pharos 4.0.0 development levels for 20,654 canonical proteins, with the evidence counts behind each"),
    Source("pharos400_full", f"{BASE}/pharos400_tdls_full.csv", 87_201_549, "2026-05-28",
           "Pharos 4.0.0 development levels for all 166,833 human UniProt rows (isoforms and alternate products)"),
    Source("pharos400_readme", f"{BASE}/pharos400_readme.html", 22_906, "2026-06-04",
           "Pharos 4.0.0 release notes: level definitions, activity cutoffs, source versions"),
]
SOURCES: dict[str, Source] = {s.key: s for s in _S}


def require(key: str) -> Path:
    """The path of a fetched source, or a SystemExit that says how to fetch it."""
    s = SOURCES[key]
    if not s.path.is_file():
        raise SystemExit(f"missing input: {s.filename} ({s.bytes / 1e9:.2f} GB, published {s.published})\n"
                         f"  {s.what}\n  fetch it with:  {s.fetch_command}")
    if s.path.stat().st_size != s.bytes:
        raise SystemExit(f"{s.filename} is {s.path.stat().st_size:,} bytes on disk, {s.bytes:,} in the registry: "
                         f"an interrupted or replaced download. Re-run:  {s.fetch_command}")
    return s.path


if __name__ == "__main__":
    for s in _S:
        have = "on disk" if s.path.is_file() and s.path.stat().st_size == s.bytes else "-"
        print(f"{s.key:22s} {s.bytes / 1e9:6.2f} GB  {s.published}  {have:8s} {s.what}")
