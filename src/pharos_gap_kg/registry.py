"""External data this package reads: where it is published, how big it is, and how to fetch it.

A missing file yields an instruction, not a traceback: `require(key)` names the fetch command.

Sizes and dates are the ones the publisher's directory listing showed on 2026-10-06. `fetch` refuses a
download whose size differs, and a digest that differs from a pinned one: every "this link appeared
later" label is a difference between two of these files, so a snapshot replaced in place would move
the labels without any code changing.

Run:  python -m pharos_gap_kg.registry          # the table, with what is on disk
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from pharos_gap_kg import paths

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
        return f"python -m pharos_gap_kg.fetch {self.key}"


_S = [
    Source("tcrd6124", f"{BASE}/archive/pharos_tcrd6124m.sql.gz", 8_251_845_338, "2022-07-07",
           "TCRD 6.12.4, the full database behind Pharos 3.12 (MySQL dump)"),
    Source("pharos316", f"{BASE}/archive/pharos316.sql.gz", 2_195_610_389, "2023-04-12",
           "Pharos 3.16 database (MySQL dump). Data version 6.13.4, the same as 3.19: its link tables are identical",
           sha256="ed42237f2e4d0f4f85a31d9435e422afc623cca92a859d5bbf1ff08ba80e0eda"),
    Source("pharos317", f"{BASE}/archive/pharos317.sql.gz", 2_198_542_435, "2023-08-11",
           "Pharos 3.17 database (MySQL dump)"),
    Source("pharos318", f"{BASE}/archive/pharos318.sql.gz", 2_206_028_018, "2023-10-24",
           "Pharos 3.18 database (MySQL dump)"),
    Source("pharos319", f"{BASE}/archive/pharos319.sql.gz", 2_257_847_063, "2024-04-29",
           "Pharos 3.19 database (MySQL dump): the last full legacy database, link evidence to early 2022",
           sha256="86ebfffdd5fd03bcc23ef43a31251d40572ccd7feec62dadd3215f276f8b0953"),
    Source("pharos400_canonicals", f"{BASE}/pharos400_tdls_canonicals.csv", 44_614_999, "2026-05-28",
           "Pharos 4.0.0 development levels for 20,654 canonical proteins, with the evidence counts behind each",
           sha256="ee8fdc79893235a451ab7a8a4b6099da390ae9ccd341179297d8f8730f801459"),
    Source("pharos400_full", f"{BASE}/pharos400_tdls_full.csv", 87_201_549, "2026-05-28",
           "Pharos 4.0.0 development levels for all 166,833 human UniProt rows (isoforms and alternate products)"),
    Source("pharos400_readme", f"{BASE}/pharos400_readme.html", 22_906, "2026-06-04",
           "Pharos 4.0.0 release notes: level definitions, activity cutoffs, source versions",
           sha256="40c03514cc43e8709cf513b1ba1f26308644287fc31ab15c1ee087b44f586039"),
]
# Reference files from outside Pharos. Neither is a model input: they define controls on the test.
_S += [
    Source("efo_obo", "https://github.com/EBISPOT/efo/releases/download/v3.94.0/efo.obo", 88_303_602, "2026-09-15",
           "Experimental Factor Ontology v3.94.0 (OBO): which GWAS traits are ancestors or descendants of which",
           sha256="1f7bdc48b91c17cd98742d8c0fca77d576d4d74998efd643e51b1b49602ca49c"),
    Source("ensembl_gtf_113", "https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz",
           64_141_785, "2024-08-18", "Ensembl release 113 gene annotation, GRCh38: where each gene sits on its chromosome",
           sha256="62f1709b40e083ce9d4cdc64a86b5ffec2c5d5371434bb7095c74dc89079c466"),
]
# The GWAS Catalog, which Pharos's GWAS source aggregates, as curated after Pharos 3.19 was built: the
# links of the prospective test (docs/PROTOCOL.md). A dated release directory, so the file does not move.
_S += [
    Source("gwas_catalog_20261005", "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/10/05/gwas-catalog-associations_ontology-annotated-full.zip",
           72_826_184, "2026-10-05", "GWAS Catalog release of 5 Oct 2026: every curated association with its mapped genes, ontology terms, study and date",
           sha256="c3b18bcde213728abe7c0a93e58f86328584c9b1dd599b7acaae456215c9cee1"),
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
