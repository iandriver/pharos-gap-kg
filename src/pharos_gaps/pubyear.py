"""The publication year of each paper a link cites, from NCBI.

A database dump has one date, the day it was written. The evidence inside it has many: a compound
activity cites the paper that measured it, a GWAS association cites its study. The year of the first
paper behind a link is when the link became knowable, which is what a split by date needs.

Years are asked of NCBI E-utilities (`esummary`, 200 ids a request, under three requests a second) and
kept in `data/tables/pubmed_year.parquet`, so an id is asked once. An id NCBI does not return is kept
with no year; it is not asked again and not guessed. NCBI asks callers to name a contact address:
pass `--email` or set `NCBI_EMAIL`. None is stored in this repository.

Run:  python -m pharos_gaps.pubyear pharos319
      python -m pharos_gaps.pubyear --selftest
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd
import pyarrow.parquet as pq

from pharos_gaps import paths

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CACHE = paths.TABLES / "pubmed_year.parquet"
YEAR = re.compile(r"^\s*((?:18|19|20)\d\d)\b")
SPLIT = re.compile(r"[|,;\s]+")

# table -> (column, whether the column holds several ids as text)
CITING = {"ncats_ligand_activity": ("pubmed_ids", True), "tiga_provenance": ("pubmedid", False), "gwas": ("pmid", False)}


def ids_in_text(s: pd.Series) -> pd.Series:
    """Every PubMed id in a column of delimited text, one per row, keeping the original index."""
    x = s.dropna().astype(str).str.strip()
    x = x[x != ""].str.split(SPLIT).explode()
    return pd.to_numeric(x, errors="coerce").dropna().astype("int64")


def cited(release: str) -> set[int]:
    out: set[int] = set()
    for table, (col, multi) in CITING.items():
        p = paths.tables(release) / f"{table}.parquet"
        if not p.is_file():
            raise SystemExit(f"{release} has no extracted {table}; run:  python -m pharos_gaps.sqldump extract {release} --tables {table}")
        s = pq.read_table(p, columns=[col]).to_pandas()[col]
        out |= set(ids_in_text(s)) if multi else set(s.dropna().astype("int64"))
    return {i for i in out if i > 0}


def parse_year(pubdate: str | None) -> int | None:
    m = YEAR.match(pubdate or "")
    return int(m.group(1)) if m else None


def parse(response: dict, asked: list[int]) -> list[dict]:
    """One record per id asked; an id NCBI did not return, or returned as an error, has no year."""
    res = response.get("result", {})
    rows = []
    for i in asked:
        rec = res.get(str(i), {})
        date = None if "error" in rec else rec.get("pubdate") or rec.get("epubdate") or None
        rows.append({"pmid": i, "year": parse_year(date), "pubdate": date})
    return rows


def lookup(ids: list[int], email: str | None, pause: float = 0.4) -> pd.DataFrame:
    rows: list[dict] = []
    for k in range(0, len(ids), 200):
        batch = ids[k:k + 200]
        q = {"db": "pubmed", "id": ",".join(map(str, batch)), "retmode": "json", "tool": "pharos-gaps"}
        if email:
            q["email"] = email
        body = urllib.parse.urlencode(q).encode()
        for attempt in range(6):
            try:
                with urllib.request.urlopen(urllib.request.Request(ESUMMARY, data=body), timeout=90) as r:
                    rows += parse(json.load(r), batch)
                break
            except Exception as e:
                last = e
                time.sleep(2 ** attempt)
        else:
            raise SystemExit(f"NCBI did not answer for ids {batch[0]}..{batch[-1]} after 6 attempts ({last}); "
                             f"what was resolved before this is not saved, re-run to retry")
        if (k // 200) % 20 == 0:
            print(f"  {min(k + 200, len(ids)):,} / {len(ids):,}", flush=True)
        time.sleep(pause)
    return pd.DataFrame(rows, columns=["pmid", "year", "pubdate"]).astype({"pmid": "int64", "year": "Int64"})


def load() -> pd.DataFrame:
    if not CACHE.is_file():
        return pd.DataFrame({"pmid": pd.Series(dtype="int64"), "year": pd.Series(dtype="Int64"), "pubdate": pd.Series(dtype="object")})
    return pd.read_parquet(CACHE)


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("a year leads the date in every form NCBI writes", [parse_year(s) for s in ("2011 Mar 3", "2019", "2020 Spring", "1998 Dec-1999 Jan")] == [2011, 2019, 2020, 1998])
    ck("no date, an empty date and a non-date give no year", [parse_year(s) for s in (None, "", "Mar 2011", "3011")] == [None] * 4)
    resp = {"result": {"uids": ["1", "2"], "1": {"pubdate": "2001 Jan"}, "2": {"pubdate": "", "epubdate": "2015 Jun 1"}, "4": {"error": "cannot get document summary"}}}
    rows = parse(resp, [1, 2, 3, 4])
    ck("each id asked gets one record", [r["pmid"] for r in rows] == [1, 2, 3, 4])
    ck("an electronic date stands in for a missing print date", rows[1]["year"] == 2015)
    ck("an id NCBI omits or errors on has no year", rows[2]["year"] is None and rows[3]["year"] is None)
    s = pd.Series(["123|456", " 789 ", None, "", "12a, 34"], index=[10, 11, 12, 13, 14])
    x = ids_in_text(s)
    ck("delimited ids are split, blanks and non-numbers dropped, and the row index kept", x.tolist() == [123, 456, 789, 34] and x.index.tolist() == [10, 10, 11, 14])
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("release")
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    a = ap.parse_args(argv)
    want = cited(a.release)
    have = load()
    todo = sorted(want - set(have.pmid))
    print(f"{a.release}: {len(want):,} distinct PubMed ids cited; {len(want) - len(todo):,} already resolved, {len(todo):,} to ask")
    if todo:
        got = lookup(todo, a.email)
        have = pd.concat([have, got], ignore_index=True).drop_duplicates("pmid").sort_values("pmid")
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        have.to_parquet(CACHE, index=False)
    mine = have[have.pmid.isin(want)]
    print(f"  with a year {int(mine.year.notna().sum()):,}; without {int(mine.year.isna().sum()):,}; "
          f"years {int(mine.year.min())}-{int(mine.year.max())} -> {CACHE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
