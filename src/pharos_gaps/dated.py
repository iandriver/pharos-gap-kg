"""Links with the year they became knowable, and how many fall after a cutoff.

Two kinds of link in a Pharos database cite the paper behind them: a compound activity against a
target, and a GWAS trait association of a gene. The first year among the papers a link cites is the
year the link became knowable. A split at year C trains on links first reported up to C and scores
on links first reported after it. This module measures how many links each cutoff leaves to score.

WHAT COUNTS AS A TEST LINK. Its year must be known for every row that asserts it. A pair with one
dated row in 2021 and one undated row may have been known in 1995, so it is not scored. A test link
is `scorable` when both of its ends already carry a link of the same kind in training: a model that
has never seen a node can only guess, and that case is counted apart, by which end is new.

Run:  python -m pharos_gaps.dated pharos319
      python -m pharos_gaps.dated --selftest
"""
from __future__ import annotations

import json
import sys

import pandas as pd
import pyarrow.parquet as pq

from pharos_gaps import paths, pubyear

CUTOFFS = [2012, 2014, 2016, 2017, 2018, 2019, 2020]


def _years(release: str) -> pd.Series:
    y = pubyear.load()
    if y.empty:
        raise SystemExit(f"no publication years yet; run:  python -m pharos_gaps.pubyear {release}")
    return y.set_index("pmid").year


def first_year(rows: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """One row per (a, b): first year over its rows, how many rows assert it, and whether every row is dated.

    `rows` has one row per assertion, with `year` the earliest year that assertion cites (NA if none).
    """
    g = rows.groupby([a, b], sort=False).year
    out = pd.DataFrame({"year": g.min(), "n_rows": g.size(), "n_dated": g.count()}).reset_index()
    out["dated"] = out.n_dated == out.n_rows
    return out


def ligand_links(release: str, years: pd.Series) -> pd.DataFrame:
    la = pq.read_table(paths.tables(release) / "ncats_ligand_activity.parquet",
                       columns=["id", "ncats_ligand_id", "target_id", "pubmed_ids"]).to_pandas()
    ids = pubyear.ids_in_text(la.pubmed_ids)
    per_row = ids.map(years).groupby(level=0).min()          # index = row of la
    la["year"] = per_row.reindex(la.index).astype("Int64")
    return first_year(la, "target_id", "ncats_ligand_id")


def trait_links(release: str, years: pd.Series) -> pd.DataFrame:
    d = paths.tables(release)
    t = pq.read_table(d / "tiga.parquet", columns=["protein_id", "ensg", "efoid"]).to_pandas()
    p = pq.read_table(d / "tiga_provenance.parquet", columns=["ensg", "efoid", "pubmedid"]).to_pandas()
    p["year"] = p.pubmedid.map(years).astype("Int64")
    m = t.merge(p[["ensg", "efoid", "year"]], on=["ensg", "efoid"], how="left")
    return first_year(m, "protein_id", "efoid")


def split_counts(links: pd.DataFrame, a: str, b: str, cutoffs: list[int]) -> list[dict]:
    out = []
    known = links[links.year.notna()]
    for c in cutoffs:
        train = known[known.year <= c]
        test = known[(known.year > c) & known.dated]
        seen_a, seen_b = set(train[a]), set(train[b])
        in_a, in_b = test[a].isin(seen_a), test[b].isin(seen_b)
        out.append({"cutoff": c, "train": len(train), "test": len(test),
                    "scorable": int((in_a & in_b).sum()), f"new_{a}": int((~in_a).sum()), f"new_{b}": int((~in_b).sum()),
                    f"{a}_with_scorable": int(test[a][in_a & in_b].nunique()),
                    f"{a}_first_seen_after": int(test[a][~in_a].nunique())})
    return out


def summary(links: pd.DataFrame, a: str, b: str) -> dict:
    y = links.year.dropna().astype(int)
    return {"pairs": len(links), "with_a_year": int(links.year.notna().sum()),
            "every_row_dated": int(links.dated.sum()), "no_year": int(links.year.isna().sum()),
            a: int(links[a].nunique()), b: int(links[b].nunique()),
            "first_year": int(y.min()), "last_year": int(y.max()),
            "by_year_from_2010": {int(k): int(v) for k, v in y[y >= 2010].value_counts().sort_index().items()},
            "splits": split_counts(links, a, b, CUTOFFS)}


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    NA = pd.NA
    rows = pd.DataFrame({"t": [1, 1, 1, 2, 2, 3, 3, 4, 1], "l": [10, 10, 11, 10, 12, 11, 13, 10, 12],
                         "year": pd.array([2010, 2015, 2012, 2019, 2019, 2020, NA, NA, 2021], dtype="Int64")})
    # (1,10):2010 two rows  (1,11):2012  (2,10):2019  (2,12):2019  (3,11):2020  (3,13):undated  (4,10):undated  (1,12):2021
    L = first_year(rows, "t", "l")
    one = L.set_index(["t", "l"])
    ck("a pair's year is its earliest row", one.loc[(1, 10), "year"] == 2010 and one.loc[(1, 10), "n_rows"] == 2)
    ck("a pair with no dated row has no year and is not dated", pd.isna(one.loc[(3, 13), "year"]) and not one.loc[(3, 13), "dated"])
    half = first_year(pd.DataFrame({"t": [5, 5], "l": [9, 9], "year": pd.array([2021, NA], dtype="Int64")}), "t", "l")
    ck("a pair with one undated row keeps its year but is not dated", half.year[0] == 2021 and not half.dated[0])
    s = split_counts(L, "t", "l", [2015])[0]
    ck("training is every pair first reported up to the cutoff", s["train"] == 2)
    ck("testing is every later pair whose rows are all dated", s["test"] == 4)
    # test pairs: (2,10) new t; (2,12) new t and new l; (3,11) new t; (1,12) new l  -> none scorable
    ck("a pair with an unseen end is not scorable, and is counted by which end is new", s["scorable"] == 0 and s["new_t"] == 3 and s["new_l"] == 2)
    ck("three test pairs with a new t are two distinct new t", s["t_first_seen_after"] == 2)
    s2 = split_counts(L, "t", "l", [2019])[0]
    # train (1,10),(1,11),(2,10),(2,12); test (3,11): new t; (1,12): both seen
    ck("a pair between two seen ends is scorable", s2["test"] == 2 and s2["scorable"] == 1 and s2["t_with_scorable"] == 1)
    ck("a half-dated pair is never a test link", split_counts(half, "t", "l", [2015])[0]["test"] == 0)
    sm = summary(L, "t", "l")
    ck("the summary counts dated and undated pairs", (sm["pairs"], sm["with_a_year"], sm["no_year"]) == (8, 6, 2))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 1:
        raise SystemExit(__doc__)
    release = argv[0]
    years = _years(release)
    res = {"release": release,
           "sha256": json.loads((paths.tables(release) / "_manifest.json").read_text())["sha256"],
           "target_ligand": summary(ligand_links(release, years), "target_id", "ncats_ligand_id"),
           "gene_trait": summary(trait_links(release, years), "protein_id", "efoid")}
    for kind, a, b in (("target_ligand", "target_id", "ncats_ligand_id"), ("gene_trait", "protein_id", "efoid")):
        r = res[kind]
        print(f"{kind}: {r['pairs']:,} pairs ({r[a]:,} x {r[b]:,}); every row dated {r['every_row_dated']:,}, "
              f"no year {r['no_year']:,}; years {r['first_year']}-{r['last_year']}")
        print("   since 2010: " + "  ".join(f"{y}:{n:,}" for y, n in r["by_year_from_2010"].items()))
        print(f"   {'cutoff':>6s} {'train':>9s} {'test':>9s} {'scorable':>9s} {'new ' + a:>14s} {'new ' + b:>20s} {'scored ' + a:>18s} {'first seen after':>17s}")
        for s in r["splits"]:
            print(f"   {s['cutoff']:>6d} {s['train']:>9,} {s['test']:>9,} {s['scorable']:>9,} {s['new_' + a]:>14,} {s['new_' + b]:>20,} "
                  f"{s[a + '_with_scorable']:>18,} {s[a + '_first_seen_after']:>17,}")
    out = paths.DOCS_DATA / f"dated_links_{release}.json"
    out.write_text(json.dumps(res, indent=1) + "\n")
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
