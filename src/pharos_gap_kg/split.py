"""The frozen split: which target-trait links a model may see, and which it is scored on.

A link is a Pharos target (a UniProt accession) and a GWAS trait. Its year is the earliest publication
year among the studies cited for it. At cutoff C:

    train     first reported in or before C
    test      first reported after C, and both its target and its trait have a training link
    unscored  first reported after C with a target or trait the training links never mention

FROZEN MEANS CHECKED. `--freeze` writes the counts and a SHA-256 of each part to `docs/data/`, once.
Every later run recomputes the parts and refuses if a digest moved, so a changed input, a changed
year lookup or a changed rule cannot quietly move a result's denominator. The digest is over the
sorted "target<TAB>trait" lines, so it does not depend on row order or file format.

Run:  python -m pharos_gap_kg.split pharos319 --cutoff 2018 --freeze     # once
      python -m pharos_gap_kg.split pharos319 --cutoff 2018              # rebuild and check
      python -m pharos_gap_kg.split --selftest
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import paths, pubyear

PARTS = ("train", "test", "unscored")


def curie(efoid: str) -> str:
    """TIGA writes a trait id as PREFIX_LOCAL; a CURIE is PREFIX:LOCAL."""
    head, sep, tail = str(efoid).partition("_")
    if not sep or not head or not tail:
        raise ValueError(f"trait id {efoid!r} is not PREFIX_LOCAL")
    return f"{head}:{tail}"


def gene_trait(release: str, years: pd.Series) -> pd.DataFrame:
    """One row per (target, trait): year first reported, studies and papers behind it, and its name."""
    d = paths.tables(release)
    t = pq.read_table(d / "tiga.parquet", columns=["protein_id", "ensg", "efoid", "trait"]).to_pandas().drop_duplicates()
    p = pq.read_table(d / "tiga_provenance.parquet", columns=["ensg", "efoid", "study_acc", "pubmedid"]).to_pandas().drop_duplicates()
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot"]).to_pandas().rename(columns={"id": "protein_id"})
    p["year"] = p.pubmedid.map(years)
    if p.year.isna().any():
        raise SystemExit(f"REFUSED: {int(p.year.isna().sum()):,} cited studies have no publication year; "
                         f"run:  python -m pharos_gap_kg.pubyear {release}")
    m = t.merge(p, on=["ensg", "efoid"], how="left").merge(prot, on="protein_id", how="left")
    if m.study_acc.isna().any() or m.uniprot.isna().any():
        raise SystemExit("REFUSED: a link has no study behind it, or a target has no accession")
    m["target"] = "UniProtKB:" + m.uniprot
    m["trait_id"] = m.efoid.map(curie)
    name = m.sort_values("trait").drop_duplicates("trait_id").set_index("trait_id").trait
    g = m.groupby(["target", "trait_id"], sort=True)
    out = pd.DataFrame({"year": g.year.min().astype(int), "last_year": g.year.max().astype(int),
                        "n_studies": g.study_acc.nunique(), "n_papers": g.pubmedid.nunique()}).reset_index()
    out = out.rename(columns={"trait_id": "trait"})
    out["trait_name"] = out.trait.map(name)
    return out


def split(links: pd.DataFrame, cutoff: int) -> dict[str, pd.DataFrame]:
    train = links[links.year <= cutoff]
    later = links[links.year > cutoff]
    seen = later.target.isin(set(train.target)) & later.trait.isin(set(train.trait))
    return {"train": train, "test": later[seen], "unscored": later[~seen]}


def digest(part: pd.DataFrame) -> str:
    lines = sorted(part.target + "\t" + part.trait)
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def describe(parts: dict[str, pd.DataFrame], cutoff: int) -> dict:
    tr, te = parts["train"], parts["test"]
    return {"cutoff": cutoff,
            **{k: {"links": len(v), "targets": int(v.target.nunique()), "traits": int(v.trait.nunique()),
                   "sha256": digest(v)} for k, v in parts.items()},
            "candidates": int(tr.target.nunique()) * int(tr.trait.nunique()) - len(tr),
            "test_by_year": {int(y): int(n) for y, n in te.year.value_counts().sort_index().items()}}


def frozen_path(cutoff: int) -> Path:
    return paths.DOCS_DATA / f"split_gene_trait_{cutoff}.json"


def out_dir(cutoff: int) -> Path:
    return paths.EVAL / f"gene_trait_{cutoff}"


def load(cutoff: int) -> dict[str, pd.DataFrame]:
    """The parts as written by the last checked run; refuses if they are not the frozen ones."""
    fz = frozen_path(cutoff)
    if not fz.is_file():
        raise SystemExit(f"no frozen split for {cutoff}; see `python -m pharos_gap_kg.split --help`")
    want = json.loads(fz.read_text())
    parts = {}
    for k in PARTS:
        p = out_dir(cutoff) / f"{k}.parquet"
        if not p.is_file():
            raise SystemExit(f"missing {p.name}; rebuild with:  python -m pharos_gap_kg.split {want['release']} --cutoff {cutoff}")
        parts[k] = pd.read_parquet(p)
        if digest(parts[k]) != want[k]["sha256"]:
            raise SystemExit(f"REFUSED: {p} is not the frozen {k} part (digest differs)")
    return parts


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    L = pd.DataFrame({"target": ["A", "A", "B", "B", "C", "A", "D"], "trait": ["x", "y", "x", "z", "x", "z", "y"],
                      "year": [2010, 2018, 2015, 2019, 2020, 2021, 2022]})
    s = split(L, 2018)
    ck("a link first reported in the cutoff year trains", len(s["train"]) == 3 and (s["train"].year <= 2018).all())
    # (B,z): z unseen -> unscored; (C,x): C unseen -> unscored; (A,z): z unseen -> unscored; (D,y): D unseen
    ck("a later link with an unseen end is not scored", len(s["unscored"]) == 4 and len(s["test"]) == 0)
    L2 = pd.concat([L, pd.DataFrame({"target": ["B"], "trait": ["y"], "year": [2020]})], ignore_index=True)
    s2 = split(L2, 2018)
    ck("the three parts partition the links", sum(len(v) for v in s2.values()) == len(L2) and len(s2["test"]) == 1)
    ck("the digest ignores row order", digest(s2["train"]) == digest(s2["train"].iloc[::-1]))
    ck("the digest changes when one link changes", digest(s2["train"]) != digest(s2["train"].assign(trait=["x", "y", "w"])))
    d = describe(s2, 2018)
    ck("candidates are every unlinked pair of seen ends", d["candidates"] == 2 * 2 - 3 and d["test"]["links"] == 1)
    ck("a trait id becomes a CURIE", curie("EFO_0000537") == "EFO:0000537" and curie("Orphanet_98") == "Orphanet:98")
    try:
        curie("EFO0000537")
        ck("a trait id with no prefix is refused", False)
    except ValueError:
        ck("a trait id with no prefix is refused", True)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("release")
    ap.add_argument("--cutoff", type=int, required=True)
    ap.add_argument("--freeze", action="store_true")
    a = ap.parse_args(argv)
    y = pubyear.load()
    links = gene_trait(a.release, y.set_index("pmid").year)
    parts = split(links, a.cutoff)
    d = {"release": a.release, "sha256": json.loads((paths.tables(a.release) / "_manifest.json").read_text())["sha256"],
         **describe(parts, a.cutoff)}
    fz = frozen_path(a.cutoff)
    if fz.is_file():
        want = json.loads(fz.read_text())
        moved = [k for k in PARTS if want[k]["sha256"] != d[k]["sha256"]]
        if moved:
            raise SystemExit(f"REFUSED: the {moved} part(s) differ from the split frozen in {fz.name} "
                             f"({', '.join(f'{k}: {want[k]['links']:,} then, {d[k]['links']:,} now' for k in moved)}). "
                             f"An input, a year or the rule changed. Nothing was written.")
        if a.freeze:
            print(f"already frozen in {fz.name}; this run matches it")
    elif a.freeze:
        fz.write_text(json.dumps(d, indent=1) + "\n")
        print(f"frozen -> {fz}")
    else:
        raise SystemExit(f"no frozen split for cutoff {a.cutoff}; pass --freeze once to write {fz.name}")
    o = out_dir(a.cutoff)
    o.mkdir(parents=True, exist_ok=True)
    for k, v in parts.items():
        v.to_parquet(o / f"{k}.parquet", index=False)
        print(f"  {k:9s} {len(v):>9,} links  {v.target.nunique():>6,} targets  {v.trait.nunique():>5,} traits  {d[k]['sha256'][:12]}")
    print(f"  candidates (unlinked pairs of seen ends): {d['candidates']:,}; test links by year: {d['test_by_year']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
