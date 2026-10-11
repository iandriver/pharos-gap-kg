"""The GWAS Catalog as curated after Pharos 3.19 was built: the links of the prospective test
(docs/PROTOCOL.md, "A prospective test: links first reported from 2023").

Pharos's GWAS links come from a source that aggregated the Catalog up to 2022. The Catalog has gone on
being curated. This turns its associations into (target, trait) links by a rule chosen where it can be
checked, on the studies Pharos already holds, and then splits them: every pair of Pharos 3.19 is
training, a pair Pharos does not hold whose first study was published in 2023 or later is a
prospective positive, and a pair Pharos does not hold whose first study is earlier is set aside.

Nothing here runs a scorer. `build` reports counts.

    python -I -m pharos_gap_kg.gwascat extract            # the columns the rule reads, from the registered zip
    python -m pharos_gap_kg.gwascat rule pharos319        # choose the rule on the studies Pharos holds
    python -m pharos_gap_kg.gwascat build pharos319 --freeze
    python -m pharos_gap_kg.gwascat --selftest

The zip is a download, so `extract` reads it in memory and is run with `-I`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import paths, registry, split

KEY = "gwas_catalog_20261005"
PROSPECTIVE = 2022                 # the split's cutoff: Pharos 3.19 holds links to 2022
FROM = 2023                        # a prospective positive was first reported in this year or later
COLUMNS = {"PUBMEDID": "pubmedid", "DATE": "date", "STUDY ACCESSION": "study_acc", "PVALUE_MLOG": "p_mlog", "SNP_GENE_IDS": "genes_in",
           "UPSTREAM_GENE_ID": "gene_up", "DOWNSTREAM_GENE_ID": "gene_down", "MAPPED_TRAIT_URI": "trait_uris", "OR or BETA": "or_beta"}
# the four declared rules, in the order that breaks a tie: which genes, and which p-values
RULES = (("in", "5e-8"), ("in", "any"), ("near", "5e-8"), ("near", "any"))
MLOG_5E8 = -np.log10(5e-8)
PASS = 0.90


def table_dir() -> Path:
    return paths.TABLES / KEY


def rule_path(release: str) -> Path:
    return paths.DOCS_DATA / f"gwascat_rule_{release}.json"


def extract() -> dict:
    """The columns the rule reads, written as a table. The archive is read in memory and never unpacked."""
    import zipfile
    src = registry.require(KEY)
    side = json.loads(src.with_name(src.name + ".lineage.json").read_text())
    if registry.SOURCES[KEY].sha256 and side["sha256"] != registry.SOURCES[KEY].sha256:
        raise SystemExit(f"REFUSED: {src.name} is not the file the registry pins")
    z = zipfile.ZipFile(src)
    names = [i.filename for i in z.infolist()]
    if len(names) != 1 or not names[0].endswith(".tsv"):
        raise SystemExit(f"REFUSED: the archive holds {names}, and one .tsv was expected")
    with z.open(names[0]) as f:
        df = pd.read_csv(f, sep="\t", dtype=str, quoting=3, usecols=list(COLUMNS), keep_default_na=False)
    df = df.rename(columns=COLUMNS)
    df["p_mlog"] = pd.to_numeric(df.p_mlog, errors="coerce")
    df["year"] = pd.to_numeric(df.date.str[:4], errors="coerce").astype("Int64")
    out = table_dir()
    out.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out / "associations.parquet", index=False)
    man = {"source": KEY, "sha256": side["sha256"], "member": names[0], "rows": int(len(df)), "studies": int(df.study_acc.nunique()), "papers": int(df.pubmedid.nunique()),
           "without_a_year": int(df.year.isna().sum()), "years": [int(df.year.min()), int(df.year.max())]}
    (out / "_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def trait_ids(uris: pd.Series) -> pd.Series:
    """Each comma-separated ontology address as the last part of its path, which is how Pharos writes a trait."""
    return uris.str.split(",").map(lambda xs: [x.strip().rsplit("/", 1)[-1] for x in xs if x.strip()])


def gene_ids(df: pd.DataFrame, which: str) -> pd.Series:
    """The Ensembl genes of each association: those the variant lies in, and for `near` also the nearest on each side."""
    inside = df.genes_in.str.split(",").map(lambda xs: [x.strip() for x in xs if x.strip()])
    if which == "in":
        return inside
    side = lambda s: s.str.strip().map(lambda x: [x] if x else [])
    return inside + side(df.gene_up) + side(df.gene_down)


def triples(df: pd.DataFrame, rule: tuple[str, str]) -> pd.DataFrame:
    """(ensg, efoid, study_acc, year), one row for each gene and trait of each association the rule keeps."""
    which, p = rule
    d = df if p == "any" else df[df.p_mlog >= MLOG_5E8 - 1e-9]
    t = pd.DataFrame({"ensg": gene_ids(d, which), "efoid": trait_ids(d.trait_uris), "study_acc": d.study_acc, "year": d.year})
    t = t.explode("ensg").explode("efoid").dropna(subset=["ensg", "efoid"])
    return t[(t.ensg != "") & (t.efoid != "")].drop_duplicates(["ensg", "efoid", "study_acc"])


def agreement(mine: pd.DataFrame, theirs: pd.DataFrame) -> dict:
    """How well one set of (ensg, efoid, study_acc) triples reproduces another."""
    key = lambda t: set(zip(t.ensg, t.efoid, t.study_acc))
    a, b = key(mine), key(theirs)
    tp = len(a & b)
    prec, rec = (tp / len(a) if a else 0.0), (tp / len(b) if b else 0.0)
    return {"produced": len(a), "pharos": len(b), "shared": tp, "precision": prec, "recall": rec, "f1": (2 * prec * rec / (prec + rec) if prec + rec else 0.0)}


def gene_to_targets(release: str) -> pd.DataFrame:
    """ensg, target: the release's own cross-references, as used for position."""
    import pyarrow.parquet as pq
    d = paths.tables(release)
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot"]).to_pandas().rename(columns={"id": "protein_id"})
    x = pq.read_table(d / "xref.parquet", columns=["xtype", "protein_id", "value"]).to_pandas()
    x = x[x.xtype == "ENSG"][["protein_id", "value"]].rename(columns={"value": "ensg"})
    t = pq.read_table(d / "tiga.parquet", columns=["protein_id", "ensg"]).to_pandas()
    m = pd.concat([x, t]).dropna()
    m["ensg"] = m.ensg.str.split(".").str[0]
    m = m.drop_duplicates().merge(prot, on="protein_id")
    m["target"] = "UniProtKB:" + m.uniprot
    return m[["ensg", "target"]].drop_duplicates()


def choose_rule(release: str) -> dict:
    """Score the four declared rules on the studies Pharos 3.19 also cites, and choose by F1."""
    import pyarrow.parquet as pq
    df = pd.read_parquet(table_dir() / "associations.parquet")
    theirs = pq.read_table(paths.tables(release) / "tiga_provenance.parquet", columns=["ensg", "efoid", "study_acc"]).to_pandas().drop_duplicates()
    shared = set(theirs.study_acc) & set(df.study_acc)
    known = set(gene_to_targets(release).ensg)
    sub, theirs_in = df[df.study_acc.isin(shared)], theirs[theirs.study_acc.isin(shared)]
    rows = {}
    for r in RULES:
        mine = triples(sub, r)
        mine = mine[mine.ensg.isin(known)]                 # a gene with no target in the release can never be a link
        rows["/".join(r)] = agreement(mine, theirs_in)
    best = max(rows, key=lambda k: rows[k]["f1"])           # the first of equals, in the declared order
    out = {"release": release, "catalog": KEY, "catalog_sha256": json.loads((table_dir() / "_manifest.json").read_text())["sha256"],
           "pharos_studies": int(theirs.study_acc.nunique()), "studies_in_both": len(shared), "pharos_triples_in_those": int(len(theirs_in)),
           "rules": rows, "chosen": best, "passes": bool(rows[best]["precision"] >= PASS and rows[best]["recall"] >= PASS), "bar": PASS,
           "compared_as": "(Ensembl gene, trait term, study accession), for genes the release maps to a target"}
    rule_path(release).write_text(json.dumps(out, indent=1) + "\n")
    return out


def classify(catalog_pairs: pd.DataFrame, pharos_pairs: set[tuple[str, str]]) -> pd.Series:
    """For each Catalog pair: 'training' if Pharos holds it, 'prospective' if not and first reported
    from FROM on, 'aside' if not and first reported earlier."""
    held = pd.Series([(a, b) in pharos_pairs for a, b in zip(catalog_pairs.target, catalog_pairs.trait)], index=catalog_pairs.index)
    return pd.Series(np.where(held, "training", np.where(catalog_pairs.year >= FROM, "prospective", "aside")), index=catalog_pairs.index)


def build(release: str, freeze: bool) -> dict:
    """The prospective split, from the chosen rule. Counts only."""
    from pharos_gap_kg import ontology
    rp = rule_path(release)
    if not rp.is_file():
        raise SystemExit(f"no rule chosen; run:  python -m pharos_gap_kg.gwascat rule {release}")
    rule = tuple(json.loads(rp.read_text())["chosen"].split("/"))
    df = pd.read_parquet(table_dir() / "associations.parquet")
    df = df[df.year.notna()]
    t = triples(df, rule).merge(gene_to_targets(release), on="ensg")
    t["trait"] = [split.curie(e) if "_" in e else None for e in t.efoid]
    t = t.dropna(subset=["trait"])
    g = t.groupby(["target", "trait"], sort=True)
    cat = pd.DataFrame({"year": g.year.min().astype(int), "last_year": g.year.max().astype(int), "n_studies": g.study_acc.nunique()}).reset_index()
    old = split.load(2018)
    train = pd.concat([old[k] for k in split.PARTS], ignore_index=True).sort_values(["target", "trait"]).reset_index(drop=True)
    held = set(zip(train.target, train.trait))
    cat["kind"] = classify(cat, held)
    terms = ontology.read(registry.require("efo_obo"))
    label = lambda c: (terms.get(c) or {}).get("name")
    pros = cat[cat.kind == "prospective"].drop(columns="kind").assign(n_papers=0)
    pros["trait_name"] = pros.trait.map(label)
    seen = pros.target.isin(set(train.target)) & pros.trait.isin(set(train.trait))
    parts = {"train": train, "test": pros[seen].reset_index(drop=True), "unscored": pros[~seen].reset_index(drop=True)}
    aside = cat[cat.kind == "aside"].drop(columns="kind").reset_index(drop=True)
    d = {"release": release, "catalog": KEY, "catalog_sha256": json.loads((table_dir() / "_manifest.json").read_text())["sha256"], "rule": "/".join(rule), "from": FROM,
         **split.describe(parts, PROSPECTIVE),
         "aside": {"links": int(len(aside)), "sha256": split.digest(aside), "between_a_training_target_and_trait": int((aside.target.isin(set(train.target)) & aside.trait.isin(set(train.trait))).sum())},
         "catalog_pairs": int(len(cat)), "catalog_pairs_pharos_holds": int((cat.kind == "training").sum()), "pharos_pairs_the_catalog_gives": float((cat.kind == "training").sum() / len(train)),
         "unscored_by_kind": {"new_trait": int((~pros.trait.isin(set(train.trait)) & pros.target.isin(set(train.target))).sum()),
                              "new_target": int((pros.trait.isin(set(train.trait)) & ~pros.target.isin(set(train.target))).sum()),
                              "both_new": int((~pros.trait.isin(set(train.trait)) & ~pros.target.isin(set(train.target))).sum())},
         "test_traits_with_a_positive": int(parts["test"].trait.nunique())}
    fz = split.frozen_path(PROSPECTIVE)
    if fz.is_file():
        want = json.loads(fz.read_text())
        moved = [k for k in (*split.PARTS, "aside") if want[k]["sha256"] != d[k]["sha256"]]
        if moved:
            raise SystemExit(f"REFUSED: the {moved} part(s) differ from the split frozen in {fz.name}. Nothing was written.")
    elif freeze:
        fz.write_text(json.dumps(d, indent=1) + "\n")
    else:
        raise SystemExit(f"no frozen split for {PROSPECTIVE}; pass --freeze once to write {fz.name}")
    o = split.out_dir(PROSPECTIVE)
    o.mkdir(parents=True, exist_ok=True)
    for k, v in parts.items():
        v.to_parquet(o / f"{k}.parquet", index=False)
    aside.to_parquet(o / "aside.parquet", index=False)
    return d


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    df = pd.DataFrame({"genes_in": ["ENSG1, ENSG2", "", "ENSG3"], "gene_up": ["", "ENSG4", ""], "gene_down": ["", "ENSG5", ""],
                       "trait_uris": ["http://www.ebi.ac.uk/efo/EFO_1, http://purl.obolibrary.org/obo/MONDO_2", "http://www.ebi.ac.uk/efo/EFO_1", ""],
                       "study_acc": ["S1", "S2", "S3"], "p_mlog": [9.0, 6.0, 12.0], "year": [2020, 2024, 2024]})
    ck("a trait address becomes its last part, and several traits of one association are each kept", trait_ids(df.trait_uris).tolist() == [["EFO_1", "MONDO_2"], ["EFO_1"], []])
    ck("the genes a variant lies in, and with `near` the nearest on each side when it lies in none",
       gene_ids(df, "in").tolist() == [["ENSG1", "ENSG2"], [], ["ENSG3"]] and gene_ids(df, "near").tolist() == [["ENSG1", "ENSG2"], ["ENSG4", "ENSG5"], ["ENSG3"]])
    key = lambda t: sorted(zip(t.ensg, t.efoid, t.study_acc))
    ck("one triple for each gene and trait, and an association with no trait gives none",
       key(triples(df, ("in", "any"))) == [("ENSG1", "EFO_1", "S1"), ("ENSG1", "MONDO_2", "S1"), ("ENSG2", "EFO_1", "S1"), ("ENSG2", "MONDO_2", "S1")])
    ck("the p-value rule drops an association above 5e-8, and the gene rule decides whether its neighbours count",
       len(triples(df, ("near", "5e-8"))) == 4 and key(triples(df, ("near", "any")))[-2:] == [("ENSG4", "EFO_1", "S2"), ("ENSG5", "EFO_1", "S2")])
    ck("exactly 5e-8 is kept", len(triples(pd.DataFrame({"genes_in": ["ENSG1"], "gene_up": [""], "gene_down": [""], "trait_uris": ["x/EFO_1"], "study_acc": ["S"], "p_mlog": [float(MLOG_5E8)], "year": [2020]}), ("in", "5e-8"))) == 1)
    mine = pd.DataFrame({"ensg": ["a", "b", "c"], "efoid": ["t"] * 3, "study_acc": ["s"] * 3})
    theirs = pd.DataFrame({"ensg": ["a", "b", "d", "e"], "efoid": ["t"] * 4, "study_acc": ["s"] * 4})
    ag = agreement(mine, theirs)
    ck("agreement counts shared triples: precision of what was produced, recall of what Pharos holds", ag["shared"] == 2 and abs(ag["precision"] - 2 / 3) < 1e-12 and abs(ag["recall"] - 0.5) < 1e-12)
    cat = pd.DataFrame({"target": ["g1", "g2", "g3", "g4"], "trait": ["a", "a", "b", "b"], "year": [2019, 2024, 2021, 2023]})
    ck("a pair Pharos holds is training whatever its year; one it lacks is prospective from 2023 and set aside before",
       classify(cat, {("g1", "a"), ("g2", "a")}).tolist() == ["training", "training", "aside", "prospective"])
    ck("four rules are declared, the stricter first", len(RULES) == 4 and RULES[0] == ("in", "5e-8") and len(set(RULES)) == 4)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=("extract", "rule", "build"))
    ap.add_argument("release", nargs="?")
    ap.add_argument("--freeze", action="store_true")
    a = ap.parse_args(argv)
    if a.stage == "extract":
        man = extract()
        print(f"  {man['rows']:,} associations, {man['studies']:,} studies, {man['papers']:,} papers, published {man['years'][0]} to {man['years'][1]}; {man['without_a_year']} without a year")
        return 0
    if not a.release:
        raise SystemExit("rule and build need the release")
    if a.stage == "rule":
        out = choose_rule(a.release)
        print(f"  {out['studies_in_both']:,} of Pharos's {out['pharos_studies']:,} studies are in the Catalog; Pharos holds {out['pharos_triples_in_those']:,} triples for them")
        for k, r in out["rules"].items():
            print(f"  {k:10s} produces {r['produced']:>9,}  shared {r['shared']:>9,}  precision {r['precision']:.3f}  recall {r['recall']:.3f}  F1 {r['f1']:.3f}" + ("   <- chosen" if k == out["chosen"] else ""))
        print(f"  both at {PASS:.0%} or more: {out['passes']}")
    else:
        d = build(a.release, a.freeze)
        print(f"  rule {d['rule']}: {d['catalog_pairs']:,} Catalog pairs; Pharos holds {d['catalog_pairs_pharos_holds']:,} ({d['pharos_pairs_the_catalog_gives']:.1%} of Pharos's {d['train']['links']:,})")
        print(f"  prospective test: {d['test']['links']:,} links over {d['test']['traits']:,} traits and {d['test']['targets']:,} targets, by year {d['test_by_year']}")
        print(f"  cold start: {d['unscored']['links']:,} ({d['unscored_by_kind']}); set aside: {d['aside']['links']:,}, of which {d['aside']['between_a_training_target_and_trait']:,} between a training target and trait")
    return 0


if __name__ == "__main__":
    sys.exit(main())
