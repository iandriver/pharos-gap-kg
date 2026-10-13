"""How large each GWAS study was, and the smallest effect a study of that size could have found.

An absent GWAS link means no study in the source reported one. How much that says depends on how
large the studies were. The GWAS Catalog gives each study's initial sample as a line of text, such as
"1,383 Danish ancestry cases, 62,656 Danish ancestry controls". This reads that line into counts and
turns the counts into a rough floor: the smallest effect a variant carried by 20% of chromosomes
needed to be found at p = 5e-8 with 80% power. The definitions are in docs/PROTOCOL.md, "Gaps for a
disease".

    cases and controls:   odds ratio  exp( sqrt( K x N / (F x cases x controls) ) ),  N = cases + controls
    individuals only:     sqrt( K / (F x N) )  standard deviations per allele

K = (5.4513 + 0.8416)^2 = 39.60, the two normal quantiles; F = 2 x 0.2 x 0.8 = 0.32.

The floor is for one variant frequency under an additive model. It ignores covariates, staged designs
and rarer variants, which need larger effects. It is here so "no link" is read as "no effect of about
this size was found", never as a measured absence.

Run:  python -I -m pharos_gap_kg.studysize build pharos319     # reads the registered archive in memory
      python -m pharos_gap_kg.studysize --selftest
"""
from __future__ import annotations

import json
import math
import re
import sys

import pandas as pd

K = (5.4513 + 0.8416) ** 2
F = 2 * 0.2 * 0.8
NUMBER = re.compile(r"\d{1,3}(?:,\d{3})+|\d+")
KEY = "gwas_catalog_20261005"


def parse(text) -> tuple[int, int, int] | None:
    """(cases, controls, individuals) from the Catalog's line, or None where it holds no number. Each
    number is read with the words up to the next number: cases or controls by whichever of the two
    words comes first, individuals where neither does. What stands in brackets is a remark and is not
    read, nor is a PubMed id: "36,989 cases, controls (see Ripke et al 2014 PMID:25056061)" is 36,989
    cases and a number of controls the line does not give."""
    s = re.sub(r"PMID:?\s*\d+", " ", re.sub(r"\([^)]*\)", " ", str(text or "")), flags=re.I)
    ms = list(NUMBER.finditer(s))
    if not ms:
        return None
    cases = controls = people = 0
    for i, m in enumerate(ms):
        words = s[m.end(): ms[i + 1].start() if i + 1 < len(ms) else len(s)].casefold()
        v = int(m.group().replace(",", ""))
        kind = re.search(r"\b(case|control)s?\b", words)
        if kind is None:
            people += v
        elif kind.group(1) == "control":
            controls += v
        else:
            cases += v
    return cases, controls, people


def effective(cases: int, controls: int, people: int) -> float:
    """The size two studies are compared by: 4 x cases x controls / N where both are given, else everyone counted."""
    if cases > 0 and controls > 0:
        return 4.0 * cases * controls / (cases + controls)
    return float(cases + controls + people)


def smallest_effect(cases: int, controls: int, people: int) -> dict | None:
    """The floor for a study of this size, or None where the size is not known."""
    if cases > 0 and controls > 0:
        n = cases + controls
        return {"design": "cases and controls", "cases": int(cases), "controls": int(controls),
                "smallest_odds_ratio_it_could_find": round(math.exp(math.sqrt(K * n / (F * cases * controls))), 3)}
    n = cases + controls + people
    if n <= 0:
        return None
    return {"design": "individuals", "individuals": int(n), "smallest_effect_it_could_find_in_standard_deviations_per_allele": round(math.sqrt(K / (F * n)), 4)}


def say(e: dict) -> str:
    """One sentence for a floor."""
    if e["design"] == "cases and controls":
        return (f"it had {e['cases']:,} cases and {e['controls']:,} controls, so a variant carried by 20% of chromosomes needed an odds ratio of about "
                f"{e['smallest_odds_ratio_it_could_find']:.2f} to be found")
    return (f"it had {e['individuals']:,} individuals, so a variant carried by 20% of chromosomes needed an effect of about "
            f"{e['smallest_effect_it_could_find_in_standard_deviations_per_allele']:.3f} standard deviations per allele to be found")


def table(studies: pd.DataFrame) -> pd.DataFrame:
    """studies: study, text. Adds cases, controls, individuals; a study whose line holds no number keeps none."""
    got = [parse(t) for t in studies.text]
    out = studies.copy()
    for j, c in enumerate(("cases", "controls", "individuals")):
        out[c] = pd.array([g[j] if g else None for g in got], dtype="Int64")
    return out


def build(release: str) -> dict:
    import zipfile

    import pyarrow.parquet as pq

    from pharos_gap_kg import gwascat, paths, registry
    src = registry.require(KEY)
    side = json.loads(src.with_name(src.name + ".lineage.json").read_text())
    if registry.SOURCES[KEY].sha256 and side["sha256"] != registry.SOURCES[KEY].sha256:
        raise SystemExit(f"REFUSED: {src.name} is not the file the registry pins")
    z = zipfile.ZipFile(src)
    names = [i.filename for i in z.infolist()]
    if len(names) != 1 or not names[0].endswith(".tsv"):
        raise SystemExit(f"REFUSED: the archive holds {names}, and one .tsv was expected")
    with z.open(names[0]) as f:
        df = pd.read_csv(f, sep="\t", dtype=str, quoting=3, usecols=["STUDY ACCESSION", "INITIAL SAMPLE SIZE"], keep_default_na=False)
    df = df.rename(columns={"STUDY ACCESSION": "study", "INITIAL SAMPLE SIZE": "text"})
    lines = df.groupby("study").text.nunique()
    t = table(df.drop_duplicates("study").sort_values("study").reset_index(drop=True))
    out = gwascat.table_dir()
    out.mkdir(parents=True, exist_ok=True)
    t.to_parquet(out / "study_sizes.parquet", index=False)
    mine = set(pq.read_table(paths.tables(release) / "tiga_provenance.parquet", columns=["study_acc"]).to_pandas().study_acc)
    m = t[t.study.isin(mine)]
    read = m.cases.notna()
    cc = read & (m.cases > 0) & (m.controls > 0)
    man = {"release": release, "catalog": KEY, "catalog_sha256": side["sha256"], "catalog_studies": len(t), "catalog_studies_read": int(t.cases.notna().sum()),
           "catalog_studies_with_more_than_one_line": int((lines > 1).sum()),
           "release_studies": len(mine), "release_studies_in_the_catalog": len(m), "release_studies_read": int(read.sum()),
           "share_of_release_studies_read": round(float(read.sum()) / len(mine), 4), "of_them_cases_and_controls": int(cc.sum()), "of_them_individuals_only": int((read & ~cc).sum()),
           "median_size": {"cases": int(m[cc].cases.median()), "controls": int(m[cc].controls.median()), "individuals": int((m[read & ~cc].cases + m[read & ~cc].controls + m[read & ~cc].individuals).median())},
           "constants": {"K": round(K, 4), "F": F, "variant_frequency": 0.2, "p": 5e-8, "power": 0.8},
           "file_sha256": __import__("pharos_gap_kg.graph", fromlist=["file_sha"]).file_sha(out / "study_sizes.parquet")}
    (paths.DOCS_DATA / f"study_sizes_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("cases and controls are read from one line, with thousands separators", parse("1,383 Danish ancestry cases, 62,656 Danish ancestry controls") == (1383, 62656, 0))
    ck("a line of individuals, with or without a separator", parse("47,745 European ancestry individuals") == (0, 0, 47745) and parse("5477 European ancestry individuals") == (0, 0, 5477))
    ck("several ancestries of one kind are added", parse("123 Han Taiwanese ancestry cases, 2,000 Japanese ancestry cases, 214,995 Han Taiwanese ancestry controls") == (2123, 214995, 0))
    ck("children and other words count as individuals", parse("7 Asian ancestry children, 352 European ancestry children") == (0, 0, 359))
    ck("a remark in brackets and a PubMed id are not people, and the first of the two words decides",
       parse("2,366 European ancestry bipolar disorder cases, 36,989 schizophrenia cases, controls (see Ripke et al 2014 PMID:25056061)") == (39355, 0, 0)
       and parse("704 responder cases, 1,662 non-responder bipolar disorder controls, 36,989 schizophrenia cases, controls (see Ripke 2014 PMID:25056061).") == (37693, 1662, 0)
       and parse("500 cases PMID: 123456") == (500, 0, 0))
    ck("a line with no number is not read", parse("NR") is None and parse("") is None and parse(None) is None)
    ck("K and F are the protocol's constants", abs(K - 39.60) < 0.01 and abs(F - 0.32) < 1e-12)
    e = smallest_effect(10_000, 10_000, 0)
    ck("ten thousand cases and as many controls: an odds ratio of 1.170", e["smallest_odds_ratio_it_could_find"] == round(math.exp(math.sqrt(39.6006 * 20000 / (0.32 * 1e8))), 3) == 1.17)
    q = smallest_effect(0, 0, 100_000)
    ck("a hundred thousand individuals: 0.0352 standard deviations per allele", q["smallest_effect_it_could_find_in_standard_deviations_per_allele"] == 0.0352 and q["design"] == "individuals")
    ck("a larger study could find a smaller effect", smallest_effect(100_000, 100_000, 0)["smallest_odds_ratio_it_could_find"] < e["smallest_odds_ratio_it_could_find"]
       and smallest_effect(0, 0, 10 ** 6)["smallest_effect_it_could_find_in_standard_deviations_per_allele"] < 0.0352)
    ck("few cases against many controls is a small study", effective(123, 214_995, 0) < effective(10_000, 10_000, 0) and round(effective(10_000, 10_000, 0)) == 20_000)
    ck("cases with no controls are counted as individuals", smallest_effect(500, 0, 0)["design"] == "individuals" and smallest_effect(0, 0, 0) is None)
    t = table(pd.DataFrame({"study": ["a", "b"], "text": ["10 cases, 20 controls", "NR"]}))
    ck("a study that is not read keeps no counts", t.cases.tolist()[0] == 10 and t.controls.tolist()[0] == 20 and pd.isna(t.cases[1]) and pd.isna(t.individuals[1]))
    ck("the sentence gives the size and the floor", "10,000 cases and 10,000 controls" in say(e) and "1.17" in say(e) and "0.035 standard deviations" in say(q))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 2 or argv[0] != "build":
        raise SystemExit(__doc__)
    m = build(argv[1])
    print(f"{m['catalog_studies']:,} studies in the Catalog file, {m['catalog_studies_read']:,} with a size read")
    print(f"{m['release']}: {m['release_studies_read']:,} of {m['release_studies']:,} studies read ({m['share_of_release_studies_read']:.1%}); "
          f"{m['of_them_cases_and_controls']:,} with cases and controls, {m['of_them_individuals_only']:,} with individuals only; median size {m['median_size']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
