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

A GWAS by proxy counts people with an affected relative as cases. Their numbers are counted apart, and
the floor of such a study is given at two ends: proxies left out, which is what an absent link
supports, and proxies counted as cases and controls; a proxy carries less information than a case, so
the study's reach lies between (docs/PROTOCOL.md, "Proxy cases in the study-size table").

The floor is for one variant frequency under an additive model. It ignores covariates, staged designs
and rarer variants, which need larger effects. It is here so "no link" is read as "no effect of about
this size was found", never as a measured absence.

Run:  python -I -m pharos_gap_kg.studysize build pharos319     # reads the registered archive in memory
      python -I -m pharos_gap_kg.studysize compare pharos319   # what counting proxies apart changed, once the bundle has the sizes
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
ODDS_RATIOS, SDS = (1.05, 1.1, 1.2), (0.02, 0.05, 0.1)        # the effects a proposed study is sized for
RARE = {"combined_frequency": 0.001, "odds_ratio": 3.0, "p": 2.5e-6}
K_RARE = (4.7081 + 0.8416) ** 2
F_RARE = 2 * 0.001 * 0.999
NUMBER = re.compile(r"\d{1,3}(?:,\d{3})+|\d+")
PROXY = re.compile(r"\bproxy\b|family history|parental history")
PROXY_CONTROL = re.compile(r"\bcontrols?\b|\bno (?:family|parental) history")
KEY = "gwas_catalog_20261005"
COUNTS = ("cases", "controls", "individuals", "proxy_cases", "proxy_controls")


def read(text) -> tuple[int, int, int, int, int] | None:
    """(cases, controls, individuals, proxy cases, proxy controls) from the Catalog's line, or None where
    it holds no number. Each number is read with the words up to the next number. Words that hold
    "proxy", "family history" or "parental history" make it proxies: proxy controls if they also hold
    "control" or "no family history" or "no parental history", proxy cases otherwise. Any other number
    is cases or controls by whichever of the two words comes first, individuals where neither does.
    What stands in brackets is a remark and is not read, nor is a PubMed id: "36,989 cases, controls
    (see Ripke et al 2014 PMID:25056061)" is 36,989 cases and a number of controls the line does not give."""
    s = re.sub(r"PMID:?\s*\d+", " ", re.sub(r"\([^)]*\)", " ", str(text or "")), flags=re.I)
    ms = list(NUMBER.finditer(s))
    if not ms:
        return None
    cases = controls = people = pcases = pcontrols = 0
    for i, m in enumerate(ms):
        words = s[m.end(): ms[i + 1].start() if i + 1 < len(ms) else len(s)].casefold()
        v = int(m.group().replace(",", ""))
        kind = re.search(r"\b(case|control)s?\b", words)
        if PROXY.search(words):
            if PROXY_CONTROL.search(words):
                pcontrols += v
            else:
                pcases += v
        elif kind is None:
            people += v
        elif kind.group(1) == "control":
            controls += v
        else:
            cases += v
    return cases, controls, people, pcases, pcontrols


def parse(text) -> tuple[int, int, int] | None:
    """(cases, controls, individuals) from the Catalog's line, proxies left out (`read`)."""
    r = read(text)
    return None if r is None else r[:3]


def effective(cases: int, controls: int, people: int) -> float:
    """The size two studies are compared by: 4 x cases x controls / N where both are given, else everyone counted."""
    if cases > 0 and controls > 0:
        return 4.0 * cases * controls / (cases + controls)
    return float(cases + controls + people)


def odds_floor(cases: int, controls: int) -> float:
    return round(math.exp(math.sqrt(K * (cases + controls) / (F * cases * controls))), 3)


def smallest_effect(cases: int, controls: int, people: int, proxy_cases: int = 0, proxy_controls: int = 0) -> dict | None:
    """The floor for a study of this size, or None where the size is not known. With proxies, the floor
    is the end with them left out, and the end with them counted is added where that gives both cases
    and controls; a study of proxies alone has only that end."""
    out = None
    if cases > 0 and controls > 0:
        out = {"design": "cases and controls", "cases": int(cases), "controls": int(controls), "smallest_odds_ratio_it_could_find": odds_floor(cases, controls)}
    elif cases + controls + people > 0:
        n = cases + controls + people
        out = {"design": "individuals", "individuals": int(n), "smallest_effect_it_could_find_in_standard_deviations_per_allele": round(math.sqrt(K / (F * n)), 4)}
    if proxy_cases + proxy_controls <= 0:
        return out
    c2, k2 = cases + proxy_cases, controls + proxy_controls
    pr = {"proxy_cases": int(proxy_cases), "proxy_controls": int(proxy_controls)}
    if c2 > 0 and k2 > 0:
        pr["smallest_odds_ratio_with_proxies_as_cases"] = odds_floor(c2, k2)
    if out is None:
        return {"design": "proxies only", **pr} if "smallest_odds_ratio_with_proxies_as_cases" in pr else None
    return {**out, **pr}


def _proxies(e: dict) -> str:
    n = f"{e['proxy_cases']:,} proxy cases" + (f" and {e['proxy_controls']:,} proxy controls" if e["proxy_controls"] else "")
    return n + " (people counted by an affected relative, or by none)" if e["proxy_controls"] else n + " (people counted by an affected relative)"


def say(e: dict) -> str:
    """One sentence for a floor."""
    if e["design"] == "proxies only":
        return (f"it counted proxies only, {_proxies(e)}; counted as cases and controls, a variant carried by 20% of chromosomes needed an odds ratio of about "
                f"{e['smallest_odds_ratio_with_proxies_as_cases']:.2f} to be found, which overstates its reach, since a proxy carries less information than a case")
    if e["design"] == "cases and controls":
        s = (f"it had {e['cases']:,} cases and {e['controls']:,} controls, so a variant carried by 20% of chromosomes needed an odds ratio of about "
             f"{e['smallest_odds_ratio_it_could_find']:.2f} to be found")
    else:
        s = (f"it had {e['individuals']:,} individuals, so a variant carried by 20% of chromosomes needed an effect of about "
             f"{e['smallest_effect_it_could_find_in_standard_deviations_per_allele']:.3f} standard deviations per allele to be found")
    if "smallest_odds_ratio_with_proxies_as_cases" in e and e["design"] == "cases and controls":
        s += (f". It also had {_proxies(e)}; counted as cases, about {e['smallest_odds_ratio_with_proxies_as_cases']:.2f}. A proxy carries less information than a case, "
              f"so the smallest effect it could find lies between {e['smallest_odds_ratio_with_proxies_as_cases']:.2f} and {e['smallest_odds_ratio_it_could_find']:.2f}")
    elif "proxy_cases" in e:
        s += f". It also had {_proxies(e)}, not counted here"
    return s


def needed(cases: int, controls: int, people: int, proxy_cases: int = 0) -> dict | None:
    """The study that would close an absent link, sized from the largest study so far: for a common
    variant, the cases (with controls at that study's ratio) or the individuals needed for each effect
    in ODDS_RATIOS or SDS, and whether that study already had them, and had them once its proxy cases
    are counted; for rare variants, one gene-based test (RARE). The floor's formula solved for the
    cases; docs/PROTOCOL.md, "The study each gap calls for". Sizes are from cases and controls without proxies."""
    if cases > 0 and controls > 0:
        r = controls / cases
        size = lambda k, f, o: math.ceil(k * (1 + 1 / r) / (f * math.log(o) ** 2))
        common = [{"odds_ratio": o, "cases": size(K, F, o), "controls": math.ceil(r * size(K, F, o)), "the_largest_study_had_them": cases >= size(K, F, o)}
                  for o in ODDS_RATIOS]
        if proxy_cases > 0:
            for c in common:
                c["the_largest_study_had_them_with_its_proxy_cases_counted"] = cases + proxy_cases >= c["cases"]
        n = size(K_RARE, F_RARE, RARE["odds_ratio"])
        return {"design": "cases and controls", "controls_per_case": round(r, 2), "largest_study_cases": int(cases), "common_variant": common,
                **({"largest_study_proxy_cases": int(proxy_cases)} if proxy_cases > 0 else {}),
                "rare_variants": {**RARE, "cases": n, "controls": math.ceil(r * n), "data": "sequence, which a genotyped GWAS is not"}}
    n0 = cases + controls + people
    if n0 <= 0:
        return None
    common = [{"sd_per_allele": s, "individuals": math.ceil(K / (F * s * s)), "the_largest_study_had_them": n0 >= math.ceil(K / (F * s * s))} for s in SDS]
    return {"design": "individuals", "largest_study_individuals": int(n0), "common_variant": common}


def two(n: int) -> str:
    """A size to two significant figures, for prose."""
    return f"{float(f'{n:.2g}'):,.0f}"


def _and(xs: list[str]) -> str:
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]


def say_needed(e: dict) -> str:
    """One or two sentences for the study an absent link calls for."""
    cv = e["common_variant"]
    unit, eff = ("cases", "odds_ratio") if e["design"] == "cases and controls" else ("individuals", "sd_per_allele")
    sizes = [f"{two(c[unit])} for {c[eff]:g}" for c in cv]
    sizes[0] = f"about {two(cv[0][unit])} {unit} for {'an odds ratio of' if unit == 'cases' else 'an effect of'} {cv[0][eff]:g}" + ("" if unit == "cases" else " standard deviations per allele")
    enough = lambda had: ("enough for none of them" if not had else "enough for all three" if len(had) == len(cv) else f"enough for {_and(had)}")
    had_say = enough([f"{c[eff]:g}" for c in cv if c["the_largest_study_had_them"]])
    if unit == "cases":
        rv = e["rare_variants"]
        px = ""
        if e.get("largest_study_proxy_cases"):
            px = f", and with its {e['largest_study_proxy_cases']:,} proxy cases counted as cases, " + enough([f"{c[eff]:g}" for c in cv if c["the_largest_study_had_them_with_its_proxy_cases_counted"]])
        return (f"The study an absent link calls for: for a common variant, {_and(sizes)}, with the largest study's {e['controls_per_case']:.1f} controls per case; the largest "
                f"study had {had_say}{px}. For rare variants (1 in {round(1 / rv['combined_frequency']):,} chromosomes combined, odds ratio {rv['odds_ratio']:g}), a gene-based "
                f"test on sequence data would need about {two(rv['cases'])} cases.")
    return f"The study an absent link calls for: for a common variant, {_and(sizes)}; the largest study had {had_say}."


def table(studies: pd.DataFrame) -> pd.DataFrame:
    """studies: study, text. Adds cases, controls, individuals, proxy_cases, proxy_controls; a study whose line holds no number keeps none."""
    got = [read(t) for t in studies.text]
    out = studies.copy()
    for j, c in enumerate(COUNTS):
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
    got = m.cases.notna()
    cc = got & (m.cases > 0) & (m.controls > 0)
    ind = got & ~cc & (m.cases + m.controls + m.individuals > 0)
    px = lambda d: d[(d.proxy_cases > 0) | (d.proxy_controls > 0)]
    man = {"release": release, "catalog": KEY, "catalog_sha256": side["sha256"], "catalog_studies": len(t), "catalog_studies_read": int(t.cases.notna().sum()),
           "catalog_studies_with_more_than_one_line": int((lines > 1).sum()),
           "release_studies": len(mine), "release_studies_in_the_catalog": len(m), "release_studies_read": int(got.sum()),
           "share_of_release_studies_read": round(float(got.sum()) / len(mine), 4), "of_them_cases_and_controls": int(cc.sum()), "of_them_individuals_only": int(ind.sum()),
           "of_them_proxies_only": int((got & ~cc & ~ind).sum()),
           "median_size": {"cases": int(m[cc].cases.median()), "controls": int(m[cc].controls.median()), "individuals": int((m[ind].cases + m[ind].controls + m[ind].individuals).median())},
           "proxies": {"catalog_studies_with_proxies": len(px(t)),
                       "release_studies_with_proxies": {r.study: {c: int(getattr(r, c)) for c in COUNTS} for r in px(m).sort_values("study").itertuples()}},
           "constants": {"K": round(K, 4), "F": F, "variant_frequency": 0.2, "p": 5e-8, "power": 0.8},
           "file_sha256": __import__("pharos_gap_kg.graph", fromlist=["file_sha"]).file_sha(out / "study_sizes.parquet")}
    (paths.DOCS_DATA / f"study_sizes_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
    return man


def read_before_29_oct(text) -> tuple[int, int, int] | None:
    """The reading before proxies were counted apart: a proxy number whose words hold "case" was cases,
    one whose words hold neither word was individuals. Kept to say what the rule of 29 Oct changed."""
    s = re.sub(r"PMID:?\s*\d+", " ", re.sub(r"\([^)]*\)", " ", str(text or "")), flags=re.I)
    ms = list(NUMBER.finditer(s))
    if not ms:
        return None
    out = [0, 0, 0]
    for i, m in enumerate(ms):
        words = s[m.end(): ms[i + 1].start() if i + 1 < len(ms) else len(s)].casefold()
        kind = re.search(r"\b(case|control)s?\b", words)
        out[2 if kind is None else 1 if kind.group(1) == "control" else 0] += int(m.group().replace(",", ""))
    return tuple(out)


def compare(release: str) -> dict:
    """Every trait's largest study and floor under the reading before 29 Oct and under the rule of 29 Oct,
    as the served answers give them; writes docs/data/proxy_floors_<release>.json. Refuses if a study whose
    line holds no proxy words is read differently by the two."""
    from pharos_gap_kg import ask, diseasegaps, gwascat, paths
    full = pd.read_parquet(gwascat.table_dir() / "study_sizes.parquet")
    plain = full[~full.text.str.casefold().str.contains(PROXY.pattern, na=False)]
    differ = [s for s, t, *n in zip(plain.study, plain.text, plain.cases, plain.controls, plain.individuals)
              if read_before_29_oct(t) != (None if pd.isna(n[0]) else tuple(int(v) for v in n))]
    if differ:
        raise SystemExit(f"REFUSED: {len(differ)} studies with no proxy words are read differently, e.g. {differ[:5]}")
    G = ask.Graph(paths.bundle(release))
    if not G.proxies:
        raise SystemExit(f"REFUSED: the bundle counts no proxies; run:  python -I -m pharos_gap_kg.studysize build {release}  then  python -m pharos_gap_kg.bundle {release} --sizes-only")
    text = dict(zip(full.study, full.text))
    before = {s: read_before_29_oct(text[s]) for s in G.sizes}
    G._evidence(0, 0)
    ev = G._ev[G._ev.year <= G.latest].drop_duplicates(["trait", "study"])
    by = {k: list(zip(g.study, g.year)) for k, g in ev.groupby("trait")}
    keep = ("study", "design", "cases", "controls", "individuals", "smallest_odds_ratio_it_could_find", "smallest_effect_it_could_find_in_standard_deviations_per_allele", "smallest_odds_ratio_with_proxies_as_cases",
            "proxy_cases", "proxy_controls")
    rows = []
    for i in range(len(G.d)):
        tid = G.d.id[i]
        sized = [(s, int(y), before[s]) for s, y in by.get(tid, []) if s in G.sizes]
        if not sized:
            continue
        s0, _, n0 = max(sized, key=lambda r: (effective(*r[2]), r[0]))
        b = {"study": s0, **(smallest_effect(*n0) or {})}
        a = (diseasegaps.study_floor(G, tid, G.latest) or {}).get("largest") or {}
        if b == {"study": s0, **{k: v for k, v in a.items() if k not in ("year",)}} and "proxy_cases" not in a:
            continue
        kind = ("the largest study changed" if a.get("study") != s0 else
                "proxies read as cases before: the floor rose" if a.get("smallest_odds_ratio_it_could_find", 0) > b.get("smallest_odds_ratio_it_could_find", 9) else
                "proxies read as individuals before: a second end added" if a.get("smallest_odds_ratio_it_could_find") == b.get("smallest_odds_ratio_it_could_find") is not None else
                "the design changed")
        rows.append({"trait": tid, "name": G.d["name"][i], "disease_like": bool(G.d.disease_like[i]), "kind": kind,
                     "before": {k: b[k] for k in keep if k in b}, "after": {k: a[k] for k in keep if k in a}})
        if a.get("study") != s0:                                   # the earlier largest study, read by the rule
            rows[-1]["earlier_largest_study_under_the_rule"] = {"study": s0, **(smallest_effect(*G.sizes[s0], *G.proxies.get(s0, (0, 0))) or {})}
    kinds = sorted({r["kind"] for r in rows})
    px = full[(full.proxy_cases > 0) | (full.proxy_controls > 0)].sort_values("study")
    man = {"release": release, "rule": "docs/PROTOCOL.md, Proxy cases in the study-size table (29 Oct 2026)",
           "catalog_studies_with_no_proxy_words_read_as_before": len(plain), "release_studies_with_proxies": len(G.proxies),
           "catalog_studies_with_proxies": {r.study: {"in_the_release": r.study in G.sizes, "text": r.text, **{c: int(getattr(r, c)) for c in COUNTS}} for r in px.itertuples()},
           "traits_with_a_sized_study": sum(1 for i in range(len(G.d)) if any(s in G.sizes for s, _ in by.get(G.d.id[i], []))),
           "traits_changed": len(rows), "by_kind": {k: sum(r["kind"] == k for r in rows) for k in kinds}, "traits": rows}
    (paths.DOCS_DATA / f"proxy_floors_{release}.json").write_text(json.dumps(man, indent=1) + "\n")
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
    ad = "24,087 European ancestry late-onset Alzheimer's disease cases, 47,793 European ancestry individuals with family history of Alzheimer's disease, 383,378 European ancestry controls"
    ck("proxies are counted apart, whether the line calls them proxy cases or individuals with a family history",
       read("39,106 European ancestry clinically diagnosed cases, 46,828 European ancestry proxy cases, 401,577 European ancestry controls") == (39106, 401577, 0, 46828, 0)
       and read(ad) == (24087, 383378, 0, 47793, 0) and parse(ad) == (24087, 383378, 0)
       and read("71,880 European ancestry Alzheimer's disease or AD-by-proxy cases, 16,666 European ancestry peptic ulcer disease cases") == (16666, 0, 0, 71880, 0))
    ck("no parental history makes proxy controls, and a remark in brackets is still not read",
       read("up to 42,034 British ancestry individuals with parental history of Alzheimer's disease, at least 272,244 British ancestry individuals with no parental "
            "history of Alzheimer's disease, 25,580 Alzheimer's disease cases, 48,466 controls") == (25580, 48466, 0, 42034, 272244)
       and read("111,263 European ancestry individuals (MTAG effective sample size boosetd by family history of depression samples)") == (0, 0, 111263, 0, 0))
    e7 = smallest_effect(*read(ad))
    ck("GCST007320's two ends are 1.077 with its proxies left out and 1.046 with them counted as cases",
       (e7["smallest_odds_ratio_it_could_find"], e7["smallest_odds_ratio_with_proxies_as_cases"], e7["proxy_cases"]) == (1.077, 1.046, 47793)
       and "between 1.05 and 1.08" in say(e7) and "47,793 proxy cases" in say(e7))
    po = smallest_effect(0, 0, 0, 42034, 272244)
    ck("a study of proxies alone has only the end with them counted, and says it overstates", po["design"] == "proxies only" and "smallest_odds_ratio_it_could_find" not in po
       and po["smallest_odds_ratio_with_proxies_as_cases"] == odds_floor(42034, 272244) and "overstates" in say(po) and smallest_effect(0, 0, 0, 500, 0) is None)
    n7 = needed(24087, 383378, 0, 47793)
    ck("the cases a study needs come from its cases without proxies, and whether it had them is given both ways",
       n7["common_variant"] == [{**c, "the_largest_study_had_them_with_its_proxy_cases_counted": 24087 + 47793 >= c["cases"]} for c in needed(24087, 383378, 0)["common_variant"]]
       and "with its 47,793 proxy cases counted as cases, enough for" in say_needed(n7) and "proxy" not in say_needed(needed(24087, 383378, 0)))
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
    nd = needed(10_000, 10_000, 0)
    back = [smallest_effect(c["cases"], c["controls"], 0)["smallest_odds_ratio_it_could_find"] for c in nd["common_variant"]]
    ck("the cases needed, solved back through the floor, give the odds ratio they were sized for", all(abs(b - o) <= 0.0015 for b, o in zip(back, ODDS_RATIOS)))
    ck("the cases needed at 1:1 for an odds ratio of 1.1 are 2K/(F (ln 1.1)^2)", nd["common_variant"][1]["cases"] == math.ceil(2 * K / (F * math.log(1.1) ** 2)) and nd["controls_per_case"] == 1)
    ck("a study whose floor is 1.17 had the cases for 1.2 and not for 1.1 or 1.05", [c["the_largest_study_had_them"] for c in nd["common_variant"]] == [False, False, True])
    rv = nd["rare_variants"]
    ck("the rare-variant test uses its own threshold and frequency, and is not marked as within a GWAS's reach",
       rv["cases"] == math.ceil(K_RARE * 2 / (F_RARE * math.log(3) ** 2)) and abs(K_RARE - 30.80) < 0.01 and abs(F_RARE - 0.001998) < 1e-12 and "the_largest_study_had_them" not in rv)
    many = needed(20_000, 400_000, 0)
    ck("more controls per case need fewer cases, and the controls follow the ratio", many["common_variant"][0]["cases"] < nd["common_variant"][0]["cases"] and many["common_variant"][0]["controls"] >= 20 * many["common_variant"][0]["cases"])
    ni = needed(0, 0, 100_000)
    ck("for individuals, K/(F sd^2) each, and whether the study had them", [c["individuals"] for c in ni["common_variant"]] == [math.ceil(K / (F * s * s)) for s in SDS]
       and [c["the_largest_study_had_them"] for c in ni["common_variant"]] == [False, True, True] and needed(0, 0, 0) is None)
    ck("sizes are said to two significant figures, with what the largest study had enough for",
       two(54_321) == "54,000" and two(3_712) == "3,700" and two(95) == "95"
       and say_needed(nd) == ("The study an absent link calls for: for a common variant, about 100,000 cases for an odds ratio of 1.05, 27,000 for 1.1 and 7,400 for 1.2, with "
                              "the largest study's 1.0 controls per case; the largest study had enough for 1.2. For rare variants (1 in 1,000 chromosomes combined, odds "
                              "ratio 3), a gene-based test on sequence data would need about 26,000 cases.")
       and say_needed(ni) == ("The study an absent link calls for: for a common variant, about 310,000 individuals for an effect of 0.02 standard deviations per allele, "
                              "50,000 for 0.05 and 12,000 for 0.1; the largest study had enough for 0.05 and 0.1."))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 2 or argv[0] not in ("build", "compare"):
        raise SystemExit(__doc__)
    if argv[0] == "compare":
        c = compare(argv[1])
        print(f"{c['release']}: {c['traits_changed']} of {c['traits_with_a_sized_study']:,} traits' floors change; {c['by_kind']}")
        for r in c["traits"]:
            print(f"  {r['trait']:14s} {r['name'][:48]:48s} {r['kind']}: {r['before']} -> {r['after']}")
        return 0
    m = build(argv[1])
    print(f"{m['catalog_studies']:,} studies in the Catalog file, {m['catalog_studies_read']:,} with a size read")
    print(f"{m['release']}: {m['release_studies_read']:,} of {m['release_studies']:,} studies read ({m['share_of_release_studies_read']:.1%}); "
          f"{m['of_them_cases_and_controls']:,} with cases and controls, {m['of_them_individuals_only']:,} with individuals only; median size {m['median_size']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
