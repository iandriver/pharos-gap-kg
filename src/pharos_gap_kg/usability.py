"""Can a small model answer questions about the graph through the typed questions? (docs/PROTOCOL.md)

Three parts, each a subcommand:

    questions   write the 35 questions and their expected answers from the bundle, and freeze them
    run         put the questions to one model under one condition (typed tools, raw SQL, or no tool)
    report      grade every stored answer by code and tabulate

EXPECTED ANSWERS DO NOT COME FROM THE QUESTION LAYER. They are computed here from the bundle's tables
with pandas, so a fault in `ask` cannot make the typed condition agree with a wrong answer. `ask` is
used only to choose names that it resolves without ambiguity, and names it finds ambiguous.

THE MODEL IS TOLD ONLY WHAT THE SERVER SAYS. In the typed condition the system prompt is the server's
own instructions and the tools are its own tool list, both fetched over MCP from a running server,
and every tool call goes through that server.

Run:  python -m pharos_gap_kg.usability questions pharos319 --freeze
      OPENROUTER_API_KEY=... python -m pharos_gap_kg.usability run pharos319 --model qwen/qwen3.7-flash --condition typed
      python -m pharos_gap_kg.usability report
      python -m pharos_gap_kg.usability --selftest
"""
from __future__ import annotations

import argparse
import asyncio
import hashlib
import json
import os
import random
import re
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import ask, paths

SEED = 20261007
# The first set is the development set: its results shaped the second round (docs/PROTOCOL.md).
SETS = {
    "dev": {"seed": 20261007, "schema": 1, "fakes": ("ZQXP1", "VRTK9", "QJMB4"), "words": ("alzheimer", "diabetes", "cholesterol"),
            "file": "usability_questions.json", "results": "usability_results.json"},
    "eval": {"seed": 20261008, "schema": 2, "fakes": ("XQVT3", "ZJKP8", "WQZR5"),
             "words": ("lymphoma", "arthritis", "leukemia", "hepatitis", "sclerosis", "carcinoma", "anemia", "dermatitis"),
             "file": "usability_questions_eval.json", "results": "usability_results_eval.json"},
    # the gaps question (docs/PROTOCOL.md, "Usability of the gaps question"): a development set, then the set the verdicts are read on
    "gaps_dev": {"seed": 20261010, "schema": 3, "fakes": ("QZVK2", "XJWP7", "VQXT4"), "words": ("ZNF", "TMEM", "SLC25"),
                 "file": "usability_questions_gaps_dev.json", "results": "usability_results_gaps_dev.json"},
    "gaps": {"seed": 20261011, "schema": 3, "fakes": ("ZWQK6", "JXVP3", "QKZT9"), "words": ("KRTAP", "CYP2", "HLA"),
             "file": "usability_questions_gaps.json", "results": "usability_results_gaps.json"},
}
STATUSES = ["answered", "no_link_reported", "not_in_graph", "ambiguous", "unknown"]
CONDITIONS = ["typed", "sql", "none"]
MAX_TOOL_CALLS = 8
PROMPT_VERSION = 1
GAPS_CAP = 0.75          # dollars, over both sets of the gaps test (docs/PROTOCOL.md)
PRICES = {"inception/mercury-2.5": (0.04, 0.15), "deepseek/deepseek-v4-flash": (0.09, 0.18), "z-ai/glm-5.3-flash": (0.15, 0.50)}

RULES = """You answer questions about a graph of dated GWAS links between human protein targets and traits, built from NIH Pharos.
- A link exists from the year of the first study behind it. "As of YEAR" means using only what had been reported by the end of that year.
- If the target or the trait is not in the graph, say so: status "not_in_graph". Do not answer from your own knowledge.
- If a name could mean several things in the graph, do not choose one: status "ambiguous".
- "No link" means no study in this source had reported one by that year. It is not evidence of no association: status "no_link_reported".
- If you cannot determine the answer, status "unknown".
Explain briefly, then end with one last line of the form  FINAL: {...}  holding the JSON object the question asks for."""

RULES_V2 = """You answer questions about a graph of dated GWAS links between human protein targets and traits, built from NIH Pharos.
- A link exists from the year of the first study behind it. "As of YEAR" means using only what had been reported by the end of that year.
- "No link" means no study in this source had reported one by that year. It is not evidence of no association.
Explain briefly, then end with one last line of the form  FINAL: {...}  holding the JSON object the question asks for.
Every FINAL object has a "problem" field:
  null             you found what was asked. This includes finding that there is no link.
  "not_in_graph"   the target or the trait is not in the graph. Do not answer from your own knowledge.
  "ambiguous"      a name in the question matches nothing in the graph exactly and could mean more than one thing. Do not choose one.
  "unknown"        you could not determine the answer.
When "problem" is not null, set the other fields to null."""

RULES_GAPS = """You answer questions about human protein targets from NIH Pharos: what each lacks for the next development level, and the dated GWAS links between targets and traits.
- Pharos grades a target Tdark, Tbio, Tchem or Tclin, lowest first. Tclin is the top level and has no gap.
- A Tdark target meets at least two of three dark criteria: a PubMed score under 5, at most 3 GeneRIFs, at most 50 antibodies. A criterion clears at a PubMed score of 5, at 4 GeneRIFs, or at 51 antibodies.
- "Criteria to clear" is the number of dark criteria the target meets, less one: 1 or 2. The nearest criterion is the one whose shortfall is the smallest share of its threshold.
- Tdark targets are grouped in bands by how far they are from leaving Tdark. For each band a back-test gives the share that was no longer Tdark in Pharos 4.0.0.
- A target with no qualifying compound lacks one at its family's potency cutoff: 30 nM for kinases, 100 nM for GPCRs and nuclear receptors, 10 uM for ion channels, 1 uM for every other target.
Explain briefly, then end with one last line of the form  FINAL: {...}  holding the JSON object the question asks for.
Every FINAL object has a "problem" field:
  null             you found what was asked.
  "not_in_graph"   the target is not in the graph. Do not answer from your own knowledge.
  "ambiguous"      a name in the question matches nothing in the graph exactly and could mean more than one thing. Do not choose one.
  "unknown"        you could not determine the answer.
When "problem" is not null, set the other fields to null."""

SQL_SCHEMA_GAPS = """
  target_gaps(target_id, level, family, lacks_compound, potency_cutoff, pubmed_score, generifs, antibodies, criteria_to_clear, route_cost, nearest_route, proposal, band)
      -- one row per target with a gap. The columns from pubmed_score to band are null unless the target is Tdark. nearest_route names the
      -- nearest criteria, joined by '|'; proposal says in words what clearing them takes; a lower route_cost is nearer to leaving Tdark
  gap_band_rates(band, promoted_share, promoted_low, promoted_high, targets)   -- of the Tdark targets in each band, the share no longer Tdark in Pharos 4.0.0"""

RULES_BY_SCHEMA = {1: RULES, 2: RULES_V2, 3: RULES_GAPS}

SQL_SCHEMA = """You have one tool, sql(query), which runs one read-only SQLite SELECT and returns up to 50 rows.
Tables:
  targets(id, symbol, name, pharos_level, family, chr, start_bp)   -- one row per Pharos protein target; id is a UniProtKB CURIE
  traits(id, name, kind, current_ontology_term)                      -- one row per GWAS trait; kind is 'disease or phenotype' or 'measurement or other'
  links(edge_id, target_id, trait_id, first_year, last_year, n_studies)  -- one row per target-trait link; it exists from first_year
  studies(target_id, trait_id, study, pmid, year)                    -- one row per study behind a link
  trait_lineage(trait_a, trait_b, relation)                          -- trait_b is the 'same term', a 'broader' or a 'narrower' form of trait_a
  neighbours(target_a, target_b, bases_apart)                        -- targets whose genes lie within 500,000 bases"""


# ----------------------------------------------------------------------------- the questions
def _tables(bundle: Path) -> dict[str, pd.DataFrame]:
    return {n: pd.read_parquet(bundle / f"{n}.parquet") for n in ("nodes", "links", "evidence", "lineage")}


def _final(fields: str) -> str:
    return (' End with a last line  FINAL: {"status": S, ' + fields + '}  where S is one of "answered", "no_link_reported", '
            '"not_in_graph", "ambiguous", "unknown". Use null for a field that does not apply.')


def _final2(fields: str) -> str:
    return ' End with a last line  FINAL: {"problem": P, ' + fields + '}  where P is null, "not_in_graph", "ambiguous" or "unknown".'


def make_questions(bundle: Path, seed: int = SEED, schema: int = 1, fakes: tuple = SETS["dev"]["fakes"], words: tuple = SETS["dev"]["words"]) -> list[dict]:
    v2 = schema == 2
    T = _tables(bundle)
    G = ask.Graph(bundle)
    rng = random.Random(seed)
    nodes, links, ev, lin = T["nodes"], T["links"], T["evidence"], T["lineage"]
    tg = nodes[nodes.kind == "target"]
    tr = nodes[nodes.kind == "trait"]
    n_sym = tg.groupby("name").id.transform("size")
    mhc = (tg.chr == ask.MHC[0]) & (tg.start >= ask.MHC[1]) & (tg.start <= ask.MHC[2])
    deg_t = links.target.value_counts()
    deg_d = links.trait.value_counts()
    ok_t = [r.id for r in tg[(n_sym == 1) & ~mhc].itertuples() if deg_t.get(r.id, 0) >= 5
            and (lambda x: x["ok"] and not x["note"] and G.t.id[x["i"]] == r.id)(G.resolve(r.name, "target"))]
    ok_d = [r.id for r in tr[tr.disease_like].itertuples() if deg_d.get(r.id, 0) >= 40
            and (lambda x: x["ok"] and not x["note"] and G.d.id[x["i"]] == r.id)(G.resolve(r.name, "trait"))]
    sym = dict(zip(tg.id, tg.name))
    nm = dict(zip(tr.id, tr.name))
    good = links[links.target.isin(set(ok_t)) & links.trait.isin(set(ok_d))].sort_values("edge_id").reset_index(drop=True)
    pair = {(a, b) for a, b in zip(links.target, links.trait)}
    qs: list[dict] = []

    def add(fam: str, text: str, gold: dict) -> None:
        qs.append({"id": f"{fam}-{sum(1 for q in qs if q['family'] == fam) + 1}", "family": fam, "text": text, "gold": gold})

    def pick(df: pd.DataFrame, n: int) -> list:
        return [df.iloc[i] for i in rng.sample(range(len(df)), n)]

    for r in pick(good[good.n_studies >= 2], 4):
        add("first_year", f"In what year was a link between {sym[r.target]} and {nm[r.trait]} first reported?"
            + (_final2('"year": YEAR') if v2 else _final('"year": YEAR')),
            {"problem": None, "year": int(r.first_year)} if v2 else {"status": "answered", "year": int(r.first_year)})

    span = ev.groupby(["target", "trait"]).year.agg(["nunique", "median"]).reset_index()
    multi = good.merge(span[span["nunique"] >= 3], on=["target", "trait"])
    for r in pick(multi, 4):
        y = int(r["median"])
        n = int(ev[(ev.target == r.target) & (ev.trait == r.trait) & (ev.year <= y)].study.nunique())
        add("studies_as_of", f"As of {y}, how many GWAS studies supported the link between {sym[r.target]} and {nm[r.trait]}?"
            + (_final2('"count": N') if v2 else _final('"count": N')),
            {"problem": None, "count": n} if v2 else {"status": "answered", "count": n})

    for r in pick(good[good.first_year >= 2012], 4):
        y = int(r.first_year) - rng.choice([1, 2, 3])
        add("not_yet", f"As of {y}, was {sym[r.target]} linked to {nm[r.trait]}? If not, was a link reported later, and in what year?"
            + (_final2('"linked": true or false, "reported_later": YEAR or null') if v2 else _final('"reported_later": YEAR or null')),
            {"problem": None, "linked": False, "reported_later": int(r.first_year)} if v2 else {"status": "no_link_reported", "reported_later": int(r.first_year)})

    n = 0
    while n < 4:
        a, b = rng.choice(ok_t), rng.choice(ok_d)
        if (a, b) not in pair:
            add("never", f"As of 2020, was {sym[a]} linked to {nm[b]}? If not, was a link reported later, and in what year?"
                + (_final2('"linked": true or false, "reported_later": YEAR or null') if v2 else _final('"reported_later": YEAR or null')),
                {"problem": None, "linked": False, "reported_later": None} if v2 else {"status": "no_link_reported", "reported_later": None})
            n += 1

    for fake in fakes:
        if G.resolve(fake, "target").get("problem") != "not in the graph":
            raise SystemExit(f"REFUSED: the invented symbol {fake} matches something in the graph")
        add("not_in_graph", f"In what year was a link between {fake} and {nm[rng.choice(ok_d)]} first reported?"
            + (_final2('"year": YEAR or null') if v2 else _final('"year": YEAR or null')),
            {"problem": "not_in_graph"} if v2 else {"status": "not_in_graph"})

    if v2:   # the first three listed words that the graph holds under several names and under none exactly
        words = tuple(w for w in words if G.resolve(w, "trait").get("problem") == "ambiguous")[:3]
        if len(words) < 3:
            raise SystemExit("REFUSED: fewer than three of the listed words are ambiguous in this graph")
    for word in words:
        if G.resolve(word, "trait").get("problem") != "ambiguous":
            raise SystemExit(f"REFUSED: {word!r} is not ambiguous in this graph")
        add("ambiguous", f"How many targets is {word} linked to in this graph?" + (_final2('"count": N or null') if v2 else _final('"count": N or null')),
            {"problem": "ambiguous"} if v2 else {"status": "ambiguous"})

    n = 0
    for d in rng.sample(ok_d, len(ok_d)):
        top = links[links.trait == d].sort_values("n_studies", ascending=False)
        if len(top) >= 2 and top.iloc[0].n_studies > top.iloc[1].n_studies and top.iloc[0].target in set(ok_t):
            t0 = top.iloc[0].target
            add("top_then_count", f"Which target is linked to {nm[d]} by the most GWAS studies in this graph, and to how many traits in total "
                f"is that target linked?" + (_final2('"symbol": GENE SYMBOL, "count": N') if v2 else _final('"symbol": GENE SYMBOL, "count": N')),
                {"problem": None, "symbol": sym[t0], "count": int(deg_t[t0])} if v2 else {"status": "answered", "symbol": sym[t0], "count": int(deg_t[t0])})
            n += 1
            if n == 4:
                break

    by_trait = links.groupby("trait").target.agg(set)
    n = 0
    while n < 3:
        a, b = rng.sample(ok_d, 2)
        k = len(by_trait[a] & by_trait[b])
        if k >= 5:
            add("shared", f"How many targets are linked to both {nm[a]} and {nm[b]} in this graph?" + (_final2('"count": N') if v2 else _final('"count": N')),
                {"problem": None, "count": k} if v2 else {"status": "answered", "count": k})
            n += 1

    for d, y in zip(rng.sample(ok_d, 3), (2015, 2017, 2019)):
        k = int(((links.trait == d) & (links.first_year > y)).sum())
        add("new_since", f"How many links to {nm[d]} were first reported after {y}?" + (_final2('"count": N') if v2 else _final('"count": N')),
            {"problem": None, "count": k} if v2 else {"status": "answered", "count": k})

    first = {(a, b): y for a, b, y in zip(links.target, links.trait, links.first_year)}
    rel = lin[lin.a.isin(set(ok_d))]
    cands = []
    for a in ok_t:
        mine = links[links.target == a]
        have = dict(zip(mine.trait, mine.first_year))
        for r in rel[rel.b.isin(set(have))].itertuples():
            y = int(have[r.b])
            if first.get((a, r.a), 9999) > y:
                cands.append((a, r.a, y))
    cands = sorted(set(cands))
    for a, d, y in rng.sample(cands, 3):
        mine = links[(links.target == a) & (links.first_year <= y)]
        okay = sorted(set(lin[(lin.a == d) & lin.b.isin(set(mine.trait))].b))
        add("related_wording", f"As of {y}, {sym[a]} was not linked to {nm[d]}. Was it already linked to the same trait under another id, or to a "
            f"broader or narrower form of it? Name one such trait."
            + (_final2('"related_trait": TRAIT NAME OR ID or null') if v2 else _final('"related_trait": TRAIT NAME OR ID or null')),
            {"problem": None, "related_trait_any_of": okay + [nm[x] for x in okay]} if v2 else {"status": "no_link_reported", "related_trait_any_of": okay + [nm[x] for x in okay]})
    return qs


CRITERIA = {"pubmed_score": "PubMed score", "generifs": "GeneRIFs", "antibodies": "antibodies"}


def make_gap_questions(bundle: Path, seed: int, fakes: tuple, words: tuple) -> list[dict]:
    """35 questions about what a target lacks, in ten kinds. Expected answers are read from the bundle's
    tables here; `ask` is used only to choose names it resolves without a note, and names it finds ambiguous."""
    G = ask.Graph(bundle)
    rng = random.Random(seed)
    nodes = pd.read_parquet(bundle / "nodes.parquet")
    links = pd.read_parquet(bundle / "links.parquet")
    tg = pd.read_parquet(bundle / "target_gaps.parquet")
    rates = json.loads((bundle / "target_gaps.json").read_text())["dark"]["bands"]
    t = nodes[nodes.kind == "target"]
    sym = dict(zip(t.id, t.name))
    level = dict(zip(t.id, t.level))
    n_sym = t.groupby("name").id.transform("size")
    clean = {r.id for r in t[n_sym == 1].itertuples()
             if (lambda x: x["ok"] and not x["note"] and G.t.id[x["i"]] == r.id)(G.resolve(r.name, "target"))}
    dark = tg[tg.flips.notna() & tg.id.isin(clean)].sort_values("id").reset_index(drop=True)
    used: set[str] = set()
    qs: list[dict] = []

    def add(fam: str, text: str, gold: dict) -> None:
        qs.append({"id": f"{fam}-{sum(1 for q in qs if q['family'] == fam) + 1}", "family": fam, "text": text, "gold": gold})

    def pick(df: pd.DataFrame, n: int) -> list:
        df = df[~df.id.isin(used)]
        rows = [df.iloc[i] for i in rng.sample(range(len(df)), n)]
        used.update(r.id for r in rows)
        return rows

    for r in pick(dark[dark.flips == 1], 2) + pick(dark[dark.flips == 2], 2):
        add("flips", f"{sym[r.id]} is a Tdark target. How many criteria does it still have to clear to leave Tdark?" + _final2('"criteria_to_clear": N'),
            {"problem": None, "criteria_to_clear": int(r.flips)})

    one = dark[dark.flips == 1]
    for route in ("pubmed_score", "generifs", "antibodies", rng.choice(list(CRITERIA))):
        r = pick(one[one.route == route], 1)[0]
        add("nearest", f"For the Tdark target {sym[r.id]}, which one criterion is nearest to clearing: the PubMed score, GeneRIFs or antibodies?"
            + _final2('"criterion": "pubmed_score" or "generifs" or "antibodies"'), {"problem": None, "criterion": route})

    for route, clears in (("generifs", 4), ("generifs", 4), ("antibodies", 51), ("antibodies", 51)):
        r = pick(one[one.route == route], 1)[0]
        add("shortfall", f"How many more {CRITERIA[route]} does the Tdark target {sym[r.id]} need for that criterion to clear?" + _final2('"more": N'),
            {"problem": None, "more": int(round(clears - r[route]))})

    for band in rng.sample(sorted(rates), 4):
        r = pick(dark[dark.band == band], 1)[0]
        add("rate", f"Of the Tdark targets in the same band as {sym[r.id]}, what percentage were no longer Tdark in Pharos 4.0.0?" + _final2('"percent": NUMBER'),
            {"problem": None, "percent": round(100 * rates[band]["promoted"]["share"], 1)})

    bio = tg[(tg.level == "Tbio") & tg.id.isin(clean)].sort_values("id")
    for fam, cut in (("Kinase", "30 nM"), ("GPCR", "100 nM"), ("IC", "10 uM"), (None, "1 uM")):
        pool = bio[bio.family.isna()] if fam is None else bio[bio.family == fam]
        r = pick(pool, 1)[0]
        if r.cutoff != cut:
            raise SystemExit(f"REFUSED: {sym[r.id]} is in the family {fam} and its cutoff is {r.cutoff}, not {cut}")
        add("ligand_cutoff", f"{sym[r.id]} has no qualifying compound. What potency would a compound have to reach for it to count, given the target's family?"
            + _final2('"cutoff": "30 nM" or "100 nM" or "10 uM" or "1 uM"'), {"problem": None, "cutoff": cut})

    clin = t[(t.level == "Tclin") & t.id.isin(clean)].sort_values("id")
    for r in pick(clin, 3):
        add("no_gap", f"What does {sym[r.id]} lack for the next Pharos level? Give its level, and whether it has a gap of this kind."
            + _final2('"level": LEVEL, "has_gap": true or false'), {"problem": None, "level": "Tclin", "has_gap": False})

    for fake in fakes:
        if G.resolve(fake, "target").get("problem") != "not in the graph":
            raise SystemExit(f"REFUSED: the invented symbol {fake} matches something in the graph")
        add("not_in_graph", f"How many criteria does {fake} still have to clear to leave Tdark?" + _final2('"criteria_to_clear": N or null'), {"problem": "not_in_graph"})

    for word in words:
        if G.resolve(word, "target").get("problem") != "ambiguous":
            raise SystemExit(f"REFUSED: {word!r} is not ambiguous among the targets of this graph")
        add("ambiguous", f"What potency cutoff does {word} lack a compound at?" + _final2('"cutoff": CUTOFF or null'), {"problem": "ambiguous"})

    alld = tg[tg.flips.notna()].sort_values(["cost", "id"])
    fams = [f for f in sorted(set(alld.family.dropna())) if (alld.family == f).sum() >= 5]
    deg = links.target.value_counts()
    usable = [f for f in fams if (lambda x: x.cost.iloc[0] < x.cost.iloc[1] and x.id.iloc[0] in clean)(alld[alld.family == f])]
    chosen = rng.sample(usable, 3)
    for f in chosen:
        top = alld[alld.family == f].iloc[0]
        add("nearest_in_family", f"Which Tdark target in the Pharos family '{f}' is nearest to leaving Tdark?" + _final2('"symbol": GENE SYMBOL'),
            {"problem": None, "symbol": sym[top.id]})
    rest = [f for f in fams if f not in chosen and (lambda x: x.cost.iloc[2] < x.cost.iloc[3])(alld[alld.family == f])]
    for f in rng.sample(rest, 3):
        top3 = alld[alld.family == f].head(3)
        add("two_step", f"Of the three Tdark targets in the Pharos family '{f}' nearest to leaving Tdark, how many have at least one GWAS link in this graph?"
            + _final2('"count": N'), {"problem": None, "count": int(sum(deg.get(i, 0) > 0 for i in top3.id))})
    if len(qs) != 35 or any(level[i] != "Tdark" for i in dark.id):
        raise SystemExit(f"REFUSED: {len(qs)} questions were made, not 35, or a dark-gap target is not Tdark in the node table")
    return qs


def via_ask(q: dict, G) -> dict | None:
    """What the question layer itself would put in FINAL for a gaps question. Used only to cross-check
    the expected answers, which are computed without it: where the two disagree, one of them is wrong."""
    text, fam = q["text"], q["family"]
    name = re.search(r"'([^']+)'", text).group(1) if fam in ("nearest_in_family", "two_step") else None
    if fam == "nearest_in_family":
        return {"problem": None, "symbol": G.gaps(family=name, top=1)["targets"][0]["symbol"]}
    if fam == "two_step":
        deg = np.bincount(G.lt, minlength=len(G.t))
        return {"problem": None, "count": sum(int(deg[G.ti[t["id"]]] > 0) for t in G.gaps(family=name, top=3)["targets"])}
    sym = {"flips": r"^(\S+) is a Tdark", "nearest": r"Tdark target (\S+),", "shortfall": r"Tdark target (\S+) need", "rate": r"same band as (\S+),",
           "ligand_cutoff": r"^(\S+) has no", "no_gap": r"What does (\S+) lack", "not_in_graph": r"does (\S+) still", "ambiguous": r"does (\S+) lack"}[fam]
    a = G.gaps(re.search(sym, text).group(1))
    if not a["found"]:
        return {"problem": {"not in the graph": "not_in_graph", "ambiguous": "ambiguous"}.get(a["problem"], a["problem"])}
    if fam == "no_gap":
        return {"problem": None, "level": a["target"]["pharos_level"], "has_gap": bool(a["gaps"])}
    if fam == "ligand_cutoff":
        return {"problem": None, "cutoff": re.search(r"compound at (.+?) or better", a["gaps"][-1]["lacks"]).group(1)}
    d = a["gaps"][0]
    if fam == "flips":
        return {"problem": None, "criteria_to_clear": d["criteria_to_clear"]}
    if fam == "nearest":
        return {"problem": None, "criterion": d["nearest_route"]}
    if fam == "shortfall":
        return {"problem": None, "more": re.match(r"(\d+) more", d["nearest_route"]).group(1)}
    return {"problem": None, "percent": re.search(r", (\d+)% \(", d["how_often_this_closed"]).group(1)}


def digest(qs: list[dict]) -> str:
    return hashlib.sha256(json.dumps(qs, sort_keys=True).encode()).hexdigest()


def frozen_path(set_: str = "dev") -> Path:
    return paths.DOCS_DATA / SETS[set_]["file"]


def load_questions(set_: str = "dev") -> list[dict]:
    p = frozen_path(set_)
    if not p.is_file():
        raise SystemExit(f"no frozen questions; run:  python -m pharos_gap_kg.usability questions <release> --set {set_} --freeze")
    d = json.loads(p.read_text())
    if digest(d["questions"]) != d["sha256"]:
        raise SystemExit(f"REFUSED: {p.name} does not match its own digest")
    return d["questions"]


# ----------------------------------------------------------------------------- grading
def parse_final(text: str) -> dict | None:
    hits = re.findall(r"FINAL:\s*(\{.*\})", text or "")
    for h in reversed(hits):
        try:
            return json.loads(h)
        except json.JSONDecodeError:
            m = re.match(r"\{[^{}]*\}", h)
            if m:
                try:
                    return json.loads(m.group(0))
                except json.JSONDecodeError:
                    pass
    return None


def _int(v) -> int | None:
    try:
        return int(str(v).replace(",", "").strip())
    except (TypeError, ValueError):
        return None


def _problem(final: dict):
    v = final.get("problem")
    return None if v in (None, "null", "", "none", "None") else v


def grade(q: dict, final: dict | None) -> bool:
    g = q["gold"]
    if not isinstance(final, dict):
        return False
    if "problem" in g:                      # second schema: a problem, or the facts
        if _problem(final) != g["problem"]:
            return False
        if g["problem"] is not None:
            return True
    elif final.get("status") != g["status"]:
        return False
    for k, want in g.items():
        if k in ("status", "problem"):
            continue
        if k == "linked":
            if final.get("linked") not in (want, str(want).lower()):
                return False
            continue
        if k == "related_trait_any_of":
            got = str(final.get("related_trait") or "").strip().casefold()
            if got not in {str(x).casefold() for x in want}:
                return False
        elif k == "symbol":
            if str(final.get("symbol") or "").strip().casefold() != want.casefold():
                return False
        elif k == "criterion":
            if _criterion(final.get("criterion")) != want:
                return False
        elif k == "cutoff":
            if _cutoff(final.get("cutoff")) != _cutoff(want):
                return False
        elif k == "percent":
            try:
                if abs(float(str(final.get("percent")).replace("%", "").strip()) - want) > 1.0:
                    return False
            except (TypeError, ValueError):
                return False
        elif k == "has_gap":
            if final.get("has_gap") not in (want, str(want).lower()):
                return False
        elif k == "level":
            if str(final.get("level") or "").strip().casefold() != want.casefold():
                return False
        elif want is None:
            if final.get(k) not in (None, "null", ""):
                return False
        elif _int(final.get(k)) != want:
            return False
    return True


def _criterion(v) -> str | None:
    s = str(v or "").casefold()
    hits = [k for k, word in (("pubmed_score", "pubmed"), ("generifs", "generif"), ("antibodies", "antibod")) if word in s]
    return hits[0] if len(hits) == 1 else None


def _cutoff(v) -> str:
    return re.sub(r"\s+", "", str(v or "")).casefold().replace("µ", "u").replace("μ", "u").rstrip(".")


def fabricated(q: dict, final: dict | None) -> bool:
    """On a not-in-graph question: the model answered as if the gene were known."""
    if q["family"] != "not_in_graph" or not isinstance(final, dict):
        return False
    if "problem" in q["gold"]:
        return _problem(final) is None
    return final.get("status") in ("answered", "no_link_reported")


# ----------------------------------------------------------------------------- the two tool sets
def sql_db(bundle: Path) -> sqlite3.Connection:
    T = _tables(bundle)
    nb = pd.read_parquet(bundle / "neighbours.parquet")
    n = T["nodes"]
    con = sqlite3.connect(":memory:", check_same_thread=False)
    t = n[n.kind == "target"]
    pd.DataFrame({"id": t.id, "symbol": t.name, "name": t.full_name, "pharos_level": t.level, "family": t.family, "chr": t.chr, "start_bp": t.start}).to_sql("targets", con, index=False)
    d = n[n.kind == "trait"]
    pd.DataFrame({"id": d.id, "name": d.name, "kind": np.where(d.disease_like, "disease or phenotype", "measurement or other"),
                  "current_ontology_term": d.live_term}).to_sql("traits", con, index=False)
    T["links"].rename(columns={"target": "target_id", "trait": "trait_id"})[["edge_id", "target_id", "trait_id", "first_year", "last_year", "n_studies"]].to_sql("links", con, index=False)
    T["evidence"].rename(columns={"target": "target_id", "trait": "trait_id"}).to_sql("studies", con, index=False)
    T["lineage"].rename(columns={"a": "trait_a", "b": "trait_b"}).to_sql("trait_lineage", con, index=False)
    nb.rename(columns={"a": "target_a", "b": "target_b", "gap": "bases_apart"}).to_sql("neighbours", con, index=False)
    if (bundle / "target_gaps.parquet").is_file():
        g = pd.read_parquet(bundle / "target_gaps.parquet")
        pd.DataFrame({"target_id": g.id, "level": g.level, "family": g.family, "lacks_compound": g.lacks_ligand, "potency_cutoff": g.cutoff,
                      "pubmed_score": g.pubmed_score, "generifs": g.generifs, "antibodies": g.antibodies, "criteria_to_clear": g.flips.astype("float"),
                      "route_cost": g.cost, "nearest_route": g.route, "proposal": g.proposal, "band": g.band}).to_sql("target_gaps", con, index=False)
        b = json.loads((bundle / "target_gaps.json").read_text())["dark"]["bands"]
        pd.DataFrame([{"band": k, "promoted_share": v["promoted"]["share"], "promoted_low": v["promoted"]["low"], "promoted_high": v["promoted"]["high"],
                       "targets": v["promoted"]["n"]} for k, v in sorted(b.items())]).to_sql("gap_band_rates", con, index=False)
        con.execute("CREATE INDEX i6 ON target_gaps(target_id)")
    for stmt in ("CREATE INDEX i1 ON links(target_id)", "CREATE INDEX i2 ON links(trait_id)", "CREATE INDEX i3 ON studies(target_id, trait_id)",
                 "CREATE INDEX i4 ON targets(symbol)", "CREATE INDEX i5 ON traits(name)", "PRAGMA query_only = ON"):
        con.execute(stmt)
    return con


def run_sql(con: sqlite3.Connection, query: str) -> str:
    q = (query or "").strip().rstrip(";")
    if not re.match(r"(?is)^\s*(select|with)\b", q) or ";" in q:
        return "error: one SELECT statement only"
    try:
        cur = con.execute(q)
        cols = [c[0] for c in cur.description]
        rows = cur.fetchmany(51)
    except sqlite3.Error as e:
        return f"error: {e}"
    out = json.dumps({"columns": cols, "rows": rows[:50], "more_rows": len(rows) > 50}, default=str)
    return out if len(out) <= 6000 else out[:6000] + ' ...[truncated; narrow the query]'


SQL_TOOL = [{"type": "function", "function": {"name": "sql", "description": "Run one read-only SQLite SELECT over the graph's tables and return up to 50 rows.",
                                              "parameters": {"type": "object", "properties": {"query": {"type": "string"}}, "required": ["query"]}}}]


# ----------------------------------------------------------------------------- one conversation
def post(base_url: str, key: str | None, body: dict, timeout: int = 180) -> dict:
    headers = {"Content-Type": "application/json"}
    if key:
        headers["Authorization"] = f"Bearer {key}"
    if "openrouter.ai" in base_url:
        headers.update({"X-Title": "pharos-gap-kg usability"})
    req = urllib.request.Request(base_url.rstrip("/") + "/chat/completions", data=json.dumps(body).encode(), headers=headers)
    last = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}: {e.read()[:300].decode('utf-8', 'replace')}"
            if e.code in (400, 401, 402, 403, 404):
                break
        except Exception as e:  # timeouts and resets are retried
            last = repr(e)
        time.sleep(2 * (attempt + 1))
    return {"error": last}


async def converse(q: dict, system: str, tools: list[dict], call_tool, model: str, base_url: str, key: str | None, max_tokens: int) -> dict:
    msgs = [{"role": "system", "content": system}, {"role": "user", "content": q["text"]}]
    usage = {"in": 0, "out": 0, "cost": 0.0}
    used: list[str] = []
    text, error = "", None
    for turn in range(MAX_TOOL_CALLS + 2):
        body = {"model": model, "messages": msgs, "max_tokens": max_tokens, "temperature": 0}
        if "openrouter.ai" in base_url:
            body["usage"] = {"include": True}
        if tools and len(used) < MAX_TOOL_CALLS:
            body["tools"] = tools
        r = await asyncio.to_thread(post, base_url, key, body)
        if "error" in r or not r.get("choices"):
            error = str(r.get("error"))[:300]
            break
        u = r.get("usage") or {}
        usage["in"] += int(u.get("prompt_tokens") or 0)
        usage["out"] += int(u.get("completion_tokens") or 0)
        usage["cost"] += float(u.get("cost") or 0.0)
        m = r["choices"][0]["message"]
        text = m.get("content") or ""
        calls = m.get("tool_calls") or []
        if not calls:
            if parse_final(text) is None and turn < MAX_TOOL_CALLS + 1:
                msgs += [{"role": "assistant", "content": text or "(no text)"},
                         {"role": "user", "content": "Give your last line now, in the form  FINAL: {...}"}]
                tools = []
                continue
            break
        msgs.append({"role": "assistant", "content": text, "tool_calls": calls})
        for c in calls:
            name = c["function"]["name"]
            try:
                args = json.loads(c["function"].get("arguments") or "{}")
                out = await call_tool(name, args)
            except Exception as e:
                out = f"error: {e}"
            used.append(name)
            msgs.append({"role": "tool", "tool_call_id": c.get("id", name), "content": out[:8000]})
    final = parse_final(text)
    return {"final": final, "text": text[-1500:], "tools_used": used, "error": error, **usage}


def _key(model: str, condition: str, qid: str, fp: str, set_: str = "dev") -> str:
    head = f"{PROMPT_VERSION}" if set_ == "dev" else f"{set_}"
    return hashlib.sha256(f"{head}|{model}|{condition}|{qid}|{fp}".encode()).hexdigest()[:24]


def calls_dir() -> Path:
    return paths.EVAL / "usability" / "calls"


async def run(bundle: Path, model: str, condition: str, base_url: str, key: str | None, max_cost: float, workers: int, max_tokens: int,
              only: set[str] | None = None, set_: str = "dev", record_failures: bool = False) -> list[dict]:
    qs = [q for q in load_questions(set_) if not only or q["id"] in only]
    schema = SETS[set_]["schema"]
    v2 = schema >= 2
    rules = RULES_BY_SCHEMA[schema]
    out_dir = calls_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    if schema == 3:      # the gaps test has one cap over its two sets, whatever a single run is allowed
        gone = sum(r.get("cost") or 0.0 for r in (json.loads(p.read_text()) for p in out_dir.glob("*.json")) if SETS.get(r.get("set", "dev"), {}).get("schema") == 3)
        max_cost = min(max_cost, GAPS_CAP - gone)
        if max_cost <= 0:
            raise SystemExit(f"REFUSED: the gaps test has spent {gone:.2f} of its {GAPS_CAP:.2f} dollars")
    spent = 0.0
    sem = asyncio.Semaphore(workers)

    async def body(session, system: str, tools: list[dict], call_tool, fp: str) -> list[dict]:
        nonlocal spent

        async def one(q: dict) -> dict:
            nonlocal spent
            f = out_dir / f"{_key(model, condition, q['id'], fp, set_)}.json"
            if f.is_file():
                return json.loads(f.read_text())
            async with sem:
                if spent >= max_cost:
                    return {"model": model, "condition": condition, "id": q["id"], "skipped": "cost cap reached"}
                r = await converse(q, system, tools, call_tool, model, base_url, key, max_tokens)
                if not r["cost"] and model in PRICES:
                    r["cost"] = (r["in"] * PRICES[model][0] + r["out"] * PRICES[model][1]) / 1e6
                spent += r["cost"]
                row = {"set": set_, "model": model, "condition": condition, "id": q["id"], "family": q["family"], "tool_fingerprint": fp, **r}
                if r["error"] and record_failures:      # the provider failed again: kept as a question that got no answer
                    row.update({"final": None, "gave_up": True})
                if not r["error"] or record_failures:
                    f.write_text(json.dumps(row, indent=1) + "\n")
                return row

        return list(await asyncio.gather(*(one(q) for q in qs)))

    if condition == "typed":
        from mcp import ClientSession, StdioServerParameters, stdio_client
        params = StdioServerParameters(command=sys.executable, args=["-m", "pharos_gap_kg.server", "--bundle", str(bundle)])
        async with stdio_client(params) as streams:
            async with ClientSession(streams[0], streams[1]) as s:
                init = await s.initialize()
                listed = (await s.list_tools()).tools
                specs = [{"type": "function", "function": {"name": t.name, "description": t.description,
                                                           "parameters": getattr(t, "input_schema", None) or getattr(t, "inputSchema")}} for t in listed]
                # first round: the server's instructions and the format line only. Second round: the same
                # rules text as every other condition, after the server's instructions.
                system = (init.instructions or "") + "\n\n" + (rules if v2 else RULES.split("\n")[-1])
                lock = asyncio.Lock()

                async def call_tool(name: str, args: dict) -> str:
                    async with lock:
                        r = await s.call_tool(name, args)
                    sc = getattr(r, "structured_content", None) or getattr(r, "structuredContent", None)
                    return json.dumps(sc) if sc else (r.content[0].text if r.content else "")

                fp = hashlib.sha256(json.dumps([system, specs], sort_keys=True).encode()).hexdigest()[:12]
                rows = await body(s, system, specs, call_tool, fp)
    elif condition == "sql":
        con = sql_db(bundle)

        async def call_tool(name: str, args: dict) -> str:
            return run_sql(con, args.get("query", "")) if name == "sql" else f"unknown tool {name}; the only tool is sql"

        system = rules + "\n\n" + SQL_SCHEMA + (SQL_SCHEMA_GAPS if schema == 3 else "")
        fp = hashlib.sha256(json.dumps([system, SQL_TOOL], sort_keys=True).encode()).hexdigest()[:12]
        rows = await body(None, system, SQL_TOOL, call_tool, fp)
    else:
        system = rules + "\n\nYou have no tool, and the graph is not available to you."
        fp = hashlib.sha256(system.encode()).hexdigest()[:12]

        async def call_tool(name: str, args: dict) -> str:
            return "no tool is available"

        rows = await body(None, system, [], call_tool, fp)
    return rows


# ----------------------------------------------------------------------------- the report
def report(set_: str = "dev") -> dict:
    qs = {q["id"]: q for q in load_questions(set_)}
    rows = [json.loads(p.read_text()) for p in sorted(calls_dir().glob("*.json"))]
    rows = [r for r in rows if r["id"] in qs and r.get("set", "dev") == set_]
    latest: dict[tuple, dict] = {}
    for r in rows:
        latest[(r["model"], r["condition"], r["id"])] = r
    cells: dict[tuple, dict] = {}
    for (model, cond, qid), r in latest.items():
        c = cells.setdefault((model, cond), {"n": 0, "correct": 0, "fabricated": 0, "no_final": 0, "cost": 0.0, "calls": 0, "by_family": {}, "right": {}})
        ok = grade(qs[qid], r.get("final"))
        c["n"] += 1
        c["correct"] += ok
        c["fabricated"] += fabricated(qs[qid], r.get("final"))
        c["no_final"] += r.get("final") is None
        c["cost"] += r.get("cost") or 0.0
        c["calls"] += len(r.get("tools_used") or [])
        f = c["by_family"].setdefault(qs[qid]["family"], [0, 0])
        f[0] += ok
        f[1] += 1
        c["right"][qid] = bool(ok)
    out = {"set": set_, "questions": len(qs), "sha256": json.loads(frozen_path(set_).read_text())["sha256"],
           "cells": {f"{m}|{c}": v for (m, c), v in sorted(cells.items())}}
    out.update(compare(out["cells"]))
    return out


HOSTED = list(PRICES)
PASS_MARK = 28


def compare(cells: dict) -> dict:
    """Typed against sql for each model, the pooled interval over hosted models, and the two verdicts
    the protocol stated in advance."""
    models = sorted({k.split("|")[0] for k in cells})
    per = {}
    for m in models:
        a, b = cells.get(f"{m}|typed"), cells.get(f"{m}|sql")
        if a and b:
            per[m] = paired(a["right"], b["right"])
    hosted = [m for m in HOSTED if f"{m}|typed" in cells and f"{m}|sql" in cells]
    pooled = None
    if hosted:
        ids = sorted(set.intersection(*(set(cells[f"{m}|{c}"]["right"]) for m in hosted for c in ("typed", "sql"))))
        d = np.array([np.mean([cells[f"{m}|typed"]["right"][i] for m in hosted]) - np.mean([cells[f"{m}|sql"]["right"][i] for m in hosted]) for i in ids])
        rng = np.random.default_rng(0)
        boot = d[rng.integers(0, len(d), (2000, len(d)))].mean(axis=1)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        pooled = {"models": hosted, "questions": len(ids), "difference": float(d.mean()), "low": float(lo), "high": float(hi)}
    typed = {m: cells[f"{m}|typed"] for m in HOSTED if f"{m}|typed" in cells}
    complete = len(typed) == len(HOSTED) and all(v["n"] == len(v["right"]) and v["n"] >= 35 for v in typed.values())
    usable = complete and sum(v["correct"] >= PASS_MARK for v in typed.values()) >= 2 and sum(v["fabricated"] for v in typed.values()) == 0
    return {"typed_minus_sql": per, "pooled_typed_minus_sql": pooled,
            "verdicts": {"all_hosted_cells_complete": complete,
                         "usable_by_a_small_model": bool(usable) if complete else None,
                         "typed_better_than_sql": (bool(pooled["low"] > 0) if pooled and complete else None)}}


def paired(a: dict[str, bool], b: dict[str, bool]) -> dict:
    """a minus b over the questions both answered: discordant counts and an exact two-sided sign test."""
    from math import comb
    ids = sorted(set(a) & set(b))
    plus = sum(a[i] and not b[i] for i in ids)
    minus = sum(b[i] and not a[i] for i in ids)
    n = plus + minus
    p = min(1.0, 2 * sum(comb(n, k) for k in range(0, min(plus, minus) + 1)) / 2 ** n) if n else 1.0
    return {"questions": len(ids), "a_only": plus, "b_only": minus, "p_sign": p}


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("the last FINAL line is the answer, and text around it is ignored",
       parse_final('thinking FINAL: {"status": "unknown"}\nmore\nFINAL: {"status": "answered", "year": 2007} done') == {"status": "answered", "year": 2007})
    ck("no FINAL line, or one that is not JSON, is no answer", parse_final("the year is 2007") is None and parse_final("FINAL: {year 2007}") is None)
    q = {"family": "first_year", "gold": {"status": "answered", "year": 2007}}
    ck("a year is right as a number or as text, and wrong by one is wrong",
       grade(q, {"status": "answered", "year": 2007}) and grade(q, {"status": "answered", "year": "2007"}) and not grade(q, {"status": "answered", "year": 2008}))
    ck("the right number under the wrong status is wrong", not grade(q, {"status": "no_link_reported", "year": 2007}) and not grade(q, None))
    q = {"family": "never", "gold": {"status": "no_link_reported", "reported_later": None}}
    ck("a field that should be null must be null", grade(q, {"status": "no_link_reported", "reported_later": None}) and not grade(q, {"status": "no_link_reported", "reported_later": 2021}))
    q = {"family": "top_then_count", "gold": {"status": "answered", "symbol": "APOE", "count": 12}}
    ck("a symbol is matched without regard to case, and both fields must be right",
       grade(q, {"status": "answered", "symbol": "apoe", "count": "12"}) and not grade(q, {"status": "answered", "symbol": "APOE", "count": 11}))
    q = {"family": "related_wording", "gold": {"status": "no_link_reported", "related_trait_any_of": ["MONDO:1", "heart disease"]}}
    ck("a related trait is right by id or by name, and another trait is wrong",
       grade(q, {"status": "no_link_reported", "related_trait": "Heart Disease"}) and grade(q, {"status": "no_link_reported", "related_trait": "MONDO:1"})
       and not grade(q, {"status": "no_link_reported", "related_trait": "lung disease"}))
    q = {"family": "not_in_graph", "gold": {"status": "not_in_graph"}}
    ck("answering about a gene that is not in the graph is counted as invention; declining is not",
       fabricated(q, {"status": "answered", "year": 2010}) and not fabricated(q, {"status": "not_in_graph"}) and not fabricated(q, {"status": "unknown"}))
    q2 = {"family": "not_yet", "gold": {"problem": None, "linked": False, "reported_later": 2020}}
    ck("second schema: a found absence is no problem, with linked false and the year",
       grade(q2, {"problem": None, "linked": False, "reported_later": 2020}) and grade(q2, {"problem": "null", "linked": "false", "reported_later": "2020"})
       and not grade(q2, {"problem": None, "linked": True, "reported_later": 2020}) and not grade(q2, {"problem": "unknown", "linked": None, "reported_later": None}))
    q2 = {"family": "not_in_graph", "gold": {"problem": "not_in_graph"}}
    ck("second schema: a problem is right by itself, and answering a not-in-graph question is invention",
       grade(q2, {"problem": "not_in_graph", "year": None}) and not grade(q2, {"problem": None, "year": 2010})
       and fabricated(q2, {"problem": None, "year": 2010}) and not fabricated(q2, {"problem": "unknown"}))
    q3 = {"family": "nearest", "gold": {"problem": None, "criterion": "generifs"}}
    ck("a criterion is right however it is written, and naming two is wrong",
       grade(q3, {"problem": None, "criterion": "GeneRIFs"}) and grade(q3, {"problem": None, "criterion": "generif count"})
       and not grade(q3, {"problem": None, "criterion": "antibodies"}) and not grade(q3, {"problem": None, "criterion": "GeneRIFs or antibodies"}))
    q3 = {"family": "ligand_cutoff", "gold": {"problem": None, "cutoff": "10 uM"}}
    ck("a cutoff is right whatever its spacing or its micro sign, and another cutoff is wrong",
       grade(q3, {"problem": None, "cutoff": "10uM"}) and grade(q3, {"problem": None, "cutoff": "10 µM"}) and not grade(q3, {"problem": None, "cutoff": "1 uM"}))
    q3 = {"family": "rate", "gold": {"problem": None, "percent": 53.6}}
    ck("a percentage is right within one point, as a number or as text",
       grade(q3, {"problem": None, "percent": 54}) and grade(q3, {"problem": None, "percent": "53.6%"}) and not grade(q3, {"problem": None, "percent": 51})
       and not grade(q3, {"problem": None, "percent": None}))
    q3 = {"family": "no_gap", "gold": {"problem": None, "level": "Tclin", "has_gap": False}}
    ck("a level and a yes or no are both needed", grade(q3, {"problem": None, "level": "tclin", "has_gap": "false"}) and not grade(q3, {"problem": None, "level": "Tclin", "has_gap": True})
       and not grade(q3, {"problem": None, "level": "Tchem", "has_gap": False}))
    mk = lambda right, fab=0: {"n": 35, "correct": sum(right.values()), "fabricated": fab, "right": right}
    good = {f"q{i}": i < 30 for i in range(35)}
    poor = {f"q{i}": i < 20 for i in range(35)}
    C = {f"{m}|{c}": mk(good if c == "typed" else poor) for m in HOSTED for c in ("typed", "sql")}
    v = compare(C)
    ck("the stated verdicts follow from the cells: usable at 28 or more for two models with nothing invented, and better than sql when the interval clears zero",
       v["verdicts"] == {"all_hosted_cells_complete": True, "usable_by_a_small_model": True, "typed_better_than_sql": True}
       and abs(v["pooled_typed_minus_sql"]["difference"] - 10 / 35) < 1e-9)
    C2 = dict(C)
    C2[f"{HOSTED[0]}|typed"] = mk(good, fab=1)
    ck("one invented answer under typed fails the usability verdict", compare(C2)["verdicts"]["usable_by_a_small_model"] is False)
    C3 = {k: v_ for k, v_ in C.items() if not k.startswith(HOSTED[2])}
    ck("with a hosted model missing there is no verdict", compare(C3)["verdicts"]["usable_by_a_small_model"] is None)
    p = paired({"a": True, "b": True, "c": False, "d": True}, {"a": True, "b": False, "c": False, "d": False})
    ck("the paired test counts questions only one side got right", (p["a_only"], p["b_only"], p["questions"]) == (2, 0, 4) and abs(p["p_sign"] - 0.5) < 1e-12)
    import tempfile
    from pharos_gap_kg import ask_selftest
    with tempfile.TemporaryDirectory() as tmp:
        ask_selftest.tiny(Path(tmp))
        con = sql_db(Path(tmp))
        ck("the SQL tool answers a SELECT", json.loads(run_sql(con, "select count(*) as n from links"))["rows"] == [[11]])
        ck("the SQL tool holds the gap table and the band rates where the bundle has them",
           json.loads(run_sql(con, "select criteria_to_clear, nearest_route from target_gaps where level = 'Tdark'"))["rows"] == [[1.0, "pubmed_score"]]
           and json.loads(run_sql(con, "select promoted_share, targets from gap_band_rates"))["rows"] == [[0.5, 40]])
        ck("the SQL tool refuses anything that is not one SELECT", run_sql(con, "delete from links").startswith("error") and run_sql(con, "select 1; select 2").startswith("error"))
        ck("a bad query returns its error as text", run_sql(con, "select nope from links").startswith("error"))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    qp = sub.add_parser("questions")
    qp.add_argument("release")
    qp.add_argument("--freeze", action="store_true")
    qp.add_argument("--set", dest="set_", choices=sorted(SETS), default="dev")
    rp = sub.add_parser("run")
    rp.add_argument("release")
    rp.add_argument("--model", required=True)
    rp.add_argument("--condition", choices=CONDITIONS, required=True)
    rp.add_argument("--base-url", default="https://openrouter.ai/api/v1")
    rp.add_argument("--key-env", default="OPENROUTER_API_KEY")
    rp.add_argument("--max-cost", type=float, default=0.25)
    rp.add_argument("--workers", type=int, default=6)
    rp.add_argument("--max-tokens", type=int, default=4000)
    rp.add_argument("--only", default="")
    rp.add_argument("--record-failures", action="store_true", help="keep a conversation the provider failed as a question with no answer (use on the third attempt)")
    rp.add_argument("--set", dest="set_", choices=sorted(SETS), default="dev")
    sub.add_parser("report").add_argument("--set", dest="set_", choices=sorted(SETS), default="dev")
    a = ap.parse_args(argv)
    if a.cmd == "questions":
        S = SETS[a.set_]
        qs = (make_gap_questions(paths.bundle(a.release), S["seed"], S["fakes"], S["words"]) if S["schema"] == 3
              else make_questions(paths.bundle(a.release), S["seed"], S["schema"], S["fakes"], S["words"]))
        d = {"release": a.release, "seed": S["seed"], "sha256": digest(qs), "questions": qs}
        if S["schema"] == 3:
            G = ask.Graph(paths.bundle(a.release))
            differ = [q["id"] for q in qs if not grade(q, via_ask(q, G))]
            if differ:
                raise SystemExit(f"REFUSED: the expected answer and the question layer disagree on {differ}; one of them is wrong. Nothing was written.")
            print(f"the question layer gives the expected answer to all {len(qs)}")
        fz = frozen_path(a.set_)
        if fz.is_file():
            old = json.loads(fz.read_text())
            if old["sha256"] != d["sha256"]:
                raise SystemExit(f"REFUSED: the questions made now differ from the frozen set in {fz.name}. Nothing was written.")
            print(f"the {len(qs)} questions match the frozen set {d['sha256'][:12]}")
        elif a.freeze:
            fz.write_text(json.dumps(d, indent=1) + "\n")
            print(f"froze {len(qs)} questions {d['sha256'][:12]} -> {fz}")
        else:
            raise SystemExit("no frozen set yet; pass --freeze to write it")
        for q in qs:
            print(f"  {q['id']:18s} {q['text'].split(' End with')[0][:110]}  => {json.dumps(q['gold'])[:80]}")
        return 0
    if a.cmd == "run":
        key = os.environ.get(a.key_env) if a.key_env else None
        if "openrouter.ai" in a.base_url and not key:
            raise SystemExit(f"{a.key_env} is not set")
        rows = asyncio.run(run(paths.bundle(a.release), a.model, a.condition, a.base_url, key, a.max_cost, a.workers, a.max_tokens,
                               set(a.only.split(",")) if a.only else None, a.set_, a.record_failures))
        qs = {q["id"]: q for q in load_questions(a.set_)}
        done = [r for r in rows if "final" in r]
        right = sum(grade(qs[r["id"]], r["final"]) for r in done)
        print(f"{a.model} / {a.condition}: {right} of {len(done)} correct; no FINAL line {sum(r['final'] is None for r in done)}; "
              f"errors {sum(bool(r.get('error')) for r in done)}; skipped {len(rows) - len(done)}; "
              f"cost ${sum(r.get('cost') or 0 for r in done):.4f}; tool calls {sum(len(r['tools_used']) for r in done)}")
        return 0
    r = report(a.set_)
    out = paths.DOCS_DATA / SETS[a.set_]["results"]
    out.write_text(json.dumps(r, indent=1) + "\n")
    for k, v in r["cells"].items():
        print(f"  {k:44s} {v['correct']:>2d} of {v['n']:>2d}  invented {v['fabricated']}  no FINAL {v['no_final']}  calls {v['calls']:>3d}  ${v['cost']:.4f}")
    fams = sorted({f for v in r["cells"].values() for f in v["by_family"]})
    print("  by kind (correct of asked):")
    for k, v in r["cells"].items():
        print(f"    {k:44s} " + " ".join(f"{f[:9]}:{v['by_family'].get(f, [0, 0])[0]}/{v['by_family'].get(f, [0, 0])[1]}" for f in fams))
    for m, pz in r["typed_minus_sql"].items():
        print(f"  typed vs sql, {m}: typed only {pz['a_only']}, sql only {pz['b_only']}, sign test p = {pz['p_sign']:.3f}")
    if r["pooled_typed_minus_sql"]:
        z = r["pooled_typed_minus_sql"]
        print(f"  pooled over {len(z['models'])} hosted models: typed minus sql {z['difference']:+.3f} [{z['low']:+.3f}, {z['high']:+.3f}] of a question")
    print(f"  verdicts: {r['verdicts']}")
    print(f"  total cost ${sum(v['cost'] for v in r['cells'].values()):.4f}")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
