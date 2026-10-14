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

from pharos_gap_kg import ask, gaptables, paths, ranktables

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
    # the rankings (docs/PROTOCOL.md, "Usability of the rankings"): a development set, then the set the verdicts are read on
    "rank_dev": {"seed": 20261012, "schema": 4, "fakes": ("QXJV5", "ZKWP2", "VJQT8"), "words": ("migraine", "cirrhosis"),
                 "file": "usability_questions_rank_dev.json", "results": "usability_results_rank_dev.json"},
    "rank": {"seed": 20261013, "schema": 4, "fakes": ("XZQJ4", "KVWP6", "JQZT3"), "words": ("cardiomyopathy", "colitis"),
             "file": "usability_questions_rank.json", "results": "usability_results_rank.json"},
    # the other sources and the gaps of a disease (docs/PROTOCOL.md, "Usability of the other sources and of the gaps of a disease").
    # fakes: a gene, a drug and a disease that are in no source
    "new_dev": {"seed": 20261014, "schema": 5, "fakes": ("QVZK7", "zelvotamab", "Harwick syndrome"), "words": ("neuropathy", "dystrophy"),
                "file": "usability_questions_new_dev.json", "results": "usability_results_new_dev.json"},
    "new": {"seed": 20261015, "schema": 5, "fakes": ("JXQZ8", "brivoquamide", "Tessler disease"), "words": ("retinopathy", "encephalopathy"),
            "file": "usability_questions_new.json", "results": "usability_results_new.json"},
}
STATUSES = ["answered", "no_link_reported", "not_in_graph", "ambiguous", "unknown"]
CONDITIONS = ["typed", "sql", "none"]
MAX_TOOL_CALLS = 8
PROMPT_VERSION = 1
GAPS_CAP = 0.75          # dollars, over both sets of the gaps test (docs/PROTOCOL.md)
RANK_CAP = 0.60          # dollars, over both sets of the rankings test
NEW_CAP = 0.60           # dollars, over both sets of the test of the other sources and the gaps of a disease
CAPS = {3: GAPS_CAP, 4: RANK_CAP, 5: NEW_CAP}
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

RULES_RANK = """You answer questions about a graph of dated GWAS links between human protein targets (from NIH Pharos) and traits, and about the ranked candidates it offers for links not yet reported.
- A candidate for a trait is a target that has links to other traits and none to this one. Candidates are ranked by a score with three parts: shared-trait neighbours, the number of traits the target is already linked to, and the number of genes within 500,000 bases of it that already have the link. Rank 1 is the highest.
- A candidate is "explained by an earlier link" when a gene within 500,000 bases of it already has the link, or when the target is already linked to a broader, narrower or second wording of the trait.
- Beside each candidate the graph gives the rate at which candidates at its ranks, in its own group (explained or not), were reported in a test on later links. That rate is the answer when a question asks how often a rank has held.
- A target with no link at all is not among those candidates. Such targets are ranked apart, by position alone: for a trait, by how many genes within 500,000 bases of the target already have the link; and for such a target, the traits are ranked by the same count.
- A target with no link and no linked gene nearby has nothing to rank. That is an answer: "problem" is null and the field asked for is null.
Explain briefly, then end with one last line of the form  FINAL: {...}  holding the JSON object the question asks for.
Every FINAL object has a "problem" field:
  null             you found what was asked.
  "not_in_graph"   the target or the trait is not in the graph. Do not answer from your own knowledge.
  "ambiguous"      a name in the question matches nothing in the graph exactly and could mean more than one thing. Do not choose one.
  "unknown"        you could not determine the answer.
When "problem" is not null, set the other fields to null."""

SQL_SCHEMA_RANK = """
  candidates(trait_id, rank, target_id, symbol, score, shared_trait_neighbours, traits_already_linked, genes_nearby_with_this_link, explained)
      -- the 200 highest-ranked candidate targets for each trait, rank 1 first. explained is 1 when an earlier link explains the candidate, else 0
  rank_rates(ranks, candidate_group, percent)
      -- of the candidates at those ranks for a trait, the percentage reported in the test on later links. ranks is '1-10', '11-50', '51-100'
      -- or '101-500'. candidate_group is 'explained', 'not explained, trait is a disease or phenotype' or 'not explained, other trait'
  no_link_targets(trait_id, rank, target_id, symbol, genes_nearby_with_this_link)
      -- for each trait, the targets that have no link at all and lie beside a gene that has this link, rank 1 first
  no_link_traits(target_id, symbol, rank, trait_id, trait_name, genes_nearby_with_this_link)
      -- for each target that has no link at all, the traits its neighbouring genes are linked to, rank 1 first. A target with no link and
      -- no row here has nothing to rank"""

RULES_NEW = """You answer questions about a graph built from NIH Pharos. It holds dated GWAS links between human protein targets and traits and, beside them, what nine disease sources assert, drugs with their targets and indications, and expression observations.
- A disease source ASSERTS a target with a condition. A row that a pruning rule withheld is not an assertion: it is reported apart, with its rule and the number of papers that mention both.
- An expression observation says a gene's expression was higher or lower in a condition. It is not an assertion, and it never counts as a source asserting the pair.
- A drug acts on targets and is indicated for conditions. A target may be marked as the drug's mechanism of action.
- A gene is TIED to a condition by a GWAS link, by assertions from two or more disease sources, or by being the mechanism target of a drug indicated for it. One source alone, an expression observation or a withheld row is not a tie.
- The gaps of a disease are counted over the genes tied to it: those with no drug acting on them and no compound at their family's potency cutoff; those still Tdark; the loci that carry the GWAS link on more than one gene; and the genes tied with no GWAS link. A count is of genes, not of database entries.
- For the largest study behind a trait's GWAS links, the graph gives the smallest odds ratio a variant carried by 20% of chromosomes needed in order to be found.
- A condition that the other sources name and that is not a GWAS trait has no dated link, so nothing about its GWAS links or its loci can be counted. That is an answer: "problem" is null and the field asked for is null.
- Only the GWAS links carry dates. Everything else is as of the release.
Explain briefly, then end with one last line of the form  FINAL: {...}  holding the JSON object the question asks for.
Every FINAL object has a "problem" field:
  null             you found what was asked.
  "not_in_graph"   the gene, the drug or the condition is not in the graph. Do not answer from your own knowledge.
  "ambiguous"      a name in the question matches nothing in the graph exactly and could mean more than one thing. Do not choose one.
  "unknown"        you could not determine the answer.
When "problem" is not null, set the other fields to null."""

SQL_SCHEMA_NEW = """
  conditions(id, name)                                               -- conditions the other sources name that are not GWAS traits. Two ids
                                                                     -- with one name are one condition: same_condition lists each under the other
  drugs(id, name)
  same_condition(asked_id, read_id, tie)
      -- one trait or condition is read under every read_id listed for its id. A row of the four tables below belongs to a trait or a
      -- condition when its condition_id is one of those read_ids, so join through this table
  assertions(target_id, condition_id, source, withheld, rule, papers_mentioning_both)
      -- one row per source that names the pair. withheld is 1 when a pruning rule withheld the row: then it is not an assertion
  expression(target_id, condition_id, condition_as_named, log2_fold_change, p)   -- observations; above zero is higher in the condition
  drug_targets(drug_id, target_id, action, is_mechanism)             -- is_mechanism is 1 when the target is the drug's mechanism of action
  drug_indications(drug_id, condition_id)
  disease_gap_summary(condition_id, name, is_a_gwas_trait, genes_tied, genes_with_a_gwas_link, genes_asserted_by_two_or_more_sources,
                      genes_asserted_by_one_source_only, drugs_indicated, genes_that_are_a_mechanism_target, no_compound, dark, loci_not_resolved,
                      genes_in_those_loci, loci_where_other_evidence_singles_out_one_gene, no_genetic_support, withheld_rows, observed_only,
                      largest_study, largest_study_year, cases, controls, individuals, smallest_odds_ratio_it_could_find, smallest_effect_it_could_find_in_sd)
      -- one row per disease trait and per condition, by its own id (join same_condition only for the four tables above). no_compound and dark
      -- count tied genes; loci_not_resolved counts loci that carry the GWAS link on more than one gene; no_genetic_support counts genes tied
      -- with no GWAS link. Those two are NULL when is_a_gwas_trait is 0
  disease_gap_genes(condition_id, list, symbol)                      -- list is 'tied', 'no_compound', 'dark', 'no_genetic_support' or 'observed_only'
  disease_loci(condition_id, locus, genes_with_the_link, genes, gene_singled_out_by_other_evidence)"""

RULES_BY_SCHEMA = {1: RULES, 2: RULES_V2, 3: RULES_GAPS, 4: RULES_RANK, 5: RULES_NEW}
SQL_EXTRA = {3: SQL_SCHEMA_GAPS, 4: SQL_SCHEMA_RANK, 5: SQL_SCHEMA_GAPS + SQL_SCHEMA_NEW}

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


TIE = 1e-9          # two scores closer than this are a tie, and a question is never put on one


def make_rank_questions(bundle: Path, seed: int, fakes: tuple, words: tuple) -> list[dict]:
    """35 questions about the rankings, in eleven kinds. Expected answers are read from `ranktables`,
    which computes every ranking from the bundle's tables without the question layer; `ask` is used
    only to choose names it resolves without a note, and names it finds ambiguous."""
    G = ask.Graph(bundle)
    rng = random.Random(seed)
    R = ranktables.build(bundle)
    C, NT, ND = R["candidates"], R["no_link_targets"], R["no_link_traits"]
    nodes = pd.read_parquet(bundle / "nodes.parquet")
    t, d = nodes[nodes.kind == "target"], nodes[nodes.kind == "trait"]
    n_sym = t.groupby("name").id.transform("size")
    clean_t = {r.id for r in t[n_sym == 1].itertuples() if (lambda x: x["ok"] and not x["note"] and G.t.id[x["i"]] == r.id)(G.resolve(r.name))}
    clean_d = {r.id for r in d.itertuples() if (lambda x: x["ok"] and not x["note"] and G.d.id[x["i"]] == r.id)(G.resolve(r.name))}
    dname, disease = dict(zip(d.id, d.name)), dict(zip(d.id, d.disease_like))
    rates = {(r.ranks, r.candidate_group): r.percent for r in R["rank_rates"].itertuples()}
    bin_of = lambda rank: next(f"{lo}-{hi}" for lo, hi in ask.TRAIT_BINS if lo <= rank <= hi)
    group_of = lambda row: ("explained" if row.explained else ranktables.GROUPS["disease_like_not_linked"] if disease[row.trait_id] else ranktables.GROUPS["not_linked"])
    by = {k: g.reset_index(drop=True) for k, g in C[C.trait_id.isin(clean_d)].groupby("trait_id", sort=True)}
    by = {k: g for k, g in by.items() if len(g) >= 150}
    alone = lambda g, i: (i == 0 or g.score[i - 1] - g.score[i] > TIE) and (i + 1 == len(g) or g.score[i] - g.score[i + 1] > TIE)       # not tied with a neighbour
    lead = lambda g: g.target_id[0] in clean_t and alone(g, 0)
    used_d, used_t = set(), set()
    qs: list[dict] = []

    def add(fam: str, text: str, gold: dict) -> None:
        qs.append({"id": f"{fam}-{sum(1 for q in qs if q['family'] == fam) + 1}", "family": fam, "text": text, "gold": gold})

    def traits(keep, n: int) -> list[str]:
        pool = sorted(k for k, g in by.items() if k not in used_d and keep(g))
        got = rng.sample(pool, n)
        used_d.update(got)
        return got

    for k in traits(lead, 4):
        add("top", f"Which target is the highest-ranked candidate for {dname[k]}?" + _final2('"symbol": GENE SYMBOL'), {"problem": None, "symbol": by[k].symbol[0]})

    def a_candidate(g, lo: int, hi: int):
        ok = [i for i in range(lo - 1, min(hi, len(g))) if g.target_id[i] in clean_t and g.target_id[i] not in used_t and alone(g, i)
              and g["rank"][i] == g.targets_scoring_higher_plus_one[i]]
        return g.iloc[rng.choice(ok)] if ok else None

    for lo, hi in ((2, 10), (11, 50), (51, 100), (101, 180)):
        k = traits(lambda g: a_candidate(g, lo, hi) is not None, 1)[0]
        r = a_candidate(by[k], lo, hi)
        used_t.add(r.target_id)
        add("rank", f"{r.symbol} has no link to {dname[k]}. By this graph's ranking, at what rank is {r.symbol} among the candidates for {dname[k]}?" + _final2('"rank": N'),
            {"problem": None, "rank": int(r["rank"])})
    for lo, hi in ((2, 10), (11, 50), (51, 100), (101, 180)):
        k = traits(lambda g: a_candidate(g, lo, hi) is not None, 1)[0]
        r = a_candidate(by[k], lo, hi)
        used_t.add(r.target_id)
        add("held", f"{r.symbol} has no link to {dname[k]}. Of candidates in its group at the ranks where it stands, what percentage were reported in the test on later links?"
            + _final2('"percent": NUMBER'), {"problem": None, "percent": float(rates[(bin_of(int(r["rank"])), group_of(r))])})

    for want in (1, 1, 0, 0):
        k = traits(lambda g: lead(g) and g.explained[0] == want, 1)[0]
        add("explained", f"Is the highest-ranked candidate target for {dname[k]} explained by an earlier link?" + _final2('"explained": true or false'),
            {"problem": None, "explained": bool(want)})
    for k in traits(lambda g: lead(g) and g.genes_nearby_with_this_link[0] >= 1, 2) + traits(lead, 1):
        add("nearby", f"For the highest-ranked candidate target for {dname[k]}, how many genes within 500,000 bases of it already have the link?" + _final2('"genes_nearby": N'),
            {"problem": None, "genes_nearby": int(by[k].genes_nearby_with_this_link[0])})

    def first_unexplained(g):
        i = next((i for i in range(len(g)) if not g.explained[i]), None)
        return i if i is not None and i >= 3 and g.target_id[i] in clean_t and alone(g, i) else None

    for k in traits(lambda g: first_unexplained(g) is not None, 3):
        add("unexplained_top", f"Leaving out the candidates that an earlier link already explains, which target is the highest-ranked candidate for {dname[k]}?"
            + _final2('"symbol": GENE SYMBOL'), {"problem": None, "symbol": by[k].symbol[first_unexplained(by[k])]})

    nt = {k: g.reset_index(drop=True) for k, g in NT[NT.trait_id.isin(clean_d)].groupby("trait_id", sort=True)}
    strict = lambda g: len(g) >= 3 and g.target_id[0] in clean_t and (g.genes_nearby_with_this_link[0], g.links_held_by_genes_nearby[0]) != (g.genes_nearby_with_this_link[1], g.links_held_by_genes_nearby[1])
    pool = sorted(k for k, g in nt.items() if k not in used_d and strict(g))
    for k in rng.sample(pool, 3):
        used_d.add(k)
        add("no_link_top", f"Among the targets that have no link at all, which is ranked highest for {dname[k]}?" + _final2('"symbol": GENE SYMBOL'),
            {"problem": None, "symbol": nt[k].symbol[0]})

    nd = {k: g.reset_index(drop=True) for k, g in ND[ND.target_id.isin(clean_t)].groupby("target_id", sort=True)}
    first = lambda g: len(g) >= 3 and (g.genes_nearby_with_this_link[0], g.targets_linked_to_this_trait[0]) != (g.genes_nearby_with_this_link[1], g.targets_linked_to_this_trait[1])
    pool = sorted(k for k, g in nd.items() if first(g))
    for k in rng.sample(pool, 3):
        used_t.add(k)
        g = nd[k]
        add("no_link_trait", f"{g.symbol[0]} has no link in this graph. Which trait is ranked highest for it?" + _final2('"trait": TRAIT NAME OR ID, or null if there is nothing to rank'),
            {"problem": None, "trait_any_of": [g.trait_id[0], g.trait_name[0]]})

    bare = R["targets_with_no_link"]
    bare = bare[bare.id.isin(clean_t) & ~bare.id.isin(set(ND.target_id))].sort_values("id")
    for i in rng.sample(range(len(bare)), 2):
        add("nothing", f"{bare.name.iloc[i]} has no link in this graph. Which trait is ranked highest for it?" + _final2('"trait": TRAIT NAME OR ID, or null if there is nothing to rank'),
            {"problem": None, "trait": None})

    for fake in fakes:
        if G.resolve(fake).get("problem") != "not in the graph":
            raise SystemExit(f"REFUSED: the invented symbol {fake} matches something in the graph")
        add("not_in_graph", f"Which trait is ranked highest for {fake}?" + _final2('"trait": TRAIT NAME OR ID, or null'), {"problem": "not_in_graph"})
    for word in words:
        if G.resolve(word).get("problem") != "ambiguous":
            raise SystemExit(f"REFUSED: {word!r} is not ambiguous in this graph")
        add("ambiguous", f"Which target is the highest-ranked candidate for {word}?" + _final2('"symbol": GENE SYMBOL or null'), {"problem": "ambiguous"})
    if len(qs) != 35:
        raise SystemExit(f"REFUSED: {len(qs)} questions were made, not 35")
    return qs


def via_ask_rank(q: dict, G) -> dict:
    """What the question layer itself would put in FINAL for a rankings question. Used only to cross-check
    the expected answers, which are computed without it: where the two disagree, one of them is wrong."""
    text, fam = q["text"].split(" End with")[0], q["family"]
    prob = lambda a: {"problem": {"not in the graph": "not_in_graph", "ambiguous": "ambiguous"}.get(a.get("problem"), a.get("problem"))}
    if fam in ("rank", "held"):
        sym, trait = re.match(r"^(\S+) has no link to (.+?)\. ", text).groups()
        w = G.why(sym, trait)
        if fam == "rank":
            return {"problem": None, "rank": w["rank_among_candidates"]}
        return {"problem": None, "percent": re.search(r", (\d+)% of candidates", w["calibration"]).group(1)}
    if fam in ("no_link_trait", "nothing", "not_in_graph"):
        a = G.predict(re.search(r"^(\S+) has no link|ranked highest for (\S+)\?", text).group(1) or re.search(r"ranked highest for (\S+)\?", text).group(1), top=1)
        if not a["found"]:
            return prob(a)
        return {"problem": None, "trait": a["candidates"][0]["id"] if a["candidates"] else None}
    trait = re.search(r"(?:candidate target for|candidate for|ranked highest for) (.+?)(?: explained by an earlier link\?|, how many genes|\?)", text).group(1)
    a = G.predict(trait, top=1, unexplained_only=fam == "unexplained_top", unlinked=fam == "no_link_top")
    if not a["found"]:
        return prob(a)
    c = a["candidates"][0]
    if fam == "explained":
        return {"problem": None, "explained": bool(c["explained_by"])}
    if fam == "nearby":
        return {"problem": None, "genes_nearby": c["score_parts"]["genes_nearby_with_this_link"]}
    return {"problem": None, "symbol": c["symbol"]}


def make_new_questions(bundle: Path, seed: int, fakes: tuple, words: tuple) -> list[dict]:
    """35 questions about the other sources and the gaps of a disease, in thirteen kinds. Expected
    answers are read from `gaptables`, which computes them from the bundle's tables without the
    question layer; `ask` is used only to choose names it resolves without a note, and names it finds
    ambiguous or absent."""
    G = ask.Graph(bundle)
    rng = random.Random(seed)
    T = gaptables.read(bundle)
    R = gaptables.build(T)
    S, summ = R["same_condition"], R["disease_gap_summary"].set_index("condition_id")
    nodes = T["nodes"]
    t, d = nodes[nodes.kind == "target"], nodes[nodes.kind == "trait"]
    cn = T["context_nodes"]
    n_sym = t.groupby("name").id.transform("size")
    plain = lambda r: r["ok"] and not r["note"]
    clean_t = {r.id for r in t[n_sym == 1].itertuples() if (lambda x: plain(x) and x["kind"] == "target" and G.t.id[x["i"]] == r.id)(G.resolve(r.name))}
    clean_d = {r.id for r in d[d.disease_like].itertuples() if (lambda x: plain(x) and x["kind"] == "trait" and G.d.id[x["i"]] == r.id)(G.resolve(r.name))}
    cond = cn[cn.kind == "condition"]
    one_name = cond.groupby(cond.name.map(gaptables.norm)).id.transform("size") == 1
    clean_c = {r.id for r in cond[one_name].itertuples() if r.id in summ.index and (lambda x: plain(x) and x["kind"] == "condition" and x["key"] == r.id)(G.resolve(r.name))}
    drugs_ = cn[cn.kind == "drug"]
    clean_g = {r.id for r in drugs_.itertuples() if (lambda x: plain(x) and x["kind"] == "drug" and x["key"] == r.id)(G.resolve(r.name))}
    sym = dict(zip(t.id, t.name))
    name = dict(zip(d.id, d.name)) | dict(zip(cn.id, cn.name))
    asked_of = S[S.asked_id.isin(clean_d | clean_c)].groupby("read_id").asked_id.apply(list).to_dict()       # an id a row sits under -> the clean things read under it
    used_c, used_t, used_g = set(), set(), set()
    qs: list[dict] = []

    def add(fam: str, text: str, gold: dict) -> None:
        qs.append({"id": f"{fam}-{sum(1 for q in qs if q['family'] == fam) + 1}", "family": fam, "text": text, "gold": gold})

    def pairs(rows: pd.DataFrame, keep, n: int) -> list[tuple[str, str, dict]]:
        """n (target, asked condition, what gaptables says of the pair), from rows of target and condition, each target and condition used once."""
        got = []
        for r in rows.sample(frac=1.0, random_state=rng.randrange(2 ** 31)).itertuples():
            if r.target not in clean_t or r.target in used_t:
                continue
            for c in asked_of.get(r.condition, []):
                if c in used_c:
                    continue
                p = gaptables.pair(T, S, r.target, c)
                if keep(p):
                    used_t.add(r.target), used_c.add(c)
                    got.append((r.target, c, p))
                    break
            if len(got) == n:
                return got
        raise SystemExit(f"REFUSED: fewer than {n} pairs of the kind wanted")

    A, X, DT, DI = T["assertions"], T["expression"], T["drug_targets"], T["drug_indications"]
    kept = A[~A.withheld]
    many = kept.groupby(["target", "condition"]).source.nunique().reset_index()
    how_many = lambda s_, c: f"How many of the disease sources assert {s_} with {c}?" + _final2('"sources": N')
    for tg_, c, p in pairs(many[many.source >= 2], lambda p: p["sources"] >= 2, 3):
        add("sources", how_many(sym[tg_], name[c]), {"problem": None, "sources": p["sources"]})
    tg_, c, p = pairs(X, lambda p: p["sources"] == 0 and p["withheld_rows"] == 0 and p["observations"] >= 1, 1)[0]
    add("sources", how_many(sym[tg_], name[c]), {"problem": None, "sources": 0})

    was_held = lambda s_, c: (f"Has a pruning rule withheld a row for {s_} and {c}? If one has, how many papers mention both?" + _final2('"withheld": true or false, "papers": N or null'))
    lit = A[A.withheld & A.rule.str.contains("literature_lift")]
    for tg_, c, p in pairs(lit, lambda p: p["withheld_rows"] == 1 and len(p["papers"]) == 1, 2):
        add("withheld", was_held(sym[tg_], name[c]), {"problem": None, "withheld": True, "papers": p["papers"][0]})
    tg_, c, p = pairs(many, lambda p: p["sources"] >= 1 and p["withheld_rows"] == 0, 1)[0]
    add("withheld", was_held(sym[tg_], name[c]), {"problem": None, "withheld": False, "papers": None})

    a_drug = lambda s_, c: (f"Does any drug act on {s_} and treat {c}? If one does, is {s_} the mechanism target of at least one such drug?" + _final2('"drug": true or false, "mechanism": true, false or null'))
    path = DT.merge(DI, on="drug").rename(columns={"condition": "condition"})
    for want in (True, False):
        tg_, c, p = pairs(path[path.mechanism == want][["target", "condition"]], lambda p: p["drugs"] >= 1 and p["mechanism"] == want, 1)[0]
        add("drug_path", a_drug(sym[tg_], name[c]), {"problem": None, "drug": True, "mechanism": want})
    tg_, c, p = pairs(many[many.source >= 2], lambda p: p["drugs"] == 0, 1)[0]
    add("drug_path", a_drug(sym[tg_], name[c]), {"problem": None, "drug": False, "mechanism": None})

    facts = {g: gaptables.drug(T, g) for g in sorted(clean_g)}
    both = sorted(g for g, f in facts.items() if 2 <= f["targets"] <= 12 and 1 <= f["mechanism_targets"] < f["targets"])
    acts_on = lambda g_: f"How many targets does the drug {g_} act on, and for how many of them is that its mechanism of action?" + _final2('"targets": N, "mechanism_targets": N')
    for g in rng.sample(both, 3):
        used_g.add(g)
        add("drug_targets", acts_on(name[g]), {"problem": None, "targets": facts[g]["targets"], "mechanism_targets": facts[g]["mechanism_targets"]})

    for tg_, c, p in pairs(X, lambda p: p["observations"] == 1, 2):
        add("observed", f"The graph holds one expression observation for {sym[tg_]} in {name[c]}. Was its expression higher or lower?" + _final2('"direction": "higher" or "lower"'),
            {"problem": None, "direction": p["directions"][0]})

    def diseases(keep, n: int, pool=clean_d) -> list[str]:
        ok = sorted(c for c in pool if c not in used_c and c in summ.index and keep(summ.loc[c]))
        got = rng.sample(ok, n)
        used_c.update(got)
        return got

    no_compound = lambda c: f"Of the genes tied to {c}, how many have no drug acting on them and no compound at their family's potency cutoff?" + _final2('"genes": N')
    still_dark = lambda c: f"Of the genes tied to {c}, how many are still Tdark?" + _final2('"genes": N')
    n_loci = lambda c: f"How many loci of {c} carry the GWAS link on more than one gene?" + _final2('"loci": N or null')
    for c in diseases(lambda r: r.no_compound >= 5 and r.no_compound < r.genes_tied, 2):
        add("no_compound", no_compound(name[c]), {"problem": None, "genes": int(summ.loc[c].no_compound)})
    for c in diseases(lambda r: r.dark >= 2, 2):
        add("dark", still_dark(name[c]), {"problem": None, "genes": int(summ.loc[c].dark)})
    for c in diseases(lambda r: r.loci_not_resolved >= 3, 3):
        add("loci", n_loci(name[c]), {"problem": None, "loci": int(summ.loc[c].loci_not_resolved)})
    for c in diseases(lambda r: r.no_genetic_support >= 3, 3):
        add("no_genetics", f"How many genes are tied to {name[c]} with no GWAS link?" + _final2('"genes": N'), {"problem": None, "genes": int(summ.loc[c].no_genetic_support)})
    for c in diseases(lambda r: pd.notna(r.smallest_odds_ratio_it_could_find), 3):
        add("floor", f"For the largest study behind the GWAS links of {name[c]}, what odds ratio did a variant carried by 20% of chromosomes need in order to be found?"
            + _final2('"odds_ratio": NUMBER'), {"problem": None, "odds_ratio": float(summ.loc[c].smallest_odds_ratio_it_could_find)})
    for c in diseases(lambda r: r.is_a_gwas_trait == 0 and r.genes_tied >= 5, 2, pool=clean_c):
        add("not_gwas", n_loci(name[c]), {"problem": None, "loci": None})

    for f in fakes:
        if G.resolve(f)["ok"] or G.resolve(f)["problem"] != "not in the graph":
            raise SystemExit(f"REFUSED: {f!r} is in the graph, or looks like something in it")
    c = diseases(lambda r: r.genes_tied >= 5, 1)[0]
    add("not_in_graph", how_many(fakes[0], name[c]), {"problem": "not_in_graph"})
    add("not_in_graph", acts_on(fakes[1]), {"problem": "not_in_graph"})
    add("not_in_graph", no_compound(fakes[2]), {"problem": "not_in_graph"})
    for w in words:
        if (lambda r: r["ok"] or r["problem"] != "ambiguous")(G.resolve(w, "trait")):
            raise SystemExit(f"REFUSED: {w!r} does not mean several conditions in this graph")
        add("ambiguous", still_dark(w), {"problem": "ambiguous"})
    if len(qs) != 35:
        raise SystemExit(f"REFUSED: {len(qs)} questions were written, 35 are wanted")
    return qs


def via_ask_new(q: dict, G) -> dict:
    """What the question layer itself would put in FINAL for one of these questions. Used only to
    cross-check the expected answers, which are computed without it."""
    text, fam = q["text"].split(" End with")[0], q["family"]
    prob = lambda a: {"problem": {"not in the graph": "not_in_graph", "ambiguous": "ambiguous"}.get(a.get("problem"), a.get("problem"))}
    m = (re.match(r"^How many of the disease sources assert (\S+) with (.+)\?$", text) or re.match(r"^Has a pruning rule withheld a row for (\S+) and (.+?)\? If one has", text)
         or re.match(r"^Does any drug act on (\S+) and treat (.+?)\? If one does", text) or re.match(r"^The graph holds one expression observation for (\S+) in (.+?)\. Was", text))
    if m:
        w = G.why(m.group(1), m.group(2))
        if "other_sources" not in w:
            return prob(w)
        o = w["other_sources"]
        if text.startswith("How many"):
            return {"problem": None, "sources": o["sources_asserting"]}
        if text.startswith("Has a pruning"):
            return {"problem": None, "withheld": o["withheld_rows"] > 0, "papers": o["withheld"][0].get("papers_mentioning_both") if o["withheld"] else None}
        if text.startswith("Does any drug"):
            return {"problem": None, "drug": o["drug_paths"] > 0, "mechanism": any(p["target_is_the_drugs_mechanism"] for p in o["through_a_drug"]) if o["drug_paths"] else None}
        return {"problem": None, "direction": o["expression_observations"][0]["direction"]}
    m = re.match(r"^How many targets does the drug (.+?) act on, and", text)
    if m:
        k = G.known(m.group(1), limit=50)
        if not k["found"]:
            return prob(k)
        return {"problem": None, "targets": len(k["targets"]) + k["targets_more"], "mechanism_targets": sum(x["is_the_drugs_mechanism"] for x in k["targets"])}
    c = (re.match(r"^Of the genes tied to (.+?), how many", text) or re.match(r"^How many loci of (.+?) carry the GWAS link", text) or re.match(r"^How many genes are tied to (.+?) with no GWAS link\?$", text)
         or re.match(r"^For the largest study behind the GWAS links of (.+?), what odds ratio", text)).group(1)
    o = G.gaps(disease=c)
    if not o["found"]:
        return prob(o)
    n = o["counts"]
    if fam == "floor":
        f = next(g for g in o["gaps"] if g["kind"] == "no_genetic_support")["what_an_absent_link_could_have_seen"]["largest"]
        return {"problem": None, "odds_ratio": f["smallest_odds_ratio_it_could_find"]}
    if fam in ("loci", "not_gwas"):
        return {"problem": None, "loci": n.get("locus_not_resolved")}
    return {"problem": None, "genes": n["no_compound" if "no compound" in text else "dark" if "Tdark" in text else "no_genetic_support"]}


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
        elif k in ("has_gap", "explained") or (k in ("withheld", "drug", "mechanism") and want is not None):
            if final.get(k) not in (want, str(want).lower()):
                return False
        elif k == "direction":
            if str(final.get(k) or "").strip().casefold() != want:
                return False
        elif k == "odds_ratio":
            try:
                if abs(float(str(final.get(k)).strip()) - want) > 0.011:
                    return False
            except (TypeError, ValueError):
                return False
        elif k == "trait_any_of":
            got = str(final.get("trait") or "").strip().casefold()
            if got not in {str(x).casefold() for x in want}:
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
def sql_db(bundle: Path, rankings: bool = False, others_and_gaps: bool = False) -> sqlite3.Connection:
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
    if rankings:
        R = ranktables.build(bundle)
        R["candidates"].drop(columns="targets_scoring_higher_plus_one").to_sql("candidates", con, index=False)
        R["rank_rates"].to_sql("rank_rates", con, index=False)
        R["no_link_targets"].drop(columns="links_held_by_genes_nearby").to_sql("no_link_targets", con, index=False)
        R["no_link_traits"].drop(columns="targets_linked_to_this_trait").to_sql("no_link_traits", con, index=False)
        for stmt in ("CREATE INDEX r1 ON candidates(trait_id, rank)", "CREATE INDEX r2 ON candidates(symbol)", "CREATE INDEX r3 ON no_link_targets(trait_id, rank)",
                     "CREATE INDEX r4 ON no_link_traits(symbol)"):
            con.execute(stmt)
    if others_and_gaps:
        X = gaptables.read(bundle)
        R = gaptables.build(X)
        cn = X["context_nodes"]
        d_ids = set(n[n.kind == "trait"].id)
        askable = set(R["same_condition"].asked_id) - d_ids                 # a condition read as part of a GWAS trait is that trait, and is not listed a second time
        cn[(cn.kind == "condition") & cn.id.isin(askable)][["id", "name"]].to_sql("conditions", con, index=False)
        cn[cn.kind == "drug"][["id", "name"]].to_sql("drugs", con, index=False)
        R["same_condition"].to_sql("same_condition", con, index=False)
        a = X["assertions"]
        pd.DataFrame({"target_id": a.target, "condition_id": a.condition, "source": a.source, "withheld": a.withheld.astype(int), "rule": a.rule,
                      "papers_mentioning_both": a.papers}).to_sql("assertions", con, index=False)
        x = X["expression"]
        pd.DataFrame({"target_id": x.target, "condition_id": x.condition, "condition_as_named": x.label, "log2_fold_change": x.log2fc, "p": x.p}).to_sql("expression", con, index=False)
        dt = X["drug_targets"]
        pd.DataFrame({"drug_id": dt.drug, "target_id": dt.target, "action": dt.action, "is_mechanism": dt.mechanism.astype(int)}).to_sql("drug_targets", con, index=False)
        X["drug_indications"].rename(columns={"drug": "drug_id", "condition": "condition_id"})[["drug_id", "condition_id"]].to_sql("drug_indications", con, index=False)
        for k in ("disease_gap_summary", "disease_gap_genes", "disease_loci"):
            R[k].to_sql(k, con, index=False)
        for stmt in ("CREATE INDEX n1 ON same_condition(asked_id)", "CREATE INDEX n2 ON assertions(target_id)", "CREATE INDEX n3 ON assertions(condition_id)", "CREATE INDEX n4 ON expression(target_id)",
                     "CREATE INDEX n5 ON drug_targets(drug_id)", "CREATE INDEX n6 ON drug_targets(target_id)", "CREATE INDEX n7 ON drug_indications(condition_id)",
                     "CREATE INDEX n8 ON disease_gap_summary(condition_id)", "CREATE INDEX n9 ON disease_gap_genes(condition_id, list)", "CREATE INDEX n10 ON conditions(name)", "CREATE INDEX n11 ON drugs(name)"):
            con.execute(stmt)
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
    if schema in CAPS:      # a later test has one cap over its two sets, whatever a single run is allowed
        gone = sum(r.get("cost") or 0.0 for r in (json.loads(p.read_text()) for p in out_dir.glob("*.json")) if SETS.get(r.get("set", "dev"), {}).get("schema") == schema)
        max_cost = min(max_cost, CAPS[schema] - gone)
        if max_cost <= 0:
            raise SystemExit(f"REFUSED: this test has spent {gone:.2f} of its {CAPS[schema]:.2f} dollars")
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
        con = sql_db(bundle, rankings=schema == 4, others_and_gaps=schema == 5)

        async def call_tool(name: str, args: dict) -> str:
            return run_sql(con, args.get("query", "")) if name == "sql" else f"unknown tool {name}; the only tool is sql"

        system = rules + "\n\n" + SQL_SCHEMA + SQL_EXTRA.get(schema, "")
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
        conr = sql_db(Path(tmp), rankings=True)
        got = json.loads(run_sql(conr, "select symbol, rank, explained from candidates where trait_id = 'MONDO:2' order by rank"))["rows"]
        ck("for the rankings test the SQL tool also holds the candidates, the rates and the two lists for targets with no link",
           [r[0] for r in got] == ["BBB", "AAA"] and json.loads(run_sql(conr, "select symbol, trait_id from no_link_targets order by trait_id"))["rows"] == [["DDD", "EFO:7"], ["DDD", "EFO:9"], ["DDD", "MONDO:2"]]
           and json.loads(run_sql(conr, "select count(*) from rank_rates"))["rows"] == [[12]]
           and json.loads(run_sql(conr, "select trait_name from no_link_traits where symbol = 'DDD' and rank = 1"))["rows"] == [["body height"]])
        G = ask.Graph(Path(tmp))
        p = G.predict("arrhythmia")
        ck("those tables, computed without the question layer, give the candidates the question layer gives",
           [r[0] for r in got] == [c["symbol"] for c in p["candidates"]] and [bool(r[2]) for r in got] == [bool(c["explained_by"]) for c in p["candidates"]])
    q = {"family": "explained", "gold": {"problem": None, "explained": True}}
    ck("a yes or no is right as a boolean or as text", grade(q, {"problem": None, "explained": True}) and grade(q, {"problem": None, "explained": "true"}) and not grade(q, {"problem": None, "explained": False}))
    q = {"family": "no_link_trait", "gold": {"problem": None, "trait_any_of": ["EFO:7", "body height"]}}
    ck("a trait is right by its id or its name, in any letter case", grade(q, {"problem": None, "trait": "efo:7"}) and grade(q, {"problem": None, "trait": "Body Height"}) and not grade(q, {"problem": None, "trait": "height"}))
    q = {"family": "nothing", "gold": {"problem": None, "trait": None}}
    ck("nothing to rank is right as a null with no problem, and wrong as a problem or a guess",
       grade(q, {"problem": None, "trait": None}) and not grade(q, {"problem": "unknown", "trait": None}) and not grade(q, {"problem": None, "trait": "body height"}))
    ck("each later test has its own cap, and the rankings test its own rules and tables", CAPS == {3: 0.75, 4: 0.60, 5: 0.60} and "ranked apart" in RULES_BY_SCHEMA[4] and "no_link_traits" in SQL_EXTRA[4]
       and SETS["rank"]["seed"] != SETS["rank_dev"]["seed"] and not set(SETS["rank"]["fakes"]) & set(SETS["rank_dev"]["fakes"]))
    # the other sources and the gaps of a disease
    ck("the newest test has its own rules and tables, and its two sets share no seed, invented name or word",
       "is not an assertion" in RULES_BY_SCHEMA[5] and "disease_gap_summary" in SQL_EXTRA[5] and "target_gaps" in SQL_EXTRA[5] and SETS["new"]["seed"] != SETS["new_dev"]["seed"]
       and not set(SETS["new"]["fakes"]) & set(SETS["new_dev"]["fakes"]) and not set(SETS["new"]["words"]) & set(SETS["new_dev"]["words"]))
    q = {"family": "withheld", "gold": {"problem": None, "withheld": True, "papers": 3}}
    ck("a withheld row is right with its count of papers, and wrong without it or with another", grade(q, {"problem": None, "withheld": True, "papers": 3}) and grade(q, {"problem": None, "withheld": "true", "papers": "3"})
       and not grade(q, {"problem": None, "withheld": True, "papers": None}) and not grade(q, {"problem": None, "withheld": False, "papers": 3}))
    q = {"family": "drug_path", "gold": {"problem": None, "drug": False, "mechanism": None}}
    ck("no drug between them is right with the mechanism left null, and wrong with one claimed", grade(q, {"problem": None, "drug": False, "mechanism": None}) and not grade(q, {"problem": None, "drug": False, "mechanism": False})
       and not grade(q, {"problem": None, "drug": True, "mechanism": None}))
    q = {"family": "observed", "gold": {"problem": None, "direction": "lower"}}
    ck("a direction is right in any letter case", grade(q, {"problem": None, "direction": "Lower"}) and not grade(q, {"problem": None, "direction": "higher"}) and not grade(q, {"problem": None}))
    q = {"family": "floor", "gold": {"problem": None, "odds_ratio": 1.084}}
    ck("an odds ratio is right within 0.011", grade(q, {"problem": None, "odds_ratio": 1.08}) and grade(q, {"problem": None, "odds_ratio": "1.09"}) and not grade(q, {"problem": None, "odds_ratio": 1.10})
       and not grade(q, {"problem": None, "odds_ratio": None}))
    q = {"family": "not_gwas", "gold": {"problem": None, "loci": None}}
    ck("a condition with no dated link is right as a null with no problem, and wrong as a count or a problem", grade(q, {"problem": None, "loci": None}) and not grade(q, {"problem": None, "loci": 0})
       and not grade(q, {"problem": "not_in_graph", "loci": None}))
    q = {"family": "sources", "gold": {"problem": None, "sources": 0}}
    ck("no source asserting a pair is right as zero, and wrong as one", grade(q, {"problem": None, "sources": 0}) and not grade(q, {"problem": None, "sources": 1}) and not grade(q, {"problem": None, "sources": None}))
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
        qs = (make_new_questions(paths.bundle(a.release), S["seed"], S["fakes"], S["words"]) if S["schema"] == 5
              else make_rank_questions(paths.bundle(a.release), S["seed"], S["fakes"], S["words"]) if S["schema"] == 4
              else make_gap_questions(paths.bundle(a.release), S["seed"], S["fakes"], S["words"]) if S["schema"] == 3
              else make_questions(paths.bundle(a.release), S["seed"], S["schema"], S["fakes"], S["words"]))
        d = {"release": a.release, "seed": S["seed"], "sha256": digest(qs), "questions": qs}
        if S["schema"] >= 3:
            G = ask.Graph(paths.bundle(a.release))
            differ = [q["id"] for q in qs if not grade(q, {5: via_ask_new, 4: via_ask_rank}.get(S["schema"], via_ask)(q, G))]
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
