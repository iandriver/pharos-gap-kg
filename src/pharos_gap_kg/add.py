"""Adding a table: a model proposes how a scientist's table joins the graph, code checks it, the person joins it.

A scientist brings a file from a paper the release does not hold (CSV, TSV, XLSX or ODS) and a
description in a fixed form. In two stages a model proposes, in a fixed JSON form, first which sheet
and which row of column names to read, then how the rows become Biolink edges. Code accepts or refuses
each proposal with its reasons, and a refused proposal goes back to the model once. The mapping has the
form of `ingest`, with a resolver for gene symbols, a separator for a cell that holds several values,
and a constant object chosen from the graph's traits; the validator is `ingest`'s nine checks. An
accepted mapping is shown to the person, and it joins, as a dated source beside the release like
`later`, only when the person runs the join. No model joins anything (docs/PROTOCOL.md, "Adding a table").

Run:  python -I -m pharos_gap_kg.add fetch kim2023_pd_supp            # the registered file of study B
      python -I -m pharos_gap_kg.add control pharos319                # the hand mapping of study A through the same checks
      python -I -m pharos_gap_kg.add run pharos319 A                  # both stages with the model (needs OPENROUTER_API_KEY)
      python -I -m pharos_gap_kg.add replay pharos319 B               # check a run's logged replies again, calling no model
      python -I -m pharos_gap_kg.add rerun pharos319 B                # the logged replies again, and the model for a step the log lacks
      python -I -m pharos_gap_kg.add show pharos319 B --rerun         # a rerun's mapping (join takes --rerun too)
      python -I -m pharos_gap_kg.add show pharos319 B                 # the accepted mapping, field by field, for the person
      python -I -m pharos_gap_kg.add join pharos319 B                 # the person's join
      python -I -m pharos_gap_kg.add compare pharos319                # A against later.py's table, B against the GWAS Catalog
      python -m pharos_gap_kg.add replay-import BUNDLE RECEIPT SOURCE  # offline saved-recipe check, no join
      python -m pharos_gap_kg.add --selftest
"""
from __future__ import annotations

import csv
import io
import json
import os
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import date, datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import adstudy, ingest, paths, registry

SEED = 20261103
CAP_USD = 0.50
MODEL, FALLBACK = "anthropic/claude-haiku-5.5", "anthropic/claude-sonnet-5.5"
PREVIEW_ROWS, PREVIEW_CUT, MAX_TRAITS = 12, 40, 8
RESOLVERS = {
    "gene_symbol": "the value is a gene symbol, such as APOE; gives every gene of the graph with that symbol (its Pharos symbol, or the Ensembl symbol "
                   "of a gene the Pharos target is cross-referenced to). Graph genes are proteins here: use the category biolink:Protein",
    **{k: v for k, v in ingest.RESOLVERS.items() if not k.startswith("pharos_")},
}
SIDE_KEYS = ingest.SIDE_KEYS | {"separator"}

# The two studies of the protocol, with the descriptions it fixes word for word.
STUDIES = {
    "A": {"key": "uffelmann2026_supp", "short": "Uffelmann et al., Nature (2026)",
          "about": {"citation": "Uffelmann et al. Alzheimer's disease genetics converges on microglial and neuronal dysregulation", "journal": "Nature",
                    "published": "2026-10-28", "doi": "10.1038/s41586-026-11206-y", "licence": "CC BY-NC-ND 4.0",
                    "what": "A genome-wide association study of Alzheimer's disease", "disease": "Alzheimer's disease",
                    "add": "Add the genes the study names at each associated locus."}},
    "B": {"key": "kim2023_pd_supp", "short": "Kim et al., Nature Genetics (2023)",
          "about": {"citation": "Kim et al. Multi-ancestry genome-wide association meta-analysis of Parkinson's disease", "journal": "Nature Genetics",
                    "published": "2023-12-28", "doi": "10.1038/s41588-023-01584-8", "licence": "CC BY 4.0",
                    "what": "A genome-wide association study of Parkinson's disease", "disease": "Parkinson's disease",
                    "add": "Add the genes the study names at each associated locus."}},
}
STUDIES["C"] = {"key": "mcgrail2026_t1d_supp", "short": "McGrail et al., Nature Genetics (2026)",
                "about": {"citation": "McGrail et al. Genetic association and machine learning improve the prediction of type 1 diabetes risk",
                          "journal": "Nature Genetics", "published": "2026-04-30", "doi": "10.1038/s41588-026-02578-y", "licence": "CC BY 4.0",
                          "what": "A genome-wide association study of type 1 diabetes", "disease": "type 1 diabetes",
                          "add": "Add the genes the study names at each associated locus."}}
STUDIES["D"] = {"key": "cheng2026_re_supp", "short": "Cheng et al., Nature Genetics (2026)",
                "about": {"citation": "Cheng et al. Multi-ancestry genome-wide association analyses of refractive error augment genetic discovery and "
                                      "polygenic prediction", "journal": "Nature Genetics", "published": "2026-04-20", "doi": "10.1038/s41588-026-02576-0",
                          "licence": "CC BY-NC-ND 4.0", "what": "A genome-wide association study of refractive error", "disease": "refractive error",
                          "add": "Add the genes the study names at each associated locus."}}
CATALOG_STUDY = {"B": "GCST90308590", "C": "GCST90824163", "D": "GCST90841196"}   # how the GWAS Catalog release of 5 Oct 2026 curates studies B to D
ABSTRACT_LOCI = {"B": 78, "C": 160, "D": 932}   # the loci, signals or variants each abstract reports (the amendments of 4 Nov 2026)

# Check 8's first figure. A run's record names the one it was made under; a record that names none was made under "resolved".
CHECK8 = {"resolved": "a row counts when its sides give identifiers of any kind, and 95% of rows must (ingest's check 8)",
          "graph": "a row of the table counts when it gives an edge between a target and a trait of the graph, and 80% of each rule's rows must "
                   "(the second amendment of 4 Nov 2026)"}
RULE, BAR = "graph", 0.80

BRIEF1 = """A scientist gives you a file from a scientific paper and a description of the paper. You will be shown every sheet of the file: its name, its numbers of rows and columns, and its first rows, numbered as read (empty rows are skipped), with each cell cut to 40 characters.

Choose the one sheet that holds what the description asks to add, and the row that holds that sheet's column names.

Return ONE JSON object and nothing else: {"sheet": the sheet's name exactly as shown, "header_row": the number of the row that holds the column names, "notes": up to 300 characters}"""

BRIEF2 = """You map one table from a scientific paper into a knowledge graph built from Pharos (the NIH Illuminating the Druggable Genome database) that follows the Biolink Model, version 4.4.3. A scientist gave the table with a description of the paper. You will be given the description, a card that describes the table (each column with its counts and its values or examples, and eight example rows), and the traits of the graph whose names match the disease the description names.

Propose how each row becomes one or more edges (associations) between a subject node and an object node. Code will check your proposal against the Biolink Model and against the data, and may send it back once with the reasons it was refused. Code sets each edge's date (the year the paper was published) and its knowledge source (the paper); do not give them.

Return ONE JSON object and nothing else, with these keys:

  "map": true, or false if this table should not become edges (then give only "map" and "reason")
  "subject": {"columns": [column, ...], "resolver": one of the resolvers below, "category": a Biolink class CURIE such as "biolink:Gene", "name_column": optional column that names the node, "separator": optional text that separates several values in one cell, each of which becomes its own edge}
  "object":  the same keys, or {"constant": one trait id from the list of traits, "category": a Biolink class CURIE} when every row's object is that one trait
  "predicate": a Biolink predicate CURIE, such as "biolink:related_to"
  "association_category": a Biolink association class CURIE, such as "biolink:Association"
  "knowledge_level": one of KNOWLEDGE_LEVELS
  "agent_type": one of AGENT_TYPES
  "qualifiers": {} or {Biolink qualifier slot name: a constant value, or {"sign_of": column} for a direction qualifier}
  "edge_properties": [columns whose values the edge should carry]
  "keep_if": null, or {"column": c, "at_least": number, "reason": why a row below it is not asserted}, or {"column": c, "in": [values], "reason": ...}, or {"column": c, "matches": a regular expression, "reason": ...}. A row that fails it is kept as a recorded gap of its subject, not dropped, and needs no object.
  "split_by": null, or a column whose values need different mappings
  "by_value": only with split_by: {value: {any of subject, object, predicate, association_category, knowledge_level, agent_type, qualifiers, edge_properties, keep_if, given as overrides}, or {"skip": reason} to leave that value's rows out}
  "notes": up to 600 characters on the choices a reader might question

"columns" lists the columns to try in order; the first that gives an identifier is used.

RESOLVERS (code turns a cell into an identifier; you choose which):
RESOLVER_MENU

Rules:
- Use only Biolink classes, predicates and qualifiers that exist in Biolink 4.4.3. The subject and object categories must fit the predicate's domain and range and the association class.
- Choose the most specific predicate the source supports, and no more specific than that.
- The knowledge level and agent type describe how the paper made the claim.
- Do not invent columns or values."""


# ---- the description -------------------------------------------------------------------------------

def text(a: dict) -> str:
    """The description as the model reads it, from its fixed fields."""
    d = date.fromisoformat(a["published"])
    return f"{a['citation']}. {a['journal']}, published {d.day} {d.strftime('%B %Y')}. doi:{a['doi']}. Licence: {a['licence']}. {a['what']}. {a['add']}"


def infores(a: dict) -> str:
    """The knowledge source code gives the paper's edges."""
    return "infores:" + re.sub(r"[^a-z0-9]+", "-", f"doi {a['doi']}".casefold()).strip("-")


# ---- reading a file --------------------------------------------------------------------------------

_S = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def _col(ref: str) -> int:
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1


def xlsx_sheets(blob: bytes) -> dict[str, list[list[str]]]:
    """Every sheet of an Office Open XML workbook as rows of text; empty rows skipped, trailing empty cells dropped."""
    z = zipfile.ZipFile(io.BytesIO(blob))
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        shared = ["".join(t.text or "" for t in si.iter(f"{_S}t")) for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(f"{_S}si")]
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    out = {}
    for s in ET.fromstring(z.read("xl/workbook.xml")).iter(f"{_S}sheet"):
        t = rels[s.get(_RID)]
        rows = []
        for r in ET.fromstring(z.read(t.lstrip("/") if t.startswith("/") else "xl/" + t)).iter(f"{_S}row"):
            cells, last = {}, -1
            for c in r.iter(f"{_S}c"):
                j = _col(c.get("r")) if c.get("r") else last + 1
                last = j
                v, kind = c.find(f"{_S}v"), c.get("t")
                if kind == "s":
                    val = shared[int(v.text)] if v is not None and v.text else ""
                elif kind == "inlineStr":
                    val = "".join(x.text or "" for x in c.iter(f"{_S}t"))
                elif kind == "b":
                    val = "TRUE" if v is not None and v.text == "1" else "FALSE"
                else:
                    val = v.text if v is not None and v.text is not None else ""
                cells[j] = val
            row = [cells.get(j, "") for j in range(max(cells) + 1)] if cells else []
            while row and row[-1].strip() == "":
                row.pop()
            if row:
                rows.append(row)
        out[s.get("name")] = rows
    return out


def ods_sheets(blob: bytes) -> dict[str, list[list[str]]]:
    root = ET.fromstring(zipfile.ZipFile(io.BytesIO(blob)).read("content.xml"))
    names = [t.get(f"{{{adstudy._T}}}name") for t in root.iter(f"{{{adstudy._T}}}table")]
    return {n: [r for r in adstudy.sheet(blob, n) if any(c.strip() for c in r)] for n in names}


def read_sheets(path: Path) -> dict[str, list[list[str]]]:
    """A CSV or TSV (one sheet, named by the file), an XLSX or an ODS file, as sheets of rows of text."""
    path = Path(path)
    suf = path.suffix.casefold()
    if suf in (".csv", ".tsv", ".txt"):
        with open(path, newline="", encoding="utf-8-sig") as f:
            rows = [r for r in csv.reader(f, delimiter="\t" if suf != ".csv" else ",")]
        rows = [r[:max(i for i, c in enumerate(r) if c.strip()) + 1] for r in rows if any(c.strip() for c in r)]
        return {path.name: rows}
    blob = path.read_bytes()
    if suf == ".xlsx":
        return xlsx_sheets(blob)
    if suf == ".ods":
        return ods_sheets(blob)
    raise SystemExit(f"REFUSED: {path.name} is not a CSV, TSV, XLSX or ODS file")


# ---- stage 1: the sheet ----------------------------------------------------------------------------

def card1(sheets: dict) -> dict:
    return {"sheets": [{"name": n, "rows": len(rows), "columns": max((len(r) for r in rows), default=0),
                        "first_rows": [{"row": i + 1, "cells": [ingest.cut(c, PREVIEW_CUT) for c in r]} for i, r in enumerate(rows[:PREVIEW_ROWS])]}
                       for n, rows in sheets.items()]}


def frame(rows: list[list[str]], header_row) -> tuple[pd.DataFrame | None, list[str]]:
    """The table below the header row, or the reasons it cannot be read. Blank names after the last named
    column are dropped, and so is a column with no name and no value below the header (padding); an empty or
    repeated name among the others refuses the header row."""
    if not isinstance(header_row, int) or isinstance(header_row, bool) or not 1 <= header_row <= len(rows):
        return None, [f"header_row {header_row!r} is not a row of the sheet (1 to {len(rows)})"]
    head, body = rows[header_row - 1], rows[header_row:]
    names = [c.strip() for c in head]
    if not any(names):
        return None, [f"row {header_row} holds no column names"]
    last = max(i for i, c in enumerate(names) if c)
    has = lambda j: any(j < len(r) and r[j].strip() for r in body)
    keep, problems = [], []
    for j in range(last + 1):
        if names[j]:
            keep.append(j)
        elif has(j):
            problems.append(f"column {j + 1} of row {header_row} has no name and holds values")
    dup = sorted({names[j] for j in keep if [names[k] for k in keep].count(names[j]) > 1})
    if dup:
        problems.append(f"row {header_row} repeats the names {dup}")
    data = [[(r[j].strip() if j < len(r) and r[j].strip() else None) for j in keep] for r in body]
    data = [r for r in data if any(v is not None for v in r)]
    if not data:
        problems.append(f"no row below row {header_row} holds a value")
    if problems:
        return None, problems
    return pd.DataFrame(data, columns=[names[j] for j in keep]), []


def check1(p, sheets: dict) -> tuple[pd.DataFrame | None, list[str]]:
    if not isinstance(p, dict) or "sheet" not in p or "header_row" not in p:
        return None, ["the reply holds no JSON object with sheet and header_row"]
    if p["sheet"] not in sheets:
        return None, [f"the file has no sheet {p['sheet']!r}; its sheets are {list(sheets)}"]
    return frame(sheets[p["sheet"]], p["header_row"])


# ---- stage 2: the mapping --------------------------------------------------------------------------

def traits(G, disease: str) -> list[dict]:
    """The graph's traits that `find` returns for the disease in words, up to MAX_TRAITS."""
    ms = G.find(disease, "trait", 25).get("matches", [])
    return [{"id": m["id"], "name": m["name"], "links": m.get("links")} for m in ms if m["id"] in G.di][:MAX_TRAITS]


def card2(df: pd.DataFrame, sheet: str, header_row: int, tr: list[dict]) -> dict:
    rng = np.random.default_rng(SEED)
    cols = []
    for c in df.columns:
        s = df[c][ingest.present(df[c])]
        vc = s.map(ingest.cut).value_counts()
        col = {"name": c, "rows_with_a_value": int(len(s)), "distinct": int(len(vc))}
        if len(vc) <= 25:
            col["values"] = {k: int(v) for k, v in vc.items()}
        else:
            col["examples"] = [vc.index[i] for i in sorted(rng.choice(len(vc), size=6, replace=False))]
        cols.append(col)
    ex = df.sample(n=min(8, len(df)), random_state=SEED)
    return {"sheet": sheet, "header_row": header_row, "rows": int(len(df)), "columns": cols,
            "example_rows": [{c: (None if v is None else ingest.cut(v)) for c, v in r.items()} for r in ex.to_dict("records")], "traits": tr}


def table_profile(df: pd.DataFrame) -> dict:
    """Read-only, full-table counts before mapping. Values use the validator's
    missing-value convention; duplicate rows are reported, never removed."""
    columns = []
    for name in df.columns:
        values = df[name]
        present = ingest.present(values)
        nonmissing = values[present]
        columns.append({"name": name, "rows_with_a_value": int(present.sum()),
                        "missing_rows": int((~present).sum()),
                        "distinct_values": int(nonmissing.nunique())})
    return {"rows": len(df), "columns": columns,
            "duplicate_rows_beyond_first": int(df.duplicated().sum()),
            "reading": "Repeated rows are not removed: they may represent separate observations. "
                       "Missing values follow the mapping validator's convention. These counts do not establish scientific validity."}


def brief2() -> str:
    from kgx_gaps.mapping import AGENT_TYPES, KNOWLEDGE_LEVELS
    B = ingest.biolink()
    menu = "\n".join(f"  {k}: {v}" for k, v in RESOLVERS.items())
    return (BRIEF2.replace("KNOWLEDGE_LEVELS", ", ".join(sorted(set(B["enums"]["KnowledgeLevelEnum"]) & KNOWLEDGE_LEVELS)))
            .replace("AGENT_TYPES", ", ".join(sorted(set(B["enums"]["AgentTypeEnum"]) & AGENT_TYPES))).replace("RESOLVER_MENU", menu))


class Symbols:
    """Gene symbol -> the graph's target ids, by the matching rule of `adstudy` (the Pharos symbol, or the
    Ensembl 113 symbol of a gene the target is cross-referenced to; the second needs the release's tables)."""

    def __init__(self, G, release: str | None):
        try:
            m = adstudy.symbol_map(G, release) if release else None
        except SystemExit:
            m = None
        if m is None:
            m = {s: {s} for s in G.t["name"]}
        ids: dict[str, list[str]] = {}
        for t, s in zip(G.t.id, G.t["name"]):
            ids.setdefault(s, []).append(t)
        self.map, self.ids, self.symbol = m, ids, dict(zip(G.t.id, G.t["name"]))

    def __call__(self, v: str) -> list[str]:
        return [t for p in sorted(self.map.get(v.strip(), ())) for t in self.ids.get(p, ())]


def side_problems(name: str, s, cols: set[str], tr_ids: set[str]) -> list[str]:
    if not isinstance(s, dict):
        return [f"1_parses: {name} must be an object"]
    if "constant" in s:
        if name != "object" or set(s) - {"constant", "category"} or not s.get("category"):
            return [f"1_parses: a constant is given only for the object, with its category and nothing else"]
        return [] if s["constant"] in tr_ids else [f"1_parses: the constant {s['constant']!r} is not one of the traits on the card {sorted(tr_ids)}"]
    out = []
    if set(s) - SIDE_KEYS:
        out.append(f"1_parses: {name} has keys {sorted(set(s) - SIDE_KEYS)} that are not {sorted(SIDE_KEYS)} or a constant")
    sep = s.get("separator")
    if sep is not None and (not isinstance(sep, str) or not 1 <= len(sep) <= 5):
        out.append(f"1_parses: {name}'s separator must be text of 1 to 5 characters")
    r = s.get("resolver")
    if r != "gene_symbol" and (not isinstance(r, str) or r.startswith("pharos_") or ingest.resolver_ok(r)):
        out.append(f"1_parses: {name}: resolver {r!r} is not on the menu")
    return out


def _resolve_side(s: dict, rows: pd.DataFrame, sym: Symbols, refs) -> list[list[tuple[str, str]]]:
    """For each row, the (identifier, value) pairs of the first column that gives any; a separator splits a cell."""
    out = []
    for r in rows[s["columns"]].itertuples(index=False):
        got = []
        for v in r:
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue
            vals = [x.strip() for x in str(v).split(s["separator"])] if s.get("separator") else [str(v).strip()]
            for x in vals:
                if not x:
                    continue
                if s["resolver"] == "gene_symbol":
                    got += [(t, x) for t in sym(x)]
                else:
                    i, _ = ingest.resolve_one(x, s["resolver"], refs)
                    if i:
                        got.append((i, x))
            if got:
                break
        out.append(list(dict.fromkeys(got)))
    return out


def _needs_rows(s: dict) -> bool:
    return isinstance(s, dict) and "constant" not in s and (s.get("resolver") == "gene_symbol" or bool(s.get("separator")))


def transform(p: dict, df: pd.DataFrame, sym: Symbols, refs, ks: str) -> tuple[dict, pd.DataFrame]:
    """The proposal and table in `ingest`'s form: a constant object becomes a column, and a side resolved
    by gene symbol or split by a separator becomes one row per identifier (a row where nothing resolves
    stays once, unresolved, so check 8 still counts rows of the table)."""
    rules = list(ingest.rules_of(p))
    split = p.get("split_by")
    parts, p2 = [], {k: v for k, v in p.items() if k not in ("subject", "object", "by_value", "primary_knowledge_source")}
    p2["primary_knowledge_source"] = ks
    for value, r, skip in rules:
        mine = df if value is None else df[df[split].astype(str) == str(value)]
        if r is None:
            parts.append(mine)
            continue
        rows = mine.copy()
        rows["__row"] = rows.index
        spec = {}
        for side in ("subject", "object"):
            s = r[side]
            if "constant" in s:
                rows[f"__{side}"] = s["constant"]
                spec[side] = {"columns": [f"__{side}"], "resolver": "curie", "category": s["category"]}
            elif _needs_rows(s):
                got = _resolve_side(s, rows, sym, refs)
                rows = rows.assign(**{f"__{side}": [[g[0] for g in x] or [None] for x in got], f"__{side}_value": [[g[1] for g in x] or [None] for x in got]})
                rows = rows.explode([f"__{side}", f"__{side}_value"], ignore_index=True)
                spec[side] = {"columns": [f"__{side}"], "resolver": "curie", "category": s["category"]}
            else:
                spec[side] = {k: v for k, v in s.items() if k != "separator"}
        parts.append(rows)
        over = {k: v for k, v in r.items() if k not in ("subject", "object", "primary_knowledge_source")}
        if value is None:
            p2.update({**over, **spec})
        else:
            p2.setdefault("by_value", {})[value] = {**over, **spec, "primary_knowledge_source": ks}
    if split:
        for value, r, skip in rules:
            if r is None:
                p2.setdefault("by_value", {})[value] = {"skip": skip}
        base = next(r for _, r, _ in rules if r is not None)
        p2.setdefault("subject", p2["by_value"][next(v for v, r, _ in rules if r is not None)]["subject"])
        p2.setdefault("object", p2["by_value"][next(v for v, r, _ in rules if r is not None)]["object"])
        for k in ("predicate", "association_category", "knowledge_level", "agent_type"):
            p2.setdefault(k, base.get(k))
    rows2 = pd.concat(parts, ignore_index=True) if parts else df.copy()
    for c in ("__subject", "__object"):
        if c not in rows2:
            rows2[c] = None
    return p2, rows2


def graph_rows(G, bar: float = BAR):
    """Check 8's first figure under the rule "graph": a row of the table counts when one of its edges joins a target and a trait of
    the graph, or, for a row the condition makes a gap, when its subject is a node of the graph. A row whose cell held several values
    counts once (the rows `transform` made from one row share its `__row`)."""
    ti, di = G.ti, G.di

    def first(mine: pd.DataFrame, sub: list, obj: list, kept: list) -> dict:
        rows = mine["__row"].tolist() if "__row" in mine else list(mine.index)
        ok: dict = {}
        for r, a, b, k in zip(rows, sub, obj, kept):
            hit = ((a in ti and b in di) or (a in di and b in ti)) if k else (a in ti or a in di)
            ok[r] = ok.get(r, False) or bool(hit)
        return {"ok": int(sum(ok.values())), "rows": len(ok), "bar": bar,
                "what": "give an edge between a target and a trait of the graph (a row the condition makes a gap: a subject in the graph)"}
    return first


def check2(p, df: pd.DataFrame, tr: list[dict], sym: Symbols, refs, ks: str, view: str, G=None, rule: str = RULE) -> dict:
    """The validator's verdict on a mapping: the additions' checks, then `ingest`'s nine on the table as transformed, with check 8's
    first figure by `rule` (CHECK8; "graph" needs the graph G)."""
    if rule not in CHECK8:
        raise ValueError(f"check 8 has no rule {rule!r}")
    if not isinstance(p, dict):
        return {"accepted": False, "problems": ["1_parses: no JSON object"]}
    if p.get("map") is False:
        return {"accepted": True, "map": False, "problems": [], "reason": p.get("reason")}
    cols, tr_ids = set(df.columns), {t["id"] for t in tr}
    problems = [f"1_parses: unknown key {k!r}" for k in sorted(set(p) - ingest.KEYS)]
    if p.get("split_by") and p["split_by"] not in cols:
        problems.append(f"1_parses: split_by names no column {p['split_by']!r}")
    rules = [] if problems else list(ingest.rules_of(p))
    for value, r, skip in rules:
        if r is not None:
            for side in ("subject", "object"):
                problems += [x + ("" if value is None else f" [{value}]") for x in side_problems(side, r.get(side), cols, tr_ids)]
    if problems:
        return {"accepted": False, "problems": problems}
    p2, rows2 = transform(p, df, sym, refs, ks)
    v = ingest.validate(p2, rows2, refs, export=True, view=view, first=graph_rows(G) if rule == "graph" else None)
    v["_rows"], v["_p2"], v["check_8"] = rows2, p2, rule
    return v


def joinable(v: dict, G, sym: Symbols, a: dict, key: str) -> pd.DataFrame:
    """The rows a join writes: each exported edge between a graph target and a graph trait, with the year,
    the study, the symbols the table gave for the target, and the properties the mapping carries."""
    e = v["_frames"]["edges"]
    e = e.loc[:, ~e.columns.duplicated()]                  # the exporter writes a carried column twice when kgx-gaps also knows its name
    rows = v["_rows"]
    given: dict[str, list[str]] = {}
    if "__subject_value" in rows:
        for t, x in zip(rows["__subject"], rows["__subject_value"]):
            if isinstance(t, str) and isinstance(x, str):
                given.setdefault(t, []).append(x)
    pc = [c for c in e.columns if c.startswith(f"{ingest.PREFIX}:") and c not in (f"{ingest.PREFIX}:source_edge_type", ingest.WITHHELD)]
    out = []
    for r in e.to_dict("records"):
        t, d = (r["subject"], r["object"]) if r["subject"] in G.ti else (r["object"], r["subject"])
        if t not in G.ti or d not in G.di:
            continue
        props = {c.split(":", 1)[1]: r[c] for c in pc if isinstance(r.get(c), str) and r[c]}
        out.append({"target": t, "symbol": sym.symbol.get(t, t), "study_symbol": "|".join(sorted(set(given.get(t, [])))) or None, "trait": d,
                    "year": date.fromisoformat(a["published"]).year, "study": key, "predicate": r["predicate"],
                    "properties": json.dumps(props, sort_keys=True)})
    return pd.DataFrame(out, columns=["target", "symbol", "study_symbol", "trait", "year", "study", "predicate", "properties"])


# ---- the model -------------------------------------------------------------------------------------

def run_dir(release: str, study: str) -> Path:
    return paths.EVAL / "add" / release / study


def ledger_path(release: str) -> Path:
    return paths.EVAL / "add" / release / "ledger.jsonl"


def spent(release: str) -> float:
    """Costs recorded, plus the ceiling of any call reserved and never settled."""
    p = ledger_path(release)
    if not p.exists():
        return 0.0
    total, open_ = 0.0, {}
    for line in p.read_text().splitlines():
        r = json.loads(line)
        if r["kind"] == "reserve":
            open_[r["id"]] = float(r["usd_ceiling"])
        else:
            total += float(r["usd"])
            open_.pop(r["id"], None)
    return total + sum(open_.values())


def _ledger(release: str, rec: dict) -> None:
    p = ledger_path(release)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a") as f:
        f.write(json.dumps({"at": datetime.now(timezone.utc).isoformat(timespec="seconds"), **rec}) + "\n")


def call(release: str, model: str, messages: list[dict], tag: str) -> dict:
    key = os.environ.get("OPENROUTER_API_KEY")
    if not key:
        raise SystemExit("set OPENROUTER_API_KEY")
    pin, pout = ingest.PRICES[model]
    ceiling = (sum(len(m["content"]) for m in messages) / 3.0 * pin + ingest.MAX_TOKENS * pout) / 1e6
    if spent(release) + ceiling > CAP_USD:
        raise SystemExit(f"STOPPED: spend {spent(release):.4f} + this call's ceiling {ceiling:.4f} would pass the cap of {CAP_USD:.2f}")
    cid = f"{tag}:{time.time_ns()}"
    _ledger(release, {"kind": "reserve", "id": cid, "model": model, "usd_ceiling": round(ceiling, 6)})
    t0 = time.time()
    r = ingest._post("https://openrouter.ai/api/v1/chat/completions",
                     {"model": model, "messages": messages, "max_tokens": ingest.MAX_TOKENS, "temperature": 0, "usage": {"include": True}},
                     {"Authorization": f"Bearer {key}", "X-Title": "pharos-gap-kg add"}, timeout=600)
    u = r.get("usage") or {}
    _ledger(release, {"kind": "cost", "id": cid, "model": model, "usd": float(u.get("cost") or 0.0), "in": u.get("prompt_tokens"),
                      "out": u.get("completion_tokens"), "provider": r.get("provider"), "error": r.get("error")})
    if "error" in r or not r.get("choices"):
        return {"error": str(r.get("error"))[:300], "seconds": round(time.time() - t0, 1)}
    return {"text": r["choices"][0]["message"].get("content") or "", "usd": float(u.get("cost") or 0.0), "seconds": round(time.time() - t0, 1)}


def _public(v: dict) -> dict:
    return {k: x for k, x in v.items() if not k.startswith("_")}


def stage(release: str, study: str, name: str, system: str, user: str, check, replay: list | None = None, mode: str = "run",
          models: tuple = (MODEL, FALLBACK)) -> dict:
    """One stage: the model, one return of a refused proposal, then the fallback model from the start under the same rule.
    replay: the logged replies of an earlier run, checked again in place of calling a model; in mode 'rerun' a step the log
    does not hold calls the model and is logged."""
    log = run_dir(release, study) / "calls.jsonl"
    tries = []
    for model in models:
        msgs = [{"role": "system", "content": system}, {"role": "user", "content": user}]
        for attempt in (1, 2):
            old = next((x for x in replay if (x["stage"], x["model"], x["attempt"]) == (name, model, attempt)), None) if replay is not None else None
            if old is not None:
                r = {"text": old["reply"], "usd": old["usd"], "seconds": old["seconds"]} if old.get("reply") is not None else {"error": old["problems"][0]}
            elif replay is not None and mode != "rerun":
                return {"accepted": False, "tries": tries, "replay": "no logged reply for the next call"}
            else:
                r = call(release, model, msgs, f"{study}:{name}:{model}:{attempt}")
            if "error" in r:
                verdict = {"accepted": False, "problems": [f"the call failed: {r['error']}"]}
                p = None
            else:
                p = ingest.parse(r["text"])
                verdict = check(p)
            rec = {"stage": name, "model": model, "attempt": attempt, "proposal": p, "accepted": verdict["accepted"], "problems": verdict["problems"],
                   "usd": r.get("usd"), "seconds": r.get("seconds")}
            if old is None:
                with open(log, "a") as f:
                    f.write(json.dumps({**rec, "reply": r.get("text"), "at": datetime.now(timezone.utc).isoformat(timespec="seconds"), "check_8": RULE},
                                       default=str) + "\n")
            tries.append(rec)
            if verdict["accepted"]:
                return {"accepted": True, "model": model, "proposal": p, "verdict": verdict, "tries": tries}
            if "text" in r:
                msgs += [{"role": "assistant", "content": r["text"]}, {"role": "user", "content": ingest.refusal(verdict)}]
    return {"accepted": False, "tries": tries}


# ---- one study, end to end -------------------------------------------------------------------------

def _setup(release: str, study: str):
    from pharos_gap_kg import ask
    s = STUDIES[study]
    G = ask.Graph(paths.bundle(release))
    return s, G, Symbols(G, release), ingest.Refs(release), registry.require(s["key"])


def record_path(release: str, study: str, rerun: bool = False) -> Path:
    return paths.DOCS_DATA / f"add_{STUDIES[study]['key']}{'_rerun' if rerun else ''}_{release}.json"


def run(release: str, study: str, proposals: tuple | None = None, mode: str = "run", back: bool | None = None) -> dict:
    """Both stages for a registered study. proposals: (sheet proposal, mapping proposal) to check without a model (the control).
    mode: 'run' calls the model and logs each reply; 'replay' checks again the replies a run logged, calling no model;
    'rerun' checks the logged replies again and calls a model only for a step the log does not hold.
    back: whether a mapping that declines sends its reason back to the sheet stage once (the amendment of 4 Nov 2026);
    by default on for a rerun and for studies C and D, and off for the runs of A and B as they were registered."""
    if back is None:
        back = mode == "rerun" or study in ("C", "D")
    s, G, sym, refs, path = _setup(release, study)
    a, key = s["about"], s["key"]
    d = run_dir(release, study if proposals is None else f"{study}_control")
    d.mkdir(parents=True, exist_ok=True)
    sheets = read_sheets(path)
    c1 = card1(sheets)
    ingest.write_json(d / "card1.json", c1)
    user1 = f"DESCRIPTION\n{text(a)}\n\n" + json.dumps(c1, indent=1, ensure_ascii=False)
    logged = None
    if mode in ("replay", "rerun"):
        f = d / "calls.jsonl"
        logged = [json.loads(x) for x in f.read_text().splitlines()] if f.is_file() else []
        if mode == "replay" and not logged:
            raise SystemExit(f"REFUSED: no logged replies in {d}")
    tr, ks = traits(G, a["disease"]), infores(a)

    def chk1(p):
        _, problems = check1(p, sheets)
        return {"accepted": not problems, "problems": problems}

    def mapping(s1: dict, name: str, models: tuple) -> tuple:
        df = check1(s1["proposal"], sheets)[0]               # the frame of the accepted proposal
        c2 = card2(df, s1["proposal"]["sheet"], s1["proposal"]["header_row"], tr)
        ingest.write_json(d / f"card2{'' if name == 'mapping' else '_' + name}.json", c2)
        chk2 = lambda p: check2(p, df, tr, sym, refs, ks, key, G)
        if proposals is not None:
            v = chk2(proposals[1])
            s2 = {"accepted": v["accepted"], "model": "hand", "proposal": proposals[1], "verdict": v,
                  "tries": [{"stage": "mapping", "model": "hand", "accepted": v["accepted"], "problems": v["problems"]}]}
        else:
            s2 = stage(release, study, name, brief2(), f"DESCRIPTION\n{text(a)}\n\n" + json.dumps(c2, indent=1, ensure_ascii=False), chk2, logged, mode, models)
        return df, s2, chk2

    if proposals is None:
        s1 = stage(release, study, "sheet", BRIEF1, user1, chk1, logged, mode)
    else:
        v = chk1(proposals[0])
        s1 = {"accepted": v["accepted"], "model": "hand", "proposal": proposals[0], "verdict": v, "tries": [{"stage": "sheet", "model": "hand", **v}]}
    out = {"release": release, "study": study, "key": key, "description": text(a), "file": path.name,
           "file_sha256": registry.SOURCES[key].sha256, "sheets": len(sheets), "stage1": {k: s1[k] for k in ("accepted", "tries")},
           "mode": mode, "back_to_the_sheet_rule": back, "check_8": RULE}
    if mode == "replay":
        out["replayed"] = "the stages were checked again from the replies logged by the run; no model was called"
    if mode == "rerun":
        out["rerun"] = ("a rerun after the registered result, under the amendment" + (" of 4 Nov 2026" if study == "B" else "s of 4 Nov 2026")
                        + ": development, not a test")
    if not s1["accepted"]:
        return _record(release, study, out, proposals is not None, mode == "rerun")
    df, s2, chk2 = mapping(s1, "mapping", (MODEL, FALLBACK))
    out["stage2"] = {k: s2[k] for k in ("accepted", "tries")}
    declined = s2["accepted"] and s2["proposal"].get("map") is False
    if declined and back and proposals is None:
        reason = s2["proposal"].get("reason") or ""
        again = (user1 + f"\n\nYOUR CHOICE WAS DECLINED\nYou chose the sheet {s1['proposal']['sheet']!r} with its column names in row "
                 f"{s1['proposal']['header_row']}. Reading that table for the mapping, a model said it should not become edges: {reason}\n"
                 "Choose the sheet again. You may name another sheet or the same one.")
        first = s1["model"]
        s1b = stage(release, study, "sheet_again", BRIEF1, again, chk1, logged, mode, tuple(dict.fromkeys((first, FALLBACK))))
        out["back_to_the_sheet"] = {"reason": reason, "stage1": {k: s1b[k] for k in ("accepted", "tries")}}
        if s1b["accepted"] and s1b["proposal"]["sheet"] != s1["proposal"]["sheet"]:
            df, s2, chk2 = mapping(s1b, "mapping_again", (MODEL, FALLBACK))
            out["back_to_the_sheet"]["stage2"] = {k: s2[k] for k in ("accepted", "tries")}
            s1 = s1b
        else:
            out["back_to_the_sheet"]["stopped"] = "the same sheet again" if s1b["accepted"] else "no sheet accepted"
    out.update({"sheet": s1["proposal"]["sheet"], "header_row": s1["proposal"]["header_row"], "rows": int(len(df)), "traits_on_the_card": tr,
                "knowledge_source": ks})
    if s2["accepted"] and s2["proposal"].get("map") is False:
        out.update({"mapping": s2["proposal"], "model": {"sheet": s1["model"], "mapping": s2["model"]}, "join": None,
                    "verdict": {"accepted": True, "map": False, "reason": s2["proposal"].get("reason")}})
    elif s2["accepted"]:
        v = chk2(s2["proposal"])                           # the verdict of the accepted proposal
        out["mapping"] = s2["proposal"]
        out["model"] = {"sheet": s1["model"], "mapping": s2["model"]}
        out["verdict"] = {k: x for k, x in _public(v).items() if k in ("accepted", "problems", "rules", "export")}
        J = joinable(v, G, sym, a, key)
        rows = v["_rows"]
        un = rows[rows["__subject"].isna()] if "__subject" in rows else rows.iloc[:0]
        vals = [x for c in (s2["proposal"].get("subject") or {}).get("columns", []) if c in un for x in un[c].dropna()]
        out["join"] = {"rows": int(len(J)), "targets": int(J.target.nunique()), "symbols": int(J.symbol.nunique()), "traits": sorted(set(J.trait)),
                       "table_rows_with_no_subject": int(len(un)), "their_values": sorted(set(map(str, vals)))[:60]}
        jf = d / ("join_rerun.parquet" if mode == "rerun" else "join.parquet")
        J.to_parquet(jf, index=False)
        out["join"].update({"file": jf.name, "file_sha256": __import__("pharos_gap_kg.graph", fromlist=["file_sha"]).file_sha(jf)})
    return _record(release, study, out, proposals is not None, mode == "rerun")


def _record(release: str, study: str, out: dict, control: bool, rerun: bool = False) -> dict:
    out["spent_usd"] = round(spent(release), 4)
    p = (paths.DOCS_DATA / f"add_{STUDIES[study]['key']}_control_{release}.json" if control else record_path(release, study, rerun))
    p.write_text(json.dumps(out, indent=1, default=str) + "\n")
    return out


def hand_proposals() -> tuple[dict, dict]:
    """Study A as `later.py` built it, in the proposal's form: SD2, its header row, the predicted effector
    genes split on '/', each linked to Alzheimer disease as the graph links a GWAS gene to a trait."""
    from pharos_gap_kg import graph
    r = graph.mapping().rules["gwas_condition"]
    sheets = read_sheets(registry.require(STUDIES["A"]["key"]))
    head = next(i for i, x in enumerate(sheets["SD2"]) if len(x) > 1 and x[1] == "Predicted effector gene(s)") + 1
    return ({"sheet": "SD2", "header_row": head},
            {"map": True, "subject": {"columns": ["Predicted effector gene(s)"], "resolver": "gene_symbol", "separator": "/", "category": "biolink:Protein"},
             "object": {"constant": adstudy.TRAIT, "category": "biolink:Disease"}, "predicate": r.predicate, "association_category": r.category,
             "knowledge_level": r.knowledge_level, "agent_type": r.agent_type, "qualifiers": {},
             "edge_properties": ["locus", "rsid", "chromosome", "start", "end", "p_value", "Novel"], "keep_if": None, "split_by": None})


def control(release: str) -> dict:
    out = run(release, "A", hand_proposals())
    if not out.get("stage2", {}).get("accepted"):
        raise SystemExit(f"REFUSED: the hand mapping of study A is refused, a defect of the validator or the form: "
                         f"{(out.get('stage2') or out['stage1'])['tries'][-1].get('problems')}")
    J = pd.read_parquet(run_dir(release, "A_control") / "join.parquet")
    H = pd.read_parquet(paths.DATA / "later" / release / f"{STUDIES['A']['key']}.parquet")
    same = set(zip(J.target, J.trait)) == set(zip(H.target, H.trait))
    out["equals_later_py"] = same
    _record(release, "A", out, True)
    if not same:
        raise SystemExit(f"REFUSED: the hand mapping's join differs from later.py's table: {len(set(zip(J.target, J.trait)) ^ set(zip(H.target, H.trait)))} pairs")
    return out


# ---- the person's view and join --------------------------------------------------------------------

def show(release: str, study: str, rerun: bool = False) -> str:
    r = json.loads(record_path(release, study, rerun).read_text())
    if not r.get("stage2", {}).get("accepted"):
        return f"study {study}: no accepted mapping; see {record_path(release, study)}"
    m, j = r["mapping"], r["join"]
    lines = [f"{r['description']}", "", f"Sheet {r['sheet']!r}, column names in row {r['header_row']}, {r['rows']:,} rows.",
             f"Proposed by {r['model']['sheet']} (sheet) and {r['model']['mapping']} (mapping); accepted by the validator.", ""]
    for side in ("subject", "object"):
        lines.append(f"  {side}: {json.dumps(m[side])}")
    for k in ("predicate", "association_category", "knowledge_level", "agent_type", "qualifiers", "edge_properties", "keep_if", "split_by", "by_value", "notes"):
        if m.get(k) not in (None, {}, []):
            lines.append(f"  {k}: {json.dumps(m[k])}")
    lines += ["", f"Knowledge source (set by code): {r['knowledge_source']}",
              f"It makes {j['rows']:,} links from {j['targets']:,} graph genes to {', '.join(j['traits'])}.",
              f"{j['table_rows_with_no_subject']:,} rows give no graph gene; among their values: {', '.join(j['their_values'][:30])}",
              "", f"To join it:  python -I -m pharos_gap_kg.add join {release} {study}"]
    return "\n".join(lines)


def join(release: str, study: str, rerun: bool = False) -> dict:
    """The person's join: the accepted mapping's rows, written beside the bundle as `later` writes its tables."""
    from pharos_gap_kg import graph
    s = STUDIES[study]
    r = json.loads(record_path(release, study, rerun).read_text())
    if not r.get("stage2", {}).get("accepted"):
        raise SystemExit(f"REFUSED: study {study} has no accepted mapping")
    if not r.get("join") or not r["join"].get("rows"):
        raise SystemExit(f"REFUSED: study {study}'s mapping makes no link from a graph gene to a graph trait, which is what the answers read")
    if study == "A":
        raise SystemExit("REFUSED: study A's table is already beside the bundle, built by hand (later.py); its mapping is compared, not joined (docs/PROTOCOL.md)")
    src = run_dir(release, study) / r["join"].get("file", "join.parquet")
    if graph.file_sha(src) != r["join"]["file_sha256"]:
        raise SystemExit(f"REFUSED: {src} is not the table the record describes; run the stages again")
    d = paths.DATA / "later" / release
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{s['key']}.parquet"
    pd.read_parquet(src).to_parquet(f, index=False)
    mf = d / "manifest.json"
    man = json.loads(mf.read_text()) if mf.is_file() else {"release": release, "files": {}, "studies": {}}
    a = s["about"]
    man["files"][f.name] = graph.file_sha(f)
    man["studies"][s["key"]] = {"citation": s["short"], "doi": a["doi"], "published": a["published"], "added": date.today().isoformat(),
                                "year": date.fromisoformat(a["published"]).year, "trait": r["join"]["traits"][0] if len(r["join"]["traits"]) == 1 else r["join"]["traits"],
                                "licence": a["licence"], "table": r["sheet"], "what": f"{a['what']}, as a model mapped its table and the person joined it",
                                "mapped_by": r["model"]}
    mf.write_text(json.dumps(man, indent=1, sort_keys=True) + "\n")
    r["joined"] = {"at": datetime.now(timezone.utc).isoformat(timespec="seconds"), "file": f"data/later/{release}/{f.name}", "file_sha256": man["files"][f.name]}
    record_path(release, study, rerun).write_text(json.dumps(r, indent=1, default=str) + "\n")
    return r


# ---- the checks the protocol states ----------------------------------------------------------------

def catalog_genes(accession: str, column: str = "REPORTED GENE(S)") -> tuple[list[str], int]:
    """The genes the GWAS Catalog release of 5 Oct 2026 gives for one study in one column, over its associations,
    and the number of associations."""
    z = zipfile.ZipFile(registry.require("gwas_catalog_20261005"))
    with z.open(z.namelist()[0]) as f:
        df = pd.read_csv(f, sep="\t", dtype=str, quoting=3, keep_default_na=False, usecols=["STUDY ACCESSION", column])
    mine = df[df["STUDY ACCESSION"] == accession][column]
    gs = [g.strip() for cell in mine for g in re.split(r"[,;]| - ", cell)]
    return sorted({g for g in gs if g and g.upper() not in ("NR", "INTERGENIC")}), int(len(mine))


def candidates(release: str, since: str = "2026-01-01", journals: tuple = ("Nat Genet", "Nat Commun", "Nature"), least: int = 50,
               share: float = 0.9) -> dict:
    """The rule that chose study D (the second amendment of 4 Nov 2026): the studies of the GWAS Catalog release of 5 Oct 2026 published
    from `since` in `journals`, whose one mapped trait is a disease-like trait of the graph other than the traits `find` gives for the
    diseases of the studies before, with at least `least` associations of which at least `share` carry a positional gene, by their
    number of associations, most first. Which of them has an open article and one spreadsheet of tables is read from its article page."""
    from pharos_gap_kg import ask
    G = ask.Graph(paths.bundle(release))
    before = {t["id"] for k, s in STUDIES.items() if k in "ABC" for t in traits(G, s["about"]["disease"])}
    disease = dict(zip(G.d.id, G.d.disease_like))
    z = zipfile.ZipFile(registry.require("gwas_catalog_20261005"))
    with z.open(z.namelist()[0]) as f:
        df = pd.read_csv(f, sep="\t", dtype=str, quoting=3, keep_default_na=False,
                         usecols=["STUDY ACCESSION", "PUBMEDID", "FIRST AUTHOR", "DATE", "JOURNAL", "DISEASE/TRAIT", "MAPPED_GENE", "MAPPED_TRAIT_URI"])
    df = df[(df["DATE"] >= since) & df["JOURNAL"].isin(journals)]
    out = []
    for acc, x in df.groupby("STUDY ACCESSION", sort=False):
        uris = set(x["MAPPED_TRAIT_URI"])
        if len(uris) != 1 or "," in next(iter(uris)):
            continue
        t = next(iter(uris)).rsplit("/", 1)[-1].replace("_", ":", 1)
        n, mapped = len(x), int((x["MAPPED_GENE"].str.strip() != "").sum())
        if t in G.di and disease.get(t) is True and t not in before and n >= least and mapped / n >= share:
            r = x.iloc[0]
            out.append({"accession": acc, "pmid": r["PUBMEDID"], "first_author": r["FIRST AUTHOR"], "published": r["DATE"], "journal": r["JOURNAL"],
                        "trait_as_curated": r["DISEASE/TRAIT"], "trait": t, "trait_name": G.d.name[G.di[t]], "associations": n,
                        "share_with_a_positional_gene": round(mapped / n, 4)})
    out.sort(key=lambda r: (-r["associations"], r["accession"]))
    rec = {"release": release, "catalog": "gwas_catalog_20261005", "since": since, "journals": list(journals), "at_least_associations": least,
           "at_least_share_with_a_positional_gene": share, "traits_of_the_studies_before": sorted(before), "candidates": out}
    (paths.DOCS_DATA / f"add_candidates_{release}.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def compare(release: str) -> dict:
    from pharos_gap_kg import clingen
    out = {"release": release}
    rA = json.loads(record_path(release, "A").read_text()) if record_path(release, "A").is_file() else None
    if rA and rA.get("stage2", {}).get("accepted"):
        J = pd.read_parquet(run_dir(release, "A") / "join.parquet")
        H = pd.read_parquet(paths.DATA / "later" / release / f"{STUDIES['A']['key']}.parquet")
        got, hand = set(zip(J.target, J.trait)), set(zip(H.target, H.trait))
        sym = dict(zip(H.target, H.symbol)) | dict(zip(J.target, J.symbol))
        out["A"] = {"hand_pairs": len(hand), "model_pairs": len(got), "reproduced": got == hand,
                    "only_the_model": sorted(f"{sym[t]} {d}" for t, d in got - hand), "only_by_hand": sorted(f"{sym[t]} {d}" for t, d in hand - got)}
    else:
        out["A"] = {"reproduced": False, "why": "no accepted mapping"}
    rB = json.loads(record_path(release, "B").read_text()) if record_path(release, "B").is_file() else None
    if rB and rB.get("stage2", {}).get("accepted") is not None:
        from pharos_gap_kg import ask
        G = ask.Graph(paths.bundle(release))
        sym = Symbols(G, release)
        mine = set(pd.read_parquet(run_dir(release, "B") / "join.parquet").symbol) if rB.get("join") else set()   # no mapping joins no gene
        reported, n_assoc = catalog_genes(CATALOG_STUDY["B"])
        mapped, _ = catalog_genes(CATALOG_STUDY["B"], "MAPPED_GENE")
        cat = {sym.symbol[t] for g in reported for t in sym(g)}
        k1, k2 = len(mine & cat), len(cat & mine)
        s1, s2 = (clingen.wilson(k1, len(mine)) if mine else None), (clingen.wilson(k2, len(cat)) if cat else None)
        if not cat:
            verdict = "not computable: the Catalog reports no gene for the study's associations"
        else:
            verdict = ("disagree" if s2[0] < 0.5 or (s1 and s1[0] < 0.5) else "agree" if s1 and s1[0] >= 0.8 and s2[0] >= 0.8 else "partial")
        out["B"] = {"catalog_study": CATALOG_STUDY["B"], "sheet": rB.get("sheet"), "joined_a_table": bool(rB.get("join")), "catalog_associations": n_assoc,
                    "catalog_mapped_genes_not_the_registered_check": len(mapped), "catalog_reported_genes": len(reported), "catalog_genes_in_the_graph": len(cat),
                    "table_genes": len(mine), "table_genes_the_catalog_reports": k1, "share_of_table_genes": s1,
                    "catalog_genes_the_table_holds": k2, "share_of_catalog_genes": s2, "reading": verdict,
                    "only_the_table": sorted(mine - cat), "only_the_catalog": sorted(cat - mine)}
    else:
        out["B"] = {"reading": None, "why": "no accepted mapping"}
    for name, study, rerun in (("B_rerun", "B", True), ("C", "C", False), ("C_rerun", "C", True), ("D", "D", False)):
        p = record_path(release, study, rerun)
        out[name] = positional(release, study, json.loads(p.read_text())) if p.is_file() else {"reading": None, "why": "not run"}
    out["rechecked"] = recheck(release)
    (paths.DOCS_DATA / f"add_compare_{release}.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def recheck(release: str) -> dict:
    """Every mapping proposal the runs logged, checked again under each rule of check 8, calling no model (the second amendment of
    4 Nov 2026): the verdict logged, the verdict under each rule, and the proposals whose verdict the rule changes."""
    from pharos_gap_kg import ask
    G = ask.Graph(paths.bundle(release))
    sym, refs = Symbols(G, release), ingest.Refs(release)
    rows, changed = [], []
    for study, s in STUDIES.items():
        f = run_dir(release, study) / "calls.jsonl"
        if not f.is_file():
            continue
        log = [json.loads(x) for x in f.read_text().splitlines()]
        sheets = read_sheets(registry.require(s["key"]))
        a = s["about"]
        tr, ks = traits(G, a["disease"]), infores(a)
        chosen = {x["stage"]: x["proposal"] for x in log if x["stage"] in ("sheet", "sheet_again") and x["accepted"]}
        for x in log:
            if x["stage"] not in ("mapping", "mapping_again") or x["proposal"] is None:
                continue
            df = check1(chosen["sheet" if x["stage"] == "mapping" else "sheet_again"], sheets)[0]
            by = {r: check2(x["proposal"], df, tr, sym, refs, ks, s["key"], G, r) for r in CHECK8}
            row = {"study": study, "stage": x["stage"], "model": x["model"], "attempt": x["attempt"], "logged": x["accepted"],
                   "logged_under": x.get("check_8", "resolved"),
                   **{r: {"accepted": v["accepted"], "check_8": [p for p in v["problems"] if p.startswith("8_identifiers")]} for r, v in by.items()}}
            rows.append(row)
            if by["resolved"]["accepted"] != by["graph"]["accepted"]:
                changed.append({k: row[k] for k in ("study", "stage", "model", "attempt")} | {"resolved": by["resolved"]["accepted"], "graph": by["graph"]["accepted"]})
    return {"rules": CHECK8, "proposals": len(rows), "changed": changed, "each": rows}


def positional(release: str, study: str, r: dict) -> dict:
    """The checks of the amendment of 4 Nov 2026: whether the sheet the stages ended on is a table of the study's loci (rows between
    half and twice the abstract's count), and the joined genes against the genes the Catalog maps to the lead variants by position."""
    from pharos_gap_kg import ask, clingen
    n = ABSTRACT_LOCI[study]
    out = {"sheet": r.get("sheet"), "rows": r.get("rows"), "abstract_loci": n, "rows_band": [n // 2, 2 * n],
           "a_table_of_the_studys_loci": r.get("rows") is not None and n // 2 <= r["rows"] <= 2 * n,
           "joined_a_table": bool(r.get("join") and r["join"].get("rows"))}
    G = ask.Graph(paths.bundle(release))
    sym = Symbols(G, release)
    mine = set(pd.read_parquet(run_dir(release, study) / r["join"]["file"]).symbol) if out["joined_a_table"] else set()
    mapped, n_assoc = catalog_genes(CATALOG_STUDY[study], "MAPPED_GENE")
    cat = {sym.symbol[t] for g in mapped for t in sym(g)}
    k = len(mine & cat)
    s1, s2 = (clingen.wilson(k, len(mine)) if mine else None), (clingen.wilson(k, len(cat)) if cat else None)
    if not cat:
        reading = "not computable: the Catalog maps no graph gene to the study's associations"
    else:
        reading = "disagree" if s2[0] < 0.5 or (s1 is not None and s1[0] < 0.5) else "agree" if s1 is not None and s1[0] >= 0.8 and s2[0] >= 0.8 else "partial"
    out.update({"catalog_study": CATALOG_STUDY[study], "catalog_associations": n_assoc, "catalog_mapped_genes": len(mapped), "of_them_graph_genes": len(cat),
                "table_genes": len(mine), "shared": k, "share_of_table_genes_the_catalog_maps": s1, "share_of_catalog_genes_the_table_holds": s2,
                "reading": reading, "what_it_says": "how far the joined genes lie at the lead variants, by the Catalog's positional mapping; not what the paper named",
                "only_the_table": sorted(mine - cat), "only_the_catalog": sorted(cat - mine)})
    return out


def fetch(key: str, src: Path | None = None) -> Path:
    """A study's registered file, from the publisher or a copy already downloaded; refused unless size and digest are the registry's."""
    from pharos_gap_kg import clingen
    s = registry.SOURCES[key]
    side = s.path.with_name(s.filename + ".lineage.json")
    if s.path.is_file() and side.is_file() and adstudy.sha(s.path.read_bytes()) == s.sha256:
        return s.path
    data = Path(src).read_bytes() if src else urllib.request.urlopen(urllib.request.Request(s.url, headers=clingen.UA), timeout=120).read()
    if len(data) != s.bytes or adstudy.sha(data) != s.sha256:
        raise SystemExit(f"REFUSED {s.filename}: {len(data):,} bytes, sha256 {adstudy.sha(data)[:16]}; the registry pins {s.bytes:,} bytes, {s.sha256[:16]}")
    s.path.parent.mkdir(parents=True, exist_ok=True)
    tmp = s.path.with_name(s.filename + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, s.path)
    side.write_text(json.dumps({"key": key, "url": s.url, "bytes": len(data), "sha256": s.sha256, "published": s.published,
                                "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"), "from": str(src) if src else s.url}, indent=1) + "\n")
    return s.path


# ---- the same stages as MCP tools, for a client's own model ----------------------------------------

ABOUT_KEYS = ("citation", "journal", "published", "doi", "licence", "what", "disease", "add")


class BundleRefs:
    """The lists `ingest`'s validator reads, from the bundle alone, where the release's tables are not on the machine:
    the graph's targets and traits. A vocabulary the bundle does not hold is not checked."""

    def __init__(self, G):
        self.uniprot, self.gene, self.target = {}, {}, {}
        self.symbol = dict(zip(G.t.id, G.t["name"]))
        self.vocab = {"UniProtKB": set(G.t.id), "MONDO": {i for i in G.d.id if i.startswith("MONDO:")}, "NCBIGene": None, "DOID": None, "UBERON": None, "GO": None}
        self.xref = set()


class Session:
    """What the tools share for one bundle: the graph, the symbol rule and the validator's lists."""

    def __init__(self, G):
        self.G, self.release = G, G.manifest["release"]
        try:
            self.refs, self.lists = ingest.Refs(self.release), "the release's tables"
            self.sym = Symbols(G, self.release)
        except (SystemExit, FileNotFoundError, OSError):
            self.refs, self.lists = BundleRefs(G), "the bundle's targets and traits only"
            self.sym = Symbols(G, None)


def _about_problems(about) -> list[str]:
    if not isinstance(about, dict) or set(ABOUT_KEYS) - set(about):
        return [f"the description needs the fields {list(ABOUT_KEYS)}"]
    bad = [f"{key} must be non-empty text" for key in ABOUT_KEYS
           if not isinstance(about[key], str) or not about[key].strip()]
    if bad:
        return bad
    try:
        date.fromisoformat(about["published"])
    except ValueError:
        return ["published must be a date written YYYY-MM-DD"]
    # The DOI determines the sidecar filename and study identity. Refuse an empty
    # or unrecognised identifier before reading the supplied file or writing it.
    if re.fullmatch(r"10\.\d{4,9}/\S+", about["doi"], flags=re.IGNORECASE) is None:
        return ["doi must be a DOI written 10.<registrant>/<suffix>"]
    return []


def tool_card(S: Session, file: str, about: dict, sheet: str | None = None, header_row: int | None = None) -> dict:
    """The first card (every sheet), or with a sheet and header row the second (the table, and the traits that match the disease)."""
    bad = _about_problems(about)
    if bad:
        return {"accepted": False, "problems": bad}
    sheets = read_sheets(Path(file).expanduser())
    if sheet is None:
        return {"stage": "sheet", "brief": BRIEF1, "description": text(about), "card": card1(sheets),
                "next": "Write the proposal the brief asks for, then call add_check with stage 'sheet'."}
    df, problems = check1({"sheet": sheet, "header_row": header_row}, sheets)
    if problems:
        return {"accepted": False, "problems": problems}
    return {"stage": "mapping", "brief": brief2(), "description": text(about), "card": card2(df, sheet, header_row, traits(S.G, about["disease"])),
            "profile": table_profile(df),
            "next": "Write the mapping the brief asks for, then call add_check with stage 'mapping' and the same sheet and header_row."}


def tool_check(S: Session, file: str, about: dict, stage: str, proposal: dict, sheet: str | None = None, header_row: int | None = None,
               rule: str = RULE) -> dict:
    """The validator's verdict on a proposal for one stage, and for an accepted mapping what it would join. rule: check 8's first
    figure (CHECK8); the served tool always uses the current one."""
    bad = _about_problems(about)
    if bad:
        return {"accepted": False, "problems": bad}
    if stage not in ("sheet", "mapping"):
        return {"accepted": False, "problems": ["stage must be sheet or mapping"]}
    sheets = read_sheets(Path(file).expanduser())
    if stage == "sheet":
        df, problems = check1(proposal, sheets)
        return {"accepted": not problems, "problems": problems, **({"rows": int(len(df))} if df is not None else {})}
    df, problems = check1({"sheet": sheet, "header_row": header_row}, sheets)
    if problems:
        return {"accepted": False, "problems": problems}
    v = check2(proposal, df, traits(S.G, about["disease"]), S.sym, S.refs, infores(about), "added", S.G, rule)
    out = {"accepted": v["accepted"], "problems": v["problems"], "checked_against": S.lists,
           "report": mapping_report(v, df, proposal, S.G, S.refs, S.sym, header_row)}
    if v["accepted"] and v.get("map") is not False:
        J = joinable(v, S.G, S.sym, about, "added")
        out["would_join"] = {"links": int(len(J)), "graph_genes": int(J.target.nunique()), "traits": sorted(set(J.trait))}
    return out


def mapping_report(v: dict, df: pd.DataFrame, proposal: dict, G, refs, sym,
                   header_row: int, sample_limit: int = 10) -> dict:
    """Account for every source row without changing the validator's decision.

    Outcomes describe resolved mapping contributions before exporter merging, not
    a promise of import success. Sheet row numbers follow read_sheets (blank rows
    skipped), the same convention as the sheet card.
    """
    from collections import Counter
    labels = {i: set() for i in df.index}
    issues = []
    issue_counts = Counter()
    def issue(i, side, reason, value=None):
        issue_counts[reason] += 1
        if len(issues) < sample_limit:
            issues.append({"source_row": int(i) + 1, "sheet_row": header_row + int(i) + 1,
                           "side": side, "reason": reason, "value": value})

    if v.get("map") is False:
        for i in labels:
            labels[i].add("declined")
    elif "_rows" not in v:
        for i in labels:
            labels[i].add("not_evaluated")
    else:
        rows, p = v["_rows"], v["_p2"]
        original_rules = ingest.rules_of(proposal)
        public_rules = v.get("rules", [])
        for n, (value, rule, skip) in enumerate(ingest.rules_of(p)):
            source = df if value is None else df[df[proposal["split_by"]].astype(str) == str(value)]
            if rule is None:
                for i in source.index:
                    labels[i].add("skipped")
                    issue(i, None, "explicit_skip", skip)
                continue
            mine = rows if value is None else rows[rows[p["split_by"]].astype(str) == str(value)]
            checks = public_rules[n].get("checks", {}) if n < len(public_rules) else {}
            if n >= len(public_rules) or checks.get("1_parses"):
                for i in source.index:
                    labels[i].add("not_evaluated")
                continue
            sid, _ = ingest.resolve(mine, rule["subject"], refs)
            oid, _ = ingest.resolve(mine, rule["object"], refs)
            k = rule.get("keep_if")
            kept = [ingest.keeps(k, x) for x in mine[k["column"]]] if k else [True] * len(mine)
            for i, a, b, keep in zip(mine["__row"], sid, oid, kept):
                if not a:
                    labels[i].add("unresolved")
                    issue(i, "subject", "unresolved_subject")
                elif not keep:
                    labels[i].add("gap")
                elif not b:
                    labels[i].add("unresolved")
                    issue(i, "object", "unresolved_object")
                elif (a in G.ti and b in G.di) or (b in G.ti and a in G.di):
                    labels[i].add("graph_pair")
                else:
                    labels[i].add("outside_graph_pair")
                    issue(i, None, "outside_graph_pair")
            # Inspect each token, including unresolved tokens in a partly resolved
            # cell. Resolution itself remains the existing exact-symbol rule.
            original = original_rules[n][1]
            for side in ("subject", "object"):
                spec = original[side]
                if spec.get("resolver") != "gene_symbol":
                    continue
                for i, record in source.iterrows():
                    for col in spec["columns"]:
                        raw = record[col]
                        if raw is None or pd.isna(raw):
                            continue
                        tokens = str(raw).split(spec["separator"]) if spec.get("separator") else [str(raw)]
                        resolved = False
                        for token in tokens:
                            token = token.strip()
                            if not token:
                                continue
                            ids = sym(token)
                            if not ids:
                                issue(i, side, "unmatched_symbol", token)
                            elif len(ids) > 1:
                                issue(i, side, "multiple_symbol_matches", token)
                            resolved |= bool(ids)
                        if resolved:
                            break  # same fallback-column rule as _resolve_side
    counts = Counter()
    examples = {}
    for i, outcomes in labels.items():
        outcome = "+".join(sorted(outcomes)) if outcomes else "unmapped_split_value"
        counts[outcome] += 1
        if len(examples.setdefault(outcome, [])) < sample_limit:
            examples[outcome].append({"source_row": int(i) + 1, "sheet_row": header_row + int(i) + 1})
    frames = v.get("_frames", {})
    return {"accepted": v["accepted"], "source_rows": len(df), "accounted_source_rows": sum(counts.values()),
            "outcomes": dict(sorted(counts.items())), "outcome_examples": examples,
            "issue_counts": dict(issue_counts), "issue_examples": issues,
            "issue_examples_truncated": sum(issue_counts.values()) > len(issues),
            "exported_edges": len(frames["edges"]) if "edges" in frames else None,
            "exported_gaps": len(frames["gaps"]) if "gaps" in frames else None,
            "reading": "Outcomes count source rows once, including combined outcomes after cell expansion. "
                       "Export totals count merged edges and gaps only when the validator produced them. "
                       "Graph-pair contributions do not establish scientific validity. source_row is one-based below the header; "
                       "sheet_row follows the sheet card's numbering after blank rows are skipped. "
                       "Unmatched tokens are not automatically corrected; multiple matches are expanded by the existing resolver."}


def tool_join(S: Session, file: str, about: dict, sheet: str, header_row: int, proposal: dict) -> dict:
    """Join an accepted mapping beside the bundle. Only on the person's word."""
    from pharos_gap_kg import ask, graph
    bad = _about_problems(about)
    if bad:
        return {"joined": False, "problems": bad}
    source_path = Path(file).expanduser()
    source_sha = graph.file_sha(source_path)
    df, problems = check1({"sheet": sheet, "header_row": header_row}, read_sheets(Path(file).expanduser()))
    if problems:
        return {"joined": False, "problems": problems}
    key = re.sub(r"[^a-z0-9]+", "_", about["doi"].casefold()).strip("_")
    v = check2(proposal, df, traits(S.G, about["disease"]), S.sym, S.refs, infores(about), key, S.G)
    if not v["accepted"] or v.get("map") is False:
        return {"joined": False, "problems": v["problems"] or ["the mapping makes no edges"]}
    J = joinable(v, S.G, S.sym, about, key)
    if J.empty:
        return {"joined": False, "problems": ["the mapping makes no link from a graph gene to a graph trait, which is what the answers read"]}
    from pharos_gap_kg import importstore
    if graph.file_sha(source_path) != source_sha:
        return {"joined": False, "problems": ["the source changed during validation; inspect and retry"]}
    metadata = {"citation": about["citation"].split(".")[0] + f", {about['journal']} ({date.fromisoformat(about['published']).year})",
                "doi": about["doi"], "published": about["published"], "added": date.today().isoformat(),
                "year": date.fromisoformat(about["published"]).year,
                "trait": sorted(set(J.trait))[0] if J.trait.nunique() == 1 else sorted(set(J.trait)),
                "licence": about["licence"], "table": sheet,
                "what": f"{about['what']}, as a client's model mapped its table and the person joined it",
                "mapped_by": "the MCP client's model"}
    def stable(value):
        if isinstance(value, dict):
            return {str(k): stable(v) for k, v in sorted(value.items(), key=lambda x: str(x[0]))}
        if isinstance(value, (set, frozenset)):
            return sorted(stable(x) for x in value)
        if isinstance(value, (list, tuple)):
            return [stable(x) for x in value]
        return value
    resolver_state = {name: stable(getattr(S.refs, name)) for name in ("uniprot", "gene", "target", "vocab", "xref")}
    resolver_state["symbols"] = stable({"map": S.sym.map, "ids": S.sym.ids})
    recipe = {"source_sha256": source_sha, "source_filename": source_path.name, "about": about,
              "sheet": sheet, "header_row": header_row, "mapping": proposal, "check_8": RULE,
              "bundle_manifest_sha256": graph.file_sha(S.G.dir / "manifest.json"),
              "code_sha256": {"add": graph.file_sha(Path(__file__)), "ingest": graph.file_sha(Path(ingest.__file__))},
              "checked_against": S.lists,
              "resolver_sha256": importstore.digest(importstore.encoded(resolver_state)),
              "biolink_sha256": ingest.sha(ingest.biolink())}
    try:
        result = importstore.commit(ask.later_dir(S.G.dir), S.release, key, J, metadata, recipe)
    except (ValueError, OSError) as e:
        return {"joined": False, "problems": [str(e)]}
    return {"joined": True, "links": int(len(J)), "graph_genes": int(J.target.nunique()), **result,
            "answer": "Already present; verified receipt and table." if result["already_present"] else
                      "Joined beside the bundle with an import receipt. Restart the server for the answers to read it."}


def replay_import(S: Session, receipt_file: str, file: str) -> dict:
    """Re-evaluate a saved recipe offline, with no join or API request.

    The source and bundle must match. New code may revalidate an old recipe;
    that is reported as a new evaluation, not a reproduction of the old verdict.
    """
    from pharos_gap_kg import graph
    receipt = json.loads(Path(receipt_file).read_text())
    recipe = receipt["recipe"]
    if graph.file_sha(Path(file).expanduser()) != recipe["source_sha256"]:
        return {"accepted": False, "problems": ["the source does not match the receipt; inspect its changes first"]}
    if graph.file_sha(S.G.dir / "manifest.json") != recipe["bundle_manifest_sha256"]:
        return {"accepted": False, "problems": ["the bundle does not match the receipt; inspect its changes first"]}
    result = tool_check(S, file, recipe["about"], "mapping", recipe["mapping"], recipe["sheet"], recipe["header_row"])
    current = {"add": graph.file_sha(Path(__file__)), "ingest": graph.file_sha(Path(ingest.__file__))}
    result["replay"] = {"writes": False, "original_output_sha256": receipt["output_sha256"],
                        "same_validator_code": current == recipe["code_sha256"],
                        "reading": "Source and bundle match; mapping was checked again under the current validator. "
                                   "This does not verify equality of exported table bytes or scientific correctness."}
    return result


# ---- self-test -------------------------------------------------------------------------------------

def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("the descriptions read as the protocol fixes them",
       text(STUDIES["A"]["about"]) == ("Uffelmann et al. Alzheimer's disease genetics converges on microglial and neuronal dysregulation. Nature, published "
                                       "28 October 2026. doi:10.1038/s41586-026-11206-y. Licence: CC BY-NC-ND 4.0. A genome-wide association study of "
                                       "Alzheimer's disease. Add the genes the study names at each associated locus.")
       and "Nature Genetics, published 28 December 2023. doi:10.1038/s41588-023-01584-8. Licence: CC BY 4.0." in text(STUDIES["B"]["about"]))
    ck("the knowledge source is written infores: and a lower-case name", re.fullmatch(r"infores:[a-z0-9][a-z0-9._-]*", infores(STUDIES["B"]["about"])) is not None)
    rows = [["Supplementary Table 2"], ["", "Gene(s)", "locus", "p"], ["", "A/B", "1", "1e-9"], ["", "C", "2", "2e-8"], ["", "Note: a footnote"]]
    df, pr = frame(rows, 2)
    ck("a blank column with no value is padding, and a footnote row is a row", pr == [] and list(df.columns) == ["Gene(s)", "locus", "p"] and len(df) == 3)
    ck("a header row with a repeated name, or a named column's empty name that holds values, is refused",
       frame([["a", "a"], ["1", "2"]], 1)[1] == ["row 1 repeats the names ['a']"] and frame([["a", "", "b"], ["1", "x", "2"]], 1)[1][0].startswith("column 2"))
    ck("a header row outside the sheet, and a sheet the file does not hold, are refused",
       frame(rows, 9)[1][0].startswith("header_row 9") and check1({"sheet": "X", "header_row": 1}, {"S": rows})[1][0].startswith("the file has no sheet"))
    blob = io.BytesIO()
    with zipfile.ZipFile(blob, "w") as z:
        z.writestr("xl/workbook.xml", f'<workbook xmlns="{_S[1:-1]}" xmlns:r="{_RID[1:].split("}")[0]}"><sheets><sheet name="T1" r:id="rId1"/></sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/sharedStrings.xml", f'<sst xmlns="{_S[1:-1]}"><si><t>gene</t></si><si><r><t>AP</t></r><r><t>OE</t></r></si></sst>')
        z.writestr("xl/worksheets/sheet1.xml", f'<worksheet xmlns="{_S[1:-1]}"><sheetData><row r="1"><c r="B1" t="s"><v>0</v></c><c r="C1" t="inlineStr"><is><t>p</t></is></c></row>'
                   '<row r="2"/><row r="3"><c r="B3" t="s"><v>1</v></c><c r="C3"><v>2.17E-11</v></c><c r="D3" t="s"/></row></sheetData></worksheet>')
    ck("an XLSX sheet is read with shared, rich and inline text, numbers as written, and empty rows skipped",
       xlsx_sheets(blob.getvalue()) == {"T1": [["", "gene", "p"], ["", "APOE", "2.17E-11"]]})

    class FakeG:
        t = pd.DataFrame({"id": ["UniProtKB:P1", "UniProtKB:P2", "UniProtKB:P3"], "name": ["AAA", "BBB", "BBB"]})
    sym = Symbols(FakeG, None)
    ck("a gene symbol gives every graph target with that symbol", sym("AAA") == ["UniProtKB:P1"] and sym(" BBB ") == ["UniProtKB:P2", "UniProtKB:P3"] and sym("NA") == [])
    d = pd.DataFrame({"g": ["AAA/BBB", "NA", "ZZZ/AAA"], "p": ["1", "2", "3"]})
    p = {"map": True, "subject": {"columns": ["g"], "resolver": "gene_symbol", "separator": "/", "category": "biolink:Protein"},
         "object": {"constant": "MONDO:1", "category": "biolink:Disease"}, "predicate": "biolink:genetically_associated_with",
         "association_category": "biolink:Association", "knowledge_level": "statistical_association", "agent_type": "data_analysis_pipeline"}
    p2, r2 = transform(p, d, sym, None, "infores:x")
    ck("a separated cell becomes one row per graph gene, a row with none stays once, and the constant fills the object",
       [x if isinstance(x, str) else None for x in r2["__subject"]] == ["UniProtKB:P1", "UniProtKB:P2", "UniProtKB:P3", None, "UniProtKB:P1"] and set(r2["__object"]) == {"MONDO:1"}
       and r2["__row"].tolist() == [0, 0, 0, 1, 2] and p2["subject"]["resolver"] == "curie" and p2["primary_knowledge_source"] == "infores:x")
    class FG:
        ti, di = {"UniProtKB:P1": 0, "UniProtKB:P2": 1}, {"MONDO:1": 0}
    f = graph_rows(FG)(r2, list(r2["__subject"]), list(r2["__object"]), [True] * len(r2))
    loose = graph_rows(FG)(r2, ["pharosgapkg:gene_ab12"] * len(r2), list(r2["__object"]), [True] * len(r2))
    ck("check 8 under 'graph' counts a table row once when one of its edges joins a graph target and a graph trait, and an ungrounded node never",
       (f["ok"], f["rows"]) == (2, 3) and (loose["ok"], loose["rows"]) == (0, 3)
       and graph_rows(FG)(r2.iloc[:1], ["UniProtKB:P1"], [None], [False])["ok"] == 1 and set(CHECK8) == {"resolved", "graph"} and RULE == "graph")
    ck("an object constant must be a trait on the card, and only the object takes one",
       side_problems("object", {"constant": "MONDO:9", "category": "biolink:Disease"}, set(), {"MONDO:1"})[0].startswith("1_parses: the constant")
       and side_problems("subject", {"constant": "MONDO:1", "category": "biolink:Disease"}, set(), {"MONDO:1"}) != []
       and side_problems("subject", p["subject"], {"g"}, set()) == [] and side_problems("subject", {**p["subject"], "resolver": "pharos_protein"}, {"g"}, set()) != [])
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) == 4 and argv[0] == "replay-import":
        from pharos_gap_kg import ask
        result = replay_import(Session(ask.Graph(Path(argv[1]))), argv[2], argv[3])
        print(json.dumps(result, indent=1))
        return 0 if result["accepted"] else 1
    if len(argv) >= 2 and argv[0] == "fetch":
        i = argv.index("--from") if "--from" in argv else None
        print(fetch(argv[1], Path(argv[i + 1]) if i else None))
        return 0
    if len(argv) == 2 and argv[0] == "candidates":
        out = candidates(argv[1])
        print(json.dumps({k: v for k, v in out.items() if k != "candidates"} | {"candidates": out["candidates"][:10]}, indent=1))
        return 0
    if len(argv) == 2 and argv[0] in ("control", "compare"):
        out = control(argv[1]) if argv[0] == "control" else compare(argv[1])
        print(json.dumps({k: v for k, v in out.items() if k not in ("stage1", "stage2", "traits_on_the_card", "verdict")}, indent=1, default=str))
        return 0
    rerun = "--rerun" in argv
    argv = [x for x in argv if x != "--rerun"]
    if len(argv) == 3 and argv[0] in ("run", "replay", "rerun", "show", "join") and argv[2] in STUDIES:
        if argv[0] == "show":
            print(show(argv[1], argv[2], rerun))
        else:
            out = run(argv[1], argv[2], mode=argv[0]) if argv[0] in ("run", "replay", "rerun") else join(argv[1], argv[2], rerun)
            print(json.dumps({k: v for k, v in out.items() if k in ("study", "sheet", "header_row", "rows", "model", "join", "joined", "spent_usd")}
                             | {"stage1_accepted": out.get("stage1", {}).get("accepted"), "stage2_accepted": out.get("stage2", {}).get("accepted")}, indent=1, default=str))
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
