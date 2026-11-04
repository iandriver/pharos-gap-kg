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
      python -I -m pharos_gap_kg.add show pharos319 B                 # the accepted mapping, field by field, for the person
      python -I -m pharos_gap_kg.add join pharos319 B                 # the person's join
      python -I -m pharos_gap_kg.add compare pharos319                # A against later.py's table, B against the GWAS Catalog
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
CATALOG_STUDY = {"B": "GCST90308590"}          # how the GWAS Catalog release of 5 Oct 2026 curates study B

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


def check2(p, df: pd.DataFrame, tr: list[dict], sym: Symbols, refs, ks: str, view: str) -> dict:
    """The validator's verdict on a mapping: the additions' checks, then `ingest`'s nine on the table as transformed."""
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
    v = ingest.validate(p2, rows2, refs, export=True, view=view)
    v["_rows"], v["_p2"] = rows2, p2
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


def stage(release: str, study: str, name: str, system: str, user: str, check, replay: list | None = None) -> dict:
    """One stage: the model, one return of a refused proposal, then the fallback model from the start under the same rule.
    replay: the logged replies of an earlier run, checked again in place of calling a model."""
    log = run_dir(release, study) / "calls.jsonl"
    tries = []
    for model in (MODEL, FALLBACK):
        msgs = [{"role": "system", "content": system}, {"role": "user", "content": user}]
        for attempt in (1, 2):
            if replay is not None:
                old = next((x for x in replay if (x["stage"], x["model"], x["attempt"]) == (name, model, attempt)), None)
                if old is None:
                    return {"accepted": False, "tries": tries, "replay": "no logged reply for the next call"}
                r = {"text": old["reply"], "usd": old["usd"], "seconds": old["seconds"]} if old.get("reply") is not None else {"error": old["problems"][0]}
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
            if replay is None:
                with open(log, "a") as f:
                    f.write(json.dumps({**rec, "reply": r.get("text"), "at": datetime.now(timezone.utc).isoformat(timespec="seconds")}, default=str) + "\n")
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


def record_path(release: str, study: str) -> Path:
    return paths.DOCS_DATA / f"add_{STUDIES[study]['key']}_{release}.json"


def run(release: str, study: str, proposals: tuple | None = None, replay: bool = False) -> dict:
    """Both stages for a registered study. proposals: (sheet proposal, mapping proposal) to check without a model (the control).
    replay: check again the replies an earlier run logged, calling no model."""
    s, G, sym, refs, path = _setup(release, study)
    a, key = s["about"], s["key"]
    d = run_dir(release, study if proposals is None else f"{study}_control")
    d.mkdir(parents=True, exist_ok=True)
    sheets = read_sheets(path)
    c1 = card1(sheets)
    ingest.write_json(d / "card1.json", c1)
    user1 = f"DESCRIPTION\n{text(a)}\n\n" + json.dumps(c1, indent=1, ensure_ascii=False)
    holder = {}

    def chk1(p):
        df, problems = check1(p, sheets)
        holder["df"] = df
        return {"accepted": not problems, "problems": problems}
    logged = [json.loads(x) for x in (d / "calls.jsonl").read_text().splitlines()] if replay else None
    if replay and not logged:
        raise SystemExit(f"REFUSED: no logged replies in {d}")
    if proposals is None:
        s1 = stage(release, study, "sheet", BRIEF1, user1, chk1, logged)
    else:
        v = chk1(proposals[0])
        s1 = {"accepted": v["accepted"], "model": "hand", "proposal": proposals[0], "verdict": v, "tries": [{"stage": "sheet", "model": "hand", **v}]}
    out = {"release": release, "study": study, "key": key, "description": text(a), "file": path.name,
           "file_sha256": registry.SOURCES[key].sha256, "sheets": len(sheets), "stage1": {k: s1[k] for k in ("accepted", "tries")}}
    if not s1["accepted"]:
        return _record(release, study, out, proposals is not None)
    df = holder["df"] if proposals is None else check1(proposals[0], sheets)[0]
    if proposals is None:                                  # the frame of the accepted proposal, not of a later refused one
        df = check1(s1["proposal"], sheets)[0]
    tr = traits(G, a["disease"])
    c2 = card2(df, s1["proposal"]["sheet"], s1["proposal"]["header_row"], tr)
    ingest.write_json(d / "card2.json", c2)
    user2 = f"DESCRIPTION\n{text(a)}\n\n" + json.dumps(c2, indent=1, ensure_ascii=False)
    ks = infores(a)
    keepv = {}

    def chk2(p):
        v = check2(p, df, tr, sym, refs, ks, key)
        keepv["v"] = v
        return v
    if proposals is None:
        s2 = stage(release, study, "mapping", brief2(), user2, chk2, logged)
    else:
        v = chk2(proposals[1])
        s2 = {"accepted": v["accepted"], "model": "hand", "proposal": proposals[1], "verdict": v, "tries": [{"stage": "mapping", "model": "hand", "accepted": v["accepted"], "problems": v["problems"]}]}
    out.update({"sheet": s1["proposal"]["sheet"], "header_row": s1["proposal"]["header_row"], "rows": int(len(df)), "traits_on_the_card": tr,
                "stage2": {k: s2[k] for k in ("accepted", "tries")}, "knowledge_source": ks})
    if replay:
        out["replayed"] = "the stages were checked again from the replies logged by the run; no model was called"
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
        J.to_parquet(d / "join.parquet", index=False)
        out["join"]["file_sha256"] = __import__("pharos_gap_kg.graph", fromlist=["file_sha"]).file_sha(d / "join.parquet")
    return _record(release, study, out, proposals is not None)


def _record(release: str, study: str, out: dict, control: bool) -> dict:
    out["spent_usd"] = round(spent(release), 4)
    p = record_path(release, study) if not control else paths.DOCS_DATA / f"add_{STUDIES[study]['key']}_control_{release}.json"
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

def show(release: str, study: str) -> str:
    r = json.loads(record_path(release, study).read_text())
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


def join(release: str, study: str) -> dict:
    """The person's join: the accepted mapping's rows, written beside the bundle as `later` writes its tables."""
    from pharos_gap_kg import graph
    s = STUDIES[study]
    r = json.loads(record_path(release, study).read_text())
    if not r.get("stage2", {}).get("accepted"):
        raise SystemExit(f"REFUSED: study {study} has no accepted mapping")
    if study == "A":
        raise SystemExit("REFUSED: study A's table is already beside the bundle, built by hand (later.py); its mapping is compared, not joined (docs/PROTOCOL.md)")
    src = run_dir(release, study) / "join.parquet"
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
    record_path(release, study).write_text(json.dumps(r, indent=1, default=str) + "\n")
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
    (paths.DOCS_DATA / f"add_compare_{release}.json").write_text(json.dumps(out, indent=1) + "\n")
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
    try:
        date.fromisoformat(about["published"])
    except (TypeError, ValueError):
        return ["published must be a date written YYYY-MM-DD"]
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
            "next": "Write the mapping the brief asks for, then call add_check with stage 'mapping' and the same sheet and header_row."}


def tool_check(S: Session, file: str, about: dict, stage: str, proposal: dict, sheet: str | None = None, header_row: int | None = None) -> dict:
    """The validator's verdict on a proposal for one stage, and for an accepted mapping what it would join."""
    bad = _about_problems(about)
    if bad:
        return {"accepted": False, "problems": bad}
    sheets = read_sheets(Path(file).expanduser())
    if stage == "sheet":
        df, problems = check1(proposal, sheets)
        return {"accepted": not problems, "problems": problems, **({"rows": int(len(df))} if df is not None else {})}
    df, problems = check1({"sheet": sheet, "header_row": header_row}, sheets)
    if problems:
        return {"accepted": False, "problems": problems}
    v = check2(proposal, df, traits(S.G, about["disease"]), S.sym, S.refs, infores(about), "added")
    out = {"accepted": v["accepted"], "problems": v["problems"], "checked_against": S.lists}
    if v["accepted"] and v.get("map") is not False:
        J = joinable(v, S.G, S.sym, about, "added")
        out["would_join"] = {"links": int(len(J)), "graph_genes": int(J.target.nunique()), "traits": sorted(set(J.trait))}
    return out


def tool_join(S: Session, file: str, about: dict, sheet: str, header_row: int, proposal: dict) -> dict:
    """Join an accepted mapping beside the bundle. Only on the person's word."""
    from pharos_gap_kg import ask, graph
    bad = _about_problems(about)
    if bad:
        return {"joined": False, "problems": bad}
    df, problems = check1({"sheet": sheet, "header_row": header_row}, read_sheets(Path(file).expanduser()))
    if problems:
        return {"joined": False, "problems": problems}
    key = re.sub(r"[^a-z0-9]+", "_", about["doi"].casefold()).strip("_")
    v = check2(proposal, df, traits(S.G, about["disease"]), S.sym, S.refs, infores(about), key)
    if not v["accepted"] or v.get("map") is False:
        return {"joined": False, "problems": v["problems"] or ["the mapping makes no edges"]}
    J = joinable(v, S.G, S.sym, about, key)
    if J.empty:
        return {"joined": False, "problems": ["the mapping makes no link from a graph gene to a graph trait, which is what the answers read"]}
    d = ask.later_dir(S.G.dir)
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{key}.parquet"
    J.to_parquet(f, index=False)
    mf = d / "manifest.json"
    man = json.loads(mf.read_text()) if mf.is_file() else {"release": S.release, "files": {}, "studies": {}}
    man["files"][f.name] = graph.file_sha(f)
    man["studies"][key] = {"citation": about["citation"].split(".")[0] + f", {about['journal']} ({date.fromisoformat(about['published']).year})",
                           "doi": about["doi"], "published": about["published"], "added": date.today().isoformat(), "year": date.fromisoformat(about["published"]).year,
                           "trait": sorted(set(J.trait))[0] if J.trait.nunique() == 1 else sorted(set(J.trait)), "licence": about["licence"], "table": sheet,
                           "what": f"{about['what']}, as a client's model mapped its table and the person joined it", "mapped_by": "the MCP client's model"}
    mf.write_text(json.dumps(man, indent=1, sort_keys=True) + "\n")
    return {"joined": True, "links": int(len(J)), "graph_genes": int(J.target.nunique()), "file": str(f),
            "answer": "Joined beside the bundle. Restart the server for the answers to read it."}


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
    ck("an object constant must be a trait on the card, and only the object takes one",
       side_problems("object", {"constant": "MONDO:9", "category": "biolink:Disease"}, set(), {"MONDO:1"})[0].startswith("1_parses: the constant")
       and side_problems("subject", {"constant": "MONDO:1", "category": "biolink:Disease"}, set(), {"MONDO:1"}) != []
       and side_problems("subject", p["subject"], {"g"}, set()) == [] and side_problems("subject", {**p["subject"], "resolver": "pharos_protein"}, {"g"}, set()) != [])
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) >= 2 and argv[0] == "fetch":
        i = argv.index("--from") if "--from" in argv else None
        print(fetch(argv[1], Path(argv[i + 1]) if i else None))
        return 0
    if len(argv) == 2 and argv[0] in ("control", "compare"):
        out = control(argv[1]) if argv[0] == "control" else compare(argv[1])
        print(json.dumps({k: v for k, v in out.items() if k not in ("stage1", "stage2", "traits_on_the_card", "verdict")}, indent=1, default=str))
        return 0
    if len(argv) == 3 and argv[0] in ("run", "replay", "show", "join") and argv[2] in STUDIES:
        if argv[0] == "show":
            print(show(argv[1], argv[2]))
        else:
            out = run(argv[1], argv[2], replay=argv[0] == "replay") if argv[0] in ("run", "replay") else join(argv[1], argv[2])
            print(json.dumps({k: v for k, v in out.items() if k in ("study", "sheet", "header_row", "rows", "model", "join", "joined", "spent_usd")}
                             | {"stage1_accepted": out.get("stage1", {}).get("accepted"), "stage2_accepted": out.get("stage2", {}).get("accepted")}, indent=1, default=str))
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
