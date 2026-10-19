"""Ingest: a model proposes how a Pharos table becomes Biolink edges, and code accepts or refuses it.

THE STAGE. A card describes one Pharos provenance view (a table, restricted the way Pharos's own
`provenance` table restricts it). A model reads the card and the fixed brief and returns one proposal
in a fixed JSON form: which column is the subject and which the object, how each becomes an
identifier (chosen from a menu of resolvers; the model writes no code), their Biolink categories, the
predicate, the association's category, the knowledge level, agent type and knowledge source,
qualifiers, carried columns, and a condition a row must meet. The validator here accepts the proposal
or refuses it with its reasons; a refused proposal goes back to the model once. Nothing a model writes
is approved by a model.

THE VALIDATOR reads Biolink 4.4.3 from docs/data/biolink_4.4.3.json (written by `biolinkslice.py` in
the standards environment) and the release's own vocabularies. Its nine checks are listed in
docs/PROTOCOL.md, "Ingest". The hand mappings of the graph go through it too, as a control.

Run:  python -m pharos_gap_kg.ingest cards pharos319            # write the cards (before any model)
      python -m pharos_gap_kg.ingest control pharos319          # the hand mappings through the validator
      python -m pharos_gap_kg.ingest propose pharos319 --model local:qwen3.6:35b-a3b --set scored
      python -m pharos_gap_kg.ingest propose pharos319 --model anthropic/claude-sonnet-5.5 --set scored
      python -m pharos_gap_kg.ingest stale pharos319           # records whose verdicts the current validator would not give
      python -m pharos_gap_kg.ingest report pharos319
      python -m pharos_gap_kg.ingest export pharos319           # the accepted new views, as KGX parts
      python -m pharos_gap_kg.ingest --selftest

A hosted model needs OPENROUTER_API_KEY; the local one needs `ollama serve`.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import graph, paths

PREFIX = graph.PREFIX
SEED = 20261019
BIOLINK = paths.ROOT / "docs" / "data" / "biolink_4.4.3.json"
CAP_USD = 1.00
PRICES = {"anthropic/claude-haiku-5.5": (0.10, 0.50), "anthropic/claude-sonnet-5.5": (2.00, 10.00), "anthropic/claude-opus-5.5": (4.00, 20.00)}
LOCAL = "local:qwen3.6:35b-a3b"
MODELS = [LOCAL, "anthropic/claude-haiku-5.5", "anthropic/claude-sonnet-5.5", "anthropic/claude-opus-5.5"]
MAX_TOKENS = 6000

# view -> (table, provenance row id). The where-clause is read from that row.
VIEWS = {
    "disease.CTD": ("disease", 58), "disease.DisGeNET": ("disease", 148), "disease.eRAM": ("disease", 59),
    "disease.Monarch": ("disease", 55), "disease.UniProt_Disease": ("disease", 8), "disease.Expression_Atlas": ("disease", 57),
    "disease.JensenLab": ("disease", 53), "tiga": ("tiga", 261), "drug_activity": ("drug_activity", 281),
    "ortholog_disease": ("ortholog_disease", 56), "phenotype.OMIM": ("phenotype", 62), "phenotype.JAX_MGI": ("phenotype", 65),
    "goa": ("goa", 5),
}
SCORED = [v for v in VIEWS if v.startswith("disease.") or v in ("tiga", "drug_activity")]
NEW = ["ortholog_disease", "phenotype.OMIM", "phenotype.JAX_MGI", "goa"]
# Two clauses in the release match no row as written; the rows are those the dataset names.
CLAUSE_FIX = {("dtype", "uniprot"): "UniProt Disease", ("ptype", "JAX/MGI Human Ortholog Phenotyp"): "JAX/MGI Human Ortholog Phenotype"}
CLAUSE = re.compile(r"^\s*(\w+)\s*(=|LIKE)\s*'([^']*)'\s*$", re.I)

RESOLVERS = {
    "pharos_protein": "the column holds a Pharos protein row id (the `protein` table's `id`); gives UniProtKB:<accession>",
    "pharos_protein_gene": "the column holds a Pharos protein row id; gives NCBIGene:<gene id> of that protein's gene",
    "pharos_target": "the column holds a Pharos target row id (the `target` table's `id`); gives UniProtKB:<accession>",
    "curie": "the value is already a compact identifier such as MONDO:0005148 (EFO_0004584 is read as EFO:0004584)",
    "prefix:P": "the value is a local identifier under the prefix P, written P:<value> (for example prefix:NCBIGene)",
    "pattern:P:REGEX": "the first group of the regular expression REGEX matched in the value, under the prefix P",
    "name:KIND": "no identifier: a local node made from the text, marked as not grounded",
}
KEYS = {"map", "reason", "subject", "object", "predicate", "association_category", "knowledge_level", "agent_type",
        "primary_knowledge_source", "qualifiers", "edge_properties", "keep_if", "split_by", "by_value", "notes"}
RULE_KEYS = {"subject", "object", "predicate", "association_category", "knowledge_level", "agent_type",
             "primary_knowledge_source", "qualifiers", "edge_properties", "keep_if"}
SIDE_KEYS = {"columns", "resolver", "category", "name_column"}
# kgx-gaps writes these three itself; a proposal's value for one is carried in a column of its own and moved in after.
NATIVE_QUALIFIERS = ("subject_aspect_qualifier", "object_direction_qualifier", "anatomical_context_qualifier")
WITHHELD = f"{PREFIX}:withheld_object"     # a gap made by the condition names the object it withheld, where it has one
VOCAB_PREFIXES = ("UniProtKB", "NCBIGene", "MONDO", "DOID", "UBERON", "GO")
XREF_PREFIXES = ("OMIM", "Orphanet", "MESH", "UMLS", "EFO")


def data_dir(release: str) -> Path:
    return paths.ROOT / "docs" / "data" / "ingest"


def cards_path(release: str) -> Path:
    return data_dir(release) / f"cards_{release}.json"


def proposals_dir(release: str, model: str) -> Path:
    return data_dir(release) / f"proposals_{release}" / re.sub(r"[^A-Za-z0-9.]+", "_", model)


def ledger_path(release: str) -> Path:
    return data_dir(release) / f"ledger_{release}.jsonl"


def report_path(release: str) -> Path:
    return paths.ROOT / "docs" / "data" / f"ingest_{release}.json"


def sha(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=str) + "\n")
    os.replace(tmp, path)


# ---- views and cards -------------------------------------------------------------------------------

def table(release: str, name: str) -> pd.DataFrame:
    return pq.read_table(paths.tables(release) / f"{name}.parquet").to_pandas()


def provenance(release: str) -> pd.DataFrame:
    pv, ds = table(release, "provenance"), table(release, "dataset")
    return pv.merge(ds.rename(columns={"id": "dataset_id"}), on="dataset_id", how="left")


def clause(where) -> tuple[str, str, str] | None:
    if where is None or (isinstance(where, float) and np.isnan(where)) or not str(where).strip():
        return None
    m = CLAUSE.match(str(where))
    if not m:
        raise SystemExit(f"REFUSED: a where-clause this reader does not know: {where!r}")
    col, op, val = m.group(1), m.group(2).upper(), m.group(3)
    return col, op, CLAUSE_FIX.get((col, val), val)


def view_rows(release: str, view: str, prov: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict]:
    name, pid = VIEWS[view]
    prov = provenance(release) if prov is None else prov
    p = prov[prov.id == pid].iloc[0]
    rows = table(release, name)
    c = clause(p.where_clause)
    if c:
        col, op, val = c
        keep = rows[col].astype(str).str.startswith(val[:-1]) if op == "LIKE" and val.endswith("%") else rows[col] == val
        rows = rows[keep].reset_index(drop=True)
    if not len(rows):
        raise SystemExit(f"REFUSED: the view {view} holds no row")
    where = None if c is None else f"{c[0]} {c[1]} '{c[2]}'"
    meta = {"table": name, "where": where, "where_in_release": None if p.where_clause is None or pd.isna(p.where_clause) else str(p.where_clause),
            "dataset": {k: (None if pd.isna(p[k]) else str(p[k])) for k in ("name", "source", "url", "comments")},
            "comment": None if pd.isna(p.comment) else str(p.comment)}
    return rows, meta


def cut(v, n: int = 80) -> str:
    s = str(int(v)) if isinstance(v, float) and np.isfinite(v) and v == int(v) else str(v)
    return s if len(s) <= n else s[:n - 3] + "..."


def present(s: pd.Series) -> pd.Series:
    return s.notna() & (s.astype(str).str.strip() != "")


def card(release: str, view: str, prov: pd.DataFrame | None = None) -> dict:
    rows, meta = view_rows(release, view, prov)
    schema = pq.read_schema(paths.tables(release) / f"{meta['table']}.parquet")
    rng = np.random.default_rng(SEED)
    cols = []
    for c in rows.columns:
        s = rows[c][present(rows[c])]
        vc = s.map(cut).value_counts()
        col = {"name": c, "type": str(schema.field(c).type), "rows_with_a_value": int(len(s)), "distinct": int(len(vc))}
        if len(vc) <= 25:
            col["values"] = {cut(k): int(v) for k, v in vc.items()}
        else:
            idx = rng.choice(len(vc), size=6, replace=False)
            col["examples"] = [cut(vc.index[i]) for i in sorted(idx)]
        cols.append(col)
    ex = rows.sample(n=min(8, len(rows)), random_state=SEED)
    examples = [{c: (None if not present(pd.Series([v])).iloc[0] else cut(v)) for c, v in r.items()} for r in ex.to_dict("records")]
    return {"view": view, **meta, "rows": int(len(rows)), "columns": cols, "example_rows": examples}


def make_cards(release: str) -> dict:
    prov = provenance(release)
    cards = {v: card(release, v, prov) for v in VIEWS}
    out = {"release": release, "seed": SEED, "brief_sha256": hashlib.sha256(BRIEF.encode()).hexdigest(), "cards": cards}
    path = cards_path(release)
    if path.exists():
        old = json.loads(path.read_text())
        if old["cards"] != json.loads(json.dumps(cards, default=str)) or old["brief_sha256"] != out["brief_sha256"]:
            raise SystemExit(f"REFUSED: {path} exists and differs; the cards are frozen once a model has read them")
        return old
    write_json(path, out)
    return out


def frozen_cards(release: str) -> dict:
    path = cards_path(release)
    if not path.exists():
        raise SystemExit(f"no cards; run:  python -m pharos_gap_kg.ingest cards {release}")
    out = json.loads(path.read_text())
    if out["brief_sha256"] != hashlib.sha256(BRIEF.encode()).hexdigest():
        raise SystemExit("REFUSED: the brief changed after the cards were frozen")
    return out


# ---- the brief -------------------------------------------------------------------------------------

BRIEF = """You map one table of Pharos (the NIH Illuminating the Druggable Genome database) into a knowledge graph that follows the Biolink Model, version 4.4.3. You will be given a card that describes the table: where Pharos says it came from, each column with its type, counts and values or examples, and eight example rows.

Propose how each row becomes one edge (an association) between a subject node and an object node. Code will check your proposal against the Biolink Model and against the data, and may send it back once with the reasons it was refused.

Return ONE JSON object and nothing else, with these keys:

  "map": true, or false if this table should not become edges (then give only "map" and "reason")
  "subject": {"columns": [column, ...], "resolver": one of the resolvers below, "category": a Biolink class CURIE such as "biolink:Gene", "name_column": optional column that names the node}
  "object":  the same four keys
  "predicate": a Biolink predicate CURIE, such as "biolink:related_to"
  "association_category": a Biolink association class CURIE, such as "biolink:Association"
  "knowledge_level": one of KNOWLEDGE_LEVELS
  "agent_type": one of AGENT_TYPES
  "primary_knowledge_source": an infores identifier for the resource that made the claim, such as "infores:hgnc"
  "qualifiers": {} or {Biolink qualifier slot name: a constant value, or {"sign_of": column} for a direction qualifier}
  "edge_properties": [columns whose values the edge should carry]
  "keep_if": null, or {"column": c, "at_least": number, "reason": why a row below it is not asserted}, or {"column": c, "in": [values], "reason": ...}, or {"column": c, "matches": a regular expression, "reason": ...}. A row that fails it is kept as a recorded gap of its subject, not dropped, and needs no object.
  "split_by": null, or a column whose values need different mappings
  "by_value": only with split_by: {value: {any of subject, object, predicate, association_category, knowledge_level, agent_type, primary_knowledge_source, qualifiers, edge_properties, keep_if, given as overrides}, or {"skip": reason} to leave that value's rows out}
  "notes": up to 600 characters on the choices a reader might question

"columns" lists the columns to try in order; the first that gives an identifier is used.

RESOLVERS (code turns a cell into an identifier; you choose which):
RESOLVER_MENU

Rules:
- Use only Biolink classes, predicates and qualifiers that exist in Biolink 4.4.3. The subject and object categories must fit the predicate's domain and range and the association class.
- Choose the most specific predicate the source supports, and no more specific than that.
- The knowledge level and agent type describe how the source made the claim.
- Do not invent columns or values."""


def brief() -> str:
    from kgx_gaps.mapping import AGENT_TYPES, KNOWLEDGE_LEVELS
    B = biolink()
    levels = sorted(set(B["enums"]["KnowledgeLevelEnum"]) & KNOWLEDGE_LEVELS)
    agents = sorted(set(B["enums"]["AgentTypeEnum"]) & AGENT_TYPES)
    menu = "\n".join(f"  {k}: {v}" for k, v in RESOLVERS.items())
    return BRIEF.replace("KNOWLEDGE_LEVELS", ", ".join(levels)).replace("AGENT_TYPES", ", ".join(agents)).replace("RESOLVER_MENU", menu)


def card_text(c: dict) -> str:
    return "CARD\n" + json.dumps(c, indent=1, ensure_ascii=False)


# ---- Biolink ---------------------------------------------------------------------------------------

_B: dict | None = None


def biolink() -> dict:
    global _B
    if _B is None:
        if not BIOLINK.exists():
            raise SystemExit(f"no {BIOLINK.name}; run:  .venv-standards/bin/python -I src/pharos_gap_kg/biolinkslice.py")
        _B = json.loads(BIOLINK.read_text())
    return _B


def class_line(c: str) -> list[str]:
    """c, then its ancestors through is_a and mixins, nearest first."""
    C = biolink()["classes"]
    out, todo = [], [c]
    while todo:
        x = todo.pop(0)
        if x in out or x not in C:
            continue
        out.append(x)
        todo += [m for m in [C[x]["is_a"]] + C[x]["mixins"] if m]
    return out


def is_a_line(c: str) -> list[str]:
    C = biolink()["classes"]
    out = []
    while c and c in C:
        out.append(c)
        c = C[c]["is_a"]
    return out


def pred_line(p: str) -> list[str]:
    P = biolink()["predicates"]
    out = []
    while p and p in P:
        out.append(p)
        p = P[p]["is_a"]
    return out


def admits(allowed: str | None, category: str) -> bool:
    if not allowed or not str(allowed).startswith("biolink:"):
        return True
    return allowed in class_line(category)


def declared(p: str, key: str) -> str | None:
    for q in pred_line(p):
        v = biolink()["predicates"][q][key]
        if v:
            return v
    return None


def usage(assoc: str, slot: str) -> dict:
    C = biolink()["classes"]
    for c in class_line(assoc):
        u = C[c]["slot_usage"].get(slot)
        if u:
            return u
    return {}


def carries(assoc: str, slot: str) -> bool:
    C = biolink()["classes"]
    return any(slot in C[c]["slots"] or slot in C[c]["slot_usage"] for c in class_line(assoc))


def prefixes_for(category: str) -> set[str]:
    C = biolink()["classes"]
    own = C.get(category, {}).get("id_prefixes") or []
    if own:
        return set(own)
    return {p for k, v in C.items() if category in is_a_line(k) for p in v["id_prefixes"]}


# ---- proposals -------------------------------------------------------------------------------------

def parse(text: str) -> dict | None:
    """The proposal: the last top-level JSON object in the reply that parses. An object nested in
    another is never taken for the proposal."""
    t = re.sub(r"```(?:json)?", "", text or "")
    dec = json.JSONDecoder()
    best, i = None, 0
    while True:
        i = t.find("{", i)
        if i < 0:
            return best
        try:
            obj, end = dec.raw_decode(t, i)
        except json.JSONDecodeError:
            i += 1
            continue
        if isinstance(obj, dict):
            best = obj
        i = end


def rules_of(p: dict) -> list[tuple[str | None, dict | None, str | None]]:
    """(split value or None, the rule with overrides applied, or None and the reason a value is left out)."""
    base = {k: p.get(k) for k in RULE_KEYS}
    if not p.get("split_by"):
        return [(None, base, None)]
    out = []
    for v, o in sorted((p.get("by_value") or {}).items()):
        o = o or {}
        if "skip" in o:
            out.append((v, None, str(o["skip"])))
            continue
        r = dict(base)
        for k, val in o.items():
            if k in ("subject", "object") and isinstance(val, dict):
                r[k] = {**(base.get(k) or {}), **val}
            else:
                r[k] = val
        out.append((v, r, None))
    return out


# ---- resolving identifiers -------------------------------------------------------------------------

class Refs:
    """The release's own lists: protein rows, target rows and vocabularies."""

    def __init__(self, release: str):
        pr = table(release, "protein")[["id", "uniprot", "geneid", "sym"]]
        self.uniprot = dict(zip(pr.id.astype(int), "UniProtKB:" + pr.uniprot.astype(str)))
        self.gene = {int(i): f"NCBIGene:{int(g)}" for i, g in zip(pr.id, pr.geneid) if pd.notna(g)}
        self.symbol = {**{f"UniProtKB:{u}": s for u, s in zip(pr.uniprot, pr.sym)},
                       **{f"NCBIGene:{int(g)}": s for g, s in zip(pr.geneid, pr.sym) if pd.notna(g)}}
        tt = table(release, "t2tc")
        self.target = {int(t): self.uniprot.get(int(p)) for t, p in zip(tt.target_id, tt.protein_id) if pd.notna(p)}
        self.vocab = {"UniProtKB": set(self.uniprot.values()), "NCBIGene": set(self.gene.values()),
                      "MONDO": set(table(release, "mondo").mondoid), "DOID": set(table(release, "do").doid),
                      "UBERON": set(table(release, "uberon").uid), "GO": set(table(release, "go").go_id)}
        self.xref = set(table(release, "mondo_xref").xref.astype(str))


def _int(v) -> int | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return int(f) if np.isfinite(f) and f == int(f) else None


def resolver_ok(spec) -> str:
    """'' if the resolver is one of the menu's, else why not."""
    if not isinstance(spec, str):
        return f"resolver {spec!r} is not a string"
    if spec in ("pharos_protein", "pharos_protein_gene", "pharos_target", "curie"):
        return ""
    if spec.startswith("prefix:") and re.fullmatch(r"prefix:[A-Za-z][\w.]*", spec):
        return ""
    if spec.startswith("name:") and re.fullmatch(r"name:[a-z][a-z0-9_]*", spec):
        return ""
    if spec.startswith("pattern:"):
        parts = spec.split(":", 2)
        if len(parts) == 3 and re.fullmatch(r"[A-Za-z][\w.]*", parts[1]):
            try:
                if re.compile(parts[2]).groups >= 1:
                    return ""
                return f"the regular expression in {spec!r} has no group"
            except re.error as e:
                return f"the regular expression in {spec!r} does not compile: {e}"
    return f"resolver {spec!r} is not on the menu"


def resolve_one(v, spec: str, refs: Refs) -> tuple[str | None, bool]:
    """(identifier or None, grounded)."""
    if v is None or (isinstance(v, float) and np.isnan(v)) or str(v).strip() == "":
        return None, True
    if spec == "pharos_protein":
        i = _int(v)
        return (refs.uniprot.get(i) if i is not None else None), True
    if spec == "pharos_protein_gene":
        i = _int(v)
        return (refs.gene.get(i) if i is not None else None), True
    if spec == "pharos_target":
        i = _int(v)
        return (refs.target.get(i) if i is not None else None), True
    s = str(v).strip()
    if spec == "curie":
        if re.fullmatch(r"[A-Za-z][\w.]*:\S+", s):
            return s, True
        m = re.fullmatch(r"([A-Za-z][A-Za-z.]*)_(\S+)", s)
        return (f"{m.group(1)}:{m.group(2)}" if m else None), True
    if spec.startswith("prefix:"):
        i = _int(v)
        return f"{spec[7:]}:{i if i is not None else s}", True
    if spec.startswith("pattern:"):
        _, p, rx = spec.split(":", 2)
        m = re.search(rx, s)
        return (f"{p}:{m.group(1)}" if m and m.group(1) else None), True
    if spec.startswith("name:"):
        return f"{PREFIX}:{spec[5:]}_{hashlib.sha256(' '.join(s.casefold().split()).encode()).hexdigest()[:12]}", False
    return None, True


def resolve(rows: pd.DataFrame, side: dict, refs: Refs) -> tuple[list, list]:
    ids, grounded = [], []
    cols = side.get("columns") or []
    for r in rows[cols].itertuples(index=False):
        got, g = None, True
        for v in r:
            got, g = resolve_one(v, side["resolver"], refs)
            if got:
                break
        ids.append(got)
        grounded.append(g)
    return ids, grounded


# ---- the validator ---------------------------------------------------------------------------------

def check_rule(r: dict, cols: set[str]) -> dict[str, list[str]]:
    """Checks 1 to 7 for one rule. name -> problems (empty when it passes)."""
    from kgx_gaps.mapping import AGENT_TYPES, KNOWLEDGE_LEVELS
    B = biolink()
    C, P, Q, E = B["classes"], B["predicates"], B["qualifiers"], B["enums"]
    out = {k: [] for k in ("1_parses", "2_predicate", "3_categories", "4_domain_range", "5_association", "6_qualifiers", "7_provenance")}
    for side in ("subject", "object"):
        s = r.get(side)
        if not isinstance(s, dict) or set(s) - SIDE_KEYS or not s.get("columns") or not s.get("resolver") or not s.get("category"):
            out["1_parses"].append(f"{side}: needs columns, resolver and category, and no other keys")
            continue
        if not isinstance(s["columns"], list) or not all(isinstance(c, str) for c in s["columns"]):
            out["1_parses"].append(f"{side}: columns must be a list of column names")
            continue
        miss = [c for c in s["columns"] + ([s["name_column"]] if s.get("name_column") else []) if c not in cols]
        if miss:
            out["1_parses"].append(f"{side}: no column {miss}")
        why = resolver_ok(s["resolver"])
        if why:
            out["1_parses"].append(f"{side}: {why}")
    for c in r.get("edge_properties") or []:
        if c not in cols:
            out["1_parses"].append(f"edge_properties: no column {c!r}")
    k = r.get("keep_if")
    if k:
        if not isinstance(k, dict) or k.get("column") not in cols or sum(x in k for x in ("at_least", "in", "matches")) != 1 or not str(k.get("reason") or "").strip():
            out["1_parses"].append("keep_if needs an existing column, one of at_least, in or matches, and a reason")
        elif "matches" in k:
            try:
                re.compile(str(k["matches"]))
            except re.error as e:
                out["1_parses"].append(f"keep_if: the regular expression does not compile: {e}")
        elif "in" in k and not isinstance(k["in"], list):
            out["1_parses"].append("keep_if: in takes a list")
    pred = r.get("predicate")
    if pred not in P:
        out["2_predicate"].append(f"{pred!r} is not a Biolink 4.4.3 predicate")
    elif P[pred]["abstract"] or P[pred]["mixin"] or P[pred]["deprecated"]:
        out["2_predicate"].append(f"{pred} is abstract, a mixin or deprecated")
    assoc = r.get("association_category")
    if assoc not in C:
        out["3_categories"].append(f"association category {assoc!r} is not a Biolink 4.4.3 class")
    elif "biolink:Association" not in is_a_line(assoc) or C[assoc]["abstract"] or C[assoc]["mixin"] or C[assoc]["deprecated"]:
        out["3_categories"].append(f"{assoc} is not a concrete class below biolink:Association")
    cats = {}
    for side in ("subject", "object"):
        cat = (r.get(side) or {}).get("category") if isinstance(r.get(side), dict) else None
        cats[side] = cat
        if cat not in C:
            out["3_categories"].append(f"{side} category {cat!r} is not a Biolink 4.4.3 class")
        elif "biolink:NamedThing" not in is_a_line(cat) or C[cat]["abstract"] or C[cat]["mixin"] or C[cat]["deprecated"]:
            out["3_categories"].append(f"{side} category {cat} is not a concrete class below biolink:NamedThing")
    if pred in P and all(c in C for c in cats.values()):
        d, g = declared(pred, "domain"), declared(pred, "range")
        if not admits(d, cats["subject"]):
            out["4_domain_range"].append(f"{pred}'s domain is {d}, which does not admit {cats['subject']}")
        if not admits(g, cats["object"]):
            out["4_domain_range"].append(f"{pred}'s range is {g}, which does not admit {cats['object']}")
    if assoc in C and pred in P and all(c in C for c in cats.values()):
        for side in ("subject", "object"):
            want = usage(assoc, side).get("range")
            if want and not admits(want, cats[side]):
                out["5_association"].append(f"{assoc} takes a {side} that is a {want}, not {cats[side]}")
        u = usage(assoc, "predicate")
        if u.get("range") in E and not set(pred_line(pred)) & set(E[u["range"]]):
            out["5_association"].append(f"{assoc} takes a predicate from {u['range']} ({', '.join(E[u['range']])}); {pred} is not one of them or below one")
        if u.get("subproperty_of") and u["subproperty_of"] not in pred_line(pred):
            out["5_association"].append(f"{assoc} takes a predicate at or below {u['subproperty_of']}; {pred} is not")
    quals = r.get("qualifiers") or {}
    if not isinstance(quals, dict):
        out["6_qualifiers"].append("qualifiers must be an object")
        quals = {}
    for q, v in quals.items():
        if q not in Q:
            out["6_qualifiers"].append(f"{q!r} is not a Biolink qualifier")
            continue
        if assoc in C and not carries(assoc, q):
            out["6_qualifiers"].append(f"{assoc} does not carry {q}")
        if isinstance(v, dict):
            if set(v) != {"sign_of"} or "direction" not in q or v["sign_of"] not in cols:
                out["6_qualifiers"].append(f"{q}: a value from a column is only {{\"sign_of\": column}} on a direction qualifier")
            continue
        rng = (usage(assoc, q).get("range") if assoc in C else None) or Q[q]["range"]
        if rng in E and str(v) not in E[rng]:
            out["6_qualifiers"].append(f"{q}: {v!r} is not in {rng}")
    lv, ag = r.get("knowledge_level"), r.get("agent_type")
    if lv not in set(B["enums"]["KnowledgeLevelEnum"]) & KNOWLEDGE_LEVELS:
        out["7_provenance"].append(f"knowledge level {lv!r} is not in both Biolink's and kgx-gaps' enumerations")
    if ag not in set(B["enums"]["AgentTypeEnum"]) & AGENT_TYPES:
        out["7_provenance"].append(f"agent type {ag!r} is not in both Biolink's and kgx-gaps' enumerations")
    if not re.fullmatch(r"infores:[a-z0-9][a-z0-9._-]*", str(r.get("primary_knowledge_source") or "")):
        out["7_provenance"].append(f"knowledge source {r.get('primary_knowledge_source')!r} is not written infores:<lower-case name>")
    return out


def identifiers(ids: list, category: str, refs: Refs) -> dict:
    """Check 8's three figures for one side, on the rows of one rule."""
    got = [i for i in ids if i]
    distinct = sorted(set(got))
    pre = Counter(i.split(":", 1)[0] for i in distinct)
    listed = prefixes_for(category) if category else set()
    off = {p: n for p, n in pre.items() if p not in listed and p != PREFIX}
    vocab = {}
    for p in VOCAB_PREFIXES:
        mine = [i for i in distinct if i.startswith(p + ":")]
        if mine:
            vocab[p] = {"distinct": len(mine), "in_release": sum(i in refs.vocab[p] for i in mine)}
    xref = {}
    for p in XREF_PREFIXES:
        mine = [i for i in distinct if i.startswith(p + ":")]
        if mine:
            xref[p] = {"distinct": len(mine), "among_mondo_xrefs": sum(i in refs.xref for i in mine)}
    return {"rows": len(ids), "resolved": len(got), "distinct": len(distinct), "prefixes": dict(pre.most_common()),
            "off_prefix": off, "off_prefix_share": round(sum(off.values()) / len(distinct), 4) if distinct else 0.0,
            "vocabulary": vocab, "mondo_xrefs": xref}


def keeps(k: dict | None, v) -> bool:
    """Whether a row meets the proposal's condition. A row with no value does not."""
    if not k:
        return True
    if v is None or (isinstance(v, float) and np.isnan(v)) or str(v).strip() == "":
        return False
    if "in" in k:
        return cut(v, 10**6) in [str(x) for x in k["in"]]
    if "matches" in k:
        return re.search(str(k["matches"]), str(v)) is not None
    try:
        return float(v) >= float(k["at_least"])
    except (TypeError, ValueError):
        return False


def check_ids(sides: dict) -> list[str]:
    bad = []
    rows = sides["subject"]["rows"]
    both = sides["_both"]
    if rows and both / rows < 0.95:
        bad.append(f"{both:,} of {rows:,} rows ({both / rows:.1%}) resolve as they must (both sides, or the subject of a row "
                   f"the condition makes a gap); at least 95% must")
    for side in ("subject", "object"):
        s = sides[side]
        if s["off_prefix_share"] > 0.01:
            bad.append(f"{side}: {sum(s['off_prefix'].values()):,} of {s['distinct']:,} distinct identifiers are under prefixes "
                       f"{sorted(s['off_prefix'])} that {s['category']}'s id_prefixes do not list")
        for p, v in s["vocabulary"].items():
            if v["in_release"] / v["distinct"] < 0.90:
                bad.append(f"{side}: {v['in_release']:,} of {v['distinct']:,} {p} identifiers are in the release's {p} list; at least 90% must")
    return bad


def validate(proposal: dict, rows: pd.DataFrame, refs: Refs, export: bool = True, view: str = "view") -> dict:
    """All nine checks. {"accepted", "problems": [..], "rules": [...], "export": {...}}."""
    cols = set(rows.columns)
    if not isinstance(proposal, dict):
        return {"accepted": False, "problems": ["1_parses: no JSON object"], "rules": []}
    unknown = set(proposal) - KEYS
    if proposal.get("map") is False:
        return {"accepted": True, "map": False, "problems": [], "rules": [], "reason": proposal.get("reason")}
    problems = [f"1_parses: unknown key {k!r}" for k in sorted(unknown)]
    split = proposal.get("split_by")
    if split:
        if split not in cols:
            problems.append(f"1_parses: split_by names no column {split!r}")
        else:
            have = set(rows[split].dropna().astype(str))
            named = set(map(str, (proposal.get("by_value") or {})))
            problems += [f"1_parses: split value {v!r} is not in the view" for v in sorted(named - have)]
            problems += [f"1_parses: the view's value {v!r} has no mapping and no skip" for v in sorted(have - named)]
    elif proposal.get("by_value"):
        problems.append("1_parses: by_value without split_by")
    rules = []
    for value, r, skip in rules_of(proposal):
        if r is None:
            rules.append({"value": value, "skip": skip})
            continue
        mine = rows if value is None else rows[rows[split].astype(str) == str(value)]
        res = {"value": value, "rows": int(len(mine)), "checks": check_rule(r, cols)}
        if not any(res["checks"]["1_parses"]):
            sides = {}
            k = r.get("keep_if")
            kept = [keeps(k, v) for v in mine[k["column"]]] if k else [True] * len(mine)
            for side in ("subject", "object"):
                ids, grounded = resolve(mine, r[side], refs)
                # An object becomes a node only on a kept row: a row the condition makes a gap names its object as a property.
                judged = ids if side == "subject" else [i for i, ok in zip(ids, kept) if ok]
                sides[side] = {**identifiers(judged, r[side]["category"], refs), "category": r[side]["category"]}
                sides[side + "_ids"], sides[side + "_grounded"] = ids, grounded
            sides["_both"] = sum(bool(a and (b or not ok)) for a, b, ok in zip(sides["subject_ids"], sides["object_ids"], kept))
            res["checks"]["8_identifiers"] = check_ids(sides)
            res["identifiers"] = {s: {k2: v for k2, v in sides[s].items()} for s in ("subject", "object")}
            res["resolved_as_needed"], res["kept"], res["made_gaps_by_condition"] = sides["_both"], sum(kept), len(kept) - sum(kept)
            res["_ids"] = {s: (sides[s + "_ids"], sides[s + "_grounded"]) for s in ("subject", "object")}
            res["_kept"] = kept
        else:
            res["checks"]["8_identifiers"] = ["not reached: the proposal does not parse"]
        rules.append(res)
        for k, v in res["checks"].items():
            problems += [f"{k}{'' if value is None else f' [{value}]'}: {x}" for x in v]
    out = {"accepted": not problems, "map": True, "problems": problems,
           "rules": [{k: v for k, v in r.items() if not k.startswith("_")} for r in rules]}
    if not problems and export:
        try:
            ex = build_export(proposal, rows, rules, refs, view)
        except Exception as e:     # a mapping the exporter cannot write is refused, with what it said
            ex = {"passed": False, "error": f"{type(e).__name__}: {e}"}
        out["export"] = {k: v for k, v in ex.items() if k not in ("nodes", "edges", "gaps")}
        if not ex.get("passed"):
            out["accepted"] = False
            out["problems"].append(f"9_conformance: {ex.get('error') or ex.get('failed')}")
        out["_frames"] = ex
    return out


# ---- the export ------------------------------------------------------------------------------------

def _name(i: str, label, refs: Refs) -> str:
    for v in (label, refs.symbol.get(i)):
        if v is not None and not (isinstance(v, float) and np.isnan(v)) and str(v).strip():
            return " ".join(str(v).split())
    return ""


def build_export(proposal: dict, rows: pd.DataFrame, rules: list[dict], refs: Refs, view: str = "view") -> dict:
    from kgx_gaps import Mapping, Rule, conformance
    from kgx_gaps.export import Exporter
    split = proposal.get("split_by")
    by_ks: dict[str, list] = {}
    for (value, r, skip), res in zip(rules_of(proposal), rules):
        if r is not None:
            by_ks.setdefault(r["primary_knowledge_source"], []).append((value, r, res))
    N, E, G, n_in, merged, failed = [], [], [], 0, 0, []
    for ks, items in sorted(by_ks.items()):
        qcols = sorted({q for _, r, _ in items for q in (r.get("qualifiers") or {})})
        stash = {q: (f"{PREFIX}:q_{q}" if q in NATIVE_QUALIFIERS else q) for q in qcols}
        pcols = sorted({f"{PREFIX}:{c}" for _, r, _ in items for c in (r.get("edge_properties") or [])})
        extra = ["aggregator_knowledge_source"] + [stash[q] for q in qcols] + pcols
        m = Mapping(prefix=PREFIX, knowledge_source=ks, edge_identity_columns=[stash[q] for q in qcols], gap_identity_columns=[WITHHELD])
        for value, r, _ in items:
            et = view if value is None else f"{view}:{value}"
            k = r.get("keep_if")
            m.rules[et] = Rule(predicate=r["predicate"], category=r["association_category"], knowledge_level=r["knowledge_level"],
                               agent_type=r["agent_type"], precondition=(lambda row, why=str(k["reason"]): "" if row["_kept"] else why) if k else None)
        bad = m.validate()
        if bad:
            raise ValueError("; ".join(bad))
        cat = {}

        def resolver(label, cat=cat):
            return label, cat[label][0], cat[label][1], cat[label][2]

        def extras(row, kind):
            if kind != "edge":
                return {"aggregator_knowledge_source": graph.AGGREGATOR, WITHHELD: row.get("object") or ""}
            return {"aggregator_knowledge_source": graph.AGGREGATOR, **{stash[q]: row["_q"].get(q, "") for q in qcols},
                    **{c: row["_p"].get(c, "") for c in pcols}}
        ex = Exporter(mapping=m, resolver=resolver, extra_edge_columns=extra, extra_gap_columns=["aggregator_knowledge_source", WITHHELD], extras=extras)
        added = 0
        for value, r, res in items:
            et = view if value is None else f"{view}:{value}"
            mine = rows if value is None else rows[rows[split].astype(str) == str(value)]
            (sid, sg), (oid, og) = res["_ids"]["subject"], res["_ids"]["object"]
            sname = mine[r["subject"]["name_column"]].tolist() if r["subject"].get("name_column") else [None] * len(mine)
            oname = mine[r["object"]["name_column"]].tolist() if r["object"].get("name_column") else [None] * len(mine)
            quals = r.get("qualifiers") or {}
            props = r.get("edge_properties") or []
            kept = res["_kept"]
            recs = {}
            for j, row in enumerate(mine.to_dict("records")):
                if not sid[j] or (kept[j] and not oid[j]):
                    continue
                cat.setdefault(sid[j], (r["subject"]["category"], sg[j], _name(sid[j], sname[j], refs)))
                if oid[j]:
                    cat.setdefault(oid[j], (r["object"]["category"], og[j], _name(oid[j], oname[j], refs)))
                q = {}
                for name, v in quals.items():
                    if isinstance(v, dict):
                        x = row.get(v["sign_of"])
                        q[name] = "" if x is None or pd.isna(x) or float(x) == 0 else ("increased" if float(x) > 0 else "decreased")
                    else:
                        q[name] = str(v)
                key = (sid[j], oid[j], tuple(sorted(q.items())), kept[j])
                g = recs.setdefault(key, {"type": et, "subject": sid[j], "object": oid[j] or "", "_q": q, "_kept": kept[j], "_p": {}})
                for c in props:
                    v = row.get(c)
                    if v is not None and not (isinstance(v, float) and np.isnan(v)) and str(v).strip():
                        g["_p"].setdefault(f"{PREFIX}:{c}", set()).add(str(int(v)) if isinstance(v, float) and v == int(v) else str(v))
            for g in recs.values():
                g["_p"] = {c: "|".join(sorted(v)) for c, v in g["_p"].items()}
                ex.add(g)
                added += 1
        n_in += added
        n, e, gp = ex.frames()
        for q in qcols:
            if stash[q] != q:
                e[q] = e.pop(stash[q]) if stash[q] in e else ""
        rep = conformance.check(n, e, gp, m, n_input_rows=added, n_merged=ex.merged)
        if not rep.passed:
            failed.append(f"{ks}: " + "; ".join(f"{c.id} {c.title}: {c.detail}" for c in rep.checks if not c.passed or c.skipped))
        N.append(n), E.append(e), G.append(gp)
        merged += ex.merged
    nodes = pd.concat(N, ignore_index=True).drop_duplicates("id") if N else pd.DataFrame()
    edges = pd.concat(E, ignore_index=True) if E else pd.DataFrame()
    gaps = pd.concat(G, ignore_index=True) if G else pd.DataFrame()
    return {"passed": not failed, "failed": "; ".join(failed) or None, "groups": n_in, "nodes": nodes, "edges": edges, "gaps": gaps,
            "n_nodes": int(len(nodes)), "n_edges": int(len(edges)), "n_gaps": int(len(gaps)),
            "predicates": dict(Counter(edges.predicate)) if len(edges) else {}}


# ---- the hand mappings, in the same form -----------------------------------------------------------

def hand(view: str) -> dict:
    """The mapping the graph's own code applies to this view, written as a proposal."""
    from pharos_gap_kg import ctxgraph, kg
    protein = {"columns": ["protein_id"], "resolver": "pharos_protein", "category": "biolink:Protein"}
    disease = {"columns": ["mondoid", "did"], "resolver": "curie", "category": "biolink:Disease", "name_column": "name"}

    def source(s: str) -> dict:
        ks, pred, level, agent = ctxgraph.SOURCES[s]
        return {"predicate": pred, "association_category": ctxgraph.ASSOC, "knowledge_level": level, "agent_type": agent, "primary_knowledge_source": ks}
    base = {"map": True, "subject": protein, "object": disease, "qualifiers": {}, "edge_properties": [], "keep_if": None, "split_by": None}
    if view == "disease.Expression_Atlas":
        r = kg.expression_mapping().rules[kg.EXPRESSION]
        return {**base, "predicate": r.predicate, "association_category": r.category, "knowledge_level": r.knowledge_level,
                "agent_type": r.agent_type, "primary_knowledge_source": kg.expression_mapping().knowledge_source,
                "qualifiers": {"subject_aspect_qualifier": r.subject_aspect_qualifier, kg.DIRECTION: {"sign_of": "log2foldchange"}},
                "edge_properties": ["pvalue"]}
    if view == "disease.JensenLab":
        vals = ["JensenLab Knowledge AmyCo", "JensenLab Knowledge MedlinePlus", "JensenLab Knowledge UniProtKB-KW", "JensenLab Text Mining"]
        return {**base, **source(vals[0]), "split_by": "dtype",
                "by_value": {**{v: source(v) for v in vals}, "JensenLab Experiment TIGA": {"skip": "the GWAS data again"}}}
    if view.startswith("disease."):
        return {**base, **source({"disease.UniProt_Disease": "UniProt Disease"}.get(view, view.split(".", 1)[1]))}
    if view == "tiga":
        m = graph.mapping()
        r = m.rules["gwas_condition"]
        return {**base, "object": {"columns": ["efoid"], "resolver": "curie", "category": graph.TRAIT_CATEGORY["EFO"]},
                "predicate": r.predicate, "association_category": r.category, "knowledge_level": r.knowledge_level,
                "agent_type": r.agent_type, "primary_knowledge_source": m.knowledge_source}
    if view == "drug_activity":
        m = kg.drug_mapping()
        r = m.rules[kg.DRUG_TARGET]
        return {**base, "subject": {"columns": ["dcid"], "resolver": "prefix:DrugCentral", "category": "biolink:Drug", "name_column": "drug"},
                "object": {"columns": ["target_id"], "resolver": "pharos_target", "category": "biolink:Protein"},
                "predicate": r.predicate, "association_category": r.category, "knowledge_level": r.knowledge_level,
                "agent_type": r.agent_type, "primary_knowledge_source": m.knowledge_source}
    raise KeyError(view)


# The graph's own name for each hand rule's edges (pharosgapkg:source_edge_type), for the cross-check with the export.
EXPORT_TYPE = {"tiga": "gwas_condition", "drug_activity": "drug_target"}
FIELDS = ["subject_column", "object_column", "subject_category", "object_category", "predicate", "association_category",
          "knowledge_level", "agent_type", "primary_knowledge_source"]


def hand_rules() -> list[tuple[str, str | None, dict]]:
    out = []
    for v in SCORED:
        for value, r, skip in rules_of(hand(v)):
            if r is not None:
                out.append((v, value, r))
    return out


def fields(r: dict) -> dict:
    return {"subject_column": (r.get("subject") or {}).get("columns"), "object_column": (r.get("object") or {}).get("columns"),
            "subject_category": (r.get("subject") or {}).get("category"), "object_category": (r.get("object") or {}).get("category"),
            **{k: r.get(k) for k in FIELDS[4:]}}


def export_values(release: str) -> dict:
    """For each hand rule, the value most of its exported edges carry, read from the graph's files."""
    g = paths.GRAPH / release
    e = pd.read_csv(g / "edges.tsv", sep="\t", dtype=str, usecols=["subject", "object", "predicate", "category", "primary_knowledge_source",
                                                                   "knowledge_level", "agent_type", f"{PREFIX}:source_edge_type"])
    n = pd.read_csv(g / "nodes.tsv", sep="\t", dtype=str, usecols=["id", "category"])
    cat = dict(zip(n.id, n.category))
    out = {}
    for v, value, r in hand_rules():
        et = value or EXPORT_TYPE.get(v) or {"disease.UniProt_Disease": "UniProt Disease"}.get(v, v.split(".", 1)[1].replace("_", " "))
        mine = e[e[f"{PREFIX}:source_edge_type"] == et]
        mode = lambda s: s.value_counts().index[0] if len(s) else None
        out[f"{v}|{value or ''}"] = {"edges": int(len(mine)), "subject_category": mode(mine.subject.map(cat)), "object_category": mode(mine.object.map(cat)),
                                     **{k: mode(mine[k]) for k in ("predicate", "category", "primary_knowledge_source", "knowledge_level", "agent_type")}}
    return out


def agree(model_rule: dict | None, hand_rule: dict) -> dict:
    """Field by field. The predicate is also 'close' when one is a Biolink ancestor of the other."""
    if model_rule is None:
        return {k: False for k in FIELDS} | {"predicate_close": False}
    m, h = fields(model_rule), fields(hand_rule)
    out = {}
    for k in FIELDS:
        if k.endswith("_column"):
            out[k] = bool(m[k]) and isinstance(m[k], list) and m[k][0] in (h[k] or [])
        else:
            out[k] = m[k] == h[k]
    p, q = m["predicate"], h["predicate"]
    out["predicate_close"] = out["predicate"] or (p in biolink()["predicates"] and (q in pred_line(p) or p in pred_line(q)))
    return out


def compared(proposal: dict | None, view: str, rows: pd.DataFrame) -> list[dict]:
    """The model's rule set against each hand rule of a view."""
    hand_r = [(value, r) for value, r, skip in rules_of(hand(view)) if r is not None]
    if not proposal or proposal.get("map") is False:
        return [{"value": value, "agree": agree(None, r)} for value, r in hand_r]
    mine = rules_of(proposal)
    split = proposal.get("split_by")
    hand_split = hand(view).get("split_by")
    out = []
    for value, r in hand_r:
        if value is None or not split:
            if split:   # the model split a view the hand mapping does not: its rule over the most rows
                sizes = {v: int((rows[split].astype(str) == str(v)).sum()) for v, rr, s in mine if rr is not None}
                pick = max(sizes, key=sizes.get) if sizes else None
                rr = next((rr for v, rr, s in mine if v == pick), None)
            else:
                rr = mine[0][1]
        elif split == hand_split:
            rr = next((rr for v, rr, s in mine if str(v) == str(value)), None)
        else:
            sub = rows[rows[hand_split] == value]
            sizes = {v: int((sub[split].astype(str) == str(v)).sum()) for v, rr, s in mine if rr is not None}
            pick = max(sizes, key=sizes.get) if sizes else None
            rr = next((rr for v, rr, s in mine if v == pick), None)
        out.append({"value": value, "agree": agree(rr, r)})
    return out


# ---- asking a model --------------------------------------------------------------------------------

def spent(release: str) -> float:
    """What the cost lines say was spent, plus the ceiling of any call reserved and never settled:
    a call whose cost was not recorded may still have been billed."""
    path = ledger_path(release)
    if not path.exists():
        return 0.0
    total, open_ = 0.0, []
    for line in path.read_text().splitlines():
        r = json.loads(line)
        if r.get("kind") == "reserve":
            open_.append(r)
        elif r.get("kind") == "cost":
            total += float(r["usd"])
            i = next((k for k in range(len(open_) - 1, -1, -1) if (open_[k]["model"], open_[k]["view"]) == (r["model"], r["view"])), None)
            if i is not None:
                open_.pop(i)
    return total + sum(float(r["usd_ceiling"]) for r in open_)


def _ledger(release: str, rec: dict) -> None:
    path = ledger_path(release)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(json.dumps({**rec, "at": time.strftime("%Y-%m-%dT%H:%M:%S")}) + "\n")


def _post(url: str, body: dict, headers: dict, timeout: int) -> dict:
    req = urllib.request.Request(url, data=json.dumps(body).encode(), headers={"Content-Type": "application/json", **headers})
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}: {e.read()[:300].decode('utf-8', 'replace')}"
            if e.code in (400, 401, 402, 403, 404):
                break
        except Exception as e:     # timeouts and resets are retried
            last = repr(e)
        time.sleep(3 * (attempt + 1))
    return {"error": last}


def ask(release: str, model: str, messages: list[dict], view: str) -> dict:
    t0 = time.time()
    if model.startswith("local:"):
        body = {"model": model[6:], "messages": messages, "stream": False, "think": True,
                "options": {"temperature": 0.6, "top_p": 0.95, "seed": SEED, "num_ctx": 32768, "num_predict": 24000}}
        r = _post("http://localhost:11434/api/chat", body, {}, timeout=3600)
        if "error" in r:
            return {"error": r["error"], "seconds": round(time.time() - t0, 1)}
        msg = r.get("message") or {}
        return {"text": msg.get("content") or "", "thinking_chars": len(msg.get("thinking") or ""),
                "usage": {"in": r.get("prompt_eval_count"), "out": r.get("eval_count"), "usd": 0.0}, "seconds": round(time.time() - t0, 1)}
    key = os.environ.get("OPENROUTER_API_KEY")
    if not key:
        raise SystemExit("set OPENROUTER_API_KEY")
    pin, pout = PRICES[model]
    estimate = (sum(len(m["content"]) for m in messages) / 3.0 * pin + MAX_TOKENS * pout) / 1e6
    if spent(release) + estimate > CAP_USD:
        raise SystemExit(f"STOPPED: spend {spent(release):.4f} + this call's ceiling {estimate:.4f} would pass the cap of {CAP_USD:.2f}")
    _ledger(release, {"kind": "reserve", "model": model, "view": view, "usd_ceiling": round(estimate, 6)})
    body = {"model": model, "messages": messages, "max_tokens": MAX_TOKENS, "temperature": 0, "usage": {"include": True}}
    r = _post("https://openrouter.ai/api/v1/chat/completions", body, {"Authorization": f"Bearer {key}", "X-Title": "pharos-gap-kg ingest"}, timeout=600)
    u = r.get("usage") or {}
    usd = float(u.get("cost") or 0.0)
    _ledger(release, {"kind": "cost", "model": model, "view": view, "usd": usd, "in": u.get("prompt_tokens"), "out": u.get("completion_tokens"),
                      "provider": r.get("provider"), "error": r.get("error")})
    if "error" in r or not r.get("choices"):
        return {"error": str(r.get("error"))[:300], "seconds": round(time.time() - t0, 1)}
    return {"text": r["choices"][0]["message"].get("content") or "", "usage": {"in": u.get("prompt_tokens"), "out": u.get("completion_tokens"), "usd": usd},
            "provider": r.get("provider"), "seconds": round(time.time() - t0, 1)}


def refusal(v: dict) -> str:
    return ("Code refused this proposal for these reasons:\n- " + "\n- ".join(v["problems"][:40]) +
            "\n\nReturn a corrected proposal: ONE JSON object in the same form, and nothing else.")


def propose(release: str, model: str, views: list[str], refs: Refs | None = None) -> dict:
    cards = frozen_cards(release)
    refs = refs or Refs(release)
    system = brief()
    done = {}
    for view in views:
        path = proposals_dir(release, model) / f"{view}.json"
        if path.exists():
            done[view] = json.loads(path.read_text())
            continue
        c = cards["cards"][view]
        rows, _ = view_rows(release, view)
        msgs = [{"role": "system", "content": system}, {"role": "user", "content": card_text(c)}]
        attempts = []
        for attempt in range(2):
            got = ask(release, model, msgs, view)
            if "error" in got:
                attempts.append({**got, "parsed": None, "validation": None})
                break
            p = parse(got["text"])
            v = validate(p, rows, refs, view=view) if p is not None else {"accepted": False, "problems": ["1_parses: no JSON object in the reply"], "rules": []}
            v.pop("_frames", None)
            attempts.append({**got, "parsed": p, "validation": v})
            if v["accepted"]:
                break
            msgs += [{"role": "assistant", "content": got["text"]}, {"role": "user", "content": refusal(v)}]
        rec = {"model": model, "view": view, "release": release, "card_sha256": sha(c), "brief_sha256": hashlib.sha256(system.encode()).hexdigest(),
               "attempts": attempts, "final": attempts[-1]["parsed"] if attempts else None,
               "accepted": bool(attempts and attempts[-1].get("validation") and attempts[-1]["validation"]["accepted"]),
               "accepted_first": bool(attempts and attempts[0].get("validation") and attempts[0]["validation"]["accepted"])}
        write_json(path, rec)
        done[view] = rec
        print(f"{model} {view}: {'accepted' if rec['accepted'] else 'refused'}"
              f"{'' if rec['accepted_first'] or not rec['accepted'] else ' after one revision'}"
              f" ({sum(a.get('seconds') or 0 for a in attempts):.0f} s)", flush=True)
    return done


def stale(release: str, refs: Refs | None = None) -> list[dict]:
    """Every recorded card whose verdicts the current validator would not give, attempt by attempt. A card
    refused for a reason the validator no longer gives was sent back with wrong reasons, so it is asked again."""
    refs = refs or Refs(release)
    root = data_dir(release) / f"proposals_{release}"
    out = []
    for d in sorted(root.iterdir()) if root.exists() else []:
        for f in sorted(d.glob("*.json")):
            rec = json.loads(f.read_text())
            rows, _ = view_rows(release, rec["view"])
            for i, a in enumerate(rec["attempts"]):
                if a.get("validation") is None:
                    continue
                now = validate(a["parsed"], rows, refs, view=rec["view"]) if a["parsed"] is not None else a["validation"]
                if now["accepted"] != a["validation"]["accepted"] or sorted(now["problems"]) != sorted(a["validation"]["problems"]):
                    out.append({"model": rec["model"], "view": rec["view"], "file": str(f.relative_to(paths.ROOT)), "attempt": i + 1,
                                "recorded": a["validation"]["problems"], "now": now["problems"]})
                    break
    return out


def set_aside(release: str, items: list[dict], why: str) -> None:
    """Move stale records to discarded_<release>, with the reason, so they are asked again."""
    for it in items:
        src = paths.ROOT / it["file"]
        dest = data_dir(release) / f"discarded_{release}" / src.parent.name / src.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        rec = json.loads(src.read_text())
        write_json(dest, {**rec, "set_aside": why, "set_aside_detail": it})
        src.unlink()


# ---- control, report, export -----------------------------------------------------------------------

def control(release: str, refs: Refs | None = None) -> dict:
    refs = refs or Refs(release)
    out = {}
    for v in SCORED:
        rows, _ = view_rows(release, v)
        val = validate(hand(v), rows, refs, view=v)
        val.pop("_frames", None)
        out[v] = val
    return out


def qualify(rows: list[dict]) -> dict:
    acc = sum(r["accepted"] for r in rows if r["view"] in SCORED)
    rules = [a for r in rows for a in r["agreement"]]
    cols = sum(a["agree"]["subject_column"] and a["agree"]["object_column"] for a in rules)
    pred = sum(a["agree"]["predicate_close"] for a in rules)
    return {"accepted_cards": acc, "of_cards": len(SCORED), "columns_agree": cols, "predicate_exact_or_close": pred,
            "predicate_exact": sum(a["agree"]["predicate"] for a in rules), "of_rules": len(rules),
            "qualifies": acc >= 8 and cols >= 11 and pred >= 9}


def report(release: str) -> dict:
    cards = frozen_cards(release)
    refs = Refs(release)
    ctrl = control(release, refs)
    ev = export_values(release)
    hand_in_export = {}
    for v, value, r in hand_rules():
        x = ev[f"{v}|{value or ''}"]
        f = fields(r)
        hand_in_export[f"{v}|{value or ''}"] = {k: f[k] == x[k2] for k, k2 in (("predicate", "predicate"), ("association_category", "category"),
                                                                              ("primary_knowledge_source", "primary_knowledge_source"),
                                                                              ("knowledge_level", "knowledge_level"), ("agent_type", "agent_type"),
                                                                              ("subject_category", "subject_category"), ("object_category", "object_category"))}
        hand_in_export[f"{v}|{value or ''}"]["edges"] = x["edges"]
    models = {}
    for model in MODELS:
        d = proposals_dir(release, model)
        if not d.exists():
            continue
        rows = []
        for view in SCORED + NEW:
            p = d / f"{view}.json"
            if not p.exists():
                continue
            rec = json.loads(p.read_text())
            vrows, _ = view_rows(release, view)
            row = {"view": view, "accepted": rec["accepted"], "accepted_first": rec["accepted_first"], "attempts": len(rec["attempts"]),
                   "seconds": round(sum(a.get("seconds") or 0 for a in rec["attempts"]), 1),
                   "usd": round(sum((a.get("usage") or {}).get("usd") or 0 for a in rec["attempts"]), 6),
                   "tokens_out": sum((a.get("usage") or {}).get("out") or 0 for a in rec["attempts"]),
                   "problems": (rec["attempts"][-1].get("validation") or {}).get("problems") if rec["attempts"] else ["no reply"],
                   "error": rec["attempts"][-1].get("error") if rec["attempts"] else None}
            if view in SCORED:
                row["agreement"] = compared(rec["final"] if rec["accepted"] or rec["final"] else None, view, vrows)
            rows.append(row)
        scored = [r for r in rows if r["view"] in SCORED]
        models[model] = {"cards": rows, "qualification": qualify(scored) if len(scored) == len(SCORED) else None,
                         "usd": round(sum(r["usd"] for r in rows), 6), "seconds": round(sum(r["seconds"] for r in rows), 1)}
    out = {"release": release, "protocol": "docs/PROTOCOL.md, Ingest", "biolink": biolink()["version"], "biolink_schema_sha256": biolink()["schema_sha256"],
           "cards_sha256": sha(cards["cards"]), "brief_sha256": cards["brief_sha256"],
           "control": {v: {"accepted": c["accepted"], "problems": c["problems"],
                           "export": c.get("export")} for v, c in ctrl.items()},
           "hand_values_in_export": hand_in_export, "models": models, "spend_usd": round(spent(release), 6), "cap_usd": CAP_USD}
    write_json(report_path(release), out)
    return out


def chosen(release: str, rep: dict) -> dict:
    """For each new view, the proposal used: the cheapest qualifying model's accepted one."""
    order = [m for m in MODELS if (rep["models"].get(m) or {}).get("qualification") and rep["models"][m]["qualification"]["qualifies"]]
    out = {}
    for view in NEW:
        pick = None
        for m in order:
            p = proposals_dir(release, m) / f"{view}.json"
            if p.exists() and json.loads(p.read_text())["accepted"]:
                pick = m
                break
        out[view] = pick
    return out


def export_new(release: str) -> dict:
    rep = json.loads(report_path(release).read_text())
    refs = Refs(release)
    out = {}
    for view, model in chosen(release, rep).items():
        if model is None:
            out[view] = {"exported": False, "why": "no qualifying model's proposal was accepted"}
            continue
        rec = json.loads((proposals_dir(release, model) / f"{view}.json").read_text())
        rows, _ = view_rows(release, view)
        v = validate(rec["final"], rows, refs, view=view)
        if not v["accepted"] or v.get("map") is False:
            out[view] = {"exported": False, "model": model, "why": v.get("reason") or v["problems"]}
            continue
        ex = v["_frames"]
        d = paths.GRAPH / release / "ingested" / view
        d.mkdir(parents=True, exist_ok=True)
        files = {}
        for name in ("nodes", "edges", "gaps"):
            p = d / f"{name}.tsv"
            tmp = p.with_suffix(".tsv.tmp")
            ex[name].to_csv(tmp, sep="\t", index=False)
            os.replace(tmp, p)
            files[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        out[view] = {"exported": True, "model": model, "dir": str(d.relative_to(paths.ROOT)), "nodes": ex["n_nodes"], "edges": ex["n_edges"],
                     "gaps": ex["n_gaps"], "predicates": ex["predicates"], "sha256": files}
    return out


# ---- selftest --------------------------------------------------------------------------------------

def _selftest() -> int:
    bad = []

    def ck(name: str, cond: bool) -> None:
        print(("ok   " if cond else "FAIL ") + name)
        if not cond:
            bad.append(name)
    ck("parse takes the last JSON object, through fences and prose",
       parse('Here:\n```json\n{"a": {"b": 1}}\n```\nor rather {"map": false, "reason": "x {y}"}') == {"map": False, "reason": "x {y}"})
    ck("parse never takes an object nested in the proposal for the proposal",
       parse('{"map": true, "subject": {"columns": ["a"]}, "qualifiers": {}, "by_value": {}}') == {"map": True, "subject": {"columns": ["a"]}, "qualifiers": {}, "by_value": {}}
       and parse('{"a": 1} then {"b": {"c": {}}}') == {"b": {"c": {}}})
    ck("a predicate inherits its domain from the nearest ancestor that declares one",
       declared("biolink:gene_associated_with_condition", "domain") == "biolink:Gene" and declared("biolink:same_as", "domain") == "biolink:NamedThing")
    ck("a protein is a gene or gene product through its mixins, and not a gene",
       admits("biolink:GeneOrGeneProduct", "biolink:Protein") and not admits("biolink:Gene", "biolink:Protein"))
    ck("a category with no id prefixes of its own takes those of the classes below it",
       {"MONDO", "HP"} <= prefixes_for("biolink:DiseaseOrPhenotypicFeature"))
    r = {"subject": {"columns": ["g"], "resolver": "prefix:NCBIGene", "category": "biolink:Gene"},
         "object": {"columns": ["d"], "resolver": "curie", "category": "biolink:Disease"},
         "predicate": "biolink:gene_associated_with_condition", "association_category": "biolink:GeneToDiseaseAssociation",
         "knowledge_level": "knowledge_assertion", "agent_type": "manual_agent", "primary_knowledge_source": "infores:ctd",
         "qualifiers": {}, "edge_properties": [], "keep_if": None}
    ck("a mapping that fits Biolink passes checks 1 to 7", not any(check_rule(r, {"g", "d"}).values()))
    w = {**r, "subject": {**r["subject"], "category": "biolink:Protein"}}
    ck("a protein under a predicate whose domain is gene is refused at check 4, and only there",
       [k for k, v in check_rule(w, {"g", "d"}).items() if v] == ["4_domain_range"])
    got = check_rule({**r, "association_category": "biolink:GeneToDiseaseOrPhenotypicFeatureAssociation"}, {"g", "d"})
    ck("an association class Biolink 4.4.3 does not have is refused at check 3", bool(got["3_categories"]))
    ck("a qualifier value outside its enumeration is refused",
       bool(check_rule({**r, "qualifiers": {"object_direction_qualifier": "up"}}, {"g", "d"})["6_qualifiers"]))
    ck("a qualifier the association class does not carry is refused",
       bool(check_rule({**r, "qualifiers": {"subject_direction_qualifier": "increased"}}, {"g", "d"})["6_qualifiers"]))
    w = {**r, "predicate": "biolink:affects", "association_category": "biolink:GeneToDiseaseAssociation"}
    ck("a predicate the association class lists passes check 5", not check_rule(w, {"g", "d"})["5_association"])
    w = {**r, "predicate": "biolink:treats"}
    ck("a predicate outside the association class's list is refused at check 5", bool(check_rule(w, {"g", "d"})["5_association"]))
    ck("resolvers: curie reads EFO_ as EFO:, a pattern takes its group, a name is local and not grounded",
       resolve_one("EFO_0004584", "curie", None) == ("EFO:0004584", True)
       and resolve_one("MIM Number: 613898; Phenotype: HMG-CoA lyase deficiency, 246450 (3)", r"pattern:OMIM:,\s*(\d{6})\s*\(\d\)", None) == ("OMIM:246450", True)
       and resolve_one("Brain", "name:tissue", None)[1] is False and resolve_one(98.0, "prefix:NCBIGene", None) == ("NCBIGene:98", True))
    ck("a resolver off the menu, and a pattern with no group, are refused",
       bool(resolver_ok("python:eval")) and bool(resolver_ok("pattern:OMIM:\\d+")) and not resolver_ok("pattern:OMIM:(\\d+)"))
    p = {"split_by": "t", "subject": r["subject"], "predicate": "biolink:related_to", "by_value": {"a": {"predicate": "biolink:affects"}, "b": {"skip": "no"}}}
    rs = rules_of(p)
    ck("a split applies each value's overrides on the shared fields, and keeps a skip with its reason",
       rs[0][0] == "a" and rs[0][1]["predicate"] == "biolink:affects" and rs[0][1]["subject"] == r["subject"] and rs[1] == ("b", None, "no"))
    a = agree({**r, "object": {**r["object"], "columns": ["d"]}, "predicate": "biolink:associated_with"},
              {**r, "object": {**r["object"], "columns": ["e", "d"]}})
    ck("agreement: a column agrees when the first one read is one the hand mapping reads; an ancestor predicate is close, not exact",
       a["subject_column"] and a["object_column"] and not a["predicate"] and a["predicate_close"])
    ck("agreement: a first column the hand mapping does not read disagrees",
       not agree({**r, "object": {**r["object"], "columns": ["x", "d"]}}, r)["object_column"])
    ck("every hand rule names a field the agreement reads", all(set(fields(x)) == set(FIELDS) for _, _, x in hand_rules()))
    ck("twelve hand rules over nine scored cards", len(hand_rules()) == 12 and len(SCORED) == 9)

    class _Refs:
        vocab = {**{p: set() for p in VOCAB_PREFIXES}, "NCBIGene": {f"NCBIGene:{i}" for i in range(1, 5)}, "MONDO": {"MONDO:1", "MONDO:2"}}
        xref = set()
    rows = pd.DataFrame({"g": ["1", "2", "3", "4"], "d": ["MONDO_1", "MONDO_2", "GO_3", "PATO_4"]})
    kept = {**r, "object": {"columns": ["d"], "resolver": "curie", "category": "biolink:Disease"},
            "keep_if": {"column": "d", "matches": "^MONDO_", "reason": "not a disease"}}
    v = validate(kept, rows, _Refs(), export=False)
    ck("the objects of rows the condition makes gaps are not judged against the category's prefixes",
       not v["problems"] and v["rules"][0]["identifiers"]["object"]["distinct"] == 2 and v["rules"][0]["made_gaps_by_condition"] == 2)
    ck("without the condition the same objects are refused", bool(validate({**kept, "keep_if": None}, rows, _Refs(), export=False)["problems"]))
    print(f"\n{len(bad)} failed" if bad else "\nall passed")
    return 1 if bad else 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", choices=["cards", "control", "propose", "stale", "report", "export"])
    ap.add_argument("--set-aside", help="with stale: move the stale records to discarded_<release> with this reason")
    ap.add_argument("release")
    ap.add_argument("--model", default=LOCAL)
    ap.add_argument("--set", choices=["scored", "new", "all"], default="scored")
    ap.add_argument("--views", nargs="*")
    a = ap.parse_args(argv)
    if a.what == "cards":
        out = make_cards(a.release)
        print(f"{cards_path(a.release).relative_to(paths.ROOT)}: {len(out['cards'])} cards, sha256 {sha(out['cards'])[:12]}")
    elif a.what == "control":
        for v, c in control(a.release).items():
            print(f"{v}: {'accepted' if c['accepted'] else 'refused'}")
            for p in c["problems"]:
                print("    " + p)
    elif a.what == "propose":
        views = a.views or {"scored": SCORED, "new": NEW, "all": SCORED + NEW}[a.set]
        propose(a.release, a.model, views)
    elif a.what == "stale":
        items = stale(a.release)
        for it in items:
            print(f"{it['model']} {it['view']} attempt {it['attempt']}:\n    recorded {it['recorded']}\n    now      {it['now']}")
        if a.set_aside and items:
            set_aside(a.release, items, a.set_aside)
            print(f"{len(items)} set aside")
        print(f"{len(items)} stale")
    elif a.what == "report":
        out = report(a.release)
        print(json.dumps({m: v["qualification"] for m, v in out["models"].items()}, indent=1))
    else:
        print(json.dumps(export_new(a.release), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
