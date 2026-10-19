"""What changed between two bundles: each table row by row on its key, and the answers to a fixed set of
questions with the edge ids masked.

Used to show that a change to how the graph is exported leaves what the questions say as it was: a
table may gain a column or move an edge id, and every other value must be equal, row for row.

Run:  python -m pharos_gap_kg.bundlediff BEFORE_DIR AFTER_DIR --out docs/data/remodel_check_pharos319.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd

from pharos_gap_kg import ask

KEYS = {"assertions": ["target", "condition", "source"], "expression": ["target", "condition", "label", "edge"], "drug_targets": ["drug", "target"],
        "drug_indications": ["drug", "condition"], "links": ["edge_id"], "nodes": ["id"], "context_nodes": ["id"], "evidence": ["target", "trait", "study"],
        "lineage": ["a", "b", "relation"], "neighbours": ["a", "b"], "target_gaps": ["id"], "study_sizes": None}
IDS = {"edge", "edge_id"}
ID = re.compile(r"pharosgapkg:[eg][0-9a-f]{12}")
QUESTIONS = [("why", "APOE", "Alzheimer disease"), ("why", "ACHE", "Alzheimer disease"), ("why", "HMGCR", "hypercholesterolemia"), ("why", "MT-RNR1", "tinnitus"),
             ("known", "PCSK9"), ("known", "atorvastatin"), ("known", "Alzheimer disease"), ("connect", "LDLR", "PCSK9"),
             ("connect", "atorvastatin", "hypercholesterolemia"), ("predict", "Parkinson disease"), ("gaps", None, "Parkinson disease"), ("gaps", "POLR3E"),
             ("whats_new", "type 2 diabetes mellitus", 2020), ("status",)]


def table(a: pd.DataFrame, b: pd.DataFrame, key: list[str] | None) -> dict:
    out = {"rows": [len(a), len(b)], "columns_added": sorted(set(b.columns) - set(a.columns)), "columns_removed": sorted(set(a.columns) - set(b.columns))}
    if key is None:
        out["equal"] = a.equals(b)
        return out
    m = a.merge(b, on=key, how="outer", suffixes=("_a", "_b"), indicator=True)
    out["rows_only_before"], out["rows_only_after"] = int((m._merge == "left_only").sum()), int((m._merge == "right_only").sum())
    both = m[m._merge == "both"]
    diff = {}
    for c in sorted(set(a.columns) & set(b.columns) - set(key)):
        x, y = both[f"{c}_a"], both[f"{c}_b"]
        n = int((~((x == y) | (x.isna() & y.isna()))).sum())
        if n:
            diff[c] = n
    out["values_that_differ"] = diff
    out["equal_but_for_ids"] = not out["rows_only_before"] and not out["rows_only_after"] and not out["columns_removed"] and not (set(diff) - IDS)
    return out


def answers(bundle: Path) -> dict:
    G = ask.Graph(bundle)
    res = {}
    for q in QUESTIONS:
        k = "|".join(map(str, q))
        if q[0] == "why":
            res[k] = G.why(q[1], q[2])
        elif q[0] == "known":
            res[k] = G.known(q[1])
        elif q[0] == "connect":
            res[k] = G.connect(q[1], q[2])
        elif q[0] == "predict":
            res[k] = G.predict(q[1])
        elif q[0] == "gaps":
            res[k] = G.gaps(q[1], None, None, q[2]) if q[1] is None else G.gaps(q[1])
        elif q[0] == "whats_new":
            res[k] = G.whats_new(q[1], q[2])
        else:
            res[k] = G.status()
    return json.loads(json.dumps(res, default=str))


def compare(before: Path, after: Path) -> dict:
    tables = {}
    for f in sorted(after.glob("*.parquet")):
        name = f.stem
        if (before / f.name).exists():
            tables[name] = table(pd.read_parquet(before / f.name), pd.read_parquet(f), KEYS.get(name))
    a, b = answers(before), answers(after)
    mask = lambda x: json.loads(ID.sub("ID", json.dumps(x)))
    qs = {k: ("same" if a[k] == b[k] else "same but for edge ids" if mask(a[k]) == mask(b[k]) else "differs") for k in a}
    sha = lambda d: hashlib.sha256((d / "manifest.json").read_bytes()).hexdigest()
    return {"before_manifest_sha256": sha(before), "after_manifest_sha256": sha(after), "tables": tables, "answers": qs,
            "every_table_equal_but_for_ids": all(t.get("equal_but_for_ids", t.get("equal")) for t in tables.values()),
            "every_answer_equal_but_for_ids": all(v != "differs" for v in qs.values())}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("before", type=Path)
    ap.add_argument("after", type=Path)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    out = compare(a.before, a.after)
    text = json.dumps(out, indent=1)
    if a.out:
        a.out.write_text(text + "\n")
    print(text)
    return 0 if out["every_table_equal_but_for_ids"] and out["every_answer_equal_but_for_ids"] else 1


if __name__ == "__main__":
    sys.exit(main())
