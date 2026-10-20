"""Withheld rows that another source asserts for the same target, on a condition of the same name.

The pruning rule counts a second source by the condition's identifier: a CTD row on MESH:D014012 is
"one source" even when DisGeNET asserts the same target with UMLS:C0040264, both named tinnitus. This
counts, for every withheld row, the kept rows of the same target from another source on a condition
whose name, normalised as the questions normalise names, is the same. It says how many withheld rows a
match by name would reach and from which sources. It does not say whether such a row is independent
evidence: DisGeNET's curated set draws on CTD.

Run:  python -m pharos_gap_kg.corroboration pharos319     # writes docs/data/withheld_named_elsewhere_<release>.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter

import pandas as pd

from pharos_gap_kg import ask, paths


def count(bundle) -> dict:
    A = pd.read_parquet(bundle / "assertions.parquet")
    names = {}
    for f in ("nodes", "context_nodes"):
        n = pd.read_parquet(bundle / f"{f}.parquet", columns=["id", "name"])
        names.update(dict(zip(n.id, n.name)))
    A["named"] = [ask._norm(names.get(c, c)) for c in A.condition]
    kept, held = A[~A.withheld], A[A.withheld]
    by_name = kept.groupby(["target", "named"]).source.apply(set).to_dict()
    ids = kept.groupby(["target", "named"]).condition.apply(set).to_dict()
    reached, others, under_another_id, examples = Counter(), Counter(), 0, []
    for r in held.itertuples(index=False):
        other = by_name.get((r.target, r.named), set()) - {r.source}
        if not other:
            continue
        reached[r.source] += 1
        others.update(other)
        under_another_id += r.condition not in ids.get((r.target, r.named), set())
        if len(examples) < 5:
            examples.append({"target": r.target, "condition": r.condition, "source": r.source, "also_asserted_by": sorted(other)})
    return {"withheld": int(len(held)), "with_a_same_named_row_from_another_source": int(sum(reached.values())),
            "of_those_under_another_identifier": int(under_another_id), "by_withheld_source": dict(reached.most_common()),
            "rows_by_the_other_source": dict(others.most_common()), "examples": examples,
            "names": "the condition's name in the bundle, normalised as the questions normalise names (ask._norm)"}


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        raise SystemExit(__doc__)
    out = {"release": argv[0], **count(paths.bundle(argv[0]))}
    dest = paths.DOCS_DATA / f"withheld_named_elsewhere_{argv[0]}.json"
    dest.write_text(json.dumps(out, indent=1) + "\n")
    print(f"{out['with_a_same_named_row_from_another_source']:,} of {out['withheld']:,} withheld rows have a kept row from another source on a "
          f"condition of the same name; by the other source: {out['rows_by_the_other_source']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
