"""How two published snapshots differ, table by table, from their inventories.

A dump's publication date is not the date of its contents. Before any link is labelled "appeared
later", this says which tables of the later dump differ from the earlier one at all. It reads the two
inventories written by `sqldump inventory` (committed under docs/data), so it needs neither dump.

A table is IDENTICAL when its payload digests match: same rows, same ids, same order. A table whose
digests differ is reported with both row counts; whether its rows differ in content or only in order
is a row-level question for `sqldump extract`.

Run:  python -m pharos_gaps.snapshots pharos316 pharos319
"""
from __future__ import annotations

import json
import sys

from pharos_gaps import paths


def load(key: str) -> dict:
    p = paths.DOCS_DATA / f"{key}_tables.json"
    if not p.is_file():
        raise SystemExit(f"no inventory for {key}; write it with:  python -m pharos_gaps.sqldump inventory {key}")
    d = json.loads(p.read_text())
    if any("payload_sha256" not in t for t in d["tables"].values()):
        raise SystemExit(f"{p.name} predates table digests; re-run:  python -m pharos_gaps.sqldump inventory {key}")
    return d


def compare(a: dict, b: dict) -> dict:
    ta, tb = a["tables"], b["tables"]
    both = sorted(set(ta) & set(tb))
    same = [t for t in both if ta[t]["payload_sha256"] == tb[t]["payload_sha256"]]
    differ = [t for t in both if t not in set(same)]
    cols = [t for t in both if [c["name"] for c in ta[t]["columns"]] != [c["name"] for c in tb[t]["columns"]]]
    return {
        "identical": same,
        "differ": [{"table": t, "rows_a": ta[t]["rows"], "rows_b": tb[t]["rows"],
                    "columns_changed": t in cols} for t in differ],
        "only_a": sorted(set(ta) - set(tb)), "only_b": sorted(set(tb) - set(ta)),
        "rows_identical": sum(ta[t]["rows"] for t in same),
        "rows_a": sum(d["rows"] for d in ta.values()), "rows_b": sum(d["rows"] for d in tb.values()),
    }


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    def tab(rows, digest, cols=("id", "x")):
        return {"rows": rows, "payload_sha256": digest, "columns": [{"name": c, "type": "int"} for c in cols]}

    a = {"tables": {"same": tab(5, "aa"), "grew": tab(5, "bb"), "gone": tab(1, "cc"), "empty": tab(0, None)}}
    b = {"tables": {"same": tab(5, "aa"), "grew": tab(7, "dd", ("id", "x", "y")), "new": tab(2, "ee"), "empty": tab(0, None)}}
    r = compare(a, b)
    ck("equal digests are identical, and two empty tables are identical", r["identical"] == ["empty", "same"])
    ck("unequal digests differ, with both row counts", r["differ"] == [{"table": "grew", "rows_a": 5, "rows_b": 7, "columns_changed": True}])
    ck("a table in one dump only is named on its side", r["only_a"] == ["gone"] and r["only_b"] == ["new"])
    ck("identical rows are counted", r["rows_identical"] == 5 and r["rows_a"] == 11 and r["rows_b"] == 14)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 2:
        raise SystemExit(__doc__)
    a, b = load(argv[0]), load(argv[1])
    r = compare(a, b)
    print(f"{argv[0]} (published {a['published']})  ->  {argv[1]} (published {b['published']})")
    print(f"  tables: {len(a['tables'])} and {len(b['tables'])}; identical {len(r['identical'])}, "
          f"differ {len(r['differ'])}, only in the first {len(r['only_a'])}, only in the second {len(r['only_b'])}")
    print(f"  rows: {r['rows_a']:,} and {r['rows_b']:,}; {r['rows_identical']:,} sit in identical tables")
    for d in sorted(r["differ"], key=lambda d: -max(d["rows_a"], d["rows_b"])):
        print(f"    differs  {d['table']:34s} {d['rows_a']:>12,} -> {d['rows_b']:>12,}"
              f"{'   (columns changed)' if d['columns_changed'] else ''}")
    for t in r["only_a"]:
        print(f"    only in {argv[0]}: {t} ({a['tables'][t]['rows']:,} rows)")
    for t in r["only_b"]:
        print(f"    only in {argv[1]}: {t} ({b['tables'][t]['rows']:,} rows)")
    out = paths.DOCS_DATA / f"snapshots_{argv[0]}_{argv[1]}.json"
    out.write_text(json.dumps({"a": argv[0], "b": argv[1], "sha256_a": a["sha256"], "sha256_b": b["sha256"], **r}, indent=1) + "\n")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
