"""Replay Cheng's saved proposal with full references and with only the release bundle.

Run from the checkout with its full reference tables and registered supplement present.
The JSON on stdout contains counts and resolver differences, without study rows.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

from pharos_gap_kg import add, ask, paths, registry, release


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preview(source, record, bundle):
    graph = ask.Graph(bundle)
    session = add.Session(graph)
    about = add.STUDIES["D"]["about"]
    proposal = record["mapping"]
    choice = {"sheet": record["sheet"], "header_row": record["header_row"]}
    check = add.tool_check(session, str(source), about, "mapping", proposal, **choice)
    if not check["accepted"]:
        raise ValueError(check["problems"])
    rows, errors = add.check1(choice, add.read_sheets(source))
    if errors:
        raise ValueError(errors)
    symbols = sorted(set(rows[proposal["subject"]["columns"][0]].dropna()))
    return {"references": check["checked_against"], **check["would_join"],
            "rows": check["report"]["source_rows"], "accounted": check["report"]["accounted_source_rows"],
            "symbols": {s: session.sym(s) for s in symbols}}


def main():
    record = json.loads(add.record_path("pharos319", "D").read_text())
    if len(sys.argv) == 3 and sys.argv[1] == "--compact":
        print(json.dumps(preview(Path(sys.argv[2]), record, paths.bundle("pharos319"))))
        return
    source = registry.require(record["key"]).resolve()
    if sha(source) != record["file_sha256"]:
        raise ValueError("The source differs from the registered Cheng workbook")
    bundle = paths.bundle("pharos319")
    full = preview(source, record, bundle)
    archive = paths.RELEASE / release.name("pharos319")
    with TemporaryDirectory() as tmp:
        env = dict(os.environ, PHAROS_GAP_KG_DATA=tmp)
        env.pop("PHAROS_GAP_KG_BUNDLE", None)
        subprocess.run([sys.executable, "-m", "pharos_gap_kg.release", "unpack", "pharos319",
                        str(archive)], env=env, check=True, stdout=sys.stderr)
        compact = json.loads(subprocess.check_output([sys.executable, __file__, "--compact", str(source)], env=env))
        if sha(bundle / "manifest.json") != sha(Path(tmp) / "bundle/pharos319/manifest.json"):
            raise ValueError("The two runs did not use the same bundle")
    a, b = full.pop("symbols"), compact.pop("symbols")
    if a.keys() != b.keys() or full["references"] == compact["references"]:
        raise ValueError("The comparison requires identical source values and distinct reference sets")
    diff = [{"symbol": s, "full": a[s], "bundle": b[s]} for s in sorted(a) if a[s] != b[s]]
    result = {"study": record["key"], "doi": add.STUDIES["D"]["about"]["doi"],
              "source_sha256": sha(source), "bundle_archive_sha256": sha(archive),
              "mapping_record": add.record_path("pharos319", "D").name,
              "mapping_record_sha256": sha(add.record_path("pharos319", "D")),
              "sheet": record["sheet"], "header_row": record["header_row"],
              "full": full, "bundle": compact, "resolver_differences": diff,
              "operation": "Read-only replay of the accepted proposal. No model calls or study join."}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
