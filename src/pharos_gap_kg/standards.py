"""The merged graph checked by the reference KGX validator against the Biolink model, offline.

kgx pins packages of its own, so this runs in an environment of its own and imports nothing from this
package:

    uv venv .venv-standards
    uv pip install --python .venv-standards/bin/python kgx==2.7.0 biolink-model==4.4.3
    .venv-standards/bin/python -I src/pharos_gap_kg/standards.py pharos319
    .venv-standards/bin/python -I src/pharos_gap_kg/standards.py --dir DIR OUT.json    # any KGX nodes.tsv and edges.tsv

kgx checks every node and edge for its required properties, the types and values of its properties
(a CURIE's prefix against Biolink's prefix map), and node categories and edge predicates against the
model. The schema and the prefix map are read from the installed biolink-model package. bmt and kgx
would otherwise fetch both from GitHub, kgx's prefix context at Biolink 2.2.5; here every network
connection is refused, so nothing is fetched. The gaps file is kgx-gaps' extension, not a KGX file,
and is not checked here; kgx-gaps' own conformance suite checks it with the edges.

Writes docs/data/standards_<release>.json.
"""
from __future__ import annotations

import collections
import json
import os
import re
import socket
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _offline() -> None:
    def refuse(*a, **k):
        raise RuntimeError("REFUSED: the standards check makes no network connection")
    socket.create_connection = refuse
    socket.socket.connect = refuse


def _local_model():
    """Point bmt and kgx at the installed biolink-model package's schema and prefix map, before kgx loads."""
    import bmt
    import bmt.toolkit as tk
    from linkml_runtime.utils.schemaview import SchemaView
    base = next(Path(p) / "biolink_model" for p in sys.path if (Path(p) / "biolink_model" / "schema" / "biolink_model.yaml").is_file())

    class Local(tk.Toolkit):
        def __init__(self, schema=None, predicate_map=None):
            self.view = SchemaView(str(base / "schema" / "biolink_model.yaml"))
            self.pmap = {}

    tk.Toolkit = bmt.Toolkit = Local
    import kgx.config as kc
    kc.jsonld_context_map["biolink"] = json.loads((base / "prefixmaps" / "biolink-model-prefix-map.json").read_text())
    return base


def check(graph_dir: Path) -> dict:
    _offline()
    _local_model()
    from importlib.metadata import version

    import kgx.utils.kgx_utils as ku
    from kgx.transformer import Transformer
    from kgx.validator import Validator
    v = Validator()
    seen = collections.Counter()

    def inspect(kind, rec):
        seen[str(kind).split(".")[-1].lower()] += 1
        v(kind, rec)

    t0 = time.time()
    Transformer(stream=True).transform(input_args={"filename": [str(graph_dir / "nodes.tsv"), str(graph_dir / "edges.tsv")], "format": "tsv"}, inspector=inspect)
    rows, nodes_bad, edges_bad = collections.defaultdict(lambda: {"ids": set(), "example": None}), set(), set()
    for level, by_type in v.get_errors().items():
        for typ, by_msg in by_type.items():
            for msg, ids in by_msg.items():
                key = (level, typ, re.sub(r"(has a value|value) '[^']*'", r"\1 <value>", msg))
                ids = set(ids) if isinstance(ids, (list, set)) else {str(ids)}
                rows[key]["ids"] |= ids
                rows[key]["example"] = rows[key]["example"] or msg
                (edges_bad if typ.startswith("INVALID_EDGE") or any("->" in i for i in ids) else nodes_bad).update(ids)
    errors = [{"level": k[0], "type": k[1], "message": k[2], "records": len(r["ids"]), "example": r["example"], "prefix": (re.search(r"CURIE prefix '([^']+)'", k[2]) or [None, None])[1]}
              for k, r in sorted(rows.items(), key=lambda kv: (-len(kv[1]["ids"]), kv[0]))]
    return {"validator": f"kgx {version('kgx')}", "biolink_model": ku.get_toolkit().get_model_version(), "bmt": version("bmt"),
            "checked": {"nodes": seen.get("node", 0), "edges": seen.get("edge", 0)}, "seconds": round(time.time() - t0, 1),
            "checks": ["required properties", "property types", "property values, including CURIE prefixes", "node categories", "edge predicates"],
            "not_checked_here": "the gaps file (kgx-gaps' extension), and whether edges point at nodes; kgx-gaps' conformance suite checks both",
            "nodes_with_an_error": len(nodes_bad), "subject_object_pairs_with_an_error": len(edges_bad), "errors": errors}


def main(argv: list[str]) -> int:
    if len(argv) == 3 and argv[0] == "--dir":
        # Any directory of KGX nodes.tsv and edges.tsv, such as a subgraph or an ingested view; the record goes to the named file.
        out = check(Path(argv[1]))
        Path(argv[2]).write_text(json.dumps({"dir": argv[1], **out}, indent=1) + "\n")
        print(f"  {out['checked']['nodes']:,} nodes and {out['checked']['edges']:,} edges; {out['nodes_with_an_error']:,} nodes and "
              f"{out['subject_object_pairs_with_an_error']:,} subject-object pairs with an error")
        for e in out["errors"]:
            print(f"  {e['records']:>6,}  {e['type']}: {e['message'][:150]}")
        return 0
    if len(argv) != 1:
        raise SystemExit(__doc__)
    release = argv[0]
    data = Path(os.environ.get("PHAROS_GAP_KG_DATA", ROOT / "data"))
    out = check(data / "graph" / release)
    (ROOT / "docs" / "data" / f"standards_{release}.json").write_text(json.dumps({"release": release, **out}, indent=1) + "\n")
    print(f"  {out['validator']}, Biolink {out['biolink_model']}: {out['checked']['nodes']:,} nodes and {out['checked']['edges']:,} edges in {out['seconds']} s; "
          f"{out['nodes_with_an_error']:,} nodes and {out['subject_object_pairs_with_an_error']:,} subject-object pairs with an error")
    for e in out["errors"]:
        print(f"  {e['records']:>6,}  {e['type']}: {e['message'][:150]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
