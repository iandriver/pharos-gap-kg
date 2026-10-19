"""The parts of the Biolink model a mapping is checked against, read from the installed schema.

The ingest validator (`ingest.py`) needs Biolink's classes, predicates, qualifiers and enumerations,
with their parents, domains, ranges and identifier prefixes. biolink-model pins packages of its own,
so the schema is read here, in the standards environment, and written once as a small JSON file that
the main package reads. Nothing is fetched.

    .venv-standards/bin/python -I src/pharos_gap_kg/biolinkslice.py

Writes docs/data/biolink_<version>.json, which records the SHA-256 of the schema it was read from.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def camel(name: str) -> str:
    return "".join(w[0].upper() + w[1:] for w in name.split())


def snake(name: str) -> str:
    return name.replace(" ", "_")


def ancestors(slots: dict, name: str) -> list[str]:
    out = []
    while name:
        out.append(name)
        name = slots.get(name, {}).get("is_a")
    return out


def slice_of(y: dict) -> dict:
    S, C, E = y["slots"], y["classes"], y["enums"]
    cls = lambda n: f"biolink:{camel(n)}" if n else None
    classes = {}
    for n, c in C.items():
        usage = {}
        for s, u in (c.get("slot_usage") or {}).items():
            u = u or {}
            keep = {k: u[k] for k in ("range", "subproperty_of") if u.get(k)}
            if "range" in keep and keep["range"] in C:
                keep["range"] = cls(keep["range"])
            if "subproperty_of" in keep:
                keep["subproperty_of"] = f"biolink:{snake(keep['subproperty_of'])}"
            if keep:
                usage[snake(s)] = keep
        classes[cls(n)] = {"is_a": cls(c.get("is_a")), "mixins": [cls(m) for m in c.get("mixins") or []],
                           "mixin": bool(c.get("mixin")), "abstract": bool(c.get("abstract")),
                           "deprecated": bool(c.get("deprecated")), "id_prefixes": list(c.get("id_prefixes") or []),
                           "slots": [snake(s) for s in c.get("slots") or []], "slot_usage": usage}
    rng = lambda r: cls(r) if r in C else r
    predicates, qualifiers = {}, {}
    for n, s in S.items():
        line = ancestors(S, n)
        if "related to" in line:
            predicates[f"biolink:{snake(n)}"] = {
                "is_a": f"biolink:{snake(s['is_a'])}" if s.get("is_a") else None,
                "mixins": [f"biolink:{snake(m)}" for m in s.get("mixins") or []],
                "domain": rng(s.get("domain")), "range": rng(s.get("range")),
                "mixin": bool(s.get("mixin")), "abstract": bool(s.get("abstract")), "deprecated": bool(s.get("deprecated")),
                "symmetric": bool(s.get("symmetric")), "inverse": f"biolink:{snake(s['inverse'])}" if s.get("inverse") else None,
                "canonical": "canonical_predicate" in (s.get("annotations") or {})}
        elif "qualifier" in line:
            qualifiers[snake(n)] = {"range": rng(s.get("range")), "is_a": snake(s["is_a"]) if s.get("is_a") else None}
    enums = {n: sorted((e or {}).get("permissible_values") or {}) for n, e in E.items()}
    return {"version": y.get("version"), "classes": classes, "predicates": predicates, "qualifiers": qualifiers, "enums": enums}


def main(argv: list[str]) -> int:
    import biolink_model
    import yaml
    path = Path(list(biolink_model.__path__)[0]) / "schema" / "biolink_model.yaml"     # a namespace package: no __file__
    raw = path.read_bytes()
    out = slice_of(yaml.safe_load(raw))
    out = {"schema": "biolink_model.yaml from the biolink-model package", "schema_sha256": hashlib.sha256(raw).hexdigest(), **out}
    dest = ROOT / "docs" / "data" / f"biolink_{out['version']}.json"
    tmp = dest.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    os.replace(tmp, dest)
    print(f"{dest.relative_to(ROOT)}: {len(out['classes'])} classes, {len(out['predicates'])} predicates, "
          f"{len(out['qualifiers'])} qualifiers, {len(out['enums'])} enumerations")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
