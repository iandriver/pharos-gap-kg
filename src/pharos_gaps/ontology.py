"""Which traits are the same trait, or a broader or narrower one, read from an OBO file.

Only `is_a` is followed. An obsolete term stands for the term named in its `replaced_by`, so two ids
that the ontology treats as one trait are one trait here. Ids are compared as CURIEs: the EFO release
writes its own terms as `efo:EFO_0000249`, which is `EFO:0000249`.

Run:  python -m pharos_gaps.ontology --selftest
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable


def norm(i: str) -> str:
    i = i.strip()
    if i.startswith("efo:EFO_"):
        return "EFO:" + i[len("efo:EFO_"):]
    return i


def load(lines: Iterable[str]) -> dict[str, dict]:
    """id -> {name, parents, obsolete, replaced_by} for every [Term] stanza."""
    terms: dict[str, dict] = {}
    cur = None
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("["):
            cur = {"name": "", "parents": set(), "obsolete": False, "replaced_by": None} if line == "[Term]" else None
            continue
        if cur is None:
            continue
        if line.startswith("id: "):
            terms[norm(line[4:])] = cur
        elif line.startswith("name: "):
            cur["name"] = line[6:]
        elif line.startswith("is_a: "):
            cur["parents"].add(norm(line[6:].split(" ! ")[0].split(" {")[0]))
        elif line.startswith("is_obsolete: true"):
            cur["obsolete"] = True
        elif line.startswith("replaced_by: "):
            cur["replaced_by"] = norm(line[13:])
    return terms


def read(path: Path) -> dict[str, dict]:
    with path.open(encoding="utf-8", errors="replace") as f:
        return load(f)


def canonical(i: str, terms: dict[str, dict]) -> str | None:
    """The live term an id stands for: itself, or what replaced it. None if the ontology lacks it or
    it is obsolete with no replacement."""
    seen = set()
    while i in terms and terms[i]["obsolete"]:
        nxt = terms[i]["replaced_by"]
        if not nxt or nxt in seen:
            return None
        seen.add(i)
        i = nxt
    return i if i in terms else None


def ancestors(i: str, terms: dict[str, dict]) -> set[str]:
    out: set[str] = set()
    stack = list(terms[i]["parents"]) if i in terms else []
    while stack:
        p = stack.pop()
        if p not in out:
            out.add(p)
            stack += list(terms[p]["parents"]) if p in terms else []
    return out


def lineage(traits: list[str], terms: dict[str, dict]) -> dict[str, set[str]]:
    """For each trait, the OTHER traits of the list in its lineage: same live term, ancestor or descendant."""
    canon = {t: canonical(t, terms) for t in traits}
    anc = {t: (ancestors(c, terms) if c else set()) for t, c in canon.items()}
    out: dict[str, set[str]] = {t: set() for t in traits}
    for a in traits:
        if not canon[a]:
            continue
        for b in traits:
            if a == b or not canon[b]:
                continue
            if canon[a] == canon[b] or canon[b] in anc[a] or canon[a] in anc[b]:
                out[a].add(b)
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    obo = """format-version: 1.2

[Term]
id: efo:EFO_0000001
name: root

[Term]
id: MONDO:1
name: disease
is_a: efo:EFO_0000001 ! root

[Term]
id: MONDO:2
name: heart disease
is_a: MONDO:1 ! disease

[Term]
id: MONDO:3
name: arrhythmia
is_a: MONDO:2 {source="x"} ! heart disease

[Term]
id: MONDO:4
name: lung disease
is_a: MONDO:1 ! disease

[Term]
id: efo:EFO_0000009
name: obsolete_arrhythmia
is_obsolete: true
replaced_by: MONDO:3

[Term]
id: efo:EFO_0000010
name: obsolete_gone
is_obsolete: true

[Typedef]
id: part_of
name: part of
""".splitlines()
    T = load(obo)
    ck("an EFO id is read as a CURIE, and a Typedef is not a term", "EFO:0000001" in T and "part_of" not in T and len(T) == 7)
    ck("a parent is read without its comment or its braces", T["MONDO:3"]["parents"] == {"MONDO:2"})
    ck("ancestors go to the root", ancestors("MONDO:3", T) == {"MONDO:2", "MONDO:1", "EFO:0000001"})
    ck("an obsolete term stands for its replacement", canonical("EFO:0000009", T) == "MONDO:3")
    ck("an obsolete term with no replacement, and an unknown id, stand for nothing", canonical("EFO:0000010", T) is None and canonical("XX:1", T) is None)
    L = lineage(["MONDO:2", "MONDO:3", "MONDO:4", "EFO:0000009", "EFO:0000010", "XX:1"], T)
    ck("a broader and a narrower trait are in one lineage", "MONDO:3" in L["MONDO:2"] and "MONDO:2" in L["MONDO:3"])
    ck("two ids for one term are in one lineage", "EFO:0000009" in L["MONDO:3"] and "MONDO:2" in L["EFO:0000009"])
    ck("siblings are not", "MONDO:4" not in L["MONDO:2"] and L["MONDO:4"] == set())
    ck("a trait is not in its own lineage, and an unknown trait has none", "MONDO:2" not in L["MONDO:2"] and L["XX:1"] == set() and L["EFO:0000010"] == set())
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else 0)
