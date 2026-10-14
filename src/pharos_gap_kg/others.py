"""What the questions read beside the dated GWAS links: the merged graph's other parts.

Which disease sources assert a target with a condition, and which rows a pruning rule withheld and
on what count; the drugs that act on a target and what they treat; and where a gene's expression
differs in a condition. All of it is as of the release: none of these sources dates a row, so an
`as_of` year does not apply to them, and every answer built here says so.

A CONDITION IS READ UNDER EVERY ID THAT NAMES IT. The GWAS table names traits by one vocabulary and
the disease sources by others, so one disease can sit under several ids. Three ties are used, and an
answer lists the ids it read:

    its own id
    the current term of a retired id       the ontology's own statement (a same_as edge in the graph)
    an id the sources give the same name   a match on the name alone; not an edge in the graph

An expression observation is kept apart from an assertion. It says a gene's expression was higher or
lower in a condition; it does not say the gene is involved in it.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import ctxgraph

UNDATED = "These sources carry no dates: this is as of the release, whatever year was asked."
OBSERVATION = "An expression observation says the gene's expression differed in the condition. It does not say the gene is involved in it."
SAME_NAME = "the sources give it the same name"
NOT_MECHANISM = "The release does not mark this target as the drug's mechanism of action, so the drug's indications say less about this target."
FILES = ("context_nodes.parquet", "assertions.parquet", "expression.parquet", "drug_targets.parquet", "drug_indications.parquet")


def norm(s: str) -> str:
    return " ".join(re.sub(r"[^a-z0-9]+", " ", str(s).casefold()).split())


def plural(n: int, one: str, many: str | None = None) -> str:
    return f"{n:,} {one if n == 1 else (many or one + 's')}"


class Others:
    def __init__(self, d: Path, traits: pd.DataFrame, symbol_of):
        """traits: the bundle's GWAS traits (id, name, live_term). symbol_of: target id -> gene symbol."""
        rd = lambda f: pd.read_parquet(Path(d) / f)
        self.N, self.A, self.X, self.DT, self.DI = (rd(f) for f in FILES)
        self.symbol_of = symbol_of
        self.name = dict(zip(self.N.id, self.N.name))
        self.name.update(dict(zip(traits.id, traits.name)))
        self.grounded = dict(zip(self.N.id, self.N.grounded))
        self.a_t, self.a_c = self.A.groupby("target").indices, self.A.groupby("condition").indices
        self.x_t, self.x_c = self.X.groupby("target").indices, self.X.groupby("condition").indices
        self.dt_t, self.dt_d = self.DT.groupby("target").indices, self.DT.groupby("drug").indices
        self.di_c, self.di_d = self.DI.groupby("condition").indices, self.DI.groupby("drug").indices
        rows = lambda c: len(self.a_c.get(c, ())) + len(self.x_c.get(c, ())) + len(self.di_c.get(c, ()))
        cond = self.N[self.N.kind == "condition"]
        by_name: dict[str, list[str]] = {}
        for i, n in zip(cond.id, cond.name):
            by_name.setdefault(norm(n), []).append(i)
        # a GWAS trait is read under its own id, its current term, and any condition the sources name the same
        self.trait_ids: dict[str, list[str]] = {}
        self.folded: dict[str, str] = {}                # an id read as part of a GWAS trait -> that trait
        self.tie: dict[tuple[str, str], str] = {}       # (the id asked about, an id read for it) -> how the two are tied
        taken: set[str] = set()
        for t, n, live in zip(traits.id, traits.name, traits.live_term):
            ids = [t] + ([live] if isinstance(live, str) and live and live != t else [])
            self.tie.update({(t, t): "as published"} | ({(t, live): "its current ontology term"} if len(ids) > 1 else {}))
            for c in by_name.get(norm(n), []):
                if c not in ids:
                    ids.append(c)
                    self.tie[(t, c)] = SAME_NAME
            self.trait_ids[t] = ids
            taken.update(ids)
            for c in ids[1:]:
                self.folded.setdefault(c.casefold(), t)
        # the other conditions, one group per name, named after the id with the most rows
        self.groups: dict[str, list[str]] = {}
        self._cname: dict[str, str] = {}
        for key, ids in by_name.items():
            left = sorted((c for c in ids if c not in taken), key=lambda c: (-rows(c), c))
            if left:
                self.groups[left[0]] = left
                self._cname[key] = left[0]
                self.tie.update({(left[0], c): ("as the sources give it" if c == left[0] else SAME_NAME) for c in left})
        self._cid = {c.casefold(): g for g, ids in self.groups.items() for c in ids}
        self._ctoks = [(k, set(k.split()), g) for k, g in self._cname.items()]
        drug = self.N[self.N.kind == "drug"]
        self._dname = {norm(n): i for i, n in zip(drug.id, drug.name)}
        self._did = {i.casefold(): i for i in drug.id}
        self.rows_of = rows

    # ------------------------------------------------------------------ names
    def match(self, text: str, kinds: tuple[str, ...]) -> list[dict]:
        """Conditions and drugs that `text` may name: tier 0 an id, 1 the exact name, 2 every word, 3 contained."""
        q, raw = norm(text), str(text).strip().casefold()
        out = []
        if "condition" in kinds:
            if raw in self._cid:
                out.append({"tier": 0, "kind": "condition", "key": self._cid[raw]})
            if q in self._cname:
                out.append({"tier": 1, "kind": "condition", "key": self._cname[q]})
            toks = set(q.split())
            if toks:
                for k, ks, g in self._ctoks:
                    if k != q and toks <= ks:
                        out.append({"tier": 2, "kind": "condition", "key": g})
                    elif k != q and len(q) >= 4 and q in k:
                        out.append({"tier": 3, "kind": "condition", "key": g})
        if "drug" in kinds:
            if raw in self._did:
                out.append({"tier": 0, "kind": "drug", "key": self._did[raw]})
            if q in self._dname:
                out.append({"tier": 1, "kind": "drug", "key": self._dname[q]})
            if len(q) >= 4:
                out += [{"tier": 3, "kind": "drug", "key": i} for n, i in self._dname.items() if n != q and n.startswith(q)]
        seen, uniq = set(), []
        for m in sorted(out, key=lambda m: (m["tier"], -self.size(m["kind"], m["key"]), m["key"])):
            if (m["kind"], m["key"]) not in seen:
                seen.add((m["kind"], m["key"]))
                uniq.append({**m, "links": self.size(m["kind"], m["key"])})
        return uniq

    def size(self, kind: str, key: str) -> int:
        if kind == "drug":
            return len(self.dt_d.get(key, ())) + len(self.di_d.get(key, ()))
        return sum(self.rows_of(c) for c in self.groups.get(key, [key]))

    def label(self, kind: str, key: str) -> dict:
        if kind == "drug":
            return {"id": key, "name": self.name.get(key, key)}
        out = {"id": key, "name": self.name.get(key, key)}
        if not self.grounded.get(key, True):
            out["no_identifier_in_the_release"] = True
        also = [c for c in self.groups.get(key, [key]) if c != key]
        if also:
            out["also_under"] = also
        return out

    def ids(self, kind: str, key: str) -> list[str]:
        return self.trait_ids.get(key, [key]) if kind == "trait" else self.groups.get(key, [key])

    def read_under(self, ids: list[str]) -> tuple[list[dict], str]:
        """The ids a condition was read under with how each is tied to the first, and a sentence where a tie rests on the name alone."""
        ties = [{"id": c, "tie": self.tie.get((ids[0], c), "as published")} for c in ids]
        by_name = [t["id"] for t in ties if t["tie"] == SAME_NAME]
        say = (f" This was also read under {', '.join(by_name)}, which {'has' if len(by_name) == 1 else 'have'} the same name in the sources: a match on the name, not an ontology statement."
               if by_name else "")
        return ties, say

    # ------------------------------------------------------------------ one target and one condition
    def _withheld(self, r) -> dict:
        rules = r.rule.split("|")
        out = {"source": r.source, "gap": r.edge, "rule": r.rule, "why": "; ".join(ctxgraph.RULES[x][0] for x in rules),
               "restored_by": "; or ".join(ctxgraph.RULES[x][1] for x in rules)}
        if "literature_lift" in rules:
            out.update({"papers_mentioning_both": int(r.papers), "lift": round(float(r.lift), 3), "cut_for_this_source": round(float(r.cut), 3)})
        return out

    def _paths(self, target: str, ids: list[str]) -> list[dict]:
        """Drugs that act on the target and treat the condition: the two edges of each path."""
        want = set(ids)
        out = []
        for p in self.dt_t.get(target, ()):
            d = self.DT.iloc[p]
            for q in self.di_d.get(d.drug, ()):
                i = self.DI.iloc[q]
                if i.condition in want:
                    out.append({"drug": self.name.get(d.drug, d.drug), "drug_id": d.drug, "action_on_target": d.action or None, "activity": d.activity or None,
                                "target_is_the_drugs_mechanism": bool(d.mechanism), "edges": [d.edge, i.edge]})
        return sorted(out, key=lambda r: (not r["target_is_the_drugs_mechanism"], r["drug"]))

    def pair(self, target: str, ids: list[str], limit: int = 8) -> dict:
        sym, want = self.symbol_of(target), set(ids)
        a = self.A.iloc[[p for p in self.a_t.get(target, ()) if self.A.condition.iat[p] in want]]
        kept, held = a[~a.withheld], a[a.withheld]
        x = self.X.iloc[[p for p in self.x_t.get(target, ()) if self.X.condition.iat[p] in want]]
        paths_ = self._paths(target, ids)
        up, down = int((x.log2fc > 0).sum()), int((x.log2fc < 0).sum())
        bits = []
        if len(kept):
            bits.append(f"{plural(kept.source.nunique(), 'source')} assert{'s' if kept.source.nunique() == 1 else ''} {sym} with this condition ({', '.join(sorted(set(kept.source)))})")
        if len(held):
            h = held.iloc[0]
            count = (f"; {plural(int(h.papers), 'paper')} mention both, a lift of {float(h.lift):.3g} against a cut of {float(h.cut):.3g}" if "literature_lift" in h.rule else "")
            bits.append(f"{plural(len(held), 'row')} from {', '.join(sorted(set(held.source)))} {'is' if len(held) == 1 else 'are'} withheld by a pruning rule ({self._withheld(h)['why']}{count})")
        if paths_:
            moa = [p for p in paths_ if p["target_is_the_drugs_mechanism"]]
            names = ", ".join(dict.fromkeys(p["drug"] for p in (moa or paths_)[:4]))
            n = len({p["drug_id"] for p in paths_})
            bits.append(f"{plural(n, 'drug')} that act{'s' if n == 1 else ''} on {sym} {'is' if n == 1 else 'are'} indicated for it ({names}"
                        + (f"; for {len({p['drug_id'] for p in moa})} of them {sym} is the mechanism target)" if moa else f"; {sym} is the mechanism target of none of them)"))
        if len(x):
            bits.append(f"its expression was higher in this condition in {plural(up, 'observation')} and lower in {down}")
        ties, tied = self.read_under(ids)
        under = lambda c: {} if c == ids[0] else {"under": c}
        out = {"read_under": ties, "sources_asserting": int(kept.source.nunique()), "withheld_rows": len(held), "expression_rows": len(x), "asserted_by": [{"source": r.source, "edge": r.edge, **under(r.condition), **({"source_evidence": r.evidence} if r.evidence else {})} for r in kept.itertuples()][:limit],
               "withheld": [self._withheld(r) for r in held.itertuples()][:limit],
               "through_a_drug": paths_[:limit], "drug_paths": len(paths_),
               "expression_observations": [{"condition_as_the_source_names_it": r.label, "direction": "higher" if r.log2fc > 0 else "lower", "log2_fold_change": round(float(r.log2fc), 3),
                                            "p": float(r.p), "edge": r.edge} for r in x.sort_values("p").itertuples()][:limit],
               "dated": False}
        if bits:
            out["answer"] = "Apart from GWAS, as of the release: " + "; ".join(bits) + "." + (f" {OBSERVATION}" if len(x) else "") + tied
        else:
            out["answer"] = (f"Apart from GWAS, nothing in this graph ties {sym} to this condition: no disease source asserts the pair, no drug path runs between them, "
                             "and no expression observation is on record. That is not evidence that there is no association.")
        return out

    # ------------------------------------------------------------------ one node
    def target(self, target: str, limit: int = 10) -> dict:
        sym = self.symbol_of(target)
        a = self.A.iloc[self.a_t.get(target, [])]
        kept, held = a[~a.withheld], a[a.withheld]
        by = kept.groupby("condition").source.agg(lambda s: sorted(set(s)))
        top = sorted(by.items(), key=lambda kv: (-len(kv[1]), self.name.get(kv[0], kv[0])))[:limit]
        dt = self.DT.iloc[self.dt_t.get(target, [])]
        drugs = sorted(({"drug": self.name.get(r.drug, r.drug), "drug_id": r.drug, "action_on_target": r.action or None, "target_is_the_drugs_mechanism": bool(r.mechanism),
                         "indications": len(self.di_d.get(r.drug, ())), "edge": r.edge} for r in dt.itertuples()), key=lambda r: (not r["target_is_the_drugs_mechanism"], -r["indications"], r["drug"]))
        x = self.X.iloc[self.x_t.get(target, [])]
        strongest = x.reindex(x.log2fc.abs().sort_values(ascending=False).index)[:5]
        bits = []
        if len(kept):
            bits.append(f"{plural(kept.source.nunique(), 'disease source')} assert{'s' if kept.source.nunique() == 1 else ''} {sym} with {plural(kept.condition.nunique(), 'condition')}")
        if len(held):
            bits.append(f"{plural(len(held), 'row')} {'is' if len(held) == 1 else 'are'} withheld by a pruning rule")
        if len(dt):
            bits.append(f"{plural(len(dt), 'drug')} act{'s' if len(dt) == 1 else ''} on it ({int(dt.mechanism.sum())} with it as the mechanism target)")
        if len(x):
            bits.append(f"its expression differs in {plural(x.condition.nunique(), 'condition')} ({plural(len(x), 'observation')}, listed apart from the assertions)")
        return {"conditions_asserted": int(kept.condition.nunique()), "sources": sorted(set(kept.source)),
                "asserted_with": [{"id": c, "name": self.name.get(c, c), "sources": s} for c, s in top], "asserted_more": max(0, len(by) - limit),
                "withheld_rows": len(held), "withheld": [{"condition": {"id": r.condition, "name": self.name.get(r.condition, r.condition)}, **self._withheld(r)} for r in held.itertuples()][:limit],
                "drugs": drugs[:limit], "drugs_more": max(0, len(drugs) - limit),
                "expression_observations": {"rows": len(x), "conditions": int(x.condition.nunique()), "higher": int((x.log2fc > 0).sum()), "lower": int((x.log2fc < 0).sum()),
                                            "largest": [{"condition": {"id": r.condition, "name": self.name.get(r.condition, r.condition)}, "condition_as_the_source_names_it": r.label,
                                                         "direction": "higher" if r.log2fc > 0 else "lower", "log2_fold_change": round(float(r.log2fc), 3), "p": float(r.p), "edge": r.edge}
                                                        for r in strongest.itertuples()], "reading": OBSERVATION},
                "dated": False,
                "answer": ("Apart from GWAS, as of the release: " + "; ".join(bits) + "." if bits else
                           f"Apart from GWAS, no disease source, drug or expression observation in this graph names {sym}. That is not evidence that none exists.")}

    def condition(self, ids: list[str], name: str, limit: int = 10) -> dict:
        pos = lambda idx: [p for c in ids for p in idx.get(c, ())]
        a = self.A.iloc[pos(self.a_c)]
        kept, held = a[~a.withheld], a[a.withheld]
        by = kept.groupby("target").source.agg(lambda s: sorted(set(s)))
        top = sorted(by.items(), key=lambda kv: (-len(kv[1]), self.symbol_of(kv[0])))[:limit]
        di = self.DI.iloc[pos(self.di_c)]
        drugs = []
        for d in dict.fromkeys(di.drug):
            t = self.DT.iloc[self.dt_d.get(d, [])]
            drugs.append({"drug": self.name.get(d, d), "drug_id": d, "mechanism_targets": sorted(self.symbol_of(x) for x in t[t.mechanism].target),
                          "other_targets": int((~t.mechanism).sum()), "edge": di[di.drug == d].edge.iat[0]})
        drugs.sort(key=lambda r: r["drug"])
        x = self.X.iloc[pos(self.x_c)]
        bits = []
        if len(kept):
            bits.append(f"{plural(kept.source.nunique(), 'disease source')} assert{'s' if kept.source.nunique() == 1 else ''} {plural(kept.target.nunique(), 'target')} with {name}")
        if len(held):
            bits.append(f"{plural(len(held), 'row')} {'is' if len(held) == 1 else 'are'} withheld by a pruning rule")
        if drugs:
            bits.append(f"{plural(len(drugs), 'drug')} {'is' if len(drugs) == 1 else 'are'} indicated for it")
        if len(x):
            bits.append(f"expression of {plural(x.target.nunique(), 'gene')} differs in it ({plural(len(x), 'observation')}, listed apart from the assertions)")
        ties, tied = self.read_under(ids)
        return {"read_under": ties, "targets_asserted": int(kept.target.nunique()), "sources": sorted(set(kept.source)),
                "asserted_with": [{"id": t, "symbol": self.symbol_of(t), "sources": s} for t, s in top], "asserted_more": max(0, len(by) - limit),
                "withheld_rows": len(held), "drugs_indicated": drugs[:limit], "drugs_more": max(0, len(drugs) - limit),
                "expression_observations": {"rows": len(x), "genes": int(x.target.nunique()), "higher": int((x.log2fc > 0).sum()), "lower": int((x.log2fc < 0).sum()), "reading": OBSERVATION},
                "dated": False,
                "answer": ("Apart from GWAS, as of the release: " + "; ".join(bits) + "." + tied if bits else
                           f"Apart from GWAS, no disease source, drug or expression observation in this graph names {name}. That is not evidence that none exists.")}

    def drug(self, drug: str, limit: int = 10) -> dict:
        t, i = self.DT.iloc[self.dt_d.get(drug, [])], self.DI.iloc[self.di_d.get(drug, [])]
        name = self.name.get(drug, drug)
        targets = sorted(({"id": r.target, "symbol": self.symbol_of(r.target), "action": r.action or None, "activity": r.activity or None,
                           "is_the_drugs_mechanism": bool(r.mechanism), "reported_by": r.source, "edge": r.edge} for r in t.itertuples()), key=lambda r: (not r["is_the_drugs_mechanism"], r["symbol"]))
        ind = sorted(({"id": r.condition, "name": self.name.get(r.condition, r.condition), "edge": r.edge,
                       **({} if self.grounded.get(r.condition, True) else {"no_identifier_in_the_release": True})} for r in i.itertuples()), key=lambda r: r["name"])
        moa = [x["symbol"] for x in targets if x["is_the_drugs_mechanism"]]
        return {"targets": targets[:limit], "targets_more": max(0, len(targets) - limit), "indications": ind[:limit], "indications_more": max(0, len(ind) - limit), "dated": False,
                "answer": (f"As of the release, {name} acts on {plural(len(targets), 'target')}" + (f" ({', '.join(moa[:5])}{' and others' if len(moa) > 5 else ''} as its mechanism of action)" if moa else
                           " (none marked as its mechanism of action)") + f" and is indicated for {plural(len(ind), 'condition')}.")}

    def counts(self) -> dict:
        a = self.A
        return {"sources": int(a.source.nunique()), "assertions": int((~a.withheld).sum()), "withheld": int(a.withheld.sum()), "conditions": int((self.N.kind == "condition").sum()),
                "drugs": int((self.N.kind == "drug").sum()), "drug_targets": len(self.DT), "drug_indications": len(self.DI), "expression_observations": len(self.X),
                "gwas_traits_also_read_under_another_id": int(sum(len(v) > 1 for v in self.trait_ids.values())),
                "gwas_traits_read_with_a_condition_of_the_same_name": len({t for (t, c), how in self.tie.items() if how == SAME_NAME and t in self.trait_ids}),
                "conditions_under_more_than_one_id_with_one_name": int(sum(len(v) > 1 for v in self.groups.values()))}


def _selftest() -> int:
    import sys
    import tempfile
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        pd.DataFrame([("MONDO:1", "condition", "biolink:Disease", "Alpha disease", True, ""), ("UMLS:9", "condition", "biolink:Disease", "BETA", True, ""),
                      ("UMLS:7", "condition", "biolink:Disease", "Gamma", True, ""), ("MESH:7", "condition", "biolink:Disease", "gamma", True, ""),
                      ("x:indication_1", "condition", "biolink:DiseaseOrPhenotypicFeature", "Pain after surgery", False, ""), ("DrugCentral:1", "drug", "biolink:Drug", "druga", True, "")],
                     columns=["id", "kind", "category", "name", "grounded", "xref"]).to_parquet(d / FILES[0])
        pd.DataFrame([("T:A", "MONDO:1", "DisGeNET", "e1", False, "", np.nan, np.nan, np.nan, "score=0.3"), ("T:A", "UMLS:9", "eRAM", "e2", False, "", np.nan, np.nan, np.nan, ""),
                      ("T:A", "MONDO:1", "CTD", "g1", True, "literature_lift", 2.0, 1.0, 2.5, ""), ("T:B", "UMLS:7", "CTD", "e3", False, "", np.nan, np.nan, np.nan, ""),
                      ("T:B", "MESH:7", "Monarch", "e4", False, "", np.nan, np.nan, np.nan, "")],
                     columns=["target", "condition", "source", "edge", "withheld", "rule", "papers", "lift", "cut", "evidence"]).to_parquet(d / FILES[1])
        pd.DataFrame([("T:A", "MONDO:1", "active alpha", 2.0, 0.01, "x1"), ("T:A", "MONDO:1", "alpha", -1.5, 0.02, "x2")], columns=["target", "condition", "label", "log2fc", "p", "edge"]).to_parquet(d / FILES[2])
        pd.DataFrame([("DrugCentral:1", "T:A", "INHIBITOR", True, "Ki=8", "CHEMBL", "d1"), ("DrugCentral:1", "T:B", "", False, "", "DRUG LABEL", "d2")],
                     columns=["drug", "target", "action", "mechanism", "activity", "source", "edge"]).to_parquet(d / FILES[3])
        pd.DataFrame([("DrugCentral:1", "MONDO:1", "i1"), ("DrugCentral:1", "x:indication_1", "i2")], columns=["drug", "condition", "edge"]).to_parquet(d / FILES[4])
        traits = pd.DataFrame([("EFO:1", "alpha as published", "MONDO:1"), ("MONDO:2", "Beta", "")], columns=["id", "name", "live_term"])
        O = Others(d, traits, lambda t: {"T:A": "AAA", "T:B": "BBB"}.get(t, t))
    ck("a retired trait id is read under its current term too, on the ontology's statement",
       O.ids("trait", "EFO:1") == ["EFO:1", "MONDO:1"] and O.read_under(O.ids("trait", "EFO:1")) == ([{"id": "EFO:1", "tie": "as published"}, {"id": "MONDO:1", "tie": "its current ontology term"}], ""))
    ties, say = O.read_under(O.ids("trait", "MONDO:2"))
    ck("a condition the sources give a trait's name is read with it, and the answer says the tie is the name alone", [t["id"] for t in ties] == ["MONDO:2", "UMLS:9"] and ties[1]["tie"] == SAME_NAME and "UMLS:9" in say and "not an ontology statement" in say)
    ck("a condition read as part of a trait is not offered as a second thing, and its id leads to the trait", not O.match("beta", ("condition",)) and O.folded == {"mondo:1": "EFO:1", "umls:9": "MONDO:2"})
    m = O.match("gamma", ("condition", "drug"))
    ck("two ids the sources give one name are one condition, named after the one with more rows", len(m) == 1 and m[0]["tier"] == 1 and sorted(O.ids("condition", m[0]["key"])) == ["MESH:7", "UMLS:7"] and O.label("condition", m[0]["key"])["also_under"])
    ck("a drug is found by its name or its id", [x["key"] for x in O.match("druga", ("drug",))] == ["DrugCentral:1"] and O.match("drugcentral:1", ("drug",))[0]["tier"] == 0)
    p = O.pair("T:A", O.ids("trait", "EFO:1"))
    ck("for a pair: who asserts it, the withheld row with its count, the drug path with both edges, and the observations apart",
       [r["source"] for r in p["asserted_by"]] == ["DisGeNET"] and p["asserted_by"][0]["under"] == "MONDO:1" and p["asserted_by"][0]["source_evidence"] == "score=0.3"
       and p["withheld"][0] == {"source": "CTD", "gap": "g1", "rule": "literature_lift", "why": ctxgraph.RULES["literature_lift"][0], "restored_by": ctxgraph.RULES["literature_lift"][1],
                                "papers_mentioning_both": 2, "lift": 1.0, "cut_for_this_source": 2.5}
       and p["through_a_drug"] == [{"drug": "druga", "drug_id": "DrugCentral:1", "action_on_target": "INHIBITOR", "activity": "Ki=8", "target_is_the_drugs_mechanism": True, "edges": ["d1", "i1"]}]
       and [(x["direction"], x["p"]) for x in p["expression_observations"]] == [("higher", 0.01), ("lower", 0.02)] and p["dated"] is False)
    ck("the sentence says each, and that an observation is not an assertion", all(w in p["answer"] for w in ("1 source asserts AAA", "1 row from CTD is withheld", "2 papers mention both", "druga", "mechanism target", "higher in this condition in 1 observation and lower in 1", OBSERVATION)))
    q = O.pair("T:B", ["MONDO:1"])
    ck("a pair tied only through a target that is not the drug's mechanism says so", q["through_a_drug"][0]["target_is_the_drugs_mechanism"] is False and "mechanism target of none of them" in q["answer"])
    ck("a pair nothing ties is reported as not evidence of absence", "not evidence that there is no association" in O.pair("T:B", ["MONDO:2"])["answer"])
    t, c, g = O.target("T:A"), O.condition(["MONDO:1"], "alpha"), O.drug("DrugCentral:1")
    ck("a target, a condition and a drug each say what the other sources hold on them",
       (t["conditions_asserted"], t["withheld_rows"], len(t["drugs"]), t["expression_observations"]["rows"]) == (2, 1, 1, 2) and (c["targets_asserted"], c["withheld_rows"], len(c["drugs_indicated"])) == (1, 1, 1)
       and c["drugs_indicated"][0]["mechanism_targets"] == ["AAA"] and [x["symbol"] for x in g["targets"]] == ["AAA", "BBB"] and g["indications"][1].get("no_identifier_in_the_release") is True)
    ck("the counts are of what the tables hold", O.counts()["assertions"] == 4 and O.counts()["withheld"] == 1 and O.counts()["drugs"] == 1)
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(_selftest() if "--selftest" in sys.argv else 0)
