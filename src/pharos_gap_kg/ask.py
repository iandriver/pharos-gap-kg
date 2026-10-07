"""Seven questions a person can put to the graph, each answered from the bundle with its sources.

    find(text)                    which target or trait is meant
    known(node)                   what is linked to it, and since when
    why(target, trait)            the studies behind one link, or why there is none
    connect(a, b)                 what ties two things together when no link does
    predict(node)                 what is likely to be linked next, with how often that has held
    whats_new(node, since)        what was first reported after a year
    status()                      what this graph is, where it came from, and what it leaves out

THE ARGUMENTS ARE RESOLVED AGAINST THE GRAPH. A name that matches nothing returns "not in the graph",
and a name that matches several things returns the list and asks which. Nothing is guessed, so a
language model in front of these questions can route and reword but cannot supply a fact.

EVERY ANSWER IS AS OF A YEAR. A link exists from the year of the first study behind it. `as_of=2018`
answers from what had been reported by the end of 2018, and the default is everything in the release.

ABSENCE IS NOT A NEGATIVE. "No link" here means no study in this source reported one by that year.
Each answer that reports an absence says so, and says what nearby evidence exists.

Run:  python -m pharos_gap_kg.ask why APOE "Alzheimer disease"
      python -m pharos_gap_kg.ask predict "Parkinson disease" --as-of 2018
      python -m pharos_gap_kg.ask --selftest
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import graph, paths

NOT_ABSENCE = ("No link here means no study in this source had reported one by that year. "
               "It is not evidence that there is no association.")
# The extended MHC on chromosome 6 (GRCh38), where one association signal spans dozens of genes. This
# only adds a caution to an answer. It does not change a score or a rank.
MHC = ("6", 25_000_000, 34_000_000)
MHC_NOTE = "in the MHC region of chromosome 6, where one association signal spans many genes; treat neighbouring MHC candidates as one lead"
TRAIT_BINS = [(1, 10), (11, 50), (51, 100), (101, 500)]
TARGET_BINS = [(1, 5), (6, 20), (21, 50)]


def _norm(s: str) -> str:
    return " ".join(re.sub(r"[^a-z0-9]+", " ", str(s).casefold()).split())


def group_leads(ranked, related, label, top: int) -> list[dict]:
    """Ranked items gathered into leads. An item joins the first earlier lead that already holds one
    of its relatives, or an item with its label; otherwise it starts a lead. Stops when a new lead
    would be the (top+1)th. `ranked` yields (rank, item); `related(item)` is a set of items.

    Why: one GWAS signal is assigned to several neighbouring genes, and Pharos holds some genes as
    several entries. Listed one per row they fill a top ten with a single finding."""
    leads: list[dict] = []
    for rank, i in ranked:
        rel, lab = related(i), label(i)
        home = next((L for L in leads if (rel & L["members"]) or lab in L["labels"]), None)
        if home is not None:
            home["members"].add(i)
            if lab not in home["labels"]:
                home["labels"].add(lab)
                home["also"].append(lab)
            continue
        if len(leads) == top:
            break
        leads.append({"rank": rank, "first": i, "members": {i}, "labels": {lab}, "also": []})
    return leads


class Graph:
    def __init__(self, bundle_dir: Path, verify: bool = True):
        self.dir = Path(bundle_dir)
        man = self.dir / "manifest.json"
        if not man.is_file():
            raise SystemExit(f"no bundle at {self.dir}; build one with:  python -m pharos_gap_kg.bundle <release>")
        self.manifest = json.loads(man.read_text())
        if verify:
            for f, want in self.manifest["files"].items():
                if graph.file_sha(self.dir / f) != want:
                    raise SystemExit(f"REFUSED: {f} in {self.dir} is not the file its manifest describes; rebuild the bundle")
        nodes = pd.read_parquet(self.dir / "nodes.parquet")
        t = nodes[nodes.kind == "target"].reset_index(drop=True)
        d = nodes[nodes.kind == "trait"].reset_index(drop=True)
        self.t, self.d = t, d
        self.ti = {v: i for i, v in enumerate(t.id)}
        self.di = {v: i for i, v in enumerate(d.id)}
        L = pd.read_parquet(self.dir / "links.parquet")
        self.lt = L.target.map(self.ti).to_numpy(np.int64)
        self.ld = L.trait.map(self.di).to_numpy(np.int64)
        self.first, self.last = L.first_year.to_numpy(), L.last_year.to_numpy()
        self.n_studies, self.edge_id = L.n_studies.to_numpy(), L.edge_id.to_numpy()
        self.key = self.lt * len(d) + self.ld
        self.order = np.argsort(self.key)
        self.sorted_key = self.key[self.order]
        self.latest, self.earliest = int(self.manifest["latest_year"]), int(self.manifest["first_year"])
        lin = pd.read_parquet(self.dir / "lineage.parquet")
        self.lineage: dict[int, list[tuple[int, str]]] = {}
        for a, b, rel in zip(lin.a.map(self.di), lin.b.map(self.di), lin.relation):
            self.lineage.setdefault(int(a), []).append((int(b), rel))
        nb = pd.read_parquet(self.dir / "neighbours.parquet")
        self.neigh: dict[int, list[tuple[int, int]]] = {}
        for a, b, gap in zip(nb.a.map(self.ti), nb.b.map(self.ti), nb.gap):
            self.neigh.setdefault(int(a), []).append((int(b), int(gap)))
        bt = self.dir / "backtest.json"
        self.backtest = json.loads(bt.read_text()) if bt.is_file() else None
        self._ev = None
        self._sym = {}
        for i, s in enumerate(t.name):
            self._sym.setdefault(_norm(s), []).append(i)
        self._tname = [_norm(x) for x in d.name]
        self._acc = {v.split(":", 1)[1].casefold(): i for v, i in self.ti.items()}
        self.mhc = ((t.chr == MHC[0]) & (t.start >= MHC[1]) & (t.start <= MHC[2])).to_numpy() if "chr" in t else np.zeros(len(t), bool)
        self.entries = t.groupby("name").id.transform("size").to_numpy()
        self._live: dict[str, list[int]] = {}
        for i, c in enumerate(d.live_term):
            if c:
                self._live.setdefault(c.casefold(), []).append(i)

    # ------------------------------------------------------------------ pieces
    def _year(self, as_of: int | None) -> int:
        y = self.latest if as_of is None else int(as_of)
        if y < self.earliest:
            raise ValueError(f"the earliest study in this graph is from {self.earliest}; as_of={y} is before any link")
        return min(y, self.latest)

    def _mask(self, y: int) -> np.ndarray:
        return self.first <= y

    def _row(self, ti: int, di: int) -> int | None:
        k = ti * len(self.d) + di
        p = np.searchsorted(self.sorted_key, k)
        return int(self.order[p]) if p < len(self.sorted_key) and self.sorted_key[p] == k else None

    def _tlabel(self, i: int) -> dict:
        r = self.t.iloc[i]
        out = {"id": r.id, "symbol": r["name"], "name": r.full_name, "pharos_level": r.level, "family": r.family}
        if self.entries[i] > 1:
            out["pharos_entries_with_this_symbol"] = int(self.entries[i])
        return out

    def _one_per_symbol(self, idx) -> list[int]:
        """Target indices in the order given, keeping the first of each gene symbol. Pharos holds some
        genes (the HLA genes above all) as many entries with one symbol and the same links."""
        seen, out = set(), []
        for i in idx:
            s = self.t["name"][int(i)]
            if s not in seen:
                seen.add(s)
                out.append(int(i))
        return out

    def _dlabel(self, i: int) -> dict:
        r = self.d.iloc[i]
        out = {"id": r.id, "name": r["name"], "kind": "disease or phenotype" if r.disease_like else "measurement or other"}
        if r.live_term and r.live_term != r.id:
            out["current_ontology_term"] = r.live_term
        return out

    def _source(self, y: int) -> dict:
        m = self.manifest
        return {"release": m["release"], "file": m["source_file"], "published": m["source_published"], "sha256": m["source_sha256"],
                "tables": m["source_tables"], "as_of": y, "data_ends": self.latest}

    def _evidence(self, ti: int, di: int) -> pd.DataFrame:
        if self._ev is None:
            e = pd.read_parquet(self.dir / "evidence.parquet")
            e["k"] = e.target.map(self.ti).to_numpy(np.int64) * len(self.d) + e.trait.map(self.di).to_numpy(np.int64)
            self._ev = e.sort_values(["k", "year", "study"]).reset_index(drop=True)
            self._evk = self._ev.k.to_numpy()
        k = ti * len(self.d) + di
        lo, hi = np.searchsorted(self._evk, [k, k + 1])
        return self._ev.iloc[lo:hi]

    # ------------------------------------------------------------------ resolving names
    def _match(self, text: str, kind: str) -> list[dict]:
        q, raw = _norm(text), str(text).strip()
        out: list[dict] = []
        if kind in ("any", "target"):
            acc = raw.split(":", 1)[-1].casefold()
            if acc in self._acc:
                out.append({"tier": 0, "kind": "target", "i": self._acc[acc]})
            for i in self._sym.get(q, []):
                out.append({"tier": 1, "kind": "target", "i": i})
            if q and len(q) >= 3:
                for s, idx in self._sym.items():
                    if s != q and s.startswith(q):
                        out += [{"tier": 3, "kind": "target", "i": i} for i in idx]
        if kind in ("any", "trait"):
            cid = raw.replace("_", ":", 1) if re.match(r"^[A-Za-z]+_[A-Za-z0-9]+$", raw) else raw
            hit = [i for i, v in enumerate(self.d.id) if v.casefold() == cid.casefold()]
            out += [{"tier": 0, "kind": "trait", "i": i} for i in hit]
            out += [{"tier": 1, "kind": "trait", "i": i} for i in self._live.get(cid.casefold(), []) if i not in hit]
            toks = q.split()
            for i, n in enumerate(self._tname):
                if n == q:
                    out.append({"tier": 1, "kind": "trait", "i": i})
                elif toks and all(tk in n.split() for tk in toks):
                    out.append({"tier": 2, "kind": "trait", "i": i})
                elif len(q) >= 4 and q in n:
                    out.append({"tier": 3, "kind": "trait", "i": i})
        seen, uniq = set(), []
        deg_t = np.bincount(self.lt, minlength=len(self.t))
        deg_d = np.bincount(self.ld, minlength=len(self.d))
        for m in sorted(out, key=lambda m: (m["tier"], -(deg_t[m["i"]] if m["kind"] == "target" else deg_d[m["i"]]))):
            if (m["kind"], m["i"]) not in seen:
                seen.add((m["kind"], m["i"]))
                m["links"] = int(deg_t[m["i"]] if m["kind"] == "target" else deg_d[m["i"]])
                uniq.append(m)
        return uniq

    def _describe(self, m: dict) -> dict:
        lab = self._tlabel(m["i"]) if m["kind"] == "target" else self._dlabel(m["i"])
        return {"kind": m["kind"], **lab, "links": m["links"]}

    def resolve(self, text: str, kind: str = "any") -> dict:
        """{'ok': True, 'kind', 'i', 'note'} or {'ok': False, 'problem', 'candidates'}."""
        ms = self._match(text, kind)
        if not ms:
            what = {"any": "target or trait", "target": "target", "trait": "trait"}[kind]
            return {"ok": False, "problem": "not in the graph",
                    "say": f"No {what} in the graph matches {text!r}. Nothing is assumed about it.", "candidates": []}
        best = [m for m in ms if m["tier"] == ms[0]["tier"]]
        if ms[0]["tier"] <= 1:
            if len(best) == 1:
                return {"ok": True, **best[0], "note": ""}
            if all(m["kind"] == "target" for m in best) and len({self.t["name"][m["i"]] for m in best}) == 1:
                return {"ok": True, **best[0], "note": f"Pharos holds {len(best)} entries with the symbol {self.t['name'][best[0]['i']]}; this answer uses "
                        f"{self.t.id[best[0]['i']]}, which has the most links."}
            if all(m["kind"] == "trait" for m in best):
                live = {self.d.live_term[m["i"]] or self.d.id[m["i"]] for m in best}
                if len(live) == 1:
                    others = ", ".join(f"{self.d.id[m['i']]} ({m['links']} links)" for m in best[1:])
                    return {"ok": True, **best[0], "note": f"{text!r} is one ontology term under more than one id; this answer uses "
                            f"{self.d.id[best[0]['i']]}, which has the most links. Also: {others}."}
        return {"ok": False, "problem": "ambiguous",
                "say": f"{text!r} matches {len(ms)} things in the graph. Say which one, by id.",
                "candidates": [self._describe(m) for m in ms[:8]]}

    def _need(self, text: str, kind: str = "any") -> tuple[dict | None, dict | None]:
        r = self.resolve(text, kind)
        if r["ok"]:
            return r, None
        return None, {"found": False, "problem": r["problem"], "answer": r["say"], "candidates": r["candidates"]}

    # ------------------------------------------------------------------ the neighbours score
    def _sim_to_trait(self, di: int, mask: np.ndarray) -> np.ndarray:
        """Cosine similarity of every trait to trait di over the targets they share, as of the mask."""
        lt, ld = self.lt[mask], self.ld[mask]
        deg = np.bincount(ld, minlength=len(self.d)).astype(float)
        inside = np.zeros(len(self.t))
        inside[lt[ld == di]] = 1.0
        shared = np.bincount(ld, weights=inside[lt], minlength=len(self.d))
        with np.errstate(divide="ignore", invalid="ignore"):
            sim = np.where((deg > 0) & (deg[di] > 0), shared / np.sqrt(deg * deg[di]), 0.0)
        sim[di] = 0.0
        return sim

    def _score_trait(self, di: int, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        sim = self._sim_to_trait(di, mask)
        return np.bincount(self.lt[mask], weights=sim[self.ld[mask]], minlength=len(self.t)), sim

    def _flags(self, ti: int, di: int, mask: np.ndarray) -> dict:
        mine = set(self.ld[mask & (self.lt == ti)].tolist())
        lin = [(o, rel) for o, rel in self.lineage.get(di, []) if o in mine]
        holders = set(self.lt[mask & (self.ld == di)].tolist())
        pos = sorted([(o, gap) for o, gap in self.neigh.get(ti, []) if o in holders], key=lambda x: x[1])
        return {"lineage": lin, "position": pos}

    def _calibration(self, axis: str, rank: int, linked: bool, disease: bool) -> str:
        if not self.backtest or "calibration" not in self.backtest:
            return ""
        c = self.backtest["calibration"]
        bins = TRAIT_BINS if axis == "by_trait" else TARGET_BINS
        b = next((f"{lo}-{hi}" for lo, hi in bins if lo <= rank <= hi), None)
        if b is None:
            return f"Below rank {bins[-1][1]} the backtest was not tallied."
        grp = "linked" if linked else ("disease_like_not_linked" if disease and axis == "by_trait" else "not_linked")
        v = c[axis][b][grp]
        if v["share"] is None:
            return ""
        what = {"linked": "that an earlier link explains by trait wording or by position",
                "not_linked": "that no earlier link explains by wording or position",
                "disease_like_not_linked": "for a disease, that no earlier link explains by wording or position"}[grp]
        return (f"In the {c['cutoff']} backtest, {v['share']:.0%} of candidates at ranks {b} {what} were reported by {c['data_ends']} "
                f"(base rate {c['base_rate']:.1%}). A candidate not reported is not known to be absent.")

    # ------------------------------------------------------------------ the seven questions
    def find(self, text: str, kind: str = "any", limit: int = 8) -> dict:
        ms = self._match(text, kind)
        if not ms:
            return {"found": False, "answer": f"Nothing in the graph matches {text!r}.", "matches": []}
        exact = ms[0]["tier"] <= 1
        if exact:
            say = f"{text!r} names something in the graph exactly; it is listed first" + (f", with {len(ms) - 1} other name(s) that contain it." if len(ms) > 1 else ".")
        elif len(ms) == 1:
            say = f"Nothing is named exactly {text!r}. One name contains it; confirm it is the one meant before using it."
        else:
            say = (f"Nothing is named exactly {text!r}, and {len(ms)} names contain it. Ask which is meant. "
                   f"Do not pick one: they are different things.")
        return {"found": True, "exact_match": exact, "ask_which_is_meant": (not exact) and len(ms) > 1, "answer": say,
                "matches": [self._describe(m) for m in ms[:limit]], "more": max(0, len(ms) - limit)}

    def known(self, node: str, as_of: int | None = None, limit: int = 10) -> dict:
        r, err = self._need(node)
        if err:
            return err
        y = self._year(as_of)
        mask = self._mask(y)
        if r["kind"] == "target":
            rows = np.flatnonzero(mask & (self.lt == r["i"]))
            lab = self._tlabel(r["i"])
            other = lambda k: self._dlabel(int(self.ld[k]))
            noun, onoun = lab["symbol"], "trait"
        else:
            rows = np.flatnonzero(mask & (self.ld == r["i"]))
            lab = self._dlabel(r["i"])
            other = lambda k: {k2: v for k2, v in self._tlabel(int(self.lt[k])).items() if k2 in ("id", "symbol", "pharos_level", "pharos_entries_with_this_symbol")}
            noun, onoun = lab["name"], "target"
        rows = rows[np.lexsort((self.first[rows], -self.n_studies[rows]))]
        if r["kind"] == "trait":
            keep = set(self._one_per_symbol(self.lt[rows]))
            shown = [k for k in rows if int(self.lt[k]) in keep]
        else:
            shown = list(rows)
        later = int(((self.first > y) & ((self.lt == r["i"]) if r["kind"] == "target" else (self.ld == r["i"]))).sum())
        out = {"found": True, "node": {"kind": r["kind"], **lab}, "as_of": y, "links": int(len(rows)),
               "top": [{**other(k), "first_reported": int(self.first[k]), "studies_in_release": int(self.n_studies[k]),
                        "edge": self.edge_id[k]} for k in shown[:limit]],
               "more": max(0, len(shown) - limit), "first_reported_after_as_of": later, "source": self._source(y), "note": r["note"]}
        if len(rows) == 0:
            out["answer"] = (f"{noun} is in the graph with no GWAS {onoun} linked to it as of {y}. {NOT_ABSENCE}"
                             + (f" {later} link(s) were first reported after {y}." if later else ""))
        else:
            if r["kind"] == "target":
                dis = int(self.d.disease_like.to_numpy()[self.ld[rows]].sum())
                split_ = f" ({dis} to diseases or phenotypes, {len(rows) - dis} to measurements)"
            else:
                split_ = ""
            out["answer"] = (f"As of {y}, {noun} has {len(rows):,} GWAS {onoun} link(s){split_}. "
                             f"The list is ordered by how many studies report each."
                             + (f" {later} more were first reported after {y}." if later else ""))
        return out

    def why(self, target: str, trait: str, as_of: int | None = None, limit: int = 15) -> dict:
        a, err = self._need(target, "target")
        if err:
            return err
        b, err = self._need(trait, "trait")
        if err:
            return err
        y = self._year(as_of)
        mask = self._mask(y)
        ti, di = a["i"], b["i"]
        tl, dl = self._tlabel(ti), self._dlabel(di)
        k = self._row(ti, di)
        base = {"target": tl, "trait": dl, "as_of": y, "source": self._source(y), "note": " ".join(x for x in (a["note"], b["note"]) if x)}
        if k is not None and self.first[k] <= y:
            ev = self._evidence(ti, di)
            then, after = ev[ev.year <= y], ev[ev.year > y]
            papers = then.drop_duplicates("pmid")
            return {**base, "found": True, "edge": self.edge_id[k],
                    "answer": (f"{tl['symbol']} is linked to {dl['name']} because {then.study.nunique()} GWAS "
                               f"stud{'y' if then.study.nunique() == 1 else 'ies'} in {len(papers)} paper(s), the first in {int(self.first[k])}, "
                               f"reported a genome-wide significant association assigned to this gene."
                               + (f" {after.study.nunique()} more studies came after {y}." if len(after) else "")),
                    "claim": {"subject": tl["id"], "predicate": "biolink:genetically_associated_with", "object": dl["id"]},
                    "studies": [{"study": s, "pmid": int(p), "year": int(yr)} for s, p, yr in zip(then.study, then.pmid, then.year)][:limit],
                    "studies_more": max(0, len(then) - limit), "studies_after_as_of": int(after.study.nunique()),
                    "knowledge": {"primary_source": graph.SOURCE, "aggregator": graph.AGGREGATOR,
                                  "level": "statistical_association", "agent": "data_analysis_pipeline"},
                    "admitted_by": "a GWAS association the source lists for this target and trait, with a study dated by the as-of year",
                    "caveat": "GWAS assigns a signal to genes by position. A link says an associated variant lies at this gene, "
                              "not that this gene is the cause." + (f" This gene is {MHC_NOTE.split(';')[0]}." if self.mhc[ti] else "")}
        fl = self._flags(ti, di, mask)
        sc, sim = self._score_trait(di, mask)
        cand = np.ones(len(self.t), bool)
        cand[self.lt[mask & (self.ld == di)]] = False
        rank = int((sc[cand] > sc[ti]).sum()) + 1
        linked = bool(fl["lineage"] or fl["position"])
        out = {**base, "found": False,
               "reported_later": int(self.first[k]) if k is not None else None,
               "same_trait_under_another_wording": [{**self._dlabel(o), "relation_to_asked_trait": rel} for o, rel in fl["lineage"][:5]],
               "neighbouring_genes_with_this_link": [{**{kk: v for kk, v in self._tlabel(o).items() if kk in ("id", "symbol")}, "bases_apart": gap}
                                                     for o, gap in fl["position"][:5]],
               "rank_among_candidates": rank if sc[ti] > 0 else None, "candidates_for_this_trait": int(cand.sum()),
               "calibration": self._calibration("by_trait", rank, linked, bool(self.d.disease_like[di])) if sc[ti] > 0 else ""}
        if k is not None:
            head = f"As of {y} there was no link between {tl['symbol']} and {dl['name']}. One was first reported in {int(self.first[k])}."
        else:
            head = f"There is no link between {tl['symbol']} and {dl['name']} in this graph (data to {self.latest}). {NOT_ABSENCE}"
        bits = []
        if fl["lineage"]:
            o, rel = fl["lineage"][0]
            bits.append(f"{tl['symbol']} was already linked to {self.d['name'][o]}, a {rel} wording of the trait"
                        if rel != "same term" else f"{tl['symbol']} was already linked to the same term under the id {self.d.id[o]}")
        if fl["position"]:
            o, gap = fl["position"][0]
            bits.append(f"{self.t['name'][o]}, {gap:,} bases away, was linked to this trait")
        if sc[ti] > 0:
            bits.append(f"by shared-trait neighbours {tl['symbol']} ranks {rank:,} of {int(cand.sum()):,} candidates for this trait")
        out["answer"] = head + (" What nearby evidence exists: " + "; ".join(bits) + "." if bits else " Nothing nearby in the graph points either way.")
        if self.mhc[ti]:
            out["caution"] = f"{tl['symbol']} is {MHC_NOTE.split(';')[0]}, so a high rank here says little about this gene in particular."
            out["answer"] += " " + out["caution"]
        return out

    def connect(self, a: str, b: str, as_of: int | None = None, limit: int = 5) -> dict:
        x, err = self._need(a)
        if err:
            return err
        z, err = self._need(b)
        if err:
            return err
        y = self._year(as_of)
        mask = self._mask(y)
        if x["kind"] == "trait" and z["kind"] == "target":
            x, z = z, x
        base = {"as_of": y, "source": self._source(y)}
        if x["kind"] == "target" and z["kind"] == "trait":
            ti, di = x["i"], z["i"]
            tl, dl = self._tlabel(ti), self._dlabel(di)
            k = self._row(ti, di)
            direct = k is not None and self.first[k] <= y
            fl = self._flags(ti, di, mask)
            sc, sim = self._score_trait(di, mask)
            mine = self.ld[mask & (self.lt == ti)]
            via = sorted(((float(sim[o]), int(o)) for o in mine if sim[o] > 0), reverse=True)[:limit]
            return {**base, "found": True, "a": {"kind": "target", **tl}, "b": {"kind": "trait", **dl},
                    "direct_link": {"edge": self.edge_id[k], "first_reported": int(self.first[k])} if direct else None,
                    "through_trait_wording": [{**self._dlabel(o), "relation_to_b": rel} for o, rel in fl["lineage"][:limit]],
                    "through_neighbouring_genes": [{"symbol": self.t["name"][o], "id": self.t.id[o], "bases_apart": gap} for o, gap in fl["position"][:limit]],
                    "through_shared_traits": [{**self._dlabel(o), "similarity_to_b": round(s, 3)} for s, o in via],
                    "answer": (f"{tl['symbol']} and {dl['name']}: " + ("a direct link, first reported in "
                               f"{int(self.first[k])}. " if direct else f"no direct link as of {y}. ")
                               + f"Traits of {tl['symbol']} that are a broader, narrower or second wording of this one: {len(fl['lineage'])}. "
                               f"Genes within 500,000 bases that carry this link: {len(fl['position'])}. "
                               f"Traits of {tl['symbol']} that share targets with this one: {len([o for o in mine if sim[o] > 0])}."
                               + ("" if direct else f" {NOT_ABSENCE}"))}
        if x["kind"] == "target":
            ia, ib = x["i"], z["i"]
            da, db = set(self.ld[mask & (self.lt == ia)].tolist()), set(self.ld[mask & (self.lt == ib)].tolist())
            both = da & db
            deg = np.bincount(self.ld[mask], minlength=len(self.d))
            rare = sorted(both, key=lambda o: deg[o])[:limit]
            gap = next((g for o, g in self.neigh.get(ia, []) if o == ib), None)
            la, lb = self._tlabel(ia), self._tlabel(ib)
            cos = len(both) / np.sqrt(len(da) * len(db)) if da and db else 0.0
            return {**base, "found": True, "a": {"kind": "target", **la}, "b": {"kind": "target", **lb},
                    "shared_traits": len(both), "traits_a": len(da), "traits_b": len(db), "similarity": round(float(cos), 3),
                    "bases_apart_if_neighbours": gap,
                    "rarest_shared_traits": [{**self._dlabel(o), "targets_with_this_trait": int(deg[o])} for o in rare],
                    "answer": (f"As of {y}, {la['symbol']} and {lb['symbol']} share {len(both)} trait(s) "
                               f"(of {len(da)} and {len(db)})."
                               + (f" Their genes are {gap:,} bases apart, so shared GWAS links may be one signal assigned to both." if gap is not None else ""))}
        ia, ib = x["i"], z["i"]
        ta, tb = set(self.lt[mask & (self.ld == ia)].tolist()), set(self.lt[mask & (self.ld == ib)].tolist())
        both = ta & tb
        rel = next((r_ for o, r_ in self.lineage.get(ia, []) if o == ib), None)
        nst = {int(self.lt[k_]): int(self.n_studies[k_]) for k_ in np.flatnonzero(mask & (self.ld == ia))}
        top = sorted(both, key=lambda o: -nst.get(o, 0))[:limit]
        la, lb = self._dlabel(ia), self._dlabel(ib)
        cos = len(both) / np.sqrt(len(ta) * len(tb)) if ta and tb else 0.0
        return {**base, "found": True, "a": {"kind": "trait", **la}, "b": {"kind": "trait", **lb},
                "shared_targets": len(both), "targets_a": len(ta), "targets_b": len(tb), "similarity": round(float(cos), 3),
                "ontology_relation": (f"b is a {rel} wording of a" if rel and rel != "same term" else ("the same term under two ids" if rel else None)),
                "shared_targets_top": [{"symbol": self.t["name"][o], "id": self.t.id[o]} for o in top],
                "answer": (f"As of {y}, {la['name']} and {lb['name']} share {len(both):,} target(s) (of {len(ta):,} and {len(tb):,})."
                           + (" The ontology files one under the other, so shared targets are expected." if rel else ""))}

    def predict(self, node: str, as_of: int | None = None, top: int = 10, unexplained_only: bool = False) -> dict:
        r, err = self._need(node)
        if err:
            return err
        y = self._year(as_of)
        mask = self._mask(y)
        top = max(1, min(int(top), 50))
        head = ("Candidates are ranked by shared-trait neighbours: a target scores high for a trait when the traits it already has "
                "share many targets with that trait. No model is fitted.")
        bt = (self.backtest or {}).get("headline")
        if bt:
            head += (f" In the {bt['cutoff']} backtest this ranking scored {bt['neighbours']:.3f} (mean average precision over traits) against "
                     f"{bt['degree']:.3f} for how well-studied the target is alone and {bt['chance']:.3f} for chance. On links that trait "
                     f"wording and gene position do not explain, it was not shown to beat the well-studied ranking.")
        rows = []
        if r["kind"] == "trait":
            di = r["i"]
            sc, sim = self._score_trait(di, mask)
            have = np.zeros(len(self.t), bool)
            have[self.lt[mask & (self.ld == di)]] = True
            if not have.any():
                return {"found": True, "node": {"kind": "trait", **self._dlabel(di)}, "as_of": y, "candidates": [],
                        "answer": f"{self.d['name'][di]} has no link as of {y}, so there is nothing to rank from. {NOT_ABSENCE}"}
            deg = np.bincount(self.lt[mask], minlength=len(self.t))
            order = np.lexsort((np.arange(len(sc)), -deg, -sc))
            order = order[~have[order] & (sc[order] > 0)]
            n_cand = int((~have).sum())
            flags: dict[int, dict] = {}

            def ranked():
                for rank, ti in enumerate(order[:5000], 1):
                    fl = flags.setdefault(int(ti), self._flags(int(ti), di, mask))
                    if unexplained_only and (fl["lineage"] or fl["position"]):
                        continue
                    yield rank, int(ti)

            leads = group_leads(ranked(), lambda i: {o for o, _ in self.neigh.get(i, [])}, lambda i: self.t["name"][i], top)
            for L in leads:
                ti, rank = L["first"], L["rank"]
                fl = flags[ti]
                linked = bool(fl["lineage"] or fl["position"])
                mine = self.ld[mask & (self.lt == ti)]
                sup = sorted(((float(sim[o]), int(o)) for o in mine if sim[o] > 0), reverse=True)[:2]
                k = self._row(ti, di)
                later = [int(self.first[kk]) for kk in (self._row(m, di) for m in L["members"]) if kk is not None]
                rows.append({"rank": rank, **{kk: v for kk, v in self._tlabel(ti).items() if kk in ("id", "symbol", "pharos_level", "pharos_entries_with_this_symbol")},
                             "score": round(float(sc[ti]), 3), "traits_already_linked": int(deg[ti]),
                             "same_locus": L["also"][:8], "same_locus_more": max(0, len(L["also"]) - 8),
                             "caution": MHC_NOTE if self.mhc[ti] else "",
                             "explained_by": ([f"already linked to {self.d['name'][o]} ({rel})" for o, rel in fl["lineage"][:1]]
                                              + [f"{self.t['name'][o]} is {gap:,} bases away and has this link" for o, gap in fl["position"][:1]]),
                             "because_it_has": [self.d["name"][o] for _, o in sup],
                             "reported_after_as_of": int(self.first[k]) if k is not None else None,
                             "any_gene_at_this_locus_reported_after_as_of": min(later) if later else None,
                             "how_often_this_held": self._calibration("by_trait", rank, linked, bool(self.d.disease_like[di]))})
            lab = {"kind": "trait", **self._dlabel(di)}
            what = f"targets for {self.d['name'][di]}"
        else:
            ti = r["i"]
            mine = self.ld[mask & (self.lt == ti)]
            if len(mine) == 0:
                return {"found": True, "node": {"kind": "target", **self._tlabel(ti)}, "as_of": y, "candidates": [],
                        "answer": f"{self.t['name'][ti]} has no link as of {y}, so there is nothing to rank from. {NOT_ABSENCE}"}
            total = np.zeros(len(self.d))
            best: dict[int, tuple[float, int]] = {}
            for o in mine:
                s = self._sim_to_trait(int(o), mask)
                total += s
                for j in np.flatnonzero(s > 0):
                    if s[j] > best.get(int(j), (0.0, -1))[0]:
                        best[int(j)] = (float(s[j]), int(o))
            have = np.zeros(len(self.d), bool)
            have[mine] = True
            degd = np.bincount(self.ld[mask], minlength=len(self.d))
            order = np.lexsort((np.arange(len(total)), -degd, -total))
            order = order[~have[order] & (total[order] > 0)]
            n_cand = int((~have & (degd > 0)).sum())
            flags = {}

            def ranked():
                for rank, dj in enumerate(order, 1):
                    fl = flags.setdefault(int(dj), self._flags(ti, int(dj), mask))
                    if unexplained_only and (fl["lineage"] or fl["position"]):
                        continue
                    yield rank, int(dj)

            leads = group_leads(ranked(), lambda j: {o for o, _ in self.lineage.get(j, [])}, lambda j: self.d.id[j], top)
            for L in leads:
                dj, rank = L["first"], L["rank"]
                fl = flags[dj]
                linked = bool(fl["lineage"] or fl["position"])
                k = self._row(ti, dj)
                rows.append({"rank": rank, **self._dlabel(dj), "score": round(float(total[dj]), 3),
                             "other_wordings_of_this_trait": [self.d["name"][self.di[x]] for x in L["also"][:5]],
                             "explained_by": ([f"already linked to {self.d['name'][o]} ({rel})" for o, rel in fl["lineage"][:1]]
                                              + [f"{self.t['name'][o]} is {gap:,} bases away and has this link" for o, gap in fl["position"][:1]]),
                             "because_it_has": [self.d["name"][best[dj][1]]] if dj in best else [],
                             "reported_after_as_of": int(self.first[k]) if k is not None else None,
                             "how_often_this_held": self._calibration("by_target", rank, linked, False)})
            lab = {"kind": "target", **self._tlabel(ti)}
            what = f"traits for {self.t['name'][ti]}"
        came = sum(1 for x in rows if x["reported_after_as_of"] or x.get("any_gene_at_this_locus_reported_after_as_of"))
        ans = (f"The {len(rows)} highest-ranked candidate {what} as of {y}, out of {n_cand:,}."
               + (" Genes within 500,000 bases of each other are shown as one lead." if r["kind"] == "trait" else
                  " Traits the ontology files under one another are shown as one lead."))
        n_mhc = sum(1 for x in rows if x.get("caution"))
        if n_mhc:
            ans += f" In the MHC region of chromosome 6, where neighbouring candidates are best read as one lead: {n_mhc} of them."
        if y < self.latest and rows:
            ans += f" Reported between {y + 1} and {self.latest}: {came} of these {len(rows)}."
        return {"found": True, "node": lab, "as_of": y, "method": head, "answer": ans, "candidates": rows,
                "caveat": "These are leads for a person to check, not findings. " + NOT_ABSENCE, "source": self._source(y), "note": r["note"]}

    def whats_new(self, node: str, since: int, limit: int = 10) -> dict:
        r, err = self._need(node)
        if err:
            return err
        since = int(since)
        sel = (self.first > since) & ((self.lt == r["i"]) if r["kind"] == "target" else (self.ld == r["i"]))
        rows = np.flatnonzero(sel)
        rows = rows[np.lexsort((-self.n_studies[rows], self.first[rows]))]
        if r["kind"] == "trait":
            keep = set(self._one_per_symbol(self.lt[rows]))
            shown = [k for k in rows if int(self.lt[k]) in keep]
        else:
            shown = list(rows)
        lab = self._tlabel(r["i"]) if r["kind"] == "target" else self._dlabel(r["i"])
        other = (lambda k: self._dlabel(int(self.ld[k]))) if r["kind"] == "target" else \
                (lambda k: {kk: v for kk, v in self._tlabel(int(self.lt[k])).items() if kk in ("id", "symbol", "pharos_level", "pharos_entries_with_this_symbol")})
        by = {int(yr): int(n) for yr, n in zip(*np.unique(self.first[rows], return_counts=True))}
        name = lab.get("symbol") or lab["name"]
        return {"found": True, "node": {"kind": r["kind"], **lab}, "since": since, "data_ends": self.latest, "new_links": int(len(rows)),
                "by_year": by, "first": [{**other(k), "first_reported": int(self.first[k]), "studies_in_release": int(self.n_studies[k]),
                                          "edge": self.edge_id[k]} for k in shown[:limit]],
                "more": max(0, len(shown) - limit), "source": self._source(self.latest), "note": r["note"],
                "answer": (f"{len(rows):,} link(s) to {name} were first reported after {since} (data to {self.latest})."
                           if len(rows) else f"No link to {name} was first reported after {since} (data to {self.latest}). {NOT_ABSENCE}")}

    def status(self) -> dict:
        m = self.manifest
        out = {"found": True,
               "answer": (f"A graph of {m['links']:,} GWAS links between {m['targets_with_a_link']:,} of Pharos's {m['targets']:,} protein targets "
                          f"and {m['traits']:,} traits, each dated by the first study behind it ({m['first_year']} to {m['latest_year']}), "
                          f"from {m['studies']:,} studies in {m['papers']:,} papers."),
               "built_from": {"file": m["source_file"], "url": m["source_url"], "published": m["source_published"], "sha256": m["source_sha256"],
                              "tables": m["source_tables"]},
               "also_reads": {"trait ontology (sha256)": m["ontology_sha256"], "gene coordinates (sha256)": m["coordinates_sha256"]},
               "checks": [m["conformance"]] + m["dated_graphs"],
               "left_out": ["Pharos's other disease associations: they carry no dates, and one of the twelve sources is this GWAS data again.",
                            "Compound activities: after any year nearly every new target-compound pair involves a compound not seen before.",
                            f"{m['targets'] - m['targets_with_a_link']:,} targets have no link. They are in the graph, and asking about one says so."],
               "reading_a_link": "A GWAS link places an associated variant at a gene. It does not say the gene causes the trait.",
               "questions": ["find", "known", "why", "connect", "predict", "whats_new", "status"]}
        if self.backtest and "headline" in self.backtest:
            out["backtest"] = self.backtest["headline"]
        return out


def default_bundle() -> Path:
    cands = sorted(p for p in paths.BUNDLE.glob("*") if (p / "manifest.json").is_file()) if paths.BUNDLE.is_dir() else []
    if (paths.BUNDLE / "manifest.json").is_file():
        return paths.BUNDLE
    if not cands:
        raise SystemExit(f"no bundle under {paths.BUNDLE}; build one with:  python -m pharos_gap_kg.bundle pharos319")
    return cands[-1]


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        from pharos_gap_kg import ask_selftest
        return ask_selftest.run()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("question", choices=["find", "known", "why", "connect", "predict", "whats_new", "status"])
    ap.add_argument("args", nargs="*")
    ap.add_argument("--as-of", type=int)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--unexplained-only", action="store_true")
    ap.add_argument("--bundle", type=Path)
    a = ap.parse_args(argv)
    G = Graph(a.bundle or default_bundle())
    q = a.question
    if q == "status":
        out = G.status()
    elif q == "find":
        out = G.find(a.args[0])
    elif q == "known":
        out = G.known(a.args[0], a.as_of, a.top)
    elif q == "why":
        out = G.why(a.args[0], a.args[1], a.as_of)
    elif q == "connect":
        out = G.connect(a.args[0], a.args[1], a.as_of)
    elif q == "predict":
        out = G.predict(a.args[0], a.as_of, a.top, a.unexplained_only)
    else:
        out = G.whats_new(a.args[0], int(a.args[1]), a.top)
    print(json.dumps(out, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
