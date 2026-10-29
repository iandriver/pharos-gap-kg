"""Eight questions a person can put to the graph, each answered from the bundle with its sources.

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

CANDIDATES ARE RANKED BY THE SCORER THE BUNDLE NAMES. Since 11 Oct 2026 that is a combination of three
things, each computed from the links as of the year asked: shared-trait neighbours, how many traits
the target already has, and how many genes within 500,000 bases of it already have the link. Its two
weights were fitted on links to 2022 and it was tested on links first reported from 2023. A bundle
with no scorer recorded is ranked by shared-trait neighbours alone.

A TARGET WITH NO LINK IS RANKED BY POSITION ALONE, where the bundle records how that has held: for a
trait, the targets with no link that lie beside a gene that has it (`predict(trait, unlinked=True)`),
and for a target with no link, the traits its neighbouring genes are linked to. Nothing is fitted.

Run:  python -m pharos_gap_kg.ask why APOE "Alzheimer disease"
      python -m pharos_gap_kg.ask predict "Parkinson disease" --as-of 2018
      python -m pharos_gap_kg.ask predict "Parkinson disease" --unlinked
      python -m pharos_gap_kg.ask gaps "Parkinson disease"
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

from pharos_gap_kg import diseasegaps, graph, others, paths

NOT_ABSENCE = ("No link here means no study in this source had reported one by that year. "
               "It is not evidence that there is no association.")
# The extended MHC on chromosome 6 (GRCh38), where one association signal spans dozens of genes. This
# only adds a caution to an answer. It does not change a score or a rank.
MHC = ("6", 25_000_000, 34_000_000)
MHC_NOTE = "in the MHC region of chromosome 6, where one association signal spans many genes; treat neighbouring MHC candidates as one lead"
TRAIT_BINS = [(1, 10), (11, 50), (51, 100), (101, 500)]
TARGET_BINS = [(1, 5), (6, 20), (21, 50)]
NO_LINK_BINS = [(1, 5), (6, 20), (21, 50)]
SAME_SIGNAL = ("A target beside a linked gene is often the same association signal assigned to one more gene, "
               "so a link reported later would say little about this gene in particular.")


def _frac(x: float) -> str:
    from fractions import Fraction
    return str(Fraction(x).limit_denominator(1024))


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
        sc = (self.backtest or {}).get("scorer")
        self.weights = sc["weights"] if sc else None                 # of log(1 + part), for parts B1, B3, B5; None: neighbours alone
        self.by_target_combined = bool(sc) and sc.get("served_by_target") == "combiner"
        self.tg, self.tg_rates = None, None
        if "target_gaps.parquet" in self.manifest["files"]:
            self.tg = pd.read_parquet(self.dir / "target_gaps.parquet").set_index("id")
            self.tg_rates = json.loads((self.dir / "target_gaps.json").read_text())
        self.sizes = None                                            # study -> (cases, controls, individuals), where the bundle holds them
        self.proxies = {}                                            # study -> (proxy cases, proxy controls), where it has any
        if "study_sizes.parquet" in self.manifest["files"]:
            z = pd.read_parquet(self.dir / "study_sizes.parquet")
            self.sizes = {s_: (int(c), int(k), int(n)) for s_, c, k, n in zip(z.study, z.cases, z.controls, z.individuals)}
            if "proxy_cases" in z:
                self.proxies = {s_: (int(c), int(k)) for s_, c, k in zip(z.study, z.proxy_cases, z.proxy_controls) if c or k}
        self.has_others = all(f in self.manifest["files"] for f in others.FILES)   # the merged graph's other parts, where the bundle holds them
        self._others = None
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
    @property
    def others(self) -> "others.Others":
        """The disease sources, drugs and expression observations, read when a question first needs them."""
        if self._others is None:
            self._others = others.Others(self.dir, self.d[["id", "name", "live_term"]], lambda i: self.t["name"][self.ti[i]] if i in self.ti else i)
        return self._others

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
        if kind in ("any", "trait", "condition"):       # "trait" and "condition" are one kind to a reader: a disease may be either here
            cid = raw.replace("_", ":", 1) if re.match(r"^[A-Za-z]+_[A-Za-z0-9]+$", raw) else raw
            hit = [i for i, v in enumerate(self.d.id) if v.casefold() == cid.casefold()]
            out += [{"tier": 0, "kind": "trait", "i": i} for i in hit]
            out += [{"tier": 1, "kind": "trait", "i": i} for i in self._live.get(cid.casefold(), []) if i not in hit]
            if self.has_others and not hit and self.others.folded.get(cid.casefold()) in self.di:      # an id of the other sources that is read as this trait
                out.append({"tier": 0, "kind": "trait", "i": self.di[self.others.folded[cid.casefold()]]})
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
        if self.has_others and kind in ("any", "trait", "condition", "drug"):
            # a condition or a drug the other sources name. A condition that is a GWAS trait under another
            # id or the same name is that trait, and is not offered twice.
            kinds = {"any": ("condition", "drug"), "trait": ("condition",), "condition": ("condition",), "drug": ("drug",)}[kind]
            uniq = sorted(uniq + self.others.match(text, kinds), key=lambda m: m["tier"])
        return uniq

    def _describe(self, m: dict) -> dict:
        if m["kind"] in ("condition", "drug"):
            return {"kind": m["kind"], **self.others.label(m["kind"], m["key"]), "rows_in_other_sources": m["links"]}
        lab = self._tlabel(m["i"]) if m["kind"] == "target" else self._dlabel(m["i"])
        return {"kind": m["kind"], **lab, "links": m["links"]}

    def _no_dated_links(self, r: dict, question: str) -> dict:
        lab = self.others.label(r["kind"], r["key"])
        return {"found": False, "problem": "no dated links", "node": {"kind": r["kind"], **lab},
                "answer": (f"{lab['name']} is a {r['kind']} named by the other sources, and this graph holds no dated GWAS link for it. "
                           f"{question} reads the dated links only. Ask known for what the other sources say about it.")}

    def resolve(self, text: str, kind: str = "any") -> dict:
        """{'ok': True, 'kind', 'i', 'note'} or {'ok': False, 'problem', 'candidates'}."""
        ms = self._match(text, kind)
        if not ms:
            what = ({"any": "target, trait, condition or drug", "target": "target", "trait": "trait or condition", "condition": "trait or condition", "drug": "drug"} if self.has_others
                    else {"any": "target or trait", "target": "target", "trait": "trait"})[kind]
            return {"ok": False, "problem": "not in the graph",
                    "say": f"No {what} in the graph matches {text!r}. Nothing is assumed about it.", "candidates": []}
        best = [m for m in ms if m["tier"] == ms[0]["tier"]]
        if ms[0]["tier"] <= 1:
            if len(best) == 1:
                return {"ok": True, **best[0], "note": ""}
            if any(m["kind"] in ("condition", "drug") for m in best):
                # a gene symbol and a condition or drug of the other sources can share letters ("ACHE", "Ache").
                # The one written exactly as asked is meant; if that does not single one out, it is ambiguous.
                raw = str(text).strip()
                spelt = [m for m in best if self._describe(m).get("symbol", self._describe(m).get("name")) == raw]
                if len(spelt) == 1 or len({(m["kind"], self._describe(m).get("symbol")) for m in spelt}) == 1 and spelt and spelt[0]["kind"] == "target":
                    rest = [self._describe(m) for m in best if m not in spelt]
                    return {"ok": True, **spelt[0], "note": (f"{raw!r} is also, in other letter case, the name of " + ", ".join(f"the {d['kind']} {d.get('name', d.get('symbol'))} ({d['id']})" for d in rest[:3])
                                                              + ". This answer is about the one written exactly as asked.") if rest else ""}
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

    def _combine(self, nb: np.ndarray, deg, pos: np.ndarray) -> np.ndarray:
        """The served score from its three parts: shared-trait neighbours, traits already linked, and
        neighbouring genes that already have the link. Each enters as log(1 + part)."""
        w = self.weights
        return w["B3"] * np.log1p(nb) + w["B1"] * np.log1p(deg) + w["B5"] * np.log1p(pos)

    def _score_trait(self, di: int, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
        """Every target's score for trait di as of the mask, the similarity of every trait to di, and
        the parts the score is made of."""
        sim = self._sim_to_trait(di, mask)
        nb = np.bincount(self.lt[mask], weights=sim[self.ld[mask]], minlength=len(self.t))
        if self.weights is None:
            return nb, sim, {"neighbours": nb}
        deg = np.bincount(self.lt[mask], minlength=len(self.t)).astype(float)
        pos = self._nearby(di, mask)
        sc = self._combine(nb, deg, pos)
        sc[deg == 0] = 0.0          # a target with no link as of the year is not ranked: the scorer was tested on targets that have one
        return sc, sim, {"neighbours": nb, "degree": deg, "position": pos}

    def _nearby(self, di: int, mask: np.ndarray) -> np.ndarray:
        """For every target, how many genes within 500,000 bases of it have a link to trait di as of the mask."""
        pos = np.zeros(len(self.t))
        for h in self.lt[mask & (self.ld == di)]:
            for o, _ in self.neigh.get(int(h), []):
                pos[o] += 1
        return pos

    def _nearest(self, position: list[tuple[int, int]], n: int = 3) -> list[dict]:
        """The nearest linked genes, each gene symbol once (Pharos holds some genes as several entries)."""
        out, seen = [], set()
        for o, gap in position:
            s = self.t["name"][o]
            if s not in seen:
                seen.add(s)
                out.append({"symbol": s, "bases_apart": gap})
            if len(out) == n:
                break
        return out

    def _no_link_for_trait(self, di: int, mask: np.ndarray) -> tuple[list[int], np.ndarray, dict[int, int]]:
        """The targets with no link as of the mask that lie beside a gene linked to trait di, best first:
        by how many such genes, ties broken by how many links the genes beside the target hold in all.
        Returns (targets in order, the count for every target, the tie-breaker for each listed target)."""
        deg = np.bincount(self.lt[mask], minlength=len(self.t))
        pos = self._nearby(di, mask)
        listed = np.flatnonzero((deg == 0) & (pos > 0))
        locus = {int(i): int(sum(deg[o] for o, _ in self.neigh.get(int(i), []))) for i in listed}
        return sorted((int(i) for i in listed), key=lambda i: (-pos[i], -locus[i], i)), pos, locus

    def _no_link_calibration(self, axis: str, rank: int) -> str:
        c = ((self.backtest or {}).get("no_link") or {}).get("calibration")
        if not c:
            return ""
        b = next((f"{lo}-{hi}" for lo, hi in NO_LINK_BINS if lo <= rank <= hi), None)
        if b is None:
            return f"Below rank {NO_LINK_BINS[-1][1]} the test was not tallied."
        v = c[axis][b]["all"]
        if v["share"] is None:
            return ""
        what = "targets with no link listed at ranks {b} for a trait were linked to it" if axis == "by_trait" else "traits listed at ranks {b} for a target with no link were linked to it"
        return (f"In the test on links first reported from {c['first_year']} to {c['data_ends']}, {v['share']:.1%} of " + what.format(b=b) + f" by {c['data_ends']}. "
                f"Of all listed pairs {c['base_rate_among_listed']:.1%} were, and of all pairs of a target with no link and a trait {c['base_rate_among_all']:.2%}. "
                "A pair not reported is not known to be absent.")

    def _no_link_method(self, by: str) -> str:
        head = ("A target with no link is ranked by position alone: " +
                ("the number of genes within 500,000 bases of it that already have this link, ties broken by how many links those neighbouring genes hold in all."
                 if by == "trait" else
                 "for each trait, the number of genes within 500,000 bases of the target that are linked to it, ties broken by how many targets the trait is linked to.")
                + " Nothing is fitted.")
        h = ((self.backtest or {}).get("no_link") or {}).get("headline")
        if h:
            x = h["by_trait" if by == "trait" else "by_target"]
            d = x["position_minus_tie_breaker"]
            over = f"{x['traits']:,} traits" if by == "trait" else f"{x['targets']:,} targets"
            head += (f" On {h['later_links']:,} links first reported from {h['first_year']} to {h['data_ends']} to targets that had no link, this ranking scored "
                     f"{x['position']:.3f} (mean average precision over {over}) against {x['tie_breaker_alone']:.3f} for the tie-breaker alone "
                     f"({d[0]:+.4f}, 95% interval {d[1]:+.4f} to {d[2]:+.4f}) and {x['chance']:.3f} for chance, and it was above the tie-breaker in two earlier "
                     f"periods as well. Position reaches only a pair where a gene nearby has the link: {h['later_links_listed']:,} of those {h['later_links']:,} links.")
        return head

    def _ranked_by(self) -> str:
        return "shared-trait neighbours" if self.weights is None else "this graph's ranking (shared-trait neighbours, traits already linked, and linked genes nearby)"

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
        test = (f"In the test on links first reported from {c['first_year']} to {c['data_ends']}" if c.get("prospective") else f"In the {c['cutoff']} backtest")
        unit = "trait" if axis == "by_trait" else "target"
        over = f", counted over every {unit} ranked" if c.get("pooled_over") else f", counted over the {unit}s that had a later link, which makes it higher than for a {unit} taken without knowing that"
        return (f"{test}, {v['share']:.0%} of candidates at ranks {b} {what} were reported by {c['data_ends']}{over} "
                f"(base rate {c['base_rate']:.1%}). A candidate not reported is not known to be absent.")

    def _method(self, by: str, y: int) -> str:
        """What a list of candidates is ranked by, and how that ranking has held."""
        bt, pr = (self.backtest or {}).get("headline"), (self.backtest or {}).get("prospective")
        if self.weights is None or (by == "target" and not self.by_target_combined):
            head = ("Candidates are ranked by shared-trait neighbours: a target scores high for a trait when the traits it already has "
                    "share many targets with that trait. No model is fitted.")
            if bt:
                head += (f" In the {bt['cutoff']} backtest this ranking scored {bt['neighbours']:.3f} (mean average precision over traits) against "
                         f"{bt['degree']:.3f} for how well-studied the target is alone and {bt['chance']:.3f} for chance. On links that trait "
                         f"wording and gene position do not explain, it was not shown to beat the well-studied ranking.")
            return head
        w = self.weights
        head = ("Candidates are ranked by a score with three parts: shared-trait neighbours (a target scores high for a trait when the traits it "
                "already has share many targets with that trait), how many traits the target already has, and how many genes within 500,000 bases "
                f"of it already have this link. Each part enters as log(1 + part), weighted {_frac(w['B3'])}, {_frac(w['B1'])} and {_frac(w['B5'])}; the two "
                f"weights that are not 1 were fitted on links to {self.backtest['scorer']['weights_fitted_on_links_to']}.")
        if pr:
            d = pr["combiner_minus_neighbours_all"]
            head += (f" On {pr['test_links']:,} links first reported from {pr['first_year']} to {pr['data_ends']}, with every link in this release as training, "
                     f"this ranking scored {pr['combiner']:.3f} (mean average precision over {pr['traits_scored']:,} traits) against {pr['neighbours']:.3f} for "
                     f"shared-trait neighbours alone ({d[0]:+.4f}, 95% interval {d[1]:+.4f} to {d[2]:+.4f}) and {pr['degree']:.3f} for how well-studied the "
                     f"target is. The gain is on candidates an earlier link already explains by gene position or trait wording. On the others it was level "
                     f"with neighbours alone.")
            if by == "target":
                x = pr["by_target"]
                head += (f" Read the other way, ranking the traits for one target, it scored {x['combiner']:.3f} against {x['neighbours']:.3f} on the same "
                         f"links. That reading was checked after the test and was not declared before it.")
        if y < self.latest:
            head += (f" The weights were fitted on links later than {y}, so what this ranking shows as reported after {y} is an illustration and not a test.")
        return head

    # ------------------------------------------------------------------ the questions
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
        if r["kind"] == "drug":
            lab = self.others.label("drug", r["key"])
            return {"found": True, "node": {"kind": "drug", **lab}, **self.others.drug(r["key"], limit), "source": self._source(y), "note": others.UNDATED}
        if r["kind"] == "condition":
            lab = self.others.label("condition", r["key"])
            o = self.others.condition(self.others.ids("condition", r["key"]), lab["name"], limit)
            return {"found": True, "node": {"kind": "condition", **lab}, "as_of": y, "links": 0, "top": [], "other_sources": o, "source": self._source(y), "note": others.UNDATED,
                    "answer": f"{lab['name']} is not a GWAS trait in this graph, so no dated link is held for it. {o['answer']}"}
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
        if self.has_others:
            o = (self.others.target(lab["id"], limit) if r["kind"] == "target" else
                 self.others.condition(self.others.ids("trait", lab["id"]), lab["name"], limit))
            out["other_sources"] = o
            out["answer"] += " " + o["answer"]
        return out

    def why(self, target: str, trait: str, as_of: int | None = None, limit: int = 15) -> dict:
        a, err = self._need(target, "target")
        if err:
            return err
        b, err = self._need(trait, "trait")
        if err:
            return err
        y = self._year(as_of)
        tid = self.t.id[a["i"]]
        if b["kind"] == "condition":
            lab = self.others.label("condition", b["key"])
            o = self.others.pair(tid, self.others.ids("condition", b["key"]))
            return {"target": self._tlabel(a["i"]), "condition": {"kind": "condition", **lab}, "as_of": y, "source": self._source(y), "note": a["note"], "found": False,
                    "asserted_elsewhere": bool(o["asserted_by"]), "other_sources": o,
                    "answer": f"{lab['name']} is not a GWAS trait in this graph, so there is no dated link to ask about. {o['answer']}"}
        out = self._why_dated(a, b, y, limit)
        if self.has_others:
            o = self.others.pair(tid, self.others.ids("trait", self.d.id[b["i"]]))
            out["asserted_elsewhere"] = bool(o["asserted_by"])
            out["other_sources"] = o
            out["answer"] += " " + o["answer"]
        return out

    def _why_dated(self, a: dict, b: dict, y: int, limit: int) -> dict:
        """The dated GWAS link between a target and a trait, or what stands near it where there is none."""
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
        sc, sim, _ = self._score_trait(di, mask)
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
            bits.append(f"by {self._ranked_by()} {tl['symbol']} ranks {rank:,} of {int(cand.sum()):,} candidates for this trait")
        if fl["position"] and "no_link" in (self.backtest or {}) and not (mask & (self.lt == ti)).any():
            order, _, _ = self._no_link_for_trait(di, mask)
            r0 = order.index(ti) + 1
            out.update({"rank_among_targets_with_no_link": r0, "targets_with_no_link_listed_for_this_trait": len(order),
                        "calibration": self._no_link_calibration("by_trait", r0)})
            bits.append(f"{tl['symbol']} has no link of its own, and by position it ranks {r0:,} of the {len(order):,} targets with no link that lie beside a gene with this one")
        out["answer"] = head + (" What nearby evidence exists: " + "; ".join(bits) + "." if bits else " Nothing nearby in the graph points either way.")
        floor = diseasegaps.study_floor(self, self.d.id[di], y)
        if floor and "largest" in floor:
            out["what_an_absent_link_could_have_seen"] = floor
            out["answer"] += " " + floor["reading"] + " " + floor["proposal"]
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
        if {x["kind"], z["kind"]} & {"condition", "drug"}:
            return {**base, **self._connect_other(x, z, mask, limit)}
        if x["kind"] == "target" and z["kind"] == "trait":
            ti, di = x["i"], z["i"]
            tl, dl = self._tlabel(ti), self._dlabel(di)
            k = self._row(ti, di)
            direct = k is not None and self.first[k] <= y
            fl = self._flags(ti, di, mask)
            sim = self._sim_to_trait(di, mask)
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
                               + ("" if direct else f" {NOT_ABSENCE}") + self._connect_elsewhere(tl, dl)[1]),
                    **self._connect_elsewhere(tl, dl)[0]}
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

    def _connect_elsewhere(self, tl: dict, dl: dict) -> tuple[dict, str]:
        """What the other sources add between a target and a GWAS trait: who asserts the pair and the drugs
        between them, and the sentence that says so."""
        if not self.has_others:
            return {}, ""
        o = self.others.pair(tl["id"], self.others.ids("trait", dl["id"]))
        return {"asserted_by_other_sources": sorted({r["source"] for r in o["asserted_by"]}), "through_a_drug": o["through_a_drug"], "other_sources": o}, " " + o["answer"]

    def _reach(self, r: dict, mask: np.ndarray) -> tuple[dict, dict, str]:
        """A node's label, the targets it reaches with how it reaches each, and its name."""
        O = self.others
        if r["kind"] == "target":
            lab = self._tlabel(r["i"])
            return {"kind": "target", **lab}, {lab["id"]: "itself"}, lab["symbol"]
        if r["kind"] == "drug":
            lab = O.label("drug", r["key"])
            t = O.DT.iloc[O.dt_d.get(r["key"], [])]
            return {"kind": "drug", **lab}, {x.target: ("acts on it, as its mechanism" if x.mechanism else "acts on it") for x in t.itertuples()}, lab["name"]
        if r["kind"] == "trait":
            lab, ids = self._dlabel(r["i"]), O.ids("trait", self.d.id[r["i"]])
            reach = {self.t.id[i]: "a GWAS link" for i in self.lt[mask & (self.ld == r["i"])]}
        else:
            lab, ids = O.label("condition", r["key"]), O.ids("condition", r["key"])
            reach = {}
        a = O.A.iloc[[p for c in ids for p in O.a_c.get(c, ())]]
        for t, n in a[~a.withheld].groupby("target").source.nunique().items():
            reach[t] = (reach[t] + " and " if t in reach else "") + f"asserted by {n} source(s)"
        return {"kind": r["kind"], **lab}, reach, lab["name"]

    def _connect_other(self, x: dict, z: dict, mask: np.ndarray, limit: int) -> dict:
        """Two things where at least one is a condition or a drug of the other sources: what ties them."""
        O = self.others
        la, ra, na = self._reach(x, mask)
        lb, rb, nb = self._reach(z, mask)
        kinds = (x["kind"], z["kind"])
        out = {"found": True, "a": la, "b": lb, "dated": False}
        if "target" in kinds and ("condition" in kinds or "trait" in kinds):
            t, c = (x, z) if x["kind"] == "target" else (z, x)
            o = O.pair(self.t.id[t["i"]], O.ids(c["kind"], c["key"] if c["kind"] == "condition" else self.d.id[c["i"]]))
            return {**out, "other_sources": o, "answer": o["answer"]}
        if "drug" in kinds and "target" in kinds:
            d, t = (x, z) if x["kind"] == "drug" else (z, x)
            tid, nd, nt = self.t.id[t["i"]], (na if x["kind"] == "drug" else nb), (nb if x["kind"] == "drug" else na)
            row = next((r for r in O.DT.iloc[O.dt_d.get(d["key"], [])].itertuples() if r.target == tid), None)
            ind = O.DI.iloc[O.di_d.get(d["key"], [])]
            asserted = {O.A.condition.iat[p] for p in O.a_t.get(tid, ()) if not O.A.withheld.iat[p]}
            both = [{"id": c, "name": O.name.get(c, c)} for c in sorted(set(ind.condition) & asserted, key=lambda c: O.name.get(c, c))]
            return {**out, "direct": ({"action": row.action or None, "activity": row.activity or None, "target_is_the_drugs_mechanism": bool(row.mechanism), "edge": row.edge} if row is not None else None),
                    "indications_a_source_also_asserts_with_the_target": both[:limit], "indications_also_asserted": len(both),
                    "answer": ((f"{nd} acts on {nt}" + (" as its mechanism of action" if row.mechanism else f". {others.NOT_MECHANISM[:-1]}") + "." if row is not None
                                else f"The release records no action of {nd} on {nt}.")
                               + f" Of the {len(ind)} condition(s) {nd} is indicated for, a disease source also asserts {nt} with {len(both)}. {others.UNDATED}")}
        if "drug" in kinds and ("condition" in kinds or "trait" in kinds):
            d, c = (x, z) if x["kind"] == "drug" else (z, x)
            nd, nc = (na, nb) if x["kind"] == "drug" else (nb, na)
            rd, rc = (ra, rb) if x["kind"] == "drug" else (rb, ra)
            ids = set(O.ids(c["kind"], c["key"] if c["kind"] == "condition" else self.d.id[c["i"]]))
            direct = next((r.edge for r in O.DI.iloc[O.di_d.get(d["key"], [])].itertuples() if r.condition in ids), None)
            shared = sorted(set(rd) & set(rc), key=lambda t: (rd[t] != "acts on it, as its mechanism", O.symbol_of(t)))
            return {**out, "indicated": direct is not None, "indication_edge": direct,
                    "targets_of_the_drug_tied_to_the_condition": [{"id": t, "symbol": O.symbol_of(t), "drug": rd[t], "condition": rc[t]} for t in shared[:limit]],
                    "targets_tied": len(shared),
                    "answer": (f"{nd} is {'indicated for' if direct else 'not recorded as indicated for'} {nc}. Of the {len(rd)} target(s) it acts on, "
                               f"{len(shared)} {'is' if len(shared) == 1 else 'are'} tied to {nc} by a GWAS link or a disease source. {others.UNDATED}")}
        shared = sorted(set(ra) & set(rb), key=O.symbol_of)
        return {**out, "shared_targets": len(shared), "targets_a": len(ra), "targets_b": len(rb),
                "shared_targets_top": [{"id": t, "symbol": O.symbol_of(t), "a": ra[t], "b": rb[t]} for t in shared[:limit]],
                "answer": f"{na} and {nb} share {len(shared):,} target(s) (of {len(ra):,} and {len(rb):,}), counting GWAS links, what the disease sources assert and what a drug acts on. {others.UNDATED}"}

    def _predict_no_link_for_trait(self, r: dict, y: int, mask: np.ndarray, top: int) -> dict:
        """For a trait: the targets with no link that lie beside a gene that has this one, by position."""
        di = r["i"]
        lab = {"kind": "trait", **self._dlabel(di)}
        if "no_link" not in (self.backtest or {}):
            return {"found": True, "node": lab, "as_of": y, "candidates": [],
                    "answer": "This bundle does not record how position rankings for targets with no link have held, so none is given."}
        order, pos, locus = self._no_link_for_trait(di, mask)
        if not order:
            return {"found": True, "node": lab, "as_of": y, "candidates": [], "method": self._no_link_method("trait"),
                    "answer": f"No target without a link lies within 500,000 bases of a gene linked to {self.d['name'][di]} as of {y}. {NOT_ABSENCE}"}
        leads = group_leads(enumerate(order, 1), lambda i: {o for o, _ in self.neigh.get(i, [])}, lambda i: self.t["name"][i], top)
        rows = []
        for L in leads:
            ti, rank = L["first"], L["rank"]
            fl = self._flags(ti, di, mask)
            k = self._row(ti, di)
            later = [int(self.first[kk]) for kk in (self._row(m, di) for m in L["members"]) if kk is not None]
            rows.append({"rank": rank, **{kk: v for kk, v in self._tlabel(ti).items() if kk in ("id", "symbol", "pharos_level", "pharos_entries_with_this_symbol")},
                         "genes_nearby_with_this_link": int(pos[ti]), "links_held_by_genes_nearby": locus[ti],
                         "nearest": self._nearest(fl["position"]),
                         "same_locus": L["also"][:8], "same_locus_more": max(0, len(L["also"]) - 8), "caution": MHC_NOTE if self.mhc[ti] else "",
                         "reported_after_as_of": int(self.first[k]) if k is not None else None,
                         "any_gene_at_this_locus_reported_after_as_of": min(later) if later else None,
                         "how_often_this_held": self._no_link_calibration("by_trait", rank)})
        came = sum(1 for x in rows if x["any_gene_at_this_locus_reported_after_as_of"])
        ans = (f"The {len(rows)} highest-ranked targets with no link as of {y} that lie beside a gene linked to {self.d['name'][di]}, out of {len(order):,}. "
               "Targets within 500,000 bases of each other are shown as one lead.")
        if y < self.latest:
            ans += f" Reported between {y + 1} and {self.latest}: {came} of these {len(rows)}."
        return {"found": True, "node": lab, "as_of": y, "method": self._no_link_method("trait"), "answer": ans, "candidates": rows,
                "caveat": "These are leads for a person to check, not findings. " + SAME_SIGNAL + " " + NOT_ABSENCE, "source": self._source(y), "note": r["note"]}

    def _predict_for_target_with_no_link(self, r: dict, y: int, mask: np.ndarray, top: int) -> dict:
        """For a target with no link: the traits that genes beside it are linked to, by position."""
        ti = r["i"]
        lab = {"kind": "target", **self._tlabel(ti)}
        near = np.zeros(len(self.d))
        for o, _ in self.neigh.get(ti, []):
            np.add.at(near, self.ld[mask & (self.lt == o)], 1)
        if "no_link" not in (self.backtest or {}) or not near.any():
            why_not = ("and no gene within 500,000 bases of it has one, " if "no_link" in (self.backtest or {}) else "")
            return {"found": True, "node": lab, "as_of": y, "candidates": [],
                    "answer": f"{self.t['name'][ti]} has no link as of {y}, {why_not}so there is nothing to rank from. {NOT_ABSENCE}"}
        degd = np.bincount(self.ld[mask], minlength=len(self.d))
        order = np.lexsort((np.arange(len(near)), -degd, -near))
        order = order[near[order] > 0]
        leads = group_leads(((rank, int(dj)) for rank, dj in enumerate(order, 1)), lambda j: {o for o, _ in self.lineage.get(j, [])}, lambda j: self.d.id[j], top)
        rows = []
        for L in leads:
            dj, rank = L["first"], L["rank"]
            fl = self._flags(ti, dj, mask)
            k = self._row(ti, dj)
            rows.append({"rank": rank, **self._dlabel(dj), "genes_nearby_with_this_link": int(near[dj]), "targets_linked_to_this_trait": int(degd[dj]),
                         "nearest": self._nearest(fl["position"]),
                         "other_wordings_of_this_trait": [self.d["name"][self.di[x]] for x in L["also"][:5]],
                         "reported_after_as_of": int(self.first[k]) if k is not None else None,
                         "how_often_this_held": self._no_link_calibration("by_target", rank)})
        came = sum(1 for x in rows if x["reported_after_as_of"])
        ans = (f"{self.t['name'][ti]} has no link as of {y}. The {len(rows)} highest-ranked traits that a gene within 500,000 bases of it is linked to, out of "
               f"{int((near > 0).sum()):,}. Traits the ontology files under one another are shown as one lead.")
        if self.mhc[ti]:
            ans += f" {self.t['name'][ti]} is {MHC_NOTE.split(';')[0]}."
        if y < self.latest:
            ans += f" Reported between {y + 1} and {self.latest}: {came} of these {len(rows)}."
        return {"found": True, "node": lab, "as_of": y, "method": self._no_link_method("target"), "answer": ans, "candidates": rows,
                "caveat": "These are leads for a person to check, not findings. " + SAME_SIGNAL + " " + NOT_ABSENCE, "source": self._source(y), "note": r["note"]}

    def predict(self, node: str, as_of: int | None = None, top: int = 10, unexplained_only: bool = False, unlinked: bool = False) -> dict:
        r, err = self._need(node)
        if err:
            return err
        if r["kind"] in ("condition", "drug"):
            return self._no_dated_links(r, "predict")
        y = self._year(as_of)
        mask = self._mask(y)
        top = max(1, min(int(top), 50))
        if unlinked and r["kind"] == "trait":
            return self._predict_no_link_for_trait(r, y, mask, top)
        if r["kind"] == "target" and not (mask & (self.lt == r["i"])).any():
            return self._predict_for_target_with_no_link(r, y, mask, top)
        head = self._method(r["kind"], y)
        rows = []
        if r["kind"] == "trait":
            di = r["i"]
            sc, sim, parts = self._score_trait(di, mask)
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
                             **({"score_parts": {"shared_trait_neighbours": round(float(parts["neighbours"][ti]), 3), "traits_already_linked": int(deg[ti]),
                                                 "genes_nearby_with_this_link": int(parts["position"][ti])}} if self.weights else {}),
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
            n_no_link = (int(((deg == 0) & (self._nearby(di, mask) > 0)).sum()) if "no_link" in (self.backtest or {}) else 0)
        else:
            ti = r["i"]
            mine = self.ld[mask & (self.lt == ti)]
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
            near = np.zeros(len(self.d))
            combined = self.weights is not None and self.by_target_combined
            if combined:
                for o, _ in self.neigh.get(ti, []):
                    np.add.at(near, self.ld[mask & (self.lt == o)], 1)
            neighbours_alone, evidence = total, (total > 0) | (near > 0)
            if combined:
                total = self._combine(neighbours_alone, len(mine), near)
            order = np.lexsort((np.arange(len(total)), -degd, -total))
            order = order[~have[order] & evidence[order]]
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
                             **({"score_parts": {"shared_trait_neighbours": round(float(neighbours_alone[dj]), 3), "traits_already_linked": int(len(mine)),
                                                 "genes_nearby_with_this_link": int(near[dj])}} if combined else {}),
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
        extra = {}
        if r["kind"] == "trait":
            extra["reading"] = ("Each row is a lead under the rank of its first gene. A gene listed under same_locus has its own rank, "
                                "which why(gene, trait) gives.")
        if r["kind"] == "trait" and n_no_link:
            extra["targets_with_no_link_beside_a_linked_gene"] = n_no_link
            ans += f" A target with no link is not in this list; {n_no_link:,} of them lie beside a gene with this link, and unlinked=true ranks those."
        return {"found": True, "node": lab, "as_of": y, "method": head, "answer": ans, "candidates": rows, **extra,
                "caveat": "These are leads for a person to check, not findings. " + NOT_ABSENCE, "source": self._source(y), "note": r["note"]}

    def whats_new(self, node: str, since: int, limit: int = 10) -> dict:
        r, err = self._need(node)
        if err:
            return err
        if r["kind"] in ("condition", "drug"):
            return self._no_dated_links(r, "whats_new")
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

    def _gap_row(self, tid: str) -> dict | None:
        return self.tg.loc[tid].to_dict() if tid in self.tg.index else None

    def _edge_of(self, g: dict, kind: str) -> dict:
        """The graph edge that states this gap, where the gaps were exported before the bundle was built."""
        i = g.get(f"{kind}_gap_id")
        return {"graph_edge": i} if isinstance(i, str) and i else {}

    def _dark_gap(self, g: dict) -> dict:
        R = self.tg_rates
        b = R["dark"]["bands"][g["band"]]
        pr, fl = b["promoted"], b["fell_in_the_mirror"]
        return {**self._edge_of(g, "dark"), "lacks": "enough recorded evidence to leave Tdark", "criteria_to_clear": int(g["flips"]),
                "evidence_now": {"pubmed_score": round(float(g["pubmed_score"]), 3), "generifs": int(g["generifs"]), "antibodies": int(g["antibodies"])},
                "nearest_route": g["proposal"], "other_route": "an experimental Gene Ontology function or process annotation, which has no distance",
                "band": g["band"],
                "how_often_this_closed": (f"Of {pr['n']:,} targets in this band in the {R['earlier']} database, {pr['share']:.0%} "
                                          f"({pr['low']:.0%} to {pr['high']:.0%}) were no longer Tdark in {R['later']}. "
                                          f"In the mirror band, {fl['share']:.0%} of {fl['n']:,} targets fell the other way, which is what "
                                          f"the pipeline alone moves at this distance.")}

    def _ligand_gap(self, g: dict) -> dict:
        R = self.tg_rates["ligand"]
        return {**self._edge_of(g, "ligand"), "lacks": f"a compound at {g['cutoff']} or better, the potency cutoff for its family",
                "distance": "not recorded: the release keeps only compounds that pass the cutoff",
                "interaction_partners_with_a_ligand": int(g["partners_with_ligand"]),
                "how_often_this_closed": (f"Of {R['candidates']:,} targets with no qualifying compound in the {self.tg_rates['earlier']} database, "
                                          f"{R['gained']:,} ({R['gained'] / R['candidates']:.0%}) held one in {self.tg_rates['later']}.")}

    def gaps(self, target: str | None = None, family: str | None = None, top: int | None = None, disease: str | None = None, as_of: int | None = None,
             kind: str | None = None) -> dict:
        """What a target lacks for the next Pharos development level; with no target, the Tdark targets
        nearest to leaving it; and for a trait or a condition, the targets tied to it and what each lacks."""
        text = disease if disease is not None else (target if family is None else None)
        if text is not None:
            r = self.resolve(text, "trait" if disease is not None else "any")
            if r["ok"] and r["kind"] in ("trait", "condition"):
                if not self.has_others:
                    return {"found": False, "problem": "not built", "answer": "This bundle does not hold the merged graph's other parts, which the gaps of a condition are read from. "
                            "They are added by:  python -m pharos_gap_kg.kg build <release>  and then  python -m pharos_gap_kg.bundle <release> --context-only"}
                return {**diseasegaps.for_condition(self, r, self._year(as_of), top, kind), "note": r["note"]}
            if r["ok"] and r["kind"] == "drug":
                return {"found": False, "problem": "not a target or a condition", "answer": f"{text!r} is a drug. gaps reads a target or a condition; ask known for what the graph holds on a drug."}
            if disease is not None:
                return {"found": False, "problem": r["problem"], "answer": r["say"], "candidates": r["candidates"]}
        reading = ("A gap is what Pharos's grade lacks, read from the release's own evidence. A PubMed score or an antibody count "
                   "measures attention to a protein. Closing the gap is not the same as learning its function.")
        if self.tg is None:
            return {"found": False, "problem": "not built", "answer": "This bundle holds no target-level gaps. They are added by:  "
                    "python -m pharos_gap_kg.tgaps <release>  and then  python -m pharos_gap_kg.bundle <release> --gaps-only"}
        R = self.tg_rates
        if target is None:
            d = self.tg[self.tg.flips.notna()]
            if family:
                d = d[d.family.fillna("").str.casefold() == family.casefold()]
                if d.empty:
                    return {"found": False, "problem": "no such family", "answer": f"No Tdark target with a gap has the family {family!r}. "
                            f"Families: {', '.join(sorted(set(self.tg.family.dropna())))}."}
            d = d.sort_values(["cost", "id"]).head(max(1, min(int(top or 10), 50)))
            rows = [{**self._tlabel(self.ti[i]), "nearest_route": r.proposal, "band": r.band,
                     "promoted_in_this_band": f"{R['dark']['bands'][r.band]['promoted']['share']:.0%}"} for i, r in d.iterrows() if i in self.ti]
            return {"found": True, "answer": (f"The {len(rows)} Tdark targets{' in the family ' + family if family else ''} nearest to leaving Tdark, of "
                                             f"{R['dark']['targets']:,} with a gap. Nearest first."),
                    "targets": rows, "how_this_ranking_did": (f"In the back-test this distance ranked a promoted target above an unpromoted one "
                                                               f"{R['dark']['route_cost_auc']:.0%} of the time ({R['dark']['best_single_count']} alone: "
                                                               f"{R['dark']['best_single_count_auc']:.0%}); {R['dark']['promoted']:,} of {R['dark']['in_both_releases']:,} were promoted."),
                    "reading": reading}
        r, miss = self._need(target, "target")
        if miss:
            return miss
        i = r["i"]
        lab, level = self._tlabel(i), self.t.level[i]
        out = {"found": True, "target": lab, "note": r["note"], "gaps": [], "reading": reading}
        g = self._gap_row(self.t.id[i])
        if level == "Tclin":
            out["answer"] = f"{lab['symbol']} is Tclin: it has an approved drug with a known mechanism, the top level. It has no gap of this kind."
        elif level == "Tchem":
            out["answer"] = (f"{lab['symbol']} is Tchem: it holds a compound that passes its family's potency cutoff. What it lacks for Tclin is an "
                             f"approved drug acting through it, which this release does not grade by distance.")
        elif g is None:
            out["answer"] = (f"{lab['symbol']} is graded {level}, yet the release's ligand table holds compound activities for it. The grade and the "
                             f"table disagree, so no gap is given.")
        else:
            if pd.notna(g["flips"]):
                out["gaps"].append({"kind": "dark", **self._dark_gap(g)})
            if g["lacks_ligand"] is True:
                out["gaps"].append({"kind": "ligand", **self._ligand_gap(g)})
            if level == "Tdark" and pd.isna(g["flips"]):
                out["answer"] = (f"{lab['symbol']} is graded Tdark, but its recorded evidence does not meet the Tdark rule, so no distance is given. "
                                 f"It lacks a compound at {g['cutoff']}.")
            elif level == "Tdark":
                d = out["gaps"][0]
                out["answer"] = (f"{lab['symbol']} is Tdark. To leave it needs {d['criteria_to_clear']} criteri{'on' if d['criteria_to_clear'] == 1 else 'a'} "
                                 f"to clear; the nearest route is: {d['nearest_route']}. {d['how_often_this_closed']}")
            else:
                out["answer"] = (f"{lab['symbol']} is {level}: it has biological evidence and no qualifying compound. "
                                 f"It lacks {out['gaps'][-1]['lacks']}. {out['gaps'][-1]['how_often_this_closed']}")
        return out

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
               "left_out": [("One of Pharos's twelve disease sources, which is this GWAS data again." if self.has_others else
                             "Pharos's other disease associations: they carry no dates, and one of the twelve sources is this GWAS data again."),
                            "Compound activities: after any year nearly every new target-compound pair involves a compound not seen before.",
                            f"{m['targets'] - m['targets_with_a_link']:,} targets have no link. They are in the graph, and asking about one says so."],
               "reading_a_link": "A GWAS link places an associated variant at a gene. It does not say the gene causes the trait.",
               "questions": ["find", "known", "why", "connect", "predict", "whats_new", "gaps", "status"]}
        if self.has_others:
            c = self.others.counts()
            out["other_sources"] = {**c, "dated": False, "reading": others.OBSERVATION}
            out["answer"] += (f" Beside them, as of the release and undated: {c['assertions']:,} assertions of a target with a condition from {c['sources']} disease sources, "
                              f"{c['withheld']:,} more rows withheld by a pruning rule and kept as gaps, {c['drugs']:,} drugs with {c['drug_targets']:,} target and "
                              f"{c['drug_indications']:,} indication edges, and {c['expression_observations']:,} expression observations, which are kept apart from the assertions.")
        if self.backtest and "headline" in self.backtest:
            out["backtest"] = self.backtest["headline"]
        if self.backtest and "scorer" in self.backtest:
            out["ranking"] = {"scorer": self.backtest["scorer"], "tested_on_later_links": self.backtest["prospective"]}
            if "no_link" in self.backtest:
                out["ranking"]["targets_with_no_link"] = {"ranked_by": "position alone, nothing fitted", **self.backtest["no_link"]["headline"]}
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
    ap.add_argument("question", choices=["find", "known", "why", "connect", "predict", "whats_new", "gaps", "status"])
    ap.add_argument("--family")
    ap.add_argument("args", nargs="*")
    ap.add_argument("--as-of", type=int)
    ap.add_argument("--top", type=int)
    ap.add_argument("--kind", help="gaps, for a disease: show one list in full, such as no_genetic_support")
    ap.add_argument("--unexplained-only", action="store_true")
    ap.add_argument("--unlinked", action="store_true", help="predict, for a trait: rank the targets that have no link and lie beside a gene that has this one")
    ap.add_argument("--bundle", type=Path)
    a = ap.parse_args(argv)
    G = Graph(a.bundle or default_bundle())
    q = a.question
    if q == "status":
        out = G.status()
    elif q == "find":
        out = G.find(a.args[0])
    elif q == "known":
        out = G.known(a.args[0], a.as_of, a.top or 10)
    elif q == "why":
        out = G.why(a.args[0], a.args[1], a.as_of)
    elif q == "connect":
        out = G.connect(a.args[0], a.args[1], a.as_of)
    elif q == "predict":
        out = G.predict(a.args[0], a.as_of, a.top or 10, a.unexplained_only, a.unlinked)
    elif q == "gaps":
        out = G.gaps(a.args[0] if a.args else None, a.family, a.top, as_of=a.as_of, kind=a.kind)
    else:
        out = G.whats_new(a.args[0], int(a.args[1]), a.top or 10)
    print(json.dumps(out, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
