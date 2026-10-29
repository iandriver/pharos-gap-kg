"""A study published after the graph: what the graph said of Alzheimer disease, against the GWAS of 28 Oct 2026.

Uffelmann et al., Nature, 28 Oct 2026 (doi:10.1038/s41586-026-11206-y): 127 loci, 48 of them new, each
with its predicted effector gene(s). Everything the graph holds predates it. Four tests, defined in
docs/PROTOCOL.md, "A study published after the graph: Alzheimer disease", before any gene was matched:

    1  gaps that close     do the genes the graph ties with no GWAS link gain one more often than the rest?
    2  a locus's order     is the gene the graph's other evidence points to first one the study names?
    3  the largest study   could the largest study behind the graph's links have found the new loci?
    4  the ranking         where does the served ranking put the genes the study names?

Run:  python -m pharos_gap_kg.adstudy fetch [--from <file>]     # the pinned Supplementary Data
      python -I -m pharos_gap_kg.adstudy audit pharos319         # -> docs/data/adstudy_pharos319.json
      python -m pharos_gap_kg.adstudy --selftest

The file is an OpenDocument spreadsheet, read here with the standard library. It is licensed CC
BY-NC-ND 4.0 and stays with the raw files; the record keeps gene symbols, counts and figures.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import sys
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import clingen, diseasegaps as dg, paths, position, registry, score, studysize

KEY = "uffelmann2026_supp"
TRAIT, TRAIT_NAME = "MONDO:0004975", "Alzheimer disease"
OTHER_TRAITS = ("EFO:1001870", "EFO:0009268")        # late-onset Alzheimer disease; family history of it
LOCI, NEW_LOCI = 127, 48
STUDY_NEFF = 416_493                                  # the paper's total effective sample size ("Data overview")
SEED, DRAWS = 20261028, 2000
GROUPS = ("tied_without_a_link", "one_source", "withheld_only", "observed_only", "nothing")
READINGS = ("primary", "broad", "data_only")
EVIDENCE = ("FLAMES", "Exome", "Coloc_MR")
Z_ALPHA, Z_POWER = 1.959964, 0.841621
_T = "urn:oasis:names:tc:opendocument:xmlns:table:1.0"
_X = "urn:oasis:names:tc:opendocument:xmlns:text:1.0"
_O = "urn:oasis:names:tc:opendocument:xmlns:office:1.0"


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def fetch(src: Path | None = None) -> Path:
    """The pinned file, from the publisher or from a copy already downloaded; refused unless its size and digest are the registry's."""
    s = registry.SOURCES[KEY]
    side = s.path.with_name(s.filename + ".lineage.json")
    if s.path.is_file() and side.is_file() and sha(s.path.read_bytes()) == s.sha256:
        print(f"{s.filename}: already fetched ({s.sha256[:12]})")
        return s.path
    modified = ""
    if src is None:
        with urllib.request.urlopen(urllib.request.Request(s.url, headers=clingen.UA), timeout=120) as r:
            data, modified = r.read(), r.headers.get("Last-Modified", "")
    else:
        data = Path(src).read_bytes()
    if len(data) != s.bytes or sha(data) != s.sha256:
        raise SystemExit(f"REFUSED {s.filename}: {len(data):,} bytes, sha256 {sha(data)[:16]}; the registry pins {s.bytes:,} bytes, {s.sha256[:16]}")
    s.path.parent.mkdir(parents=True, exist_ok=True)
    tmp = s.path.with_name(s.filename + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, s.path)
    side.write_text(json.dumps({"key": s.key, "url": s.url, "bytes": len(data), "sha256": s.sha256, "server_last_modified": modified,
                                "published": s.published, "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                                "from": str(src) if src else s.url}, indent=1) + "\n")
    print(f"{s.filename}: {len(data):,} bytes, sha256 {s.sha256}")
    return s.path


# ---------------------------------------------------------------- the study's side

def sheet(blob: bytes, name: str) -> list[list[str]]:
    """One table of an OpenDocument spreadsheet as rows of text. A numeric cell gives its full value."""
    root = ET.fromstring(zipfile.ZipFile(io.BytesIO(blob)).read("content.xml"))
    for t in root.iter(f"{{{_T}}}table"):
        if t.get(f"{{{_T}}}name") != name:
            continue
        rows = []
        for r in t.iter(f"{{{_T}}}table-row"):
            cells = []
            for c in r:
                if c.tag not in (f"{{{_T}}}table-cell", f"{{{_T}}}covered-table-cell"):
                    continue
                v = c.get(f"{{{_O}}}value")
                text = v if v is not None else "\n".join("".join(p.itertext()) for p in c.findall(f"{{{_X}}}p"))
                cells += [text] * min(int(c.get(f"{{{_T}}}number-columns-repeated", "1")), 1000)
            while cells and cells[-1] == "":
                cells.pop()
            if cells:                                  # an empty row repeated a million times is padding
                rows += [cells] * int(r.get(f"{{{_T}}}number-rows-repeated", "1"))
        return rows
    raise SystemExit(f"REFUSED: the spreadsheet has no table {name!r}")


def genes(cell: str) -> list[str]:
    """Gene symbols of one cell: separated by '/', 'NA' for none; a gene cluster names no gene."""
    out = []
    for g in str(cell).split("/"):
        g = g.strip()
        if g and g != "NA" and "cluster" not in g.lower() and g not in out:
            out.append(g)
    return out


def study_loci(rows: list[list[str]]) -> pd.DataFrame:
    """SD2: one row per locus, its genes by reading, and the lead SNP's frequency and effect."""
    head = next(i for i, r in enumerate(rows) if len(r) > 1 and r[1] == "Predicted effector gene(s)")
    cols = rows[head]
    body = [(r + [""] * len(cols))[:len(cols)] for r in rows[head + 1:] if len(r) > 2 and r[2].strip()]
    t = pd.DataFrame(body, columns=cols)
    num = lambda c: pd.to_numeric(t[c], errors="coerce")
    out = pd.DataFrame({"locus": t["locus"].astype(int), "new": t["Novel"].astype(int) == 1,
                        "cell": t["Predicted effector gene(s)"], "effector": t["Predicted effector gene(s)"].map(genes),
                        "evidence": [sum((genes(r[c]) for c in EVIDENCE), []) for _, r in t.iterrows()],
                        "closest": t["ClosestGene"].map(genes),
                        "eaf": num("effect_allele_frequency"), "beta": num("beta"), "odds_ratio": num("OR"), "neff": num("neff")})
    if len(out) != LOCI or int(out.new.sum()) != NEW_LOCI or out.locus.duplicated().any():
        raise SystemExit(f"REFUSED: SD2 gives {len(out)} loci, {int(out.new.sum())} new; the paper reports {LOCI} and {NEW_LOCI}")
    if out[out.new][["eaf", "beta"]].isna().any().any():
        raise SystemExit("REFUSED: a new locus has no effect allele frequency or beta in SD2")
    return out


# ---------------------------------------------------------------- matching

def symbol_map(G, release: str) -> dict[str, set[str]]:
    """Study symbol -> the Pharos symbols it names: a target's own symbol, and the Ensembl 113 symbol of
    every Ensembl gene the target is cross-referenced to (as position.target_genes places it)."""
    m: dict[str, set[str]] = {}
    name = dict(zip(G.t.id, G.t["name"]))
    for s in name.values():
        m.setdefault(s, set()).add(s)
    gtf = position.genes(registry.require("ensembl_gtf_113"))
    tg = position.target_genes(release, gtf)
    ens = dict(zip(gtf.ensg, gtf.symbol))
    for t, e in zip(tg.target, tg.ensg):
        if t in name and ens.get(e):
            m.setdefault(ens[e], set()).add(name[t])
    return m


def named(L: pd.DataFrame, m: dict[str, set[str]]) -> dict[str, pd.Series]:
    """For each reading, each locus's named graph genes (Pharos symbols)."""
    to = lambda gs: set().union(*(m.get(g, set()) for g in gs)) if gs else set()
    eff = L.effector.map(to)
    ev = L.evidence.map(to)
    return {"primary": eff, "broad": pd.Series([e | v | to(c) for e, v, c in zip(eff, ev, L.closest)], index=L.index),
            "data_only": pd.Series([e & v for e, v in zip(eff, ev)], index=L.index)}


# ---------------------------------------------------------------- the graph's side

def graph_sets(G, r: dict, y: int) -> dict:
    """The sets diseasegaps.for_condition builds for a trait, in full, as target ids."""
    O, di = G.others, r["i"]
    ids, mask = O.ids("trait", G.d.id[di]), G._mask(y)
    gw = {G.t.id[G.lt[k]]: int(k) for k in np.flatnonzero(mask & (G.ld == di))}
    pos = lambda idx: [p for c in ids for p in idx.get(c, ())]
    a = O.A.iloc[pos(O.a_c)]
    kept, held = a[~a.withheld], a[a.withheld]
    src = {t: sorted(set(s)) for t, s in kept.groupby("target").source}
    two = {t for t, s in src.items() if len(s) >= dg.MIN_SOURCES}
    one = {t for t, s in src.items() if len(s) < dg.MIN_SOURCES}
    moa: dict[str, list[str]] = {}
    for d in dict.fromkeys(O.DI.iloc[pos(O.di_c)].drug):
        for row in O.DT.iloc[O.dt_d.get(d, [])].itertuples():
            if row.mechanism:
                moa.setdefault(row.target, []).append(d)
    x = O.X.iloc[pos(O.x_c)]
    return {"gw": gw, "two": two, "one": one, "moa": moa, "src": src, "held": set(held.target), "held_rows": len(held), "x": x, "mask": mask}


def groups(G, S: dict) -> dict[str, set[str]]:
    """Every graph gene (Pharos symbol) with no GWAS link to the trait, in the first group it qualifies for."""
    sym = lambda t: G.t["name"][G.ti[t]]
    of = lambda ts: {sym(t) for t in ts if t in G.ti}
    linked = of(S["gw"])
    out, taken = {}, set(linked)
    for g, ts in zip(GROUPS[:-1], ((S["two"] | set(S["moa"])), S["one"], S["held"], set(S["x"].target))):
        out[g] = of(ts) - taken
        taken |= out[g]
    out["nothing"] = {s for s in G.t["name"] if isinstance(s, str) and s} - taken
    out["linked"] = linked
    return out


def crosscheck(G, r: dict, y: int, S: dict, L_graph: list[dict]) -> dict:
    """The served answer's counts and first 25 rows of each list must agree with the sets computed here."""
    served = dg.for_condition(G, r, y, top=25)
    sym = lambda t: G.t["name"][G.ti[t]]
    ng = {sym(t) for t in (S["two"] | set(S["moa"])) & set(G.ti) - set(S["gw"])}
    x = S["x"]
    tied_one = (set(S["gw"]) | S["two"] | set(S["moa"]) | S["one"]) & set(G.ti)
    only = {sym(t) for t in x.target if t in G.ti and t not in tied_one}
    mine = {"no_genetic_support": len(ng), "locus_not_resolved": len(L_graph), "withheld": S["held_rows"], "observed_only": len(only)}
    by = {g["kind"]: g for g in served["gaps"]}
    bad = [k for k, v in mine.items() if served["counts"][k] != v]
    bad += [f"no_genetic_support row {row['symbol']}" for row in by["no_genetic_support"]["top"] if row["symbol"] not in ng]
    orders = [tuple(l["genes"][:8]) for l in L_graph]
    bad += [f"locus row {row['genes']}" for row in by["locus_not_resolved"]["top"] if tuple(row["genes"]) not in orders]
    if bad:
        raise SystemExit(f"REFUSED: the sets computed here disagree with the served gaps answer: {bad[:5]}")
    return {"served_counts": served["counts"], "agree": True, "rows_compared": {k: len(by[k]["top"]) for k in ("no_genetic_support", "locus_not_resolved")}}


def graph_loci(G, S: dict) -> list[dict]:
    """The locus_not_resolved loci, each with its genes in the served order, the leading level, and the genes with the most studies."""
    gw, sym = S["gw"], (lambda i: G.t["name"][i])
    studies = lambda t: int(G.n_studies[gw[t]])
    out = []
    for members in dg.loci([G.ti[t] for t in gw], G.neigh, sym):
        ids_ = [G.t.id[i] for i in members]
        hold = {t: dg.holds(t, S["two"], S["moa"], S["src"], S["x"], lambda t_: {})[0] for t in ids_}
        order = sorted(ids_, key=lambda t: (hold[t], -studies(t), sym(G.ti[t])))
        level = min(hold.values())
        top_n = max(studies(t) for t in ids_)
        other = [sym(G.ti[t]) for t in ids_ if t in S["two"] or t in S["moa"]]
        out.append({"genes": [sym(G.ti[t]) for t in order], "level": level,
                    "first": [sym(G.ti[t]) for t in order if hold[t] == level] if level < 3 else [],
                    "most_studies": sorted(sym(G.ti[t]) for t in ids_ if studies(t) == top_n),
                    "singled_out": other[0] if len(other) == 1 else None})
    return out


# ---------------------------------------------------------------- measures

def detectable(n1: int, n0: int, p0: float) -> float | None:
    """The smallest rise over p0 that a two-sided test at 0.05 detects with 80% power, n1 and n0 genes in the two groups."""
    if n1 == 0 or n0 == 0:
        return None
    lo, hi = p0, 1.0
    for _ in range(80):
        p1 = (lo + hi) / 2
        se = math.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
        lo, hi = (lo, p1) if se > 0 and (p1 - p0) / se >= Z_ALPHA + Z_POWER else (p1, hi)
    return round(hi - p0, 5)


def below(s: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """For each positive, the share of negatives scored below it, ties counted half."""
    neg = np.sort(s[~pos])
    lo = np.searchsorted(neg, s[pos], "left")
    hi = np.searchsorted(neg, s[pos], "right")
    return (lo + (hi - lo) / 2) / len(neg)


def interval(x: np.ndarray) -> list:
    x = x[~np.isnan(x)]
    return [round(float(np.percentile(x, 2.5)), 4), round(float(np.percentile(x, 97.5)), 4)]


def test_gaps(Gr: dict[str, set[str]], N: dict[str, pd.Series], L: pd.DataFrame, other_linked: set[str]) -> dict:
    out = {}
    for reading in READINGS:
        for scope, sel in (("all_loci", L.index), ("new_loci", L.index[L.new])):
            nm = set().union(*N[reading][sel])
            rows = {}
            for g in GROUPS:
                k, n = len(Gr[g] & nm), len(Gr[g])
                rows[g] = {"genes": n, "closed": k, "share": clingen.wilson(k, n),
                           "with_a_link_to_another_alzheimer_trait": len(Gr[g] & other_linked)}
            t, z = rows["tied_without_a_link"], rows["nothing"]
            out[f"{reading}:{scope}"] = {
                "groups": rows, "linked_named": len(Gr["linked"] & nm), "linked": len(Gr["linked"]),
                "tied_minus_nothing": clingen.newcombe(t["closed"], t["genes"], z["closed"], z["genes"]),
                "one_source_minus_nothing": clingen.newcombe(rows["one_source"]["closed"], rows["one_source"]["genes"], z["closed"], z["genes"]),
                "tied_minus_one_source": clingen.newcombe(t["closed"], t["genes"], rows["one_source"]["closed"], rows["one_source"]["genes"]),
                "detectable_at_80_percent_power": detectable(t["genes"], z["genes"], z["closed"] / z["genes"] if z["genes"] else 0.0),
                "tied_genes_closed": sorted(Gr["tied_without_a_link"] & nm)}
    return out


def test_order(Lg: list[dict], N: dict[str, pd.Series], rng: np.random.Generator) -> dict:
    out = {}
    for reading in READINGS:
        nm = set().union(*N[reading])
        rows = []
        for l in Lg:
            E = set(l["genes"]) & nm
            if not E:
                continue
            share = lambda gs: len(set(gs) & E) / len(gs)
            rows.append({"genes": l["genes"][:8], "n": len(l["genes"]), "named": sorted(E), "expected": len(E) / len(l["genes"]),
                         "graph": share(l["first"]) if l["first"] else None, "first": l["first"], "baseline": share(l["most_studies"]),
                         "singled_out": l["singled_out"], "singled_out_named": None if l["singled_out"] is None else l["singled_out"] in E})
        w = [r_ for r_ in rows if r_["graph"] is not None]
        res = {"loci": len(Lg), "assessed": len(rows), "with_first_genes": len(w)}
        if w:
            g_, e_, b_ = (np.array([r_[k] for r_ in w]) for k in ("graph", "expected", "baseline"))
            ix = rng.integers(0, len(w), (DRAWS, len(w)))
            res |= {"graph": round(float(g_.sum()), 3), "chance": round(float(e_.sum()), 3), "baseline": round(float(b_.sum()), 3),
                    "graph_minus_chance": [round(float(g_.sum() - e_.sum()), 3), *interval((g_[ix] - e_[ix]).sum(1))],
                    "graph_minus_baseline": [round(float(g_.sum() - b_.sum()), 3), *interval((g_[ix] - b_[ix]).sum(1))]}
        s = [r_ for r_ in rows if r_["singled_out"] is not None]
        res |= {"singled_out_assessed": len(s), "singled_out_named": sum(bool(r_["singled_out_named"]) for r_ in s), "rows": rows}
        out[reading] = res
    return out


def test_floor(G, L: pd.DataFrame, y: int) -> dict:
    floor = dg.study_floor(G, TRAIT, y)
    if not floor or "largest" not in floor or floor["largest"].get("design") != "cases and controls":
        return {"floor": floor, "assessed": False}
    c, k = floor["largest"]["cases"], floor["largest"]["controls"]
    reach = lambda c_, f: math.exp(math.sqrt(studysize.K * (c_ + k) / (2 * f * (1 - f) * c_ * k)))
    could = lambda f: reach(c, f)
    new = L[L.new]
    rows = [{"locus": int(r.locus), "genes": r.cell, "effect_allele_frequency": round(float(r.eaf), 4), "odds_ratio_above_1": round(math.exp(abs(r.beta)), 4),
             "could_have_found": round(could(r.eaf), 4)} for r in new.itertuples()]
    for r_ in rows:
        r_["should_have_been_seen"] = r_["odds_ratio_above_1"] >= r_["could_have_found"]
    need = next(e for e in floor["the_study_that_would_close_it"]["common_variant"] if e["odds_ratio"] == 1.05)
    out = {"assessed": True, "largest": floor["largest"], "check_at_20_percent": round(could(0.2), 3) == floor["largest"]["smallest_odds_ratio_it_could_find"],
           "new_loci": len(rows), "should_have_been_seen": sum(r_["should_have_been_seen"] for r_ in rows),
           "seen": [r_ for r_ in rows if r_["should_have_been_seen"]],
           "odds_ratio_above_1_range": [min(r_["odds_ratio_above_1"] for r_ in rows), max(r_["odds_ratio_above_1"] for r_ in rows)],
           "needed_at_1_05": {"cases": need["cases"], "controls": need["controls"], "effective": round(studysize.effective(need["cases"], need["controls"], 0))},
           "this_study_effective": STUDY_NEFF, "rows": rows}
    # The amendment of 28 Oct: the largest study's individuals, its proxy cases, counted as cases. A check, not the test.
    people = int(G.sizes[floor["largest"]["study"]][2])
    if people:
        c2 = c + people
        seen2 = [{"locus": r_["locus"], "genes": r_["genes"], "effect_allele_frequency": r_["effect_allele_frequency"], "odds_ratio_above_1": r_["odds_ratio_above_1"],
                  "could_have_found": round(reach(c2, r_["effect_allele_frequency"]), 4)} for r_ in rows]
        seen2 = [r_ for r_ in seen2 if r_["odds_ratio_above_1"] >= r_["could_have_found"]]
        out["check_with_individuals_as_cases"] = {"individuals": people, "cases": c2, "at_20_percent": round(reach(c2, 0.2), 3),
                                                  "should_have_been_seen": len(seen2), "seen": seen2}
    return out


def test_ranking(G, r: dict, S: dict, Gr: dict, N: dict[str, pd.Series], L: pd.DataFrame, rng: np.random.Generator) -> dict:
    sc, _, parts = G._score_trait(r["i"], S["mask"])
    deg, nb = parts["degree"], parts["neighbours"]
    have = np.zeros(len(G.t), bool)
    have[G.lt[S["mask"] & (G.ld == r["i"])]] = True
    cand = np.flatnonzero((deg > 0) & ~have)
    C = pd.DataFrame({"sym": G.t["name"].to_numpy()[cand], "served": sc[cand], "neighbours": nb[cand], "degree": deg[cand]})
    C = C[~C.sym.isin(Gr["linked"])].groupby("sym").max()
    linked_any = set(G.t["name"].to_numpy()[np.unique(G.lt[S["mask"]])])
    out = {"candidates": int(len(C))}
    for reading in READINGS:
        for scope, sel in (("new_loci", L.index[L.new]), ("all_loci", L.index)):
            per = [sorted(N[reading][i] & set(C.index)) for i in sel]
            P = sorted(set().union(*map(set, per)))
            pos = C.index.isin(P)
            nm = set().union(*N[reading][sel])
            res = {"positives": len(P), "named_with_no_link_to_any_trait": len(nm - Gr["linked"] - linked_any), "named_already_linked": len(nm & Gr["linked"])}
            if 0 < len(P) < len(C):
                u = {k: below(C[k].to_numpy(), pos) for k in ("served", "neighbours", "degree")}
                where = {p: i for i, p in enumerate(C.index[pos])}
                W = np.zeros((DRAWS, len(P)))
                for d, draw in enumerate(rng.integers(0, len(per), (DRAWS, len(per)))):
                    for j in draw:
                        for p in per[j]:
                            W[d, where[p]] += 1
                tot = W.sum(1)
                ok = tot > 0
                auc = lambda k: np.where(ok, (W * u[k]).sum(1) / np.where(ok, tot, 1), np.nan)
                rngs = [np.random.default_rng(s) for s in score.SEEDS]
                ap = lambda k: round(float(np.mean([score.ranked(C[k].to_numpy(), pos, g)[0] for g in rngs])), 4)
                res |= {"roc_area": {k: round(float(u[k].mean()), 4) for k in u}, "average_precision": {k: ap(k) for k in u},
                        "served_minus_neighbours": [round(float(u["served"].mean() - u["neighbours"].mean()), 4), *interval(auc("served") - auc("neighbours"))],
                        "served_minus_degree": [round(float(u["served"].mean() - u["degree"].mean()), 4), *interval(auc("served") - auc("degree"))],
                        "draws_with_a_positive": int(ok.sum()),
                        "positives_by_served_percentile": {p: round(float(x), 4) for p, x in sorted(zip(C.index[pos], u["served"]), key=lambda z: -z[1])}}
            out[f"{reading}:{scope}"] = res
    return out


# ---------------------------------------------------------------- the audit

def audit(release: str) -> dict:
    from pharos_gap_kg import ask
    blob = registry.require(KEY).read_bytes()
    if sha(blob) != registry.SOURCES[KEY].sha256:
        raise SystemExit(f"REFUSED: {registry.SOURCES[KEY].filename} is not the file the registry pins")
    L = study_loci(sheet(blob, "SD2"))
    G = ask.Graph(paths.bundle(release))
    r = G.resolve(TRAIT_NAME, "trait")
    if not (r["ok"] and r["kind"] == "trait" and G.d.id[r["i"]] == TRAIT):
        raise SystemExit(f"REFUSED: {TRAIT_NAME!r} does not resolve to the trait {TRAIT}")
    y = G._year(None)
    m = symbol_map(G, release)
    N = named(L, m)
    S = graph_sets(G, r, y)
    Lg = graph_loci(G, S)
    check = crosscheck(G, r, y, S, Lg)
    Gr = groups(G, S)
    other = set()
    for tid in OTHER_TRAITS:
        if tid in G.di:
            other |= {G.t["name"][G.lt[k]] for k in np.flatnonzero(S["mask"] & (G.ld == G.di[tid]))}
    rng = np.random.default_rng(SEED)
    study_syms = sorted({g for gs in L.effector for g in gs} | {g for gs in L.evidence for g in gs} | {g for gs in L.closest for g in gs})
    eff_syms = sorted({g for gs in L.effector for g in gs})
    arch = json.loads((paths.DOCS_DATA / f"bundle_archive_{release}.json").read_text())
    return {
        "release": release, "protocol": "docs/PROTOCOL.md, A study published after the graph: Alzheimer disease (28 Oct 2026)",
        "bundle_archive_sha256": arch["sha256"],
        "study": {"key": KEY, "file": registry.SOURCES[KEY].filename, "sha256": registry.SOURCES[KEY].sha256, "doi": "10.1038/s41586-026-11206-y",
                  "loci": int(len(L)), "new_loci": int(L.new.sum()), "clusters": sorted(L.cell[L.cell.str.contains("cluster", case=False)]),
                  "effective_size": STUDY_NEFF},
        "trait": {"id": TRAIT, "name": G.d["name"][r["i"]], "as_of": y, "other_alzheimer_traits": list(OTHER_TRAITS)},
        "matching": {"effector_symbols": len(eff_syms), "effector_symbols_matched": sum(g in m for g in eff_syms),
                     "effector_symbols_not_in_the_graph": [g for g in eff_syms if g not in m],
                     "matched_through_ensembl_symbol": sorted(g for g in eff_syms if g in m and g not in set(G.t["name"])),
                     "all_symbols": len(study_syms), "all_symbols_matched": sum(g in m for g in study_syms),
                     "loci_naming_a_graph_gene": {k: int((N[k].map(len) > 0).sum()) for k in READINGS}},
        "crosscheck": check,
        "groups": {g: len(Gr[g]) for g in (*GROUPS, "linked")},
        "seed": SEED, "draws": DRAWS,
        "gaps_that_close": test_gaps(Gr, N, L, other),
        "order_of_a_locus": test_order(Lg, N, rng),
        "largest_study": test_floor(G, L, y),
        "ranking": test_ranking(G, r, S, Gr, N, L, rng),
    }


def _selftest() -> int:
    fails = []

    def ck(name: str, cond: bool) -> None:
        print(("  ok    " if cond else "  FAIL  ") + name)
        if not cond:
            fails.append(name)

    cell = lambda v, n=1, num=False: (f'<table:table-cell table:number-columns-repeated="{n}" ' + (f'office:value-type="float" office:value="{v}"/>' if num
                                      else f'office:value-type="string"><text:p>{v}</text:p></table:table-cell>'))
    xml = (f'<office:document-content xmlns:office="{_O}" xmlns:table="{_T}" xmlns:text="{_X}"><office:body><office:spreadsheet>'
           '<table:table table:name="SD2"><table:table-row>' + cell("a") + cell("b", 2) + '</table:table-row>'
           '<table:table-row table:number-rows-repeated="1000000"><table:table-cell table:number-columns-repeated="1024"/></table:table-row>'
           '<table:table-row>' + cell("0.54590000000000005", num=True) + '<table:table-cell/>' + cell("c") + '</table:table-row>'
           '</table:table></office:spreadsheet></office:body></office:document-content>')
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("content.xml", xml)
    rows = sheet(buf.getvalue(), "SD2")
    ck("a sheet reads repeated cells, skips empty padding rows and keeps a number's full value", rows == [["a", "b", "b"], ["0.54590000000000005", "", "c"]])
    ck("an absent sheet is refused", str(_raises(lambda: sheet(buf.getvalue(), "SD9"))).startswith("REFUSED"))
    ck("genes split on '/', drop NA and clusters, keep order", genes(" A/B /NA/MS4A gene cluster/A") == ["A", "B"] and genes("NA") == [])
    s = np.array([5.0, 3.0, 3.0, 1.0, 0.0])
    p = np.array([True, False, True, False, False])
    ck("below: the share of negatives under each positive, ties half", np.allclose(below(s, p), [1.0, 2.5 / 3]))
    ck("below agrees with score.ranked's ROC area without ties",
       abs(below(np.array([4., 3, 2, 1]), np.array([True, False, True, False])).mean() - score.ranked(np.array([4., 3, 2, 1]), np.array([1, 0, 1, 0]), np.random.default_rng(0))[1]) < 1e-12)
    d1, d2 = detectable(100, 10_000, 0.01), detectable(400, 10_000, 0.01)
    ck("a larger group detects a smaller rise", d1 is not None and d2 is not None and 0 < d2 < d1 < 1)
    ck("no genes, no power figure", detectable(0, 100, 0.1) is None)
    c, k = 21_982, 41_944
    could = math.exp(math.sqrt(studysize.K * (c + k) / (2 * 0.2 * 0.8 * c * k)))
    ck("at a frequency of 20% the floor is studysize's", round(could, 3) == studysize.smallest_effect(c, k, 0)["smallest_odds_ratio_it_could_find"])
    print("selftest:", "PASS" if not fails else f"FAIL {fails}")
    return 0 if not fails else 1


def _raises(fn):
    try:
        fn()
    except SystemExit as e:
        return e
    return None


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if argv[:1] == ["fetch"]:
        fetch(Path(argv[argv.index("--from") + 1]) if "--from" in argv else None)
        return 0
    if len(argv) == 2 and argv[0] == "audit":
        out = audit(argv[1])
        dest = paths.DOCS_DATA / f"adstudy_{argv[1]}.json"
        tmp = dest.with_name(dest.name + ".tmp")
        tmp.write_text(json.dumps(out, indent=1) + "\n")
        os.replace(tmp, dest)
        print(f"  {dest}")
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
