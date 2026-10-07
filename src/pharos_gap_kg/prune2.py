"""Second pruning round: a literature check scored on matched planted pairs, and the sources' own fields.

Round one (findings/PRUNING.md) showed that nothing inside the graph tells a real pair from a planted
one at an acceptable cost. This round asks two things from outside it (docs/PROTOCOL.md):

  - does a PubMed check, alone or with a model reading the abstracts, withhold planted pairs and keep
    real ones;
  - does a source's own evidence field order its real pairs the way that check does.

A planted pair has no evidence value, so the fields are scored on real pairs only.

Stages, each writing what the next reads (data/eval/pruning2/<release>):

    python -m pharos_gap_kg.prune2 freeze pharos319     # plant, draw the sample, write its digest
    python -m pharos_gap_kg.prune2 ingraph pharos319    # round one's rule on the matched planting
    python -m pharos_gap_kg.prune2 search pharos319     # the code-only check (PubMed)
    python -m pharos_gap_kg.prune2 read pharos319 dev   # three models read the development part
    python -m pharos_gap_kg.prune2 choose pharos319     # the reader, by the stated rule
    python -m pharos_gap_kg.prune2 read pharos319 test  # the chosen model reads the test part
    python -m pharos_gap_kg.prune2 report pharos319
    python -m pharos_gap_kg.prune2 --selftest
"""
from __future__ import annotations

import json
import math
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import ctxtables, evidence, graph, litcheck, ontology, paths, planted as round1, registry

SEED = 20261009
FRACTION, AT_LEAST, SAMPLE, DEV, TEST = 0.05, 200, 200, 20, 80
MATCHED = (0.40, 0.60)
QUALIFY, WITHIN = 0.95, 0.02
SPEND_CAP = 2.00
N_BOOT = 2000
GRADES_AT, MIN_SIDE = 0.60, 20
PROVISIONAL_BELOW = 0.70
KEY = ["target", "disease", "source"]


def work(release: str) -> Path:
    return paths.EVAL / "pruning2" / release


def dump(obj) -> str:
    """JSON with a number that is not a number written as null."""
    def fix(x):
        if isinstance(x, dict):
            return {k: fix(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [fix(v) for v in x]
        return None if isinstance(x, float) and math.isnan(x) else x
    return json.dumps(fix(obj), indent=1, allow_nan=False) + "\n"


# ---- planting and the sample -------------------------------------------------------------------------------

def one_source(links: pd.DataFrame) -> pd.DataFrame:
    n = links.groupby(["target", "disease"]).source.transform("nunique")
    return links[n == 1].drop_duplicates(["target", "disease"]).sort_values(KEY).reset_index(drop=True)


def plant_matched(links: pd.DataFrame, fraction: float = FRACTION, at_least: int = AT_LEAST, seed: int = SEED) -> pd.DataFrame:
    """target, disease, source, order, in_graph_set: false pairs made by crossing two of a source's
    ONE-SOURCE pairs. `order` is the order drawn; the first 5% of a source's are its in-graph set."""
    rng = np.random.default_rng(seed)
    taken = set(zip(links.target, links.disease))
    out = []
    for s, g in one_source(links).groupby("source"):
        t, d = g.target.to_numpy(), g.disease.to_numpy()
        floor = math.ceil(fraction * len(g))
        want, got, tries = max(floor, at_least), 0, 0
        while got < want and tries < 200 * want + 1000:
            tries += 1
            i, j = rng.integers(0, len(g), 2)
            pair = (t[i], d[j])
            if pair not in taken:
                taken.add(pair)
                out.append((t[i], d[j], s, got, got < floor))
                got += 1
        if got < floor:
            raise SystemExit(f"REFUSED: could plant {got} of the {floor} false pairs {s} needs; its one-source pairs cross into real ones too often")
    return pd.DataFrame(out, columns=KEY + ["order", "in_graph_set"])


def draw_real(single: pd.DataFrame, n: int = SAMPLE, seed: int = SEED) -> pd.DataFrame:
    """target, disease, source, order: `n` of each source's one-source pairs, in the order drawn."""
    rng = np.random.default_rng(seed + 1)
    out = []
    for s, g in single.groupby("source"):
        take = rng.permutation(len(g))[:n]
        out.append(g.iloc[take][KEY].assign(order=np.arange(len(take))))
    return pd.concat(out, ignore_index=True)


def part(order: int) -> str:
    return "dev" if order < DEV else "test" if order < DEV + TEST else "rest"


def names(release: str) -> tuple[dict, dict]:
    """(source, disease id) -> the release's normalised name most often given; target -> gene symbol."""
    d = paths.tables(release)
    r = pq.read_table(d / "disease.parquet", columns=["dtype", "did", "mondoid", "ncats_name"]).to_pandas()
    r["disease"] = [ctxtables.disease_id(m, x) for m, x in zip(r.mondoid, r.did)]
    r = r.dropna(subset=["disease", "ncats_name"])
    top = r.groupby(["dtype", "disease", "ncats_name"]).size().reset_index(name="n").sort_values(["n", "ncats_name"], ascending=[False, True])
    top = top.drop_duplicates(["dtype", "disease"])
    p = pq.read_table(d / "protein.parquet", columns=["uniprot", "sym"]).to_pandas().dropna(subset=["sym"])
    return dict(zip(zip(top.dtype, top.disease), top.ncats_name)), dict(zip("UniProtKB:" + p.uniprot, p.sym))


def make_sample(links: pd.DataFrame, pl: pd.DataFrame, name_of: dict, symbol_of: dict) -> pd.DataFrame:
    a = pl[pl.order < SAMPLE][KEY + ["order"]].assign(planted=True)
    b = draw_real(one_source(links)).assign(planted=False)
    s = pd.concat([a, b], ignore_index=True)
    s["part"] = s.order.map(part)
    s["symbol"] = s.target.map(symbol_of)
    s["name"] = [name_of.get((src, d)) for src, d in zip(s.source, s.disease)]
    return s.sort_values(["source", "planted", "order"]).reset_index(drop=True)


def freeze(release: str) -> dict:
    C = paths.context(release)
    links = pd.read_parquet(C / "disease_links.parquet")
    pl = plant_matched(links)
    name_of, symbol_of = names(release)
    s = make_sample(links, pl, name_of, symbol_of)
    W = work(release)
    W.mkdir(parents=True, exist_ok=True)
    pl.to_parquet(W / "planted.parquet", index=False)
    s.to_parquet(W / "sample.parquet", index=False)
    single = one_source(links)
    per = {}
    for src in sorted(set(links.source)):
        x = s[s.source == src]
        per[src] = {"one_source_pairs": int((single.source == src).sum()), "planted": int((pl.source == src).sum()),
                    "in_graph_set": int((pl[pl.source == src].in_graph_set).sum()),
                    "sample_planted": int(x.planted.sum()), "sample_real": int((~x.planted).sum()),
                    "sample_without_a_symbol": int(x.symbol.isna().sum()), "sample_without_a_name": int(x.name.isna().sum())}
    out = {"release": release, "seed": SEED, "fraction": FRACTION, "at_least": AT_LEAST, "sample_per_source": SAMPLE, "dev": DEV, "test": TEST,
           "disease_links_sha256": json.loads((C / "manifest.json").read_text())["tables"]["disease_links"]["sha256"],
           "files": {f: graph.file_sha(W / f) for f in ("planted.parquet", "sample.parquet")}, "sources": per}
    (paths.DOCS_DATA / f"pruning2_sample_{release}.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def frozen(release: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The planted pairs and the sample, refused unless they are the ones whose digests were committed."""
    f = paths.DOCS_DATA / f"pruning2_sample_{release}.json"
    if not f.is_file():
        raise SystemExit(f"no frozen sample; run:  python -m pharos_gap_kg.prune2 freeze {release}")
    want = json.loads(f.read_text())["files"]
    W = work(release)
    for name, sha in want.items():
        if not (W / name).is_file() or graph.file_sha(W / name) != sha:
            raise SystemExit(f"REFUSED: {W / name} is not the file whose digest is in {f.name}; the sample is frozen and is not redrawn")
    return pd.read_parquet(W / "planted.parquet"), pd.read_parquet(W / "sample.parquet")


# ---- small statistics --------------------------------------------------------------------------------------

def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if not n:
        return float("nan"), float("nan")
    p, d = k / n, 1 + z * z / n
    c, h = p + z * z / (2 * n), z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def share(flags) -> dict:
    flags = np.asarray(flags, bool)
    lo, hi = wilson(int(flags.sum()), len(flags))
    return {"n": int(len(flags)), "share": float(flags.mean()) if len(flags) else float("nan"), "low": lo, "high": hi}


def wauc(codes: np.ndarray, label: np.ndarray, w: np.ndarray | None = None) -> float:
    """The chance a labelled row's value is above an unlabelled row's, ties a half, rows weighted.
    `codes` are the values as ranks (np.unique's inverse), so a resample only changes `w`."""
    w = np.ones(len(codes)) if w is None else w
    P = np.bincount(codes, weights=w * label)
    N = np.bincount(codes, weights=w * ~label)
    if not P.sum() or not N.sum():
        return float("nan")
    return float((P * (np.cumsum(N) - N + 0.5 * N)).sum() / (P.sum() * N.sum()))


def auc_interval(value: np.ndarray, label: np.ndarray, cluster: np.ndarray | None = None, seed: int = SEED, n_boot: int = N_BOOT) -> dict:
    """The area with a 95% interval from resampling rows, or whole clusters where `cluster` is given."""
    codes = np.unique(value, return_inverse=True)[1]
    label = np.asarray(label, bool)
    cl = np.arange(len(value)) if cluster is None else np.unique(cluster, return_inverse=True)[1]
    k = int(cl.max()) + 1 if len(cl) else 0
    rng = np.random.default_rng(seed)
    b = np.array([wauc(codes, label, rng.multinomial(k, np.full(k, 1 / k))[cl].astype(float)) for _ in range(n_boot)])
    b = b[~np.isnan(b)]
    lo, hi = (np.quantile(b, [0.025, 0.975]) if len(b) else (float("nan"), float("nan")))
    return {"auc": wauc(codes, label), "low": float(lo), "high": float(hi), "labelled": int(label.sum()), "unlabelled": int((~label).sum())}


def matched_cut(real: np.ndarray, pl: np.ndarray, real_share: float) -> tuple[float, float]:
    """A cut on a score set to withhold at least `real_share` of real pairs: the shares of planted and
    of real pairs at or below it."""
    k = math.ceil(real_share * len(real) - 1e-9)
    if k <= 0 or not len(real):
        return 0.0, 0.0
    tau = np.sort(real)[k - 1]
    return float((pl <= tau).mean()), float((real <= tau).mean())


# ---- the code-only check -----------------------------------------------------------------------------------

def classify(symbol, name, pair: dict | None, gene: dict | None, disease: dict | None) -> tuple[int | None, str | None]:
    """(papers that mention both, None) for a searchable pair, or (None, why it is unsearchable)."""
    if not isinstance(symbol, str) or not symbol:
        return None, "the target has no gene symbol"
    if not isinstance(name, str) or not litcheck.clean_name(name):
        return None, "the disease has no searchable name"
    if disease["not_found"]:
        return None, "PubMed did not find the disease name"
    if gene["count"] == 0:
        return 0, None
    if pair["count"] > min(gene["count"], disease["count"]):
        return None, "the counts do not agree"
    return pair["count"], None


def search(release: str) -> pd.DataFrame:
    _, s = frozen(release)
    ok = s.symbol.notna() & s.name.map(lambda n: isinstance(n, str) and bool(litcheck.clean_name(n)))
    q = s[ok]
    terms = [litcheck.EVERYTHING] + [litcheck.pair_term(a, b) for a, b in zip(q.symbol, q.name)] + [litcheck.gene_term(a) for a in q.symbol] + [litcheck.disease_term(b) for b in q.name]
    R = litcheck.search_many(terms)
    bad = [t for t, r in R.items() if r["error"]]
    if bad:
        raise SystemExit(f"REFUSED: {len(bad)} searches came back as errors, e.g. {bad[0]!r}: {R[bad[0]]['error']}; re-run to ask them again")
    total = R[litcheck.EVERYTHING]["count"]
    rows = []
    for sym, name, good in zip(s.symbol, s.name, ok):
        pr = R[litcheck.pair_term(sym, name)] if good else None
        g = R[litcheck.gene_term(sym)] if good else None
        d = R[litcheck.disease_term(name)] if good else None
        both, why = classify(sym, name, pr, g, d)
        rows.append({"both": both, "gene_papers": g["count"] if good else None, "disease_papers": d["count"] if good else None,
                     "lift": litcheck.lift(both, g["count"], d["count"], total) if both is not None else None,
                     "ids": "|".join(pr["ids"]) if both else "", "unsearchable": why})
    out = pd.concat([s[KEY + ["planted", "order", "part"]].reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    out["both"] = out.both.astype("Int64")
    out.attrs = {}
    W = work(release)
    out.to_parquet(W / "search.parquet", index=False)
    (W / "search.json").write_text(json.dumps({"pubmed_records": total, "searches": len(R), "sample_sha256": graph.file_sha(W / "sample.parquet")}, indent=1) + "\n")
    return out


# ---- the reader --------------------------------------------------------------------------------------------

def read(release: str, which: str, models: list[str], key: str | None, base_url: str, workers: int = 6) -> pd.DataFrame:
    from pharos_gap_kg.usability import PRICES
    _, s = frozen(release)
    W = work(release)
    se = pd.read_parquet(W / "search.parquet")
    q = s.merge(se, on=KEY + ["planted", "order", "part"])
    q = q[(q.part == which) & q.unsearchable.isna() & (q.both > 0)].reset_index(drop=True)
    recs = []
    for k, ids in enumerate(q.ids, 1):
        recs.append(litcheck.records(ids.split("|")))
        if k % 200 == 0:
            print(f"  records for {k:,} / {len(q):,} pairs", flush=True)
    lock = threading.Lock()
    total = [litcheck.spent()]
    rows = []
    for model in models:
        def one(i: int) -> dict:
            base = {**{c: q.at[i, c] for c in KEY + ["planted", "order", "part"]}, "model": model}
            if not recs[i]:
                return {**base, "verdict": None, "why_refused": "no record has a title or an abstract", "attempts": 0, "cost": 0.0, "error": None}
            with lock:
                if total[0] >= SPEND_CAP:
                    return {**base, "verdict": None, "why_refused": None, "attempts": 0, "cost": 0.0, "error": "spend cap reached"}
            r = litcheck.read(q.at[i, "symbol"], q.at[i, "name"], recs[i], model, key, base_url, PRICES.get(model))
            with lock:
                total[0] += 0.0 if r.get("cached") else r.get("new_cost", r["cost"])
            return {**base, **{c: r.get(c) for c in ("verdict", "pmid", "quote", "why_refused", "attempts", "cost", "error")}}

        before = litcheck.spent()
        with ThreadPoolExecutor(workers) as ex:
            got = list(ex.map(one, range(len(q))))
        total[0] = litcheck.spent()
        rows += got
        n_err = sum(1 for r in got if r["error"])
        print(f"  {model}: {len(got)} pairs, {sum(1 for r in got if r['verdict'])} with a verdict, {n_err} errors, spent {total[0] - before:.3f} (round total {total[0]:.3f})", flush=True)
    out = pd.DataFrame(rows)
    f = W / f"reader_{which}.parquet"
    if f.is_file():
        old = pd.read_parquet(f)
        out = pd.concat([old[~old.model.isin(models)], out], ignore_index=True)
    out.to_parquet(f, index=False)
    return out


def withheld_l1(sample_search: pd.DataFrame, reader: pd.DataFrame, model: str) -> pd.DataFrame:
    """The searchable pairs of the reader's part with `l0` and `l1`: withheld by no paper, or by that or the reader's flag."""
    r = reader[reader.model == model]
    if r.error.notna().any():
        raise SystemExit(f"REFUSED: {int(r.error.notna().sum())} pairs were not read by {model} ({r.error.dropna().iloc[0]}); read them before scoring")
    x = sample_search[sample_search.unsearchable.isna()].merge(r[KEY + ["planted", "verdict", "pmid", "quote"]], on=KEY + ["planted"], how="left")
    unread = (x.both > 0) & ~x.set_index(KEY + ["planted"]).index.isin(r.set_index(KEY + ["planted"]).index)
    if unread.any():
        raise SystemExit(f"REFUSED: {int(unread.sum())} mentioned pairs have no reading by {model}")
    x["l0"] = (x.both == 0).to_numpy()
    x["l1"] = x.l0 | (x.verdict == "mentions_only").to_numpy()
    x["no_verdict"] = ((x.both > 0) & x.verdict.isna()).to_numpy()
    return x


def choose(release: str) -> dict:
    from pharos_gap_kg.usability import PRICES
    W = work(release)
    se = pd.read_parquet(W / "search.parquet")
    rd = pd.read_parquet(W / "reader_dev.parquet")
    dev = se[se.part == "dev"]
    res = {}
    for model in sorted(set(rd.model)):
        x = withheld_l1(dev, rd, model)
        asked = x[x.both > 0]
        p, r = x[x.planted], x[~x.planted]
        res[model] = {"read": int(len(asked)), "passed_the_check": float(asked.verdict.notna().mean()),
                      "planted_withheld": float(p.l1.mean()), "real_withheld": float(r.l1.mean()),
                      "difference": float(p.l1.mean() - r.l1.mean()), "price_in_out": list(PRICES.get(model, (None, None))),
                      "cost": float(rd[rd.model == model].cost.sum())}
        res[model]["qualifies"] = bool(res[model]["passed_the_check"] >= QUALIFY)
    ok = {m: v for m, v in res.items() if v["qualifies"]}
    chosen = None
    if ok:
        best = max(v["difference"] for v in ok.values())
        near = [m for m, v in ok.items() if best - v["difference"] <= WITHIN]
        chosen = min(near, key=lambda m: (sum(PRICES.get(m, (9, 9))), m))
    out = {"release": release, "brief_sha": litcheck.BRIEF_SHA, "dev_pairs": int(len(dev)), "dev_unsearchable": int(dev.unsearchable.notna().sum()),
           "models": res, "chosen": chosen, "rule": f"passes the check on at least {QUALIFY:.0%} of answers; then the largest planted minus real withheld; within {WITHIN}, the cheaper"}
    (paths.DOCS_DATA / f"pruning2_reader_{release}.json").write_text(dump(out))
    return out


# ---- the report --------------------------------------------------------------------------------------------

def bar(planted: dict, real: dict) -> dict:
    ok = planted["share"] >= round1.KEEP_AT_LEAST and real["share"] <= round1.LOSE_AT_MOST
    return {"passes": bool(ok), "provisional": bool(ok and planted["low"] < PROVISIONAL_BELOW)}


def rule_rows(x: pd.DataFrame, col: str) -> dict:
    out = {}
    for s, g in list(x.groupby("source")) + [("all sources", x)]:
        p, r = share(g[g.planted][col]), share(g[~g.planted][col])
        out[s] = {"planted_withheld": p, "real_withheld": r, **bar(p, r)}
    return out


def field_rows(release: str, se: pd.DataFrame, lineage: pd.DataFrame) -> dict:
    ev = pd.read_parquet(paths.context(release) / "disease_evidence.parquet")
    real = se[~se.planted & se.unsearchable.isna()]
    out = {}
    for (src, field), g in ev.groupby(["source", "field"]):
        x = real[real.source == src].merge(g[["target", "disease", "value"]], on=["target", "disease"])
        label = (x.both > 0).to_numpy()
        v = x.value.to_numpy()
        row = {"independent_of_the_label": (src, field) not in evidence.NOT_INDEPENDENT, "sampled_real_pairs": int(len(x)),
               "mentioned": int(label.sum()), "not_mentioned": int((~label).sum())}
        row["measurable"] = bool(min(label.sum(), (~label).sum()) >= MIN_SIDE)
        row["auc_for_mentioned"] = auc_interval(v, label)
        a = row["auc_for_mentioned"]
        row["grades_its_source"] = bool(row["measurable"] and a["auc"] >= GRADES_AT and a["low"] > 0.5)
        med = float(np.median(v)) if len(v) else float("nan")
        row["median"] = med
        row["not_mentioned_by_band"] = {name: share(~label[m]) for name, m in (("below_median", v < med), ("at_median", v == med), ("above_median", v > med)) if m.any()}
        pl = se[se.planted & se.unsearchable.isna() & (se.source == src)]
        row["planted_not_mentioned"] = share((pl.both == 0).to_numpy())
        full = lineage[lineage.source == src].merge(g[["target", "disease", "value"]], on=["target", "disease"])
        row["auc_for_lineage_all_one_source_pairs"] = {**auc_interval(full.value.to_numpy(), full.lineage.to_numpy(), cluster=full.target.to_numpy()),
                                                       "lineage_share": float(full.lineage.mean())}
        out[f"{src} | {field}"] = row
    return out


def split_cut(ok: pd.DataFrame) -> dict:
    """Rule L0b (protocol amendment): for each source a cut on the lift, set to withhold 80% of the planted
    pairs of the part the reader does not use, and scored on the source's other pairs."""
    fit, score = ok[ok.part == "rest"], ok[ok.part != "rest"]
    cut = {s: float(np.quantile(g[g.planted].lift, round1.KEEP_AT_LEAST)) for s, g in fit.groupby("source")}
    x = score.assign(l0b=(score.lift.to_numpy(float) <= score.source.map(cut).to_numpy(float)))
    out = rule_rows(x, "l0b")
    for s in cut:
        out[s]["cut"] = cut[s]
        out[s]["planted_pairs_the_cut_was_set_on"] = int((fit[fit.source == s].planted).sum())
    return out


def reader_against_lift(t: pd.DataFrame, seed: int = SEED) -> dict:
    """On the test part: L1 against a cut on the lift set to withhold at least as many real pairs."""
    p, r = t[t.planted], t[~t.planted]
    pl_lift, re_lift, pl_l1, re_l1 = p.lift.to_numpy(float), r.lift.to_numpy(float), p.l1.to_numpy(), r.l1.to_numpy()

    def diff(pi: np.ndarray, ri: np.ndarray) -> tuple[float, float, float, float]:
        cut_p, cut_r = matched_cut(re_lift[ri], pl_lift[pi], float(re_l1[ri].mean()))
        return float(pl_l1[pi].mean() - cut_p), cut_p, cut_r, float(pl_l1[pi].mean())

    d, cut_p, cut_r, l1_p = diff(np.arange(len(p)), np.arange(len(r)))
    rng = np.random.default_rng(seed)
    b = np.array([diff(rng.integers(0, len(p), len(p)), rng.integers(0, len(r), len(r)))[0] for _ in range(N_BOOT)])
    lo, hi = np.quantile(b, [0.025, 0.975])
    return {"l1_planted_withheld": l1_p, "l1_real_withheld": float(re_l1.mean()), "lift_cut_planted_withheld": cut_p, "lift_cut_real_withheld": cut_r,
            "difference": d, "low": float(lo), "high": float(hi), "reader_adds": bool(lo > 0)}


def report(release: str) -> dict:
    C, W = paths.context(release), work(release)
    pl, s = frozen(release)
    se = pd.read_parquet(W / "search.parquet")
    links = pd.read_parquet(C / "disease_links.parquet")
    chosen = json.loads((paths.DOCS_DATA / f"pruning2_reader_{release}.json").read_text())["chosen"]
    ing = json.loads((W / "ingraph.json").read_text())
    ok = se[se.unsearchable.isna()].copy()
    ok["l0"] = (ok.both == 0).to_numpy()
    out = {"release": release, "seed": SEED, "sample_sha256": graph.file_sha(W / "sample.parquet"), "pubmed": json.loads((W / "search.json").read_text()),
           "bar": {"planted_withheld_at_least": round1.KEEP_AT_LEAST, "real_withheld_at_most": round1.LOSE_AT_MOST, "provisional_if_planted_interval_reaches_below": PROVISIONAL_BELOW},
           "unsearchable": {str(k): int(v) for k, v in se.unsearchable.value_counts().items()},
           "unsearchable_by_source": {src: {"planted": int(g[g.planted].unsearchable.notna().sum()), "real": int(g[~g.planted].unsearchable.notna().sum())} for src, g in se.groupby("source")}}

    l0 = rule_rows(ok, "l0")
    for src, g in list(ok.groupby("source")) + [("all sources", ok)]:
        p, r = g[g.planted], g[~g.planted]
        l0[src].update({"auc_papers": round1.auc(r.both.to_numpy(float), p.both.to_numpy(float)), "auc_lift": round1.auc(r.lift.to_numpy(float), p.lift.to_numpy(float)),
                        "real_withheld_by_a_lift_cut_that_withholds_80pct_of_planted": round1.cost_at(r.lift.to_numpy(float), p.lift.to_numpy(float), round1.KEEP_AT_LEAST),
                        "median_papers": {"planted": float(p.both.median()), "real": float(r.both.median())}})
    out["code_only_L0"] = l0
    out["lift_cut_L0b"] = split_cut(ok)

    out["reader"] = {"chosen": chosen}
    if chosen:
        rd = pd.read_parquet(W / "reader_test.parquet")
        t = withheld_l1(se[se.part == "test"], rd, chosen)
        l1 = rule_rows(t, "l1")
        for src, g in list(t.groupby("source")) + [("all sources", t)]:
            l1[src].update({"no_verdict": int(g.no_verdict.sum()), "read": int((g.both > 0).sum()),
                            "l0_planted_withheld_same_pairs": float(g[g.planted].l0.mean()), "l0_real_withheld_same_pairs": float(g[~g.planted].l0.mean())})
        kept = t[t.planted & (t.verdict == "asserts")].merge(s[KEY + ["planted", "symbol", "name"]], on=KEY + ["planted"])
        kept[KEY + ["symbol", "name", "pmid", "quote"]].to_parquet(W / "planted_kept_by_reader.parquet", index=False)
        out["reader"].update({"test_L1": l1, "against_lift": reader_against_lift(t), "cost_of_test_reading": float(rd[rd.model == chosen].cost.sum()),
                              "planted_kept_with_a_passage": int(len(kept)), "planted_read": int((t.planted & (t.both > 0)).sum()),
                              "planted_kept_examples": [{"symbol": a, "disease": b, "source": c, "pmid": d} for a, b, c, d in
                                                        kept.sort_values(["source", "order"]).head(12)[["symbol", "name", "source", "pmid"]].itertuples(index=False)]})
    out["spent_on_model_calls"] = litcheck.spent()
    out["spend_cap"] = SPEND_CAP

    terms = ontology.read(registry.require("efo_obo"))
    single = one_source(links)
    rel, onto = round1.lineage_of(sorted(set(links.disease)), terms)
    gwas = pd.read_parquet(paths.bundle(release) / "links.parquet", columns=["target", "trait"])
    lineage = round1.corroboration(single, links, rel, onto, gwas, pd.DataFrame({"a": [], "b": []}), terms)[KEY + ["lineage"]]
    out["fields"] = field_rows(release, se, lineage)
    out["in_graph_rule_on_matched_planting"] = ing
    # Two descriptions added after the stated measures were computed.
    out["planted_pairs_asserted_under_another_identifier"] = same_name(links, s, names(release)[0])
    out["in_graph_against_literature"] = agreement(ok[~ok.planted], pd.read_parquet(W / "ingraph_detail.parquet"))
    (paths.DOCS_DATA / f"pruning2_{release}.json").write_text(dump(out))
    return out


def show(out: dict) -> None:
    pc = lambda d: f"{d['share']:6.1%} [{d['low']:5.1%}, {d['high']:5.1%}] n={d['n']:<4d}"
    order = lambda rows: sorted(rows.items(), key=lambda kv: (kv[0] == "all sources", kv[0]))
    print(f"  unsearchable: {out['unsearchable']}")
    print("  L0, no paper mentions the pair")
    for s, r in order(out["code_only_L0"]):
        print(f"    {s:34s} planted {pc(r['planted_withheld'])}  real {pc(r['real_withheld'])}  AUC papers {r['auc_papers']:.3f}  lift {r['auc_lift']:.3f}  "
              f"lift-cut cost {r['real_withheld_by_a_lift_cut_that_withholds_80pct_of_planted']:.3f}  {'PASS' if r['passes'] else 'no'}{' (provisional)' if r['provisional'] else ''}")
    print("  L0b, a lift cut set on the part the reader does not use, scored on the rest")
    for s, r in order(out["lift_cut_L0b"]):
        print(f"    {s:34s} planted {pc(r['planted_withheld'])}  real {pc(r['real_withheld'])}  cut {r.get('cut', float('nan')):8.3f}  {'PASS' if r['passes'] else 'no'}{' (provisional)' if r['provisional'] else ''}")
    rd = out["reader"]
    if rd.get("chosen"):
        print(f"  L1 on the test part, reader {rd['chosen']}")
        for s, r in order(rd["test_L1"]):
            print(f"    {s:34s} planted {pc(r['planted_withheld'])}  real {pc(r['real_withheld'])}  read {r['read']:<4d} no verdict {r['no_verdict']:<3d} "
                  f"{'PASS' if r['passes'] else 'no'}{' (provisional)' if r['provisional'] else ''}")
        a = rd["against_lift"]
        print(f"  reader against a lift cut of the same cost: L1 {a['l1_planted_withheld']:.3f} planted at {a['l1_real_withheld']:.3f} real; lift cut {a['lift_cut_planted_withheld']:.3f} at {a['lift_cut_real_withheld']:.3f}; "
              f"difference {a['difference']:+.3f} [{a['low']:+.3f}, {a['high']:+.3f}]  adds: {a['reader_adds']}")
        print(f"  planted pairs the reader kept with a passage: {rd['planted_kept_with_a_passage']} of {rd['planted_read']} read")
    print("  fields")
    for k, r in out["fields"].items():
        a, f = r["auc_for_mentioned"], r["auc_for_lineage_all_one_source_pairs"]
        bands = "  ".join(f"{n.split('_')[0]} {d['share']:.2f} (n={d['n']})" for n, d in r["not_mentioned_by_band"].items())
        print(f"    {k:44s} AUC mentioned {a['auc'] if a['auc'] is not None else float('nan'):.3f} [{a['low']:.3f}, {a['high']:.3f}] ({r['mentioned']}/{r['not_mentioned']})  "
              f"{'GRADES' if r['grades_its_source'] else 'no' if r['measurable'] else 'not measurable'}{'' if r['independent_of_the_label'] else ' (not independent)'}  "
              f"not mentioned: {bands}; planted {r['planted_not_mentioned']['share']:.2f}  | lineage AUC {f['auc']:.3f} [{f['low']:.3f}, {f['high']:.3f}]")
    sn = out["planted_pairs_asserted_under_another_identifier"]["all sources"]
    print(f"  planted pairs a source asserts under another identifier of the same name: {sn['n']} of {sn['of']}")
    print("  share of sampled real pairs no paper mentions, by whether the in-graph rule corroborates them")
    for s, r in order(out["in_graph_against_literature"]):
        print(f"    {s:34s} uncorroborated {pc(r['uncorroborated'])}  corroborated {pc(r['corroborated'])}  flagged by both {r['flagged_by_both']}")
    print(f"  spent on model calls: {out['spent_on_model_calls']:.3f} of {out['spend_cap']:.2f}")


def same_name(links: pd.DataFrame, sample: pd.DataFrame, name_of: dict) -> dict:
    """Planted pairs that some source does assert, for the same target, under another identifier whose
    name is the same. Pairs are told apart by identifier, so these were planted as false and are not."""
    named = set(zip(links.target, (str(name_of.get(k, "")).lower() for k in zip(links.source, links.disease)))) - {(t, "") for t in set(links.target)}
    pl = sample[sample.planted]
    hit = np.array([(t, str(n).lower()) in named for t, n in zip(pl.target, pl.name)])
    return {"all sources": {"n": int(hit.sum()), "of": int(len(pl))}, **{src: {"n": int(hit[(pl.source == src).to_numpy()].sum()), "of": int((pl.source == src).sum())} for src in sorted(set(pl.source))}}


def agreement(real: pd.DataFrame, detail: pd.DataFrame) -> dict:
    """For the sampled real one-source pairs: is a pair the in-graph rule leaves uncorroborated more often
    one that no paper mentions? Two checks that share no evidence, on the same pairs."""
    x = real.merge(detail[~detail.planted][KEY + ["corroborated"]], on=KEY)
    out = {}
    for src, g in list(x.groupby("source")) + [("all sources", x)]:
        unc, cor = g[~g.corroborated], g[g.corroborated]
        out[src] = {"uncorroborated": share((unc.both == 0).to_numpy()), "corroborated": share((cor.both == 0).to_numpy()),
                    "flagged_by_both": int(((g.both == 0) & ~g.corroborated).sum()), "pairs": int(len(g))}
    return out


def ingraph(release: str) -> dict:
    C, W = paths.context(release), work(release)
    pl, _ = frozen(release)
    links = pd.read_parquet(C / "disease_links.parquet")
    partners = pd.read_parquet(C / "partners.parquet")
    gwas = pd.read_parquet(paths.bundle(release) / "links.parquet", columns=["target", "trait"])
    terms = ontology.read(registry.require("efo_obo"))
    res, detail = round1.evaluate(links, pl[pl.in_graph_set][KEY], gwas, partners, terms)
    detail.to_parquet(W / "ingraph_detail.parquet", index=False)
    for s, r in res.items():
        r["matched"] = bool(MATCHED[0] <= r["auc_degree_product"] <= MATCHED[1])
    # Added after the first run, as description. The stated check counts a pair's OTHER pairs (the
    # degrees less the pair itself). Beside it: the same area on the degrees as they stand in the
    # source, which says whether a planted pair's target and disease are as common as a real pair's.
    for s, r in commonness(links, pl[pl.in_graph_set]).items():
        res[s].update(r)
    for s, r in standardised(detail).items():
        res[s].update(r)
    m = res["Monarch"]
    out = {"sources": res, "monarch_round_one_gaps": "kept" if m["rule_applied"] else "restored",
           "rule_applied_to": sorted(s for s, r in res.items() if s != "all sources" and r["rule_applied"])}
    (W / "ingraph.json").write_text(dump(out))
    if out["monarch_round_one_gaps"] == "restored":
        restore(release)
    return out


def commonness(links: pd.DataFrame, pl: pd.DataFrame) -> dict:
    """For each source: the area, real one-source against planted, of the product of the target's and
    the disease's number of pairs in the source, and the medians of the two for real one-source pairs."""
    dt = links.groupby(["source", "target"]).disease.nunique().to_dict()
    dd = links.groupby(["source", "disease"]).target.nunique().to_dict()
    out = {}
    for src, g in one_source(links).groupby("source"):
        p = pl[pl.source == src]
        rt, rd = np.array([dt[(src, t)] for t in g.target]), np.array([dd[(src, d)] for d in g.disease])
        pt, pd_ = np.array([dt[(src, t)] for t in p.target]), np.array([dd[(src, d)] for d in p.disease])
        out[src] = {"auc_degree_counting_the_pair_itself": round1.auc((rt * rd).astype(float), (pt * pd_).astype(float)),
                    "median_diseases_per_target": float(np.median(rt)), "median_targets_per_disease": float(np.median(rd))}
    return out


def standardised(detail: pd.DataFrame) -> dict:
    """For each source: the share of planted pairs withheld when they are weighted to have the real
    one-source pairs' spread of other-pair counts (fifths of the real pairs' degree product)."""
    out = {}
    for src, g in detail.groupby("source"):
        r, p = g[~g.planted], g[g.planted]
        edges = np.unique(np.quantile(r.degree_product, np.linspace(0, 1, 6)))
        rb, pb = np.searchsorted(edges, r.degree_product, side="right"), np.searchsorted(edges, p.degree_product, side="right")
        num = den = 0.0
        for k in np.unique(rb):
            m = pb == k
            if m.any():
                w = float((rb == k).mean())
                num += w * float((~p.corroborated[m]).mean())
                den += w
        out[src] = {"planted_withheld_at_the_real_pairs_other_pair_counts": num / den if den else float("nan"), "share_of_real_pairs_in_strata_with_a_planted_pair": den}
    return out


def restore(release: str) -> int:
    """Put back the pairs round one withheld: the pruned table is the whole table and nothing is withheld."""
    C = paths.context(release)
    links = pd.read_parquet(C / "disease_links.parquet")
    w = pd.read_parquet(C / "disease_links_withheld.parquet")
    links.to_parquet(C / "disease_links_pruned.parquet", index=False)
    w.iloc[0:0].to_parquet(C / "disease_links_withheld.parquet", index=False)
    return len(w)


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    # source s: 300 one-source pairs over 30 diseases, and 40 pairs a second source also asserts
    rows = [(f"T{i}", f"D{i % 30}", "s") for i in range(300)] + [(f"M{i}", f"D{i % 5}", "s") for i in range(40)] + [(f"M{i}", f"D{i % 5}", "u") for i in range(40)]
    links = pd.DataFrame(rows, columns=KEY)
    single = one_source(links)
    ck("one-source pairs are those no second source asserts", len(single) == 300 and set(single.source) == {"s"})
    pl = plant_matched(links, at_least=50, seed=1)
    ck("planting crosses one-source pairs only: no planted target comes from a pair with two sources", set(pl.target) <= set(single.target) and not set(pl.target) & {f"M{i}" for i in range(40)})
    ck("at least the asked number are planted, the first 5% marked as the in-graph set", len(pl) == 50 and int(pl.in_graph_set.sum()) == 15 and pl.in_graph_set.tolist() == [True] * 15 + [False] * 35)
    ck("no planted pair is a real pair or a repeat", not set(zip(pl.target, pl.disease)) & set(zip(links.target, links.disease)) and not pl.duplicated(["target", "disease"]).any())
    ck("planting is the same under the same seed", plant_matched(links, at_least=50, seed=1).equals(pl) and not plant_matched(links, at_least=50, seed=2).equals(pl))
    dr = draw_real(single, n=25, seed=1)
    ck("the real sample is drawn from the one-source pairs without repeats, in a fixed order",
       len(dr) == 25 and not dr.duplicated(["target", "disease"]).any() and dr.order.tolist() == list(range(25)) and draw_real(single, 25, 1).equals(dr))
    ck("the first 20 of a draw are development, the next 80 test", [part(0), part(19), part(20), part(99), part(100)] == ["dev", "dev", "test", "test", "rest"])
    sm = make_sample(links, plant_matched(links, at_least=30, seed=1), {("s", "D1"): "disease one"}, {"T1": "SYM1"})
    ck("a sample row carries its symbol and its disease name, or neither", sm[sm.target == "T1"].symbol.iloc[0] == "SYM1" and sm.symbol.isna().sum() > 0 and set(sm.name.dropna()) == {"disease one"})

    lo, hi = wilson(8, 10)
    ck("a Wilson interval is the textbook one", abs(lo - 0.4902) < 5e-4 and abs(hi - 0.9433) < 5e-4 and math.isnan(wilson(0, 0)[0]))
    v, lab = np.array([1.0, 2.0, 2.0, 3.0, 5.0]), np.array([False, False, True, True, True])
    codes = np.unique(v, return_inverse=True)[1]
    ck("the weighted area equals the plain one at unit weights, ties a half", abs(wauc(codes, lab) - round1.auc(v[lab], v[~lab])) < 1e-12 and abs(wauc(codes, lab) - (0.5 + 1 + 1 + 1 + 1 + 1) / 6) < 1e-12)
    ck("a row weighted twice counts twice", abs(wauc(codes, lab, np.array([1.0, 1, 2, 1, 1])) - (2 * 1.5 + 2 + 2) / 8) < 1e-12)
    ai = auc_interval(np.r_[np.zeros(50), np.ones(50)], np.r_[np.zeros(50, bool), np.ones(50, bool)], n_boot=50)
    ck("a field that separates perfectly has area 1 with no width", ai["auc"] == 1.0 and ai["low"] == 1.0 and ai["labelled"] == 50)
    ck("a cut matched to a real share withholds at least that share of real pairs", matched_cut(np.array([0.0, 0, 1, 2, 3, 4, 5, 6, 7, 8]), np.array([0.0, 0, 0, 1, 9]), 0.25) == (0.8, 0.3))
    ck("a cut matched to nothing withholds nothing", matched_cut(np.array([1.0, 2]), np.array([0.0]), 0.0) == (0.0, 0.0))
    ck("a number that is not a number is written as null", json.loads(dump({"a": float("nan"), "b": [1.5, float("nan")]})) == {"a": None, "b": [1.5, None]})

    cm = commonness(links, pl)["s"]
    ck("commonness is read on the degrees as they stand in the source", 0.3 < cm["auc_degree_counting_the_pair_itself"] < 0.7 and cm["median_diseases_per_target"] == 1.0 and cm["median_targets_per_disease"] == 10.0)

    det = pd.DataFrame({"source": ["s"] * 12, "planted": [False] * 6 + [True] * 6, "degree_product": [0, 0, 0, 9, 9, 9, 0, 9, 9, 9, 9, 9],
                        "corroborated": [False] * 6 + [False, True, True, True, True, True]})
    st = standardised(det)["s"]
    # real pairs are half at 0 and half at 9; planted withheld is 1.0 at 0 and 0.0 at 9, so 0.5 once weighted (0.167 unweighted)
    ck("planted pairs are reweighted to the real pairs' other-pair counts", abs(st["planted_withheld_at_the_real_pairs_other_pair_counts"] - 0.5) < 1e-12 and st["share_of_real_pairs_in_strata_with_a_planted_pair"] == 1.0)

    lk = pd.DataFrame([("A", "MONDO:1", "s"), ("B", "DOID:9", "u")], columns=KEY)
    sp_ = pd.DataFrame({"target": ["A", "B", "B"], "disease": ["DOID:9", "MONDO:1", "MONDO:7"], "source": ["s", "s", "u"], "planted": [True, True, True], "name": ["Gout", "gout", "other"]})
    sn = same_name(lk, sp_, {("s", "MONDO:1"): "gout", ("u", "DOID:9"): "Gout"})
    ck("a planted pair that a source asserts under another identifier of the same name is counted", sn["all sources"] == {"n": 2, "of": 3} and sn["u"] == {"n": 0, "of": 1})
    ag = agreement(pd.DataFrame({"target": list("abcd"), "disease": ["x"] * 4, "source": ["s"] * 4, "both": pd.array([0, 0, 3, 0], "Int64")}),
                   pd.DataFrame({"target": list("abcd"), "disease": ["x"] * 4, "source": ["s"] * 4, "planted": [False] * 4, "corroborated": [False, False, True, True]}))["s"]
    ck("the two checks are set side by side on the same real pairs", ag["uncorroborated"]["share"] == 1.0 and ag["corroborated"]["share"] == 0.5 and ag["flagged_by_both"] == 2)

    g1, d1, dx = {"count": 10, "not_found": []}, {"count": 500, "not_found": []}, {"count": 500, "not_found": ["zz"]}
    ck("a pair is unsearchable without a symbol, without a name, or when PubMed does not find the name",
       classify(None, "x", None, None, None)[1] == "the target has no gene symbol" and classify("A", "()", None, None, None)[1] == "the disease has no searchable name"
       and classify("A", "x", {"count": 3}, g1, dx) == (None, "PubMed did not find the disease name"))
    ck("a symbol no paper uses gives zero papers, whatever the pair search returned", classify("A", "x", {"count": 77}, {"count": 0, "not_found": ['"A"[tiab]']}, d1) == (0, None))
    ck("more papers for the pair than for either alone is not believed", classify("A", "x", {"count": 11}, g1, d1) == (None, "the counts do not agree") and classify("A", "x", {"count": 4}, g1, d1) == (4, None))

    se = pd.DataFrame({"target": list("abcdefgh"), "disease": ["x"] * 8, "source": ["s"] * 8, "planted": [True] * 4 + [False] * 4, "order": [0, 1, 2, 3] * 2, "part": ["dev"] * 8,
                       "both": pd.array([0, 3, 2, 5, 0, 4, 6, 1], "Int64"), "lift": [0.0, 1, 2, 3, 0, 9, 8, 7], "unsearchable": [None] * 7 + ["the target has no gene symbol"]})
    rd = pd.DataFrame({"target": list("bcdfg"), "disease": ["x"] * 5, "source": ["s"] * 5, "planted": [True, True, True, False, False], "model": ["m"] * 5,
                       "verdict": ["mentions_only", "asserts", None, "asserts", "mentions_only"], "pmid": [None, "1", None, "2", None], "quote": [None, "q", None, "q", None], "error": [None] * 5})
    x = withheld_l1(se, rd, "m").set_index("target")
    ck("L0 withholds a pair no paper mentions; L1 also one the reader flags; no verdict keeps the pair",
       x.l0.to_dict() == {"a": True, "b": False, "c": False, "d": False, "e": True, "f": False, "g": False}
       and x.l1.to_dict() == {"a": True, "b": True, "c": False, "d": False, "e": True, "f": False, "g": True} and x.no_verdict.to_dict()["d"] and "h" not in x.index)
    try:
        withheld_l1(se, rd[rd.target != "g"], "m")
        ck("a mentioned pair the model did not read stops the scoring", False)
    except SystemExit:
        ck("a mentioned pair the model did not read stops the scoring", True)
    rr = rule_rows(x.reset_index(), "l1")["s"]
    ck("shares are of searchable pairs, with an interval", rr["planted_withheld"]["n"] == 4 and rr["planted_withheld"]["share"] == 0.5 and rr["real_withheld"]["share"] == 2 / 3 and not rr["passes"])
    ck("the bar is 80% of planted and at most 50% of real, provisional when the interval reaches under 70%",
       bar(share([True] * 9 + [False]), share([False] * 10)) == {"passes": True, "provisional": True} and bar(share([True] * 190 + [False] * 10), share([True] * 90 + [False] * 110)) == {"passes": True, "provisional": False}
       and not bar(share([True] * 10), share([True] * 6 + [False] * 4))["passes"])
    sc = pd.DataFrame({"source": ["s"] * 20, "planted": [True] * 10 + [False] * 10, "part": (["rest"] * 5 + ["test"] * 5) * 2,
                       "lift": [0.0, 0, 0, 1, 9, 0, 0, 1, 2, 9] + [0.0, 5, 5, 5, 5, 0, 1, 5, 5, 5]})
    sp = split_cut(sc)["s"]
    # the cut is the 80th percentile of the five planted lifts it is set on (0, 0, 0, 1, 9), which is 2.6
    ck("a lift cut is set on one part and scored on the other", abs(sp["cut"] - 2.6) < 1e-9 and sp["planted_withheld"]["share"] == 0.8 and sp["real_withheld"]["share"] == 0.4 and sp["planted_withheld"]["n"] == 5)
    t = pd.DataFrame({"planted": [True] * 10 + [False] * 10, "lift": [0.0] * 5 + [1.0] * 5 + [0.0] * 2 + [5.0] * 8, "l1": [True] * 9 + [False] + [True] * 2 + [False] * 8})
    ra = reader_against_lift(t)
    ck("the reader is set against a lift cut of the same cost in real pairs", ra["l1_planted_withheld"] == 0.9 and ra["lift_cut_planted_withheld"] == 0.5 and ra["lift_cut_real_withheld"] == 0.2 and abs(ra["difference"] - 0.4) < 1e-12)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) < 2:
        raise SystemExit(__doc__)
    stage, release = argv[0], argv[1]
    if stage == "freeze":
        out = freeze(release)
        for s, r in out["sources"].items():
            print(f"  {s:34s} one-source {r['one_source_pairs']:>8,}  planted {r['planted']:>6,}  in-graph set {r['in_graph_set']:>6,}  sample {r['sample_planted']}+{r['sample_real']}  no symbol {r['sample_without_a_symbol']}  no name {r['sample_without_a_name']}")
    elif stage == "ingraph":
        out = ingraph(release)
        for s, r in sorted(out["sources"].items(), key=lambda kv: (kv[0] == "all sources", -kv[1]["pairs"])):
            print(f"  {s:34s} planted {r['planted']:>6,}  planted withheld {r['planted_withheld']:.3f}  real withheld {r['real_single_source_withheld']:.3f}  AUC neighbours {r['auc_neighbours']:.3f}  "
                  f"AUC degree {r['auc_degree_product']:.3f}  matched {r['matched']}  applied {r['rule_applied']}")
        print(f"  Monarch's round-one gaps: {out['monarch_round_one_gaps']}")
    elif stage == "search":
        out = search(release)
        print(f"  {len(out):,} pairs; unsearchable {int(out.unsearchable.notna().sum())}; no paper mentions {int((out.both == 0).sum())}")
    elif stage == "read":
        which = argv[2]
        from pharos_gap_kg.usability import PRICES
        models = argv[3:] or (list(PRICES) if which == "dev" else [json.loads((paths.DOCS_DATA / f"pruning2_reader_{release}.json").read_text())["chosen"]])
        key = os.environ.get("OPENROUTER_API_KEY")
        if not key:
            raise SystemExit("set OPENROUTER_API_KEY")
        read(release, which, models, key, "https://openrouter.ai/api/v1")
    elif stage == "choose":
        out = choose(release)
        for m, r in out["models"].items():
            print(f"  {m:30s} read {r['read']}  passed the check {r['passed_the_check']:.3f}  planted withheld {r['planted_withheld']:.3f}  real withheld {r['real_withheld']:.3f}  difference {r['difference']:+.3f}  qualifies {r['qualifies']}")
        print(f"  chosen: {out['chosen']}")
    elif stage == "report":
        show(report(release))
        print(f"  -> {paths.DOCS_DATA / f'pruning2_{release}.json'}")
    else:
        raise SystemExit(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
