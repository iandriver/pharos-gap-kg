"""The prospective test: links first reported from 2023, scored against everything Pharos 3.19 holds
(docs/PROTOCOL.md, "A prospective test" and the addition of 10 Oct 2026 that fixes what is scored).

The split is built by `gwascat`. Training is Pharos's pairs by live trait term. A pair Pharos lacks
that was first reported by 2022 is set aside: it is never a positive and never a candidate.

What is scored, each once:

    B1 B3 B4 B5                 the parameter-free scorers, from all training links
    combine ... combine_lin_pos the four combiners with the weights frozen on 9 Oct
    mf, lightgcn_lin_pos        the 2018 split's models with their frozen settings and passes
    lightgcn_lin_pos:L05:steps  the rehearsal's best graph arm, tuned with a cutoff of 2020
    combine[B1+B3+B5]:F1        the rehearsal arm that carried, weights fitted with a cutoff of 2020
    combine[B1+B3+B4+B5]:F1     the rehearsal's highest combiner arm, the same way

Stages. The first three read training links only. `score` and `coldstart` read the prospective links.

    python -m pharos_gap_kg.prospect tune pharos319                       # the graph candidate's setting, at 2020
    python -m pharos_gap_kg.prospect weights pharos319                    # the two candidate combiners' weights, at 2020
    python -m pharos_gap_kg.prospect fit pharos319 --model lightgcn_lin_pos    # all training links, five seeds; once
    python -m pharos_gap_kg.prospect score pharos319
    python -m pharos_gap_kg.prospect coldstart pharos319
    python -m pharos_gap_kg.prospect --selftest
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import baselines, combine, controls, graph, model, paths, position, rehearse, score, split

CUTOFF, INNER = 2022, 2020
FROZEN = ("mf", "lightgcn_lin_pos")                               # the 2018 split's models, as frozen
CANDIDATE = {"name": "lightgcn_lin_pos:L05:steps", "model": "lightgcn_lin_pos", "draw": "L05", "carry": "steps", "rehearsal": "did not carry"}
COMBINERS = {"combine[B1+B3+B5]:F1": {"scores": ("B1", "B3", "B5"), "form": "F1", "rehearsal": "carried"},
             "combine[B1+B3+B4+B5]:F1": {"scores": ("B1", "B3", "B4", "B5"), "form": "F1", "rehearsal": "did not carry"}}
FITTED = (*FROZEN, CANDIDATE["name"])
VIEWS = ("full", "low", "neutral", "locus")


def safe(name: str) -> str:
    return name.replace(":", "_").replace("[", "_").replace("]", "").replace("+", "-")


def model_dir(name: str) -> Path:
    return split.out_dir(CUTOFF) / "models" / safe(name)


def chosen_path() -> Path:
    return paths.DOCS_DATA / f"prospect_candidate_{CUTOFF}.json"


def weights_path() -> Path:
    return paths.DOCS_DATA / f"prospect_weights_{CUTOFF}.json"


def training() -> pd.DataFrame:
    """The training part alone, checked against the frozen split. No prospective link is read."""
    tr = pd.read_parquet(split.out_dir(CUTOFF) / "train.parquet")
    if split.digest(tr) != json.loads(split.frozen_path(CUTOFF).read_text())["train"]["sha256"]:
        raise SystemExit("REFUSED: train.parquet is not the frozen training part of the prospective split")
    return tr


def train_matrix(tr: pd.DataFrame) -> dict:
    targets, traits = sorted(tr.target.unique()), sorted(tr.trait.unique())
    ti, di = {t: i for i, t in enumerate(targets)}, {d: i for i, d in enumerate(traits)}
    A = np.zeros((len(targets), len(traits)), bool)
    A[tr.target.map(ti).to_numpy(), tr.trait.map(di).to_numpy()] = True
    return {"train": A, "targets": targets, "traits": traits}


def inner(tr: pd.DataFrame) -> tuple[dict, dict]:
    """The split used to tune a candidate: fit on the links to INNER, validate on the later training links."""
    s = split.split(tr, INNER)
    M = baselines.matrices({"train": s["train"], "test": s["test"]})
    return M, {"fit": {"links": len(s["train"]), "sha256": split.digest(s["train"])}, "validate": {"links": len(s["test"]), "sha256": split.digest(s["test"])}}


def pair_matrix(pairs: pd.DataFrame, targets: list[str], traits: list[str]) -> np.ndarray:
    """[target, trait] true for the pairs whose target and trait are both listed."""
    ti, di = {t: i for i, t in enumerate(targets)}, {d: i for i, d in enumerate(traits)}
    p = pairs[pairs.target.isin(ti) & pairs.trait.isin(di)]
    X = np.zeros((len(targets), len(traits)), bool)
    X[p.target.map(ti).to_numpy(), p.trait.map(di).to_numpy()] = True
    return X


# ---- stages that read training links only ----------------------------------------------------------------------

def tune(release: str) -> dict:
    """The graph candidate's eight settings at a 2020 cutoff, the best on the links of 2021 and 2022."""
    tr = training()
    M, parts = inner(tr)
    A, Y = M["train"], M["test"]
    m, draw = CANDIDATE["model"], CANDIDATE["draw"]
    edges, n, _ = model.graph_for(release, m, M)
    base = {k: model.primary(v, A, Y) for k, v in baselines.scorers(A).items() if k in ("B1", "B3")}
    rows = []
    for i, s in enumerate(rehearse.settings(m)):
        f = model_dir(CANDIDATE["name"]) / "tune" / f"{i}.json"
        if not f.is_file():
            t0 = time.time()
            _, hist = model.fit(A, s, 0, model.device(), model.MAX_PASSES, edges, n, validate=lambda S: model.primary(S, A, Y), balance=rehearse.DRAWS[draw])
            best = max(hist, key=lambda h: h["validation"])
            model.put(f, json.dumps({**s, "validation": best["validation"], "passes": best["passes"], "passes_run": hist[-1]["passes"], "seconds": round(time.time() - t0, 1),
                                     "history": [[h["passes"], round(h["validation"], 6)] for h in hist], "parts": parts, "implementation": model.implementation()}) + "\n")
        rows.append(json.loads(f.read_text()))
        print(f"    {i + 1}/8  {s}  validation {rows[-1]['validation']:.4f} at {rows[-1]['passes']} passes  ({rows[-1]['seconds']} s)", flush=True)
    if any(r["parts"] != json.loads(json.dumps(parts)) or r["implementation"] != model.implementation() for r in rows):
        raise SystemExit("REFUSED: a kept tuning result is from other links or other fitting code")
    top = max(rows, key=lambda r: r["validation"])
    out = {**{k: CANDIDATE[k] for k in ("name", "model", "draw", "carry", "rehearsal")}, "cutoff": CUTOFF, "inner_cutoff": INNER, "parts": parts, "baselines_on_validation": base,
           "chosen": {k: top[k] for k in ("dim", "layers", "lr", "l2", "passes")}, "validation": top["validation"],
           "passes_on_all_links": rehearse.carried(top["passes"], parts["fit"]["links"], len(tr), CANDIDATE["carry"]),
           "tried": [{k: v for k, v in r.items() if k not in ("parts", "implementation")} for r in rows], "implementation": model.implementation()}
    model.put(chosen_path(), json.dumps(out, indent=1) + "\n")
    return out


def weights(release: str) -> dict:
    """The two candidate combiners' weights, fitted at a 2020 cutoff on the links of 2021 and 2022."""
    tr = training()
    M, parts = inner(tr)
    A, Y = M["train"], M["test"]
    pairs, R, used = model.load_structure(release, M["targets"], M["traits"], ("position", "lineage"))
    S = combine.raw_scores(A, pairs, R)
    cols = np.flatnonzero(Y.any(axis=0))
    b3 = score.by_axis(S["B3"][:, cols], A[:, cols], Y[:, cols], axis=1, seeds=(0,))["ap"]
    out = {"cutoff": CUTOFF, "inner_cutoff": INNER, "parts": parts, "traits": int(len(cols)), "neighbours_on_validation": float(b3.mean()), "reference_files": used,
           "implementation": combine.implementation(), "combiners": {}}
    for name, c in COMBINERS.items():
        t0 = time.time()
        X = combine.design(S, A, c["scores"], c["form"], cols)
        w, fitted = combine.fit_weights(X, ~A[:, cols], Y[:, cols], c["scores"].index("B3"))
        out["combiners"][name] = {"scores": list(c["scores"]), "form": c["form"], "rehearsal": c["rehearsal"], "terms": combine.terms(c["scores"], c["form"]), "weights": [float(x) for x in w],
                                  "fitted": fitted, "read_as": "the weights were fitted on these links, so this figure flatters the combiner"}
        print(f"    {name}: fitted {fitted:.4f} (neighbours {b3.mean():.4f})  " + " ".join(f"{x:g}" for x in w) + f"  ({time.time() - t0:.0f} s)", flush=True)
    model.put(weights_path(), json.dumps(out, indent=1) + "\n")
    return out


def recipe(name: str, n_train: int) -> dict:
    """The setting, passes and way of drawing a fitted model is to be fitted with on all training links."""
    if name in FROZEN:
        ch = json.loads(model.chosen_path(2018, name).read_text())["chosen"]
        return {"model": name, "setting": {k: ch[k] for k in ("dim", "layers", "lr", "l2")}, "passes": ch["passes"], "balance": 0.0, "from": f"model_{name}_2018.json, as frozen"}
    if name == CANDIDATE["name"]:
        if not chosen_path().is_file():
            raise SystemExit("the candidate has no chosen setting; run:  python -m pharos_gap_kg.prospect tune <release>  and commit the file it writes")
        ch = json.loads(chosen_path().read_text())
        return {"model": CANDIDATE["model"], "setting": {k: ch["chosen"][k] for k in ("dim", "layers", "lr", "l2")}, "passes": ch["passes_on_all_links"],
                "balance": rehearse.DRAWS[CANDIDATE["draw"]], "from": chosen_path().name}
    raise SystemExit(f"{name} is not a fitted model of the prospective test")


def fit(release: str, name: str) -> dict:
    """Fit on all training links under five seeds and keep the mean matrix. Once. No prospective link is read."""
    d = model_dir(name)
    if (d / "final.json").is_file():
        raise SystemExit(f"REFUSED: {name} has been fitted for the prospective test ({d / 'final.json'}); it is fitted once")
    tr = training()
    M = train_matrix(tr)
    rc = recipe(name, len(tr))
    model.claim(d, {"model": name, "cutoff": CUTOFF, "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "recipe": rc, "implementation": model.implementation()})
    edges, n, joined = model.graph_for(release, rc["model"], M)
    mats = []
    for seed in model.SEEDS:
        t0 = time.time()
        S, _ = model.fit(M["train"], rc["setting"], seed, model.device(), rc["passes"], edges, n, balance=rc["balance"])
        mats.append(S)
        print(f"    {name} seed {seed}: {rc['passes']} passes  ({time.time() - t0:.0f} s)", flush=True)
    S = np.mean(mats, axis=0).astype(np.float32)
    np.save(d / "scores.npy", S)
    np.save(d / "seeds.npy", np.stack(mats).astype(np.float32))
    out = {"model": name, "cutoff": CUTOFF, "recipe": rc, "joined": joined, "training_links": int(len(tr)), "scores_sha256": graph.file_sha(d / "scores.npy"),
           "seeds_sha256": graph.file_sha(d / "seeds.npy"), "implementation": model.implementation()}
    model.put(d / "final.json", json.dumps(out, indent=1) + "\n")
    return out


# ---- scoring ---------------------------------------------------------------------------------------------------

def load_all() -> dict:
    parts = split.load(CUTOFF)
    aside = pd.read_parquet(split.out_dir(CUTOFF) / "aside.parquet")
    if split.digest(aside) != json.loads(split.frozen_path(CUTOFF).read_text())["aside"]["sha256"]:
        raise SystemExit("REFUSED: aside.parquet is not the frozen set-aside part")
    return {**parts, "aside": aside}


def views(S: np.ndarray, taken: np.ndarray, Y: np.ndarray, linked: np.ndarray, pairs: np.ndarray, seeds=score.SEEDS) -> dict[str, np.ndarray]:
    """Per-trait average precision on the four views. `taken` marks the pairs that are not candidates:
    training links and pairs set aside."""
    new = Y & ~linked
    return {"full": score.by_axis(S, taken, Y, axis=1, seeds=seeds)["ap"], "low": score.by_axis(S, taken, new, axis=1, exclude=Y & linked, seeds=seeds)["ap"],
            "neutral": score.by_axis(S, taken, new, axis=1, exclude=linked, seeds=seeds)["ap"], "locus": controls.collapsed(S, taken, Y, pairs, seeds=seeds)["ap"]}


def paired(per: dict, a: str, b: str) -> dict:
    out = {}
    for v in VIEWS:
        x, y = per[a][v], per[b][v]
        if len(x) != len(y):          # the two new-link views drop traits with no positive left; they are the same traits for every scorer
            raise SystemExit(f"REFUSED: {a} and {b} were not scored on the same traits on the view {v}")
        m, lo, hi = score.interval(x - y)
        out[v] = {"difference": m, "low": lo, "high": hi, "half_width": (hi - lo) / 2, "traits_higher": int((x > y).sum()), "traits": int(len(x))}
    return out


def scorers(release: str, M: dict) -> tuple[dict[str, np.ndarray], dict]:
    """Every scorer's matrix, from the training links and the frozen records. No prospective link is read here."""
    A = M["train"]
    pairs, R, used = model.load_structure(release, M["targets"], M["traits"], ("position", "lineage"))
    raw = combine.raw_scores(A, pairs, R)
    S = {"B1": raw["B1"], "B3": raw["B3"], "B4": raw["B4"], "B5": raw["B5"]}
    notes = {}
    frozen = json.loads(combine.record_path(2018).read_text())
    if frozen["implementation"] != combine.implementation() or json.loads(combine.frozen_path(2018).read_text())["record_sha256"] != graph.file_sha(combine.record_path(2018)):
        raise SystemExit("REFUSED: the combiners' record is not the one frozen on 9 Oct, or the code that applies it has changed")
    for name, c in frozen["combiners"].items():
        S[name] = combine.combined(raw, A, combine.COMBINERS[name], frozen["form"], c["weights"])
        notes[name] = {"kind": "the 2018 split's combiner, weights as frozen", "form": frozen["form"], "weights": dict(zip(c["terms"], c["weights"]))}
    w = json.loads(weights_path().read_text())
    if w["implementation"] != combine.implementation():
        raise SystemExit("REFUSED: the candidate combiners' weights were fitted by other code")
    for name, c in w["combiners"].items():
        S[name] = combine.combined(raw, A, tuple(c["scores"]), c["form"], c["weights"])
        notes[name] = {"kind": "rehearsal candidate, weights fitted at 2020", "rehearsal": c["rehearsal"], "weights": dict(zip(c["terms"], c["weights"]))}
    for name in FITTED:
        d = model_dir(name)
        fin = json.loads((d / "final.json").read_text())
        if graph.file_sha(d / "scores.npy") != fin["scores_sha256"]:
            raise SystemExit(f"REFUSED: {d / 'scores.npy'} is not the matrix its final.json describes")
        S[name] = np.load(d / "scores.npy")
        notes[name] = {"kind": "the 2018 split's model, setting and passes as frozen" if name in FROZEN else "rehearsal candidate, tuned at 2020", "recipe": fin["recipe"],
                       **({"rehearsal": CANDIDATE["rehearsal"]} if name == CANDIDATE["name"] else {})}
    return S, {"pairs": pairs, "R": R, "used": used, "notes": notes}


def scored(release: str) -> dict:
    """Every scorer on the prospective test links. The matrices are fixed before this runs, so it gives
    the same answer however often it is run; the marker records when the links were first read."""
    marker = split.out_dir(CUTOFF) / "scored_started.json"
    if not marker.is_file():
        model.put(marker, json.dumps({"started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "scorers": ["B1", "B3", "B4", "B5", *combine.COMBINERS, *COMBINERS, *FITTED],
                                      "test_sha256": json.loads(split.frozen_path(CUTOFF).read_text())["test"]["sha256"]}, indent=1) + "\n")
    P = load_all()
    M = train_matrix(P["train"])
    A = M["train"]
    Y = pair_matrix(P["test"], M["targets"], M["traits"])
    if int(Y.sum()) != len(P["test"]) or (A & Y).any():
        raise SystemExit("REFUSED: a prospective link is a training link, or its target or trait is not a training one")
    taken = A | pair_matrix(P["aside"], M["targets"], M["traits"])
    if (taken & Y).any():
        raise SystemExit("REFUSED: a prospective link is among the pairs set aside")
    S, X = scorers(release, M)
    if list(S) != json.loads(marker.read_text())["scorers"]:
        raise SystemExit("REFUSED: the scorers are not the ones written down when the prospective links were first read")
    linked = (((A.astype(np.float32) @ X["R"].astype(np.float32).T) > 0) | (position.neighbour_count(X["pairs"], A) > 0)) & ~taken
    per, res, extra = {}, {}, {}
    for k, s in S.items():
        t0 = time.time()
        per[k] = views(s, taken, Y, linked, X["pairs"])
        full = score.by_axis(s, taken, Y, axis=1)
        res[k] = {v: score.interval(a) for v, a in per[k].items()}
        extra[k] = {"trait_auc": score.interval(full["auc"]), "trait_recall_at_100": score.interval(full["recall"])}
        print(f"    {k:28s} all {res[k]['full'][0]:.4f}  new, lower bound {res[k]['low'][0]:.4f}  new, neutral {res[k]['neutral'][0]:.4f}  by locus {res[k]['locus'][0]:.4f}  ({time.time() - t0:.0f} s)", flush=True)
    pr = {f"{k}-B3": paired(per, k, "B3") for k in S if k != "B3"}
    pr.update({f"{k}-B1": paired(per, k, "B1") for k in S if k not in ("B1", "B3")})
    for k in ("combine_pos", "combine_lin", "combine_lin_pos"):
        pr[f"{k}-combine"] = paired(per, k, "combine")
    pr["lightgcn_lin_pos:L05:steps-lightgcn_lin_pos"] = paired(per, CANDIDATE["name"], "lightgcn_lin_pos")
    seeds = {}
    for name in FITTED:
        mats = np.load(model_dir(name) / "seeds.npy")
        seeds[name] = [float(np.nanmean(score.by_axis(m_, taken, Y, axis=1)["ap"])) for m_ in mats]
    year = np.zeros(A.shape, np.int16)
    ti, di = {t: i for i, t in enumerate(M["targets"])}, {d: i for i, d in enumerate(M["traits"])}
    year[P["test"].target.map(ti).to_numpy(), P["test"].trait.map(di).to_numpy()] = P["test"].year.to_numpy()
    by_year = {}
    for y in sorted(int(v) for v in np.unique(year[Y])):
        this = Y & (year == y)
        by_year[str(y)] = {"links": int(this.sum()), "traits": int(this.any(axis=0).sum()),
                           **{k: float(np.nanmean(score.by_axis(S[k], taken, this, axis=1, exclude=Y & ~this, seeds=(0,))["ap"])) for k in S}}
    fz = json.loads(split.frozen_path(CUTOFF).read_text())
    out = {"cutoff": CUTOFF, "release": release, "split": {k: fz[k]["sha256"] for k in (*split.PARTS, "aside")}, "reference_files": X["used"],
           "counts": {"training_links": int(A.sum()), "targets": A.shape[0], "traits": A.shape[1], "test_links": int(Y.sum()), "traits_scored": int(Y.any(axis=0).sum()),
                      "set_aside_between_a_training_target_and_trait": int((taken & ~A).sum()), "candidates": int((~taken).sum()),
                      "positives_lineage_or_position_linked": int((Y & linked).sum()), "candidates_linked": int(linked.sum())},
           "scorers": X["notes"], "results": res, "secondary": extra, "paired": pr, "by_seed": seeds, "by_year": by_year, "scoring": json.loads(marker.read_text())}
    (paths.DOCS_DATA / f"prospect_gene_trait_{CUTOFF}.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def cold(release: str) -> dict:
    """The cold-start scorers on the prospective split's new traits and new targets."""
    import pyarrow.parquet as pq
    from pharos_gap_kg import coldstart, ontology, registry
    P = load_all()
    terms = ontology.read(registry.require("efo_obo"))
    tg = position.target_genes(release, position.genes(registry.require("ensembl_gtf_113")))
    prot = pq.read_table(paths.tables(release) / "protein.parquet", columns=["uniprot"]).to_pandas()
    all_targets = sorted(set("UniProtKB:" + prot.uniprot.dropna()))
    out = coldstart.evaluate(release, P["train"], P["unscored"], terms, tg, all_targets, log=lambda s: print(s, flush=True), aside=P["aside"])
    fz = json.loads(split.frozen_path(CUTOFF).read_text())
    rec = {"cutoff": CUTOFF, "release": release, "split": {k: fz[k]["sha256"] for k in (*split.PARTS, "aside")}, "scorers": {**coldstart.TRAIT_SCORERS, **coldstart.TARGET_SCORERS, **coldstart.BY_TARGET}, **out}
    (paths.DOCS_DATA / f"prospect_coldstart_{CUTOFF}.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("three fitted models: the two frozen on the 2018 split and the rehearsal's graph candidate", FITTED == ("mf", "lightgcn_lin_pos", "lightgcn_lin_pos:L05:steps"))
    ck("two candidate combiners, the one that carried and the one highest in value, both holding neighbours",
       {c["rehearsal"] for c in COMBINERS.values()} == {"carried", "did not carry"} and all("B3" in c["scores"] and c["form"] == "F1" for c in COMBINERS.values()))
    ck("each scorer keeps its own folder", len({safe(n) for n in (*FITTED, *COMBINERS)}) == 5 and ":" not in safe(CANDIDATE["name"]))
    pairs = pd.DataFrame({"target": ["t1", "t2", "zz"], "trait": ["a", "b", "a"]})
    X = pair_matrix(pairs, ["t1", "t2"], ["a", "b"])
    ck("a pair whose target or trait is not listed is left out of the matrix", X.tolist() == [[True, False], [False, True]])
    tr = pd.DataFrame({"target": ["t1", "t2", "t1", "t2", "t3"], "trait": ["a", "a", "b", "b", "a"], "year": [2015, 2019, 2020, 2021, 2022]})
    M, parts = inner(tr)
    ck("a candidate is tuned on the links to 2020 against the later training links between a target and trait it has seen",
       parts["fit"]["links"] == 3 and parts["validate"]["links"] == 1 and M["train"].sum() == 3 and M["test"].sum() == 1)
    # six targets, one trait: target 0 is a training link, 1 is set aside, 2 and 3 are positives, 4 and 5 are not
    taken = np.array([[1], [1], [0], [0], [0], [0]], bool)
    Y = np.array([[0], [0], [1], [1], [0], [0]], bool)
    S = np.array([[9.0], [8.0], [7.0], [1.0], [5.0], [0.5]])
    none = np.zeros_like(Y)
    v = views(S, taken, Y, none, np.zeros((0, 2), np.int64), seeds=(0,))
    ck("a pair set aside is not a candidate: it cannot be ranked above a positive", abs(v["full"][0] - (1 / 1 + 2 / 3) / 2) < 1e-12 and abs(v["locus"][0] - v["full"][0]) < 1e-12)
    linked = np.array([[0], [0], [1], [0], [1], [0]], bool)
    v2 = views(S, taken, Y, linked, np.zeros((0, 2), np.int64), seeds=(0,))
    # new positive is target 3 only. lower bound: linked positive 2 removed, linked negative 4 stays above it -> rank 2 of 3. neutral: 2 and 4 removed -> rank 1 of 2
    ck("the two ways of setting linked pairs aside: removing only the linked positive, or every linked candidate", abs(v2["low"][0] - 1 / 2) < 1e-12 and abs(v2["neutral"][0] - 1.0) < 1e-12)
    per = {"a": {k: np.array([0.5, 0.2]) for k in VIEWS}, "b": {k: np.array([0.1, 0.2]) for k in VIEWS}}
    ck("a paired difference is taken trait by trait", abs(paired(per, "a", "b")["full"]["difference"] - 0.2) < 1e-12 and paired(per, "a", "b")["full"]["traits_higher"] == 1)
    ck("carried by steps, the passes shrink with the links", rehearse.carried(100, 197553, 294954, "steps") == 67)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=("tune", "weights", "fit", "score", "coldstart"))
    ap.add_argument("release")
    ap.add_argument("--model", choices=FITTED)
    a = ap.parse_args(argv)
    if a.stage == "tune":
        out = tune(a.release)
        print(f"  {out['name']}: {out['chosen']} validation {out['validation']:.4f} (neighbours {out['baselines_on_validation']['B3']:.4f}); {out['passes_on_all_links']} passes on all links")
    elif a.stage == "weights":
        weights(a.release)
    elif a.stage == "fit":
        if not a.model:
            raise SystemExit("fit needs --model")
        out = fit(a.release, a.model)
        print(f"  {a.model}: fitted on {out['training_links']:,} links, {out['recipe']}")
    elif a.stage == "score":
        out = scored(a.release)
        for k, v in out["paired"].items():
            if k.endswith("-B3"):
                print(f"  {k:36s} " + "  ".join(f"{x} {y['difference']:+.4f} [{y['low']:+.4f}, {y['high']:+.4f}]" for x, y in v.items()))
    else:
        cold(a.release)
    return 0


if __name__ == "__main__":
    sys.exit(main())
