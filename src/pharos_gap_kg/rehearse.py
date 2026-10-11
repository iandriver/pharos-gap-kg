"""The tune, refit, score procedure rehearsed one step back in time (docs/PROTOCOL.md, "Rehearsal inside
the training years"). Only training links are read.

    tuning   fits on the links to 2014 and validates on those of 2015 and 2016
    refit    on the links to 2016, under five seeds
    target   the links of 2017 and 2018, which were the validation positives of the scored models

An arm is a whole procedure. For the three fitted models it is how a link is drawn in fitting (by
link L0, by the square root of its trait's size L05, by trait L1) and how the length of fitting is
carried from tuning to the refit (the same passes, or the same steps). For the combiners it is which
scores are combined and in which form. What an arm does on the target says whether the procedure
carries forward. It is evidence about a procedure and not a test of a model: the target has been used
to tune models before, and the arms were chosen knowing that.

    python -m pharos_gap_kg.rehearse tune pharos319 --model lightgcn --draw L1
    python -m pharos_gap_kg.rehearse refit pharos319 --model lightgcn --draw L1
    python -m pharos_gap_kg.rehearse combiners pharos319
    python -m pharos_gap_kg.rehearse report pharos319
    python -m pharos_gap_kg.rehearse --selftest
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import baselines, combine, model, paths, score, split

CUTOFF, REFIT, TUNE = 2018, 2016, 2014          # links first reported after CUTOFF are never read here
MODELS = ("mf", "lightgcn", "lightgcn_lin_pos")
DRAWS = {"L0": 0.0, "L05": 0.5, "L1": 1.0}
CARRIES = ("passes", "steps")
GRID = {"dim": (64, 128), "layers": (2, 4), "lr": (0.001, 0.005), "l2": (1e-4,)}
SETS = (("B3", "B5"), ("B1", "B3"), ("B1", "B3", "B5"), ("B1", "B3", "B4"), ("B1", "B3", "B4", "B5"))
FORMS = ("F1", "F3")


def settings(m: str) -> list[dict]:
    layers = GRID["layers"] if model.MODELS[m]["graph"] else (0,)
    return [{"dim": d, "layers": k, "lr": lr, "l2": l2} for d, k, lr, l2 in itertools.product(GRID["dim"], layers, GRID["lr"], GRID["l2"])]


def carried(passes: int, tune_links: int, refit_links: int, carry: str) -> int:
    """The passes of the refit: the passes tuning chose, or as many as make the same number of steps."""
    if carry == "passes":
        return passes
    return max(1, round(passes * tune_links / refit_links))


def model_arms() -> list[tuple[str, str, str]]:
    return [(m, d, c) for m in MODELS for d in DRAWS for c in CARRIES]


def arm_name(m: str, draw: str, carry: str) -> str:
    return f"{m}:{draw}:{carry}"


def combiner_name(names: tuple[str, ...], form: str) -> str:
    return f"combine[{'+'.join(names)}]:{form}"


def periods(links: pd.DataFrame) -> tuple[dict, dict]:
    """The two splits of the rehearsal, from links first reported up to CUTOFF: (tuning, refit)."""
    if int(links.year.max()) > CUTOFF:
        raise SystemExit("REFUSED: the rehearsal reads training links only")
    at16 = split.split(links, REFIT)
    at14 = split.split(at16["train"], TUNE)
    return at14, at16


def matrices() -> tuple[dict, dict, dict]:
    """(matrices for tuning, matrices for the refit, digests). Read from the training part alone."""
    d = split.out_dir(CUTOFF)
    tr = pd.read_parquet(d / "train.parquet")
    frozen = json.loads(split.frozen_path(CUTOFF).read_text())
    if split.digest(tr) != frozen["train"]["sha256"]:
        raise SystemExit(f"REFUSED: {d / 'train.parquet'} is not the frozen training part")
    at14, at16 = periods(tr)
    M14 = baselines.matrices({"train": at14["train"], "test": at14["test"]})
    M16 = baselines.matrices({"train": at16["train"], "test": at16["test"]})
    parts = {"tune_fit": {"links": len(at14["train"]), "sha256": split.digest(at14["train"])}, "tune_validate": {"links": len(at14["test"]), "sha256": split.digest(at14["test"])},
             "refit": {"links": len(at16["train"]), "sha256": split.digest(at16["train"])}, "target": {"links": len(at16["test"]), "sha256": split.digest(at16["test"])}}
    return M14, M16, parts


def work() -> Path:
    return split.out_dir(CUTOFF) / "rehearsal"


def out_path() -> Path:
    return paths.DOCS_DATA / f"rehearsal_gene_trait_{CUTOFF}.json"


def setting_file(m: str, draw: str, s: dict, parts: dict) -> Path:
    what = {"model": m, "draw": draw, "setting": s, "parts": {k: v["sha256"] for k, v in parts.items()}, "how": [model.BATCH, model.MAX_PASSES, model.EVERY, model.PATIENCE], "code": model.implementation()}
    return work() / "tune" / f"{m}_{draw}" / (hashlib.sha256(json.dumps(what, sort_keys=True).encode()).hexdigest()[:20] + ".json")


def tune(release: str, m: str, draw: str) -> dict:
    """Try the arm's settings at 2014 and choose the best on the links of 2015 and 2016."""
    M14, _, parts = matrices()
    A, Y = M14["train"], M14["test"]
    edges, n, _ = model.graph_for(release, m, M14)
    ap = lambda S: score.by_axis(S, A, Y, axis=1, seeds=(0,))["ap"]
    rows = []
    for i, s in enumerate(settings(m)):
        f = setting_file(m, draw, s, parts)
        if not f.is_file():
            t0 = time.time()
            S, hist = model.fit(A, s, 0, model.device(), model.MAX_PASSES, edges, n, validate=lambda S: model.primary(S, A, Y), balance=DRAWS[draw])
            best = max(hist, key=lambda h: h["validation"])
            model.put(f, json.dumps({**s, "validation": best["validation"], "passes": best["passes"], "passes_run": hist[-1]["passes"], "seconds": round(time.time() - t0, 1),
                                     "history": [[h["passes"], round(h["validation"], 6)] for h in hist], "ap": ap(S).tolist()}) + "\n")
        rows.append(json.loads(f.read_text()))
        print(f"    {m} {draw} {i + 1}/{len(settings(m))}  {s}  validation {rows[-1]['validation']:.4f} at {rows[-1]['passes']} passes  ({rows[-1]['seconds']} s)", flush=True)
    top = max(rows, key=lambda r: r["validation"])          # the first of equals, in the order tried
    out = {"model": m, "draw": draw, "balance": DRAWS[draw], "chosen": {k: top[k] for k in ("dim", "layers", "lr", "l2", "passes")}, "validation": top["validation"], "ap": top["ap"],
           "tried": [{k: v for k, v in r.items() if k != "ap"} for r in rows], "parts": parts, "implementation": model.implementation()}
    model.put(work() / "tune" / f"{m}_{draw}" / "chosen.json", json.dumps(out) + "\n")
    return out


def refit(release: str, m: str, draw: str) -> dict:
    """Refit the chosen setting on the links to 2016 under five seeds, once for each way of carrying the
    length of fitting, and keep the mean matrix of each."""
    ch = json.loads((work() / "tune" / f"{m}_{draw}" / "chosen.json").read_text())
    _, M16, parts = matrices()
    if ch["parts"] != json.loads(json.dumps(parts)) or ch["implementation"] != model.implementation():
        raise SystemExit(f"REFUSED: the tuning of {m} {draw} was on other links or by other fitting code; tune it again")
    A, Y = M16["train"], M16["test"]
    edges, n, _ = model.graph_for(release, m, M16)
    s, out = ch["chosen"], {}
    for carry in CARRIES:
        passes = carried(s["passes"], parts["tune_fit"]["links"], parts["refit"]["links"], carry)
        mats, by_seed = [], {}
        t0 = time.time()
        for seed in model.SEEDS:
            S, _ = model.fit(A, s, seed, model.device(), passes, edges, n, balance=DRAWS[draw])
            mats.append(S)
            by_seed[seed] = model.primary(S, A, Y, score.SEEDS)
        S = np.mean(mats, axis=0).astype(np.float32)
        name = arm_name(m, draw, carry)
        (work() / "scores").mkdir(parents=True, exist_ok=True)
        np.save(work() / "scores" / f"{name.replace(':', '_')}.npy", S)
        out[carry] = {"arm": name, "setting": s, "passes_in_tuning": s["passes"], "passes": passes, "by_seed": by_seed, "primary": model.primary(S, A, Y, score.SEEDS), "seconds": round(time.time() - t0, 1)}
        print(f"    {name}: {passes} passes, each seed {min(by_seed.values()):.4f} to {max(by_seed.values()):.4f}, mean of five {out[carry]['primary']:.4f}  ({out[carry]['seconds']} s)", flush=True)
    model.put(work() / "refit" / f"{m}_{draw}.json", json.dumps(out) + "\n")
    return out


def combiners(release: str) -> dict:
    """Fit each combiner's weights at 2014 on the links of 2015 and 2016, and apply them at 2016."""
    M14, M16, parts = matrices()
    scores, cols, out = {}, {}, {}
    for key, M in (("tune", M14), ("refit", M16)):
        pairs, R, _ = model.load_structure(release, M["targets"], M["traits"], ("position", "lineage"))
        scores[key] = combine.raw_scores(M["train"], pairs, R)
        cols[key] = np.flatnonzero(M["test"].any(axis=0))
    A14, Y14, A16, Y16 = M14["train"], M14["test"], M16["train"], M16["test"]
    c14 = cols["tune"]
    (work() / "scores").mkdir(parents=True, exist_ok=True)
    for names, form in itertools.product(SETS, FORMS):
        t0 = time.time()
        X = combine.design(scores["tune"], A14, names, form, c14)
        w, fitted = combine.fit_weights(X, ~A14[:, c14], Y14[:, c14], names.index("B3"))
        Z = combine.combined(scores["refit"], A16, names, form, w)
        name = combiner_name(names, form)
        np.save(work() / "scores" / f"{name.replace(':', '_')}.npy", Z)
        out[name] = {"arm": name, "scores": list(names), "form": form, "terms": combine.terms(names, form), "weights": [float(x) for x in w], "fitted_at_tuning": fitted,
                     "ap_at_tuning": score.by_axis(X @ w, A14[:, c14], Y14[:, c14], axis=1, seeds=(0,))["ap"].tolist(), "primary": model.primary(Z, A16, Y16, score.SEEDS)}
        print(f"    {name}: fitted {fitted:.4f} at 2014, {out[name]['primary']:.4f} on the target  " + " ".join(f"{x:+.3f}" for x in w) + f"  ({time.time() - t0:.0f} s)", flush=True)
    model.put(work() / "combiners.json", json.dumps({"parts": parts, "implementation": combine.implementation(), "arms": out}) + "\n")
    return out


def against(a: np.ndarray, b: np.ndarray) -> dict:
    m, lo, hi = score.interval(a - b)
    return {"difference": m, "low": lo, "high": hi, "half_width": (hi - lo) / 2, "units_higher": int((a > b).sum()), "units": int(len(a))}


def thirds(size: np.ndarray) -> np.ndarray:
    """0, 1 or 2 for each trait: the third of the traits it falls in by its number of fitting links."""
    return np.digitize(size, np.percentile(size, [100 / 3, 200 / 3]))


def report(release: str) -> dict:
    """Every arm on the rehearsal target, against the neighbours baseline computed from the links to 2016."""
    M14, M16, parts = matrices()
    A, Y = M16["train"], M16["test"]
    ap = lambda S: score.by_axis(S, A, Y, axis=1)["ap"]
    base16 = {k: ap(v) for k, v in baselines.scorers(A).items() if k in ("B1", "B3")}
    a14 = lambda S: score.by_axis(S, M14["train"], M14["test"], axis=1, seeds=(0,))["ap"]
    base14 = a14(baselines.scorers(M14["train"])["B3"])
    units = np.flatnonzero(Y.any(axis=0))
    size = A.sum(axis=0)[units]
    third = thirds(size)
    load = lambda name: np.load(work() / "scores" / f"{name.replace(':', '_')}.npy")
    per, arms, missing = {}, {}, []
    for m, draw, carry in model_arms():
        name = arm_name(m, draw, carry)
        f, c = work() / "refit" / f"{m}_{draw}.json", work() / "tune" / f"{m}_{draw}" / "chosen.json"
        if not f.is_file() or not c.is_file():
            missing.append(name)
            continue
        r, ch = json.loads(f.read_text())[carry], json.loads(c.read_text())
        per[name] = ap(load(name))
        arms[name] = {"model": m, "draw": draw, "carry": carry, "setting": r["setting"], "passes_in_tuning": r["passes_in_tuning"], "passes": r["passes"],
                      "by_seed": r["by_seed"], "in_tuning": {"validation": ch["validation"], "against_neighbours": against(np.array(ch["ap"]), base14)}}
    cf = work() / "combiners.json"
    comb = json.loads(cf.read_text())["arms"] if cf.is_file() else {}
    for names, form in itertools.product(SETS, FORMS):
        name = combiner_name(names, form)
        if name not in comb:
            missing.append(name)
            continue
        per[name] = ap(load(name))
        arms[name] = {"scores": list(names), "form": form, "weights": dict(zip(comb[name]["terms"], comb[name]["weights"])),
                      "in_tuning": {"validation": comb[name]["fitted_at_tuning"], "against_neighbours": against(np.array(comb[name]["ap_at_tuning"]), base14),
                                    "read_as": "the weights were fitted on these links, so this figure flatters the arm"}}
    for name, a in arms.items():
        d = against(per[name], base16["B3"])
        a.update({"primary": float(per[name].mean()), "against_neighbours": d, "carries": bool(d["low"] > 0),
                  "against_neighbours_by_third": [float((per[name] - base16["B3"])[third == k].mean()) for k in range(3)]})
    pairs = {}
    for m in MODELS:
        for carry in CARRIES:
            for draw in ("L05", "L1"):
                x, y = arm_name(m, draw, carry), arm_name(m, "L0", carry)
                if x in per and y in per:
                    pairs[f"{x} minus {y}"] = {"what": "draw", **against(per[x], per[y])}
        for draw in DRAWS:
            x, y = arm_name(m, draw, "steps"), arm_name(m, draw, "passes")
            if x in per and y in per:
                pairs[f"{x} minus {y}"] = {"what": "carry", **against(per[x], per[y])}
    for names in SETS:
        x, y = combiner_name(names, "F1"), combiner_name(names, "F3")
        if x in per and y in per:
            pairs[f"{x} minus {y}"] = {"what": "form", **against(per[x], per[y])}
    out = {"cutoff": CUTOFF, "refit_to": REFIT, "tuned_to": TUNE, "release": release, "parts": parts, "traits_on_target": int(len(units)),
           "thirds": [{"traits": int((third == k).sum()), "fitting_links_from": int(size[third == k].min()), "fitting_links_to": int(size[third == k].max()),
                       "share_of_fitting_links": float(size[third == k].sum() / size.sum())} for k in range(3)],
           "baselines_on_target": {k: float(v.mean()) for k, v in base16.items()}, "neighbours_in_tuning": float(base14.mean()),
           "neighbours_by_third": [float(base16["B3"][third == k].mean()) for k in range(3)],
           "arms": arms, "paired": pairs, "missing": missing, "complete": not missing,
           "implementation": {"model": model.implementation(), "combine": combine.implementation()}}
    out_path().write_text(json.dumps(out, indent=1) + "\n")
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("eighteen model arms and ten combiner arms, every name once", len(model_arms()) == 18 and len({arm_name(*a) for a in model_arms()} | {combiner_name(n, f) for n in SETS for f in FORMS}) == 28)
    ck("eight settings for a graph model and four for the one without steps", len(settings("lightgcn")) == 8 == len(settings("lightgcn_lin_pos")) and len(settings("mf")) == 4 and {s["layers"] for s in settings("mf")} == {0})
    ck("carried by passes the length is unchanged; by steps it shrinks with the links, and is never zero",
       carried(110, 26503, 57411, "passes") == 110 and carried(110, 26503, 57411, "steps") == 51 and carried(1, 10, 1000, "steps") == 1)
    ck("the three ways of drawing are by link, by the square root of the trait and by trait", DRAWS == {"L0": 0.0, "L05": 0.5, "L1": 1.0})
    ck("every combiner arm holds the neighbours score, whose weight the fit keeps at 1", all("B3" in s for s in SETS))
    # a toy history of links: years 2013 to 2018
    rows = [("t1", "a", 2013), ("t2", "a", 2014), ("t1", "b", 2014), ("t2", "b", 2015), ("t3", "a", 2016), ("t1", "c", 2016), ("t3", "b", 2017), ("t2", "c", 2018), ("t9", "a", 2018)]
    links = pd.DataFrame(rows, columns=["target", "trait", "year"])
    at14, at16 = periods(links)
    ck("tuning fits to 2014 and validates on 2015 and 2016, between a target and a trait it has seen",
       sorted(at14["train"].year) == [2013, 2014, 2014] and sorted(zip(at14["test"].target, at14["test"].trait)) == [("t2", "b")])
    ck("the refit is on the links to 2016 and the target is 2017 and 2018, between a target and a trait the refit has seen",
       len(at16["train"]) == 6 and sorted(zip(at16["test"].target, at16["test"].trait)) == [("t2", "c"), ("t3", "b")])
    try:
        periods(pd.concat([links, pd.DataFrame([("t1", "a", 2019)], columns=links.columns)]))
        ck("a link first reported after the cutoff is refused", False)
    except SystemExit:
        ck("a link first reported after the cutoff is refused", True)
    ck("traits fall into three groups by size", thirds(np.array([1, 2, 3, 10, 11, 12, 100, 200, 300])).tolist() == [0, 0, 0, 1, 1, 1, 2, 2, 2])
    d = against(np.array([0.3, 0.2, 0.5]), np.array([0.1, 0.2, 0.1]))
    ck("a paired difference counts the units that are higher", d["units_higher"] == 2 and d["units"] == 3 and abs(d["difference"] - 0.2) < 1e-12)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=("tune", "refit", "combiners", "report"))
    ap.add_argument("release")
    ap.add_argument("--model", choices=MODELS)
    ap.add_argument("--draw", choices=list(DRAWS))
    a = ap.parse_args(argv)
    if a.stage in ("tune", "refit") and not (a.model and a.draw):
        raise SystemExit("tune and refit need --model and --draw")
    if a.stage == "tune":
        out = tune(a.release, a.model, a.draw)
        print(f"  {a.model} {a.draw}: chosen {out['chosen']} validation {out['validation']:.4f}")
    elif a.stage == "refit":
        refit(a.release, a.model, a.draw)
    elif a.stage == "combiners":
        combiners(a.release)
    else:
        out = report(a.release)
        print(f"  neighbours on the target {out['baselines_on_target']['B3']:.4f}; {len(out['arms'])} arms" + (f"; missing {out['missing']}" if out["missing"] else ""))
        for name, x in sorted(out["arms"].items(), key=lambda kv: -kv[1]["primary"]):
            d = x["against_neighbours"]
            print(f"  {name:36s} {x['primary']:.4f}  {d['difference']:+.4f} [{d['low']:+.4f}, {d['high']:+.4f}]{'  carries' if x['carries'] else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
