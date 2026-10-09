"""Parameter-free scores combined with fitted weights (docs/PROTOCOL.md, "Lineage and position as inputs").

    C0 combine          B1 target degree and B3 neighbours
    C1 combine_pos      C0 and B5 position: how many of the target's chromosome neighbours are linked to the trait
    C2 combine_lin      C0 and B4 lineage: how many of the target's traits are in the candidate trait's lineage
    C3 combine_lin_pos  C0, B4 and B5

A candidate's score is a weighted sum of terms, one or two for each score. The weights are fitted inside
the training years: the scores are computed from the links of 2016 and earlier, and the weights
maximise the primary measure on the links of 2017 and 2018. The weight on B3 is fixed at 1 and the
rest start at 0, so the search starts at the neighbours baseline and is a search one weight at a time
over a fixed grid, with nothing random. The form in which a score enters (F1, F2 or F3) is chosen by
cross-validation over traits.

The weights were first fitted by a likelihood (commit 86b57d5, `docs/data/combine_gene_trait_2018_first.json`).
That fit put the combiner of B1 and B3 below B3 alone, and the protocol was amended before any of these
was scored on the test links.

All four were added after `mf` and `lightgcn` had been scored on the test links. C0 uses neither
lineage nor position: it shows what fitting two weights does without them.

Stages, in order. `final` scores on the test links, so it runs once, after the file `fit` wrote has been
committed:

    python -m pharos_gap_kg.combine fit pharos319 --cutoff 2018
    python -m pharos_gap_kg.combine freeze pharos319 --cutoff 2018
    python -m pharos_gap_kg.combine final pharos319 --cutoff 2018
    python -m pharos_gap_kg.combine score pharos319 --cutoff 2018
    python -m pharos_gap_kg.combine versus pharos319 --cutoff 2018     # the declared comparisons between scored models
    python -m pharos_gap_kg.combine --selftest
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from pharos_gap_kg import baselines, graph, model, paths, score, split

COMBINERS = {"combine": ("B1", "B3"), "combine_pos": ("B1", "B3", "B5"), "combine_lin": ("B1", "B3", "B4"), "combine_lin_pos": ("B1", "B3", "B4", "B5")}
FULL, BASE = "combine_lin_pos", "combine"
FORMS = ("F1", "F2", "F3")
FOLDS = 5
# the values a weight may take: 0, and plus or minus 2 to the power -6 up to 4
GRID = (0.0,) + tuple(s * 2.0 ** k for k in range(-6, 5) for s in (1.0, -1.0))
MAX_SWEEPS = 10
FITTING = ("raw_scores", "terms", "design", "combined", "tie_lists", "mean_ap", "fit_weights", "cross_validated")
VIEWS = ("full", "low", "neutral", "locus")
# the comparisons between two scored models that the protocol declares: (model, the one it is read against)
VERSUS = (("lightgcn_pos", "lightgcn"), ("lightgcn_lin_pos", "lightgcn_pos"), ("lightgcn_pos", "combine_pos"), ("lightgcn_lin_pos", "combine_lin_pos"))


def implementation(text: str | None = None) -> str:
    """An identity for the code that fits and applies the weights: the source of those functions, hashed."""
    import ast
    text = Path(__file__).read_text() if text is None else text
    found = {n.name: ast.get_source_segment(text, n) for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name in FITTING}
    missing = [f for f in FITTING if f not in found]
    if missing:
        raise SystemExit(f"REFUSED: the fitting code has no function {missing}")
    return hashlib.sha256("\n".join(found[f] for f in FITTING).encode()).hexdigest()[:16]


# ---- the scores and the weights ------------------------------------------------------------------------------

def raw_scores(A: np.ndarray, pairs: np.ndarray, R: np.ndarray) -> dict[str, np.ndarray]:
    """The four parameter-free scores of every (target, trait), from the links in A alone. B4 and B5 are
    computed as the controls compute them."""
    from pharos_gap_kg import position
    b = baselines.scorers(A)
    return {"B1": b["B1"], "B3": b["B3"], "B4": A.astype(np.float32) @ R.astype(np.float32).T, "B5": position.neighbour_count(pairs, A)}


def terms(names: tuple[str, ...], form: str) -> list[str]:
    """What each weight multiplies, in order."""
    if form == "F1":
        return [f"log(1 + {n})" for n in names]
    if form == "F2":
        return [f"percentile of {n}" for n in names]
    if form == "F3":
        return [f"log(1 + {n})" for n in names] + [f"log(1 + {n}) x log(1 + trait links)" for n in names]
    raise ValueError(form)


def design(S: dict[str, np.ndarray], A: np.ndarray, names: tuple[str, ...], form: str, cols: np.ndarray) -> np.ndarray:
    """The terms of `form` for the traits in `cols`, as [targets, len(cols), terms]."""
    if form == "F2":
        from pharos_gap_kg import enrich
        return np.stack([enrich.percentile(S[n][:, cols], ~A[:, cols]) for n in names], axis=2).astype(np.float64)
    X = np.stack([np.log1p(S[n][:, cols]) for n in names], axis=2).astype(np.float64)
    if form == "F1":
        return X
    if form == "F3":
        size = np.log1p(A.sum(axis=0)[cols]).astype(np.float64)
        return np.concatenate([X, X * size[None, :, None]], axis=2)
    raise ValueError(form)


def combined(S: dict[str, np.ndarray], A: np.ndarray, names: tuple[str, ...], form: str, w) -> np.ndarray:
    """The combiner's score for every (target, trait)."""
    out = np.zeros(A.shape, np.float64)
    for k in range(0, A.shape[1], 200):
        cols = np.arange(k, min(k + 200, A.shape[1]))
        out[:, cols] = design(S, A, names, form, cols) @ np.asarray(w, np.float64)
    return out.astype(np.float32)


def tie_lists(cand: np.ndarray, pos: np.ndarray) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """For each trait: its candidates, which of them are positives, and the random numbers that break
    ties, drawn as `score.by_axis` draws them under seed 0."""
    rng = np.random.default_rng(0)
    out = []
    for j in range(cand.shape[1]):
        idx = np.flatnonzero(cand[:, j])
        out.append((idx, pos[idx, j], rng.random(len(idx))))
    return out


def mean_ap(Z: np.ndarray, lists: list[tuple[np.ndarray, np.ndarray, np.ndarray]]) -> float:
    """The primary measure: the mean over traits of average precision among the trait's candidates."""
    total = 0.0
    for j, (idx, y, r) in enumerate(lists):
        ranked = y[np.lexsort((r, -Z[idx, j]))]
        at = np.flatnonzero(ranked)
        total += float((np.cumsum(ranked)[at] / (at + 1)).mean())
    return total / len(lists)


def fit_weights(X: np.ndarray, cand: np.ndarray, pos: np.ndarray, fixed: int) -> tuple[np.ndarray, float]:
    """The weights that maximise `mean_ap`, found one weight at a time over GRID. The term `fixed` keeps
    weight 1 and the others start at 0, so the search starts at that term alone and never ends below it.
    A weight keeps its value unless another value scores strictly higher."""
    if (pos.sum(axis=0) == 0).any() or (pos & ~cand).any():
        raise SystemExit("REFUSED: a trait has no positive, or a positive is not a candidate")
    lists = tie_lists(cand, pos)
    w = np.zeros(X.shape[2])
    w[fixed] = 1.0
    best = mean_ap(X @ w, lists)
    for _ in range(MAX_SWEEPS):
        changed = False
        for i in range(len(w)):
            if i == fixed:
                continue
            rest = X @ np.where(np.arange(len(w)) == i, 0.0, w)
            for g in GRID:
                if g == w[i]:
                    continue
                v = mean_ap(rest + g * X[..., i], lists)
                if v > best:
                    best, changed = v, True
                    w[i] = g
        if not changed:
            break
    return w, best


def cross_validated(X: np.ndarray, cand: np.ndarray, pos: np.ndarray, fixed: int, folds: int = FOLDS) -> tuple[np.ndarray, list[list[int]]]:
    """Each trait scored by weights fitted on the other folds: the out-of-fold score matrix, and the folds."""
    order = np.random.default_rng(0).permutation(X.shape[1])
    parts = [np.sort(p) for p in np.array_split(order, folds)]
    Z = np.zeros(X.shape[:2])
    for held in parts:
        rest = np.setdiff1d(np.arange(X.shape[1]), held)
        w, _ = fit_weights(X[:, rest], cand[:, rest], pos[:, rest], fixed)
        Z[:, held] = X[:, held] @ w
    return Z, [p.tolist() for p in parts]


# ---- the stages ----------------------------------------------------------------------------------------------

def record_path(cutoff: int) -> Path:
    return paths.DOCS_DATA / f"combine_gene_trait_{cutoff}.json"


def frozen_path(cutoff: int) -> Path:
    return paths.DOCS_DATA / f"combine_gene_trait_{cutoff}_frozen.json"


def out_dir(cutoff: int) -> Path:
    return split.out_dir(cutoff) / "models" / "combiners"


def reference_files() -> dict:
    """The digests of the two reference files lineage and position are read from, and the window."""
    from pharos_gap_kg import position, registry
    side = lambda k: json.loads(registry.SOURCES[k].path.with_name(registry.SOURCES[k].filename + ".lineage.json").read_text())["sha256"]
    return {"gtf_sha256": side("ensembl_gtf_113"), "window": position.WINDOW, "efo_sha256": side("efo_obo")}


def against(a: np.ndarray, b: np.ndarray) -> dict:
    m, lo, hi = score.interval(a - b)
    return {"difference": m, "low": lo, "high": hi, "half_width": (hi - lo) / 2, "traits_higher": int((a > b).sum()), "traits": int(len(a))}


def fit_stage(release: str, cutoff: int) -> dict:
    """Choose the form and fit the four combiners' weights on the inner split. No test link is read."""
    import scipy
    M, parts = model.inner_matrices(cutoff)
    A, Y = M["train"], M["test"]
    pairs, R, used = model.load_structure(release, M["targets"], M["traits"], ("position", "lineage"))
    S = raw_scores(A, pairs, R)
    cols = np.flatnonzero(Y.any(axis=0))
    cand, pos = ~A[:, cols], Y[:, cols]
    ap = lambda Z: score.by_axis(Z, A[:, cols], pos, axis=1, seeds=(0,))["ap"]
    base = {k: ap(S[k][:, cols]) for k in ("B1", "B3", "B4", "B5")}
    lin, pon = (S["B4"][:, cols] > 0) & cand, (S["B5"][:, cols] > 0) & cand
    cover = {"traits": int(len(cols)), "candidates": int(cand.sum()), "positives": int(pos.sum()),
             "candidates_lineage_linked": int(lin.sum()), "candidates_position_linked": int(pon.sum()),
             "positives_lineage_linked": int((pos & lin).sum()), "positives_position_linked": int((pos & pon).sum()), "positives_either": int((pos & (lin | pon)).sum()),
             "neighbour_pairs": int(len(pairs) // 2), "lineage_pairs": int(R.sum() // 2)}
    print(f"  inner split: fit {int(A.sum()):,} links, {cover['positives']:,} validation positives over {cover['traits']} traits; "
          f"lineage-linked {cover['positives_lineage_linked']:,}, position-linked {cover['positives_position_linked']:,}, either {cover['positives_either']:,}", flush=True)
    by_form, folds = {}, None
    for f in FORMS:
        t0 = time.time()
        Z, folds = cross_validated(design(S, A, COMBINERS[FULL], f, cols), cand, pos, COMBINERS[FULL].index("B3"))
        by_form[f] = ap(Z)
        print(f"  {FULL} under {f}: cross-validated {by_form[f].mean():.4f}  ({time.time() - t0:.0f} s)", flush=True)
    form = max(FORMS, key=lambda f: float(by_form[f].mean()))          # the first of equals, so the earlier form on a tie
    cv, out = {}, {}
    for name, names in COMBINERS.items():
        X = design(S, A, names, form, cols)
        cv[name] = by_form[form] if name == FULL else ap(cross_validated(X, cand, pos, names.index("B3"))[0])
        w, value = fit_weights(X, cand, pos, names.index("B3"))
        if abs(value - float(ap(X @ w).mean())) > 1e-12:
            raise SystemExit("REFUSED: the measure the weights were fitted to is not the primary measure as `score` computes it")
        out[name] = {"scores": list(names), "terms": terms(names, form), "weights": [float(x) for x in w],
                     "cross_validated": float(cv[name].mean()), "fitted_on_all_traits": value}
    for name in COMBINERS:
        out[name]["cross_validated_minus_B3"] = against(cv[name], base["B3"])
        if name != BASE:
            out[name]["cross_validated_minus_" + BASE] = against(cv[name], cv[BASE])
        d = out[name]["cross_validated_minus_B3"]
        print(f"  {name:16s} {form}: cross-validated {out[name]['cross_validated']:.4f}, against neighbours {d['difference']:+.4f} [{d['low']:+.4f}, {d['high']:+.4f}]  "
              + "  ".join(f"{t} {x:+.3f}" for t, x in zip(out[name]["terms"], out[name]["weights"])), flush=True)
    rec = {"release": release, "cutoff": cutoff, "inner_cutoff": model.INNER, "inner_split": parts, "reference_files": used, "coverage_on_validation": cover,
           "baselines_on_validation": {k: float(v.mean()) for k, v in base.items()}, "folds": FOLDS, "fold_traits": folds,
           "forms": {f: {"cross_validated": float(by_form[f].mean()), "terms": terms(COMBINERS[FULL], f)} for f in FORMS}, "form": form,
           "search": {"fixed": "the weight on B3's own term is 1", "grid": list(GRID), "max_sweeps": MAX_SWEEPS, "start": "every other weight 0"},
           "combiners": out, "implementation": implementation(), "numpy": np.__version__, "scipy": scipy.__version__,
           "read_as": "cross-validated: each trait is scored by weights fitted without it. The form was chosen on these same scores, so the chosen form's figure flatters it a little."}
    model.put(record_path(cutoff), json.dumps(rec, indent=1) + "\n")
    return rec


def preflight(rec: dict) -> list[str]:
    """What stands between a fitted record and a scoring on the test links; empty when nothing does."""
    out = []
    if rec.get("implementation") != implementation():
        out.append("the record was made by other fitting code than this")
    if set(rec.get("combiners", {})) != set(COMBINERS) or any(tuple(c["scores"]) != COMBINERS[k] for k, c in rec.get("combiners", {}).items() if k in COMBINERS):
        out.append("the record does not hold the four declared combiners")
    if any(c.get("weights", [None] * 9)[COMBINERS[k].index("B3")] != 1.0 or not set(c.get("weights", [9])) <= set(GRID)
           for k, c in rec.get("combiners", {}).items() if k in COMBINERS):
        out.append("a weight is not a value of the grid, or the weight on B3 is not 1")
    if rec.get("form") not in FORMS or set(rec.get("forms", {})) != set(FORMS):
        out.append("the record does not say how each declared form scored")
    elif rec["form"] != max(FORMS, key=lambda f: rec["forms"][f]["cross_validated"]):
        out.append("the form is not the best of the three by cross-validation")
    if rec.get("reference_files") != reference_files():
        out.append("the ontology or the gene annotation is not the file the weights were fitted with")
    return out


def freeze(release: str, cutoff: int) -> dict:
    """Record which committed file holds the weights to be scored. Runs in the checkout, where git can
    say the file is committed."""
    import subprocess
    rp = record_path(cutoff)
    if paths.ROOT is None or not (paths.ROOT / ".git").exists():
        raise SystemExit("REFUSED: freeze runs in the git checkout, where a file can be shown to be committed")
    if not rp.is_file():
        raise SystemExit(f"no fitted weights; run:  python -m pharos_gap_kg.combine fit {release} --cutoff {cutoff}")
    rec = json.loads(rp.read_text())
    problems = preflight(rec)
    git = lambda *a: subprocess.run(["git", "-C", str(paths.ROOT), *a], capture_output=True, text=True)
    rel = str(rp.relative_to(paths.ROOT))
    if git("ls-files", "--error-unmatch", rel).returncode != 0:
        problems.append(f"{rel} is not tracked by git")
    elif git("diff", "--quiet", "HEAD", "--", rel).returncode != 0:
        problems.append(f"{rel} differs from the committed copy")
    if problems:
        raise SystemExit("REFUSED: the combiners cannot be frozen: " + "; ".join(problems))
    out = {"cutoff": cutoff, "form": rec["form"], "weights": {k: c["weights"] for k, c in rec["combiners"].items()}, "implementation": rec["implementation"],
           "record_file": rel, "record_sha256": graph.file_sha(rp), "commit": git("log", "-1", "--format=%H", "--", rel).stdout.strip()}
    model.put(frozen_path(cutoff), json.dumps(out, indent=1) + "\n")
    return out


def final(release: str, cutoff: int) -> dict:
    """The four combiners' matrices on all training links, each scored on the test links. Once."""
    rp, d, fz = record_path(cutoff), out_dir(cutoff), frozen_path(cutoff)
    if not rp.is_file():
        raise SystemExit(f"no fitted weights; run:  python -m pharos_gap_kg.combine fit {release} --cutoff {cutoff}  and commit {rp.name}")
    if (d / "final.json").is_file():
        raise SystemExit(f"REFUSED: the combiners have been scored on the test links ({d / 'final.json'}); the protocol scores them once")
    rec = json.loads(rp.read_text())
    problems = preflight(rec)
    if not fz.is_file():
        problems.append(f"the weights are not frozen; commit {rp.name}, then run:  python -m pharos_gap_kg.combine freeze {release} --cutoff {cutoff}")
    elif json.loads(fz.read_text())["record_sha256"] != graph.file_sha(rp):
        problems.append(f"{rp.name} is not the file that was frozen ({fz.name})")
    if problems:
        raise SystemExit("REFUSED: the combiners are not ready to be scored on the test links: " + "; ".join(problems))
    model.claim(d, {"models": list(COMBINERS), "cutoff": cutoff, "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "record_sha256": graph.file_sha(rp),
                    "frozen_at": json.loads(fz.read_text())["commit"], "implementation": implementation()})
    M = baselines.matrices(split.load(cutoff))
    A, Y = M["train"], M["test"]
    pairs, R, _ = model.load_structure(release, M["targets"], M["traits"], ("position", "lineage"))
    S = raw_scores(A, pairs, R)
    out = {"release": release, "cutoff": cutoff, "form": rec["form"], "record_sha256": graph.file_sha(rp), "models": {}}
    for name, c in rec["combiners"].items():
        Z = combined(S, A, COMBINERS[name], rec["form"], c["weights"])
        dd = model.out_dir(cutoff, name)
        dd.mkdir(parents=True, exist_ok=True)
        np.save(dd / "scores.npy", Z)
        one = {"model": name, "release": release, "cutoff": cutoff, "form": rec["form"], "terms": c["terms"], "weights": c["weights"],
               "primary_of_the_matrix": model.primary(Z, A, Y, score.SEEDS), "scores_sha256": graph.file_sha(dd / "scores.npy"), "record_sha256": graph.file_sha(rp)}
        model.put(dd / "final.json", json.dumps(one, indent=1) + "\n")
        out["models"][name] = {k: one[k] for k in ("primary_of_the_matrix", "scores_sha256")}
        print(f"    {name}: primary {one['primary_of_the_matrix']:.4f}", flush=True)
    model.put(d / "final.json", json.dumps(out, indent=1) + "\n")
    return out


def load_scored(cutoff: int, name: str) -> np.ndarray:
    """A scored model's matrix, checked to be the one its scoring wrote."""
    d = model.out_dir(cutoff, name)
    if not (d / "final.json").is_file() or not (d / "scores.npy").is_file():
        raise FileNotFoundError(name)
    if graph.file_sha(d / "scores.npy") != json.loads((d / "final.json").read_text())["scores_sha256"]:
        raise SystemExit(f"REFUSED: {d / 'scores.npy'} is not the matrix its final.json describes")
    return np.load(d / "scores.npy")


def views(release: str, cutoff: int, mats: dict[str, np.ndarray]) -> tuple[dict, dict]:
    """Per-trait average precision of each matrix on the four views of the Controls, beside B1 and B3
    computed here the same way and checked against the committed controls. Returns (per trait, the split's matrices)."""
    from pharos_gap_kg import controls, position
    M = baselines.matrices(split.load(cutoff))
    A, Y = M["train"], M["test"]
    pairs, R, _ = model.load_structure(release, M["targets"], M["traits"], ("position", "lineage"))
    linked = (((A.astype(np.float32) @ R.astype(np.float32).T) > 0) | (position.neighbour_count(pairs, A) > 0)) & ~A
    new = Y & ~linked
    per = {}
    for k, s in {**{b: v for b, v in baselines.scorers(A).items() if b in ("B1", "B3")}, **mats}.items():
        per[k] = {"full": score.by_axis(s, A, Y, axis=1)["ap"], "low": score.by_axis(s, A, new, axis=1, exclude=Y & linked)["ap"],
                  "neutral": score.by_axis(s, A, new, axis=1, exclude=linked)["ap"], "locus": controls.collapsed(s, A, Y, pairs)["ap"]}
    stored = json.loads((paths.DOCS_DATA / f"controls_gene_trait_{cutoff}.json").read_text())["results"]
    for k in ("B1", "B3"):
        was = [stored[k][v][0] for v in ("primary", "new_links_lower_bound", "new_links_neutral", "locus_collapsed")]
        now = [float(np.nanmean(per[k][v])) for v in VIEWS]
        if not np.allclose(was, now, atol=1e-9):
            raise SystemExit(f"REFUSED: {k} recomputed here is {now}, and the committed controls say {was}; the models would be compared with a different baseline")
    return per, M


def paired(per: dict, a: str, b: str) -> dict:
    out = {}
    for v in VIEWS:
        m, lo, hi = score.interval(per[a][v] - per[b][v])
        out[v] = {"difference": m, "low": lo, "high": hi, "half_width": (hi - lo) / 2}
    return out


def scored(release: str, cutoff: int) -> dict:
    """Every measure and control for the four combiners, and the paired differences the protocol names."""
    d = out_dir(cutoff)
    fin = json.loads((d / "final.json").read_text())
    mats = {name: load_scored(cutoff, name) for name in COMBINERS}
    per, M = views(release, cutoff, mats)
    res = {k: {v: score.interval(a) for v, a in per[k].items()} for k in per}
    pd_ = {f"{name}-{other}": paired(per, name, other) for name in COMBINERS for other in ("B3", "B1")}
    pd_.update({f"{name}-{BASE}": paired(per, name, BASE) for name in COMBINERS if name != BASE})
    years = sorted(int(y) for y in np.unique(M["year"][M["test"]]) if (M["year"][M["test"]] == y).sum() >= 100)
    second = {}
    for name in COMBINERS:
        full = baselines.evaluate(mats[name], M, years)
        second[name] = {**{k: full[k] for k in ("trait_auc", "trait_recall_at_100", "target_ap", "target_auc", "overall_ap_by_seed")},
                        "by_year": {str(y): v for y, v in full["by_year"].items()}}
    rec = json.loads(record_path(cutoff).read_text())
    out = {"cutoff": cutoff, "form": fin["form"], "weights": {k: dict(zip(c["terms"], c["weights"])) for k, c in rec["combiners"].items()},
           "results": res, "paired": pd_, "secondary": second, "scores_sha256": {k: v["scores_sha256"] for k, v in fin["models"].items()},
           "scoring": {**json.loads((d / "final_started.json").read_text()), "primary_of_the_matrix": {k: v["primary_of_the_matrix"] for k, v in fin["models"].items()}}}
    (paths.DOCS_DATA / f"combine_gene_trait_{cutoff}_scored.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def versus(release: str, cutoff: int) -> dict:
    """The declared comparisons between two scored models, for every pair where both have been scored."""
    have, waiting = {}, []
    for a, b in VERSUS:
        try:
            for k in (a, b):
                if k not in have:
                    have[k] = load_scored(cutoff, k)
        except FileNotFoundError as e:
            waiting.append({"comparison": f"{a}-{b}", "not_scored": str(e)})
    per, _ = views(release, cutoff, have) if have else ({}, None)
    done = {f"{a}-{b}": paired(per, a, b) for a, b in VERSUS if a in have and b in have}
    out = {"cutoff": cutoff, "declared": [f"{a}-{b}" for a, b in VERSUS], "paired": done, "waiting": waiting,
           "results": {k: {v: score.interval(per[k][v]) for v in VIEWS} for k in have}}
    (paths.DOCS_DATA / f"versus_gene_trait_{cutoff}.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    here = Path(__file__).read_text()
    ck("the fitting code has one identity, and it changes when a fitting function changes and not when another does",
       implementation() == implementation(here) and implementation(here.replace("if v > best:", "if v >= best:")) != implementation()
       and implementation(here.replace("def record_path(", "def record_path_(")) == implementation())
    ck("C0 uses neither input, and each other combiner adds lineage, position or both to it",
       COMBINERS[BASE] == ("B1", "B3") and all(set(COMBINERS[BASE]) < set(v) for k, v in COMBINERS.items() if k != BASE) and set(COMBINERS[FULL]) == {"B1", "B3", "B4", "B5"})
    ck("the third form has two terms for each score and the others one", [len(terms(COMBINERS[FULL], f)) for f in FORMS] == [4, 4, 8])

    ck("the grid holds zero and both signs of every power of two from -6 to 4", len(GRID) == 23 and GRID[0] == 0.0 and min(GRID) == -16.0 and 2.0 ** -6 in GRID and 1.0 in GRID)

    # 40 traits, 300 targets: positives are the top of a sum of two terms; a third term is noise
    rng = np.random.default_rng(1)
    T, U = 300, 40
    X = rng.normal(size=(T, U, 3))
    cand = rng.random((T, U)) > 0.1
    z = np.where(cand, X[..., 0] + 2.0 * X[..., 1] + rng.normal(scale=0.5, size=(T, U)), -np.inf)
    pos = z >= np.sort(z, axis=0)[-30]
    lists = tie_lists(cand, pos)
    flat = np.round(X[..., 0])          # many ties, so the tie-breaking is exercised
    ck("the measure the weights are fitted to is the primary measure as `score` computes it, ties and all",
       abs(mean_ap(flat, lists) - float(score.by_axis(flat, ~cand, pos, axis=1, seeds=(0,))["ap"].mean())) < 1e-12)
    w, value = fit_weights(X, cand, pos, 0)
    alone = mean_ap(X[..., 0], lists)
    ck(f"the search starts at the fixed term alone and ends above it ({alone:.3f} -> {value:.3f}, weights {w.tolist()})", w[0] == 1.0 and value > alone + 0.1 and abs(value - mean_ap(X @ w, lists)) < 1e-12)
    ck("the term that helps gets its weight, the noise term gets little, and every weight is a value of the grid", w[1] in (1.0, 2.0, 4.0) and abs(w[2]) <= 0.25 and set(w.tolist()) <= set(GRID))
    w2, _ = fit_weights(X.copy(), cand.copy(), pos.copy(), 0)
    ck("fitting twice gives the same weights", np.array_equal(w, w2))
    wn, vn = fit_weights(X[..., [0, 2]], cand, pos, 0)
    ck("a term that does not help cannot pull the fit below the fixed term alone", vn >= alone)
    Z, folds = cross_validated(X, cand, pos, 0)
    ck("every trait is held out once, and is scored by weights fitted without it",
       sorted(i for f in folds for i in f) == list(range(U)) and len(folds) == FOLDS
       and np.allclose(Z[:, folds[0]], X[:, folds[0]] @ fit_weights(*(a[:, np.setdiff1d(np.arange(U), folds[0])] for a in (X, cand, pos)), 0)[0]))
    try:
        fit_weights(X, cand, pos & False, 0)
        ck("a trait with no positive is refused", False)
    except SystemExit:
        ck("a trait with no positive is refused", True)

    # the scores, on a toy: 4 targets, 3 traits; targets 0 and 1 are neighbours; traits 0 and 2 are in one lineage
    A = np.array([[1, 0, 0], [0, 1, 0], [0, 1, 1], [1, 1, 0]], bool)
    pairs = np.array([[0, 1], [1, 0]])
    R = np.zeros((3, 3), bool)
    R[0, 2] = R[2, 0] = True
    S = raw_scores(A, pairs, R)
    ck("B4 counts the target's traits in the candidate trait's lineage, and B5 its neighbours linked to the trait",
       S["B4"].tolist() == [[0, 0, 1], [0, 0, 0], [1, 0, 0], [0, 0, 1]] and S["B5"].tolist() == [[0, 1, 0], [1, 0, 0], [0, 0, 0], [0, 0, 0]])
    cols = np.arange(3)
    X1, X3 = design(S, A, COMBINERS[FULL], "F1", cols), design(S, A, COMBINERS[FULL], "F3", cols)
    ck("the first form is log(1 + score), and the third repeats it times log(1 + the trait's links)",
       X1.shape == (4, 3, 4) and np.allclose(X1[..., 0], np.log1p(S["B1"])) and X3.shape == (4, 3, 8)
       and np.allclose(X3[..., 4:], X1 * np.log1p(A.sum(axis=0))[None, :, None]))
    X2 = design(S, A, COMBINERS[FULL], "F2", cols)
    ck("the second form is a rank among the trait's candidates, and tied scores share it", X2.shape == (4, 3, 4) and X2[1, 0, 2] == X2[2, 0, 3] and 0 < X2[:, :, 0][~A].min() <= X2.max() <= 1)
    wts = np.array([0.5, 1.0, 2.0, 3.0])
    C = combined(S, A, COMBINERS[FULL], "F1", wts)
    ck("a combiner's matrix is the weighted sum of its terms", np.allclose(C, X1 @ wts, atol=1e-6))
    unlinked = (S["B4"] == 0) & (S["B5"] == 0)
    ck("where a candidate is neither lineage- nor position-linked, the two inputs add nothing to its score",
       np.allclose(C[unlinked], combined(S, A, COMBINERS[BASE], "F1", wts[:2])[unlinked], atol=1e-6))
    one = combined(S, A, ("B3",), "F1", [2.0])
    ck("with one score the combiner ranks as that score does", all((np.argsort(-one[:, j], kind="stable") == np.argsort(-S["B3"][:, j], kind="stable")).all() for j in range(3)))

    good = {"implementation": implementation(), "form": "F3", "forms": {"F1": {"cross_validated": 0.1}, "F2": {"cross_validated": 0.05}, "F3": {"cross_validated": 0.2}},
            "combiners": {k: {"scores": list(v), "weights": [0.25, 1.0] + [0.0] * (2 * len(v) - 2)} for k, v in COMBINERS.items()}}
    ref = reference_files
    try:
        globals()["reference_files"] = lambda: {"x": 1}
        good["reference_files"] = {"x": 1}
        ck("a record of the four combiners, the best form and this code's fitting is ready to be scored", preflight(good) == [])
        ck("a record made by other code, with another form than the best, or with another ontology is each refused",
           len(preflight({**good, "implementation": "0" * 16})) == 1 and preflight({**good, "form": "F1"}) == ["the form is not the best of the three by cross-validation"]
           and len(preflight({**good, "reference_files": {"x": 2}})) == 1)
        ck("a record that lacks a combiner is refused", len(preflight({**good, "combiners": {k: v for k, v in good["combiners"].items() if k != BASE}})) == 1)
        off = {**good, "combiners": {**good["combiners"], BASE: {"scores": ["B1", "B3"], "weights": [0.3, 1.0, 0.0, 0.0]}}}
        ck("a weight that is not a value of the grid is refused", preflight(off) == ["a weight is not a value of the grid, or the weight on B3 is not 1"])
    finally:
        globals()["reference_files"] = ref
    ck("every declared comparison is between a graph model and the model it is read against", all(a in model.MODELS and (b in model.MODELS or b in COMBINERS) for a, b in VERSUS))
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=("fit", "freeze", "final", "score", "versus"))
    ap.add_argument("release")
    ap.add_argument("--cutoff", type=int, required=True)
    a = ap.parse_args(argv)
    if a.stage == "fit":
        rec = fit_stage(a.release, a.cutoff)
        print(f"  form {rec['form']} of {[(f, round(v['cross_validated'], 4)) for f, v in rec['forms'].items()]}; neighbours on validation {rec['baselines_on_validation']['B3']:.4f}  -> {record_path(a.cutoff)}")
    elif a.stage == "freeze":
        out = freeze(a.release, a.cutoff)
        print(f"  the combiners, form {out['form']}, frozen from {out['record_file']} at {out['commit'][:7]}")
    elif a.stage == "final":
        final(a.release, a.cutoff)
    elif a.stage == "score":
        out = scored(a.release, a.cutoff)
        for k, v in out["paired"].items():
            print(f"  {k}: " + "  ".join(f"{x} {y['difference']:+.4f} [{y['low']:+.4f}, {y['high']:+.4f}]" for x, y in v.items()))
    else:
        out = versus(a.release, a.cutoff)
        for k, v in out["paired"].items():
            print(f"  {k}: " + "  ".join(f"{x} {y['difference']:+.4f} [{y['low']:+.4f}, {y['high']:+.4f}]" for x, y in v.items()))
        for w in out["waiting"]:
            print(f"  {w['comparison']}: waiting, {w['not_scored']} has not been scored")
    return 0


if __name__ == "__main__":
    sys.exit(main())
