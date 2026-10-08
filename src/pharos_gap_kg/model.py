"""Fitted scorers on the frozen split, tuned inside the training years (docs/PROTOCOL.md, "Fitted models").

    M1 mf                   two tables of vectors and a number per target and per trait
    M2 lightgcn             the training links as a graph, vectors smoothed over one to three steps
    M3 lightgcn_ctx         M2 with the undated context joined on
    M4 lightgcn_ctx_pruned  M3 without the disease pairs the pruning rules withheld

Each scores a (target, trait) pair by the dot product of two learned vectors and is trained to rank a
linked target above one drawn at random for the same trait. Tuning never sees a test link: the
training links are split again at 2016 and every setting is chosen on the links of 2017 and 2018.

IN THE CONTEXT GRAPH (M3, M4) a node is a training target, a training trait, or a pathway, GO term,
compound or disease that at least two training targets share. A disease with the same identifier as a
trait is that trait's node. Interaction partners join two targets directly. Proteins outside the
training targets are left out.

Stages, in order. `final` scores on the test links, so it runs once per model, after the setting
`tune` chose has been committed:

    python -m pharos_gap_kg.model tune  pharos319 --cutoff 2018 --model lightgcn
    python -m pharos_gap_kg.model final pharos319 --cutoff 2018 --model lightgcn
    python -m pharos_gap_kg.model score pharos319 --cutoff 2018 --model lightgcn
    python -m pharos_gap_kg.model --selftest

Needs the [model] extra (torch). A GPU is used when there is one.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from pharos_gap_kg import baselines, graph, paths, score, split

INNER = 2016
SEEDS = (0, 1, 2, 3, 4)
GRID = {"dim": (64, 128), "layers": (1, 2, 3), "lr": (0.001, 0.005), "l2": (1e-6, 1e-4)}
BATCH, MAX_PASSES, EVERY, PATIENCE = 4096, 300, 10, 50
MODELS = {"mf": {"graph": False, "context": None}, "lightgcn": {"graph": True, "context": None},
          "lightgcn_ctx": {"graph": True, "context": "disease_links"}, "lightgcn_ctx_pruned": {"graph": True, "context": "disease_links_pruned"}}
MEMBERSHIPS = ("pathways", "go_terms", "compounds")


def settings(model: str) -> list[dict]:
    """Every setting the protocol says is tried for this model, in a fixed order."""
    layers = GRID["layers"] if MODELS[model]["graph"] else (0,)
    return [{"dim": d, "layers": k, "lr": lr, "l2": l2} for d, k, lr, l2 in itertools.product(GRID["dim"], layers, GRID["lr"], GRID["l2"])]


def out_dir(cutoff: int, model: str) -> Path:
    return split.out_dir(cutoff) / "models" / model


def chosen_path(cutoff: int, model: str) -> Path:
    return paths.DOCS_DATA / f"model_{model}_{cutoff}.json"


def device() -> str:
    import torch
    return "cuda" if torch.cuda.is_available() else "cpu"


# ---- the graph a model is fitted on ------------------------------------------------------------------------

def link_edges(A: np.ndarray) -> np.ndarray:
    """The links as edges between target nodes [0, T) and trait nodes [T, T + D)."""
    t, d = np.nonzero(A)
    return np.stack([t, d + A.shape[0]], axis=1)


def context_edges(tables: dict[str, pd.DataFrame], targets: list[str], traits: list[str]) -> tuple[np.ndarray, int, dict]:
    """Edges that join context onto the target and trait nodes, the number of nodes in all, and counts.

    `tables`: partners (a, b), disease (target, disease), and memberships (target, grp) by name."""
    ti = {t: i for i, t in enumerate(targets)}
    di = {d: len(targets) + i for i, d in enumerate(traits)}
    n = len(targets) + len(traits)
    edges, counts = [], {}
    p = tables["partners"]
    p = p[p.a.isin(ti) & p.b.isin(ti)]
    edges.append(np.stack([p.a.map(ti).to_numpy(), p.b.map(ti).to_numpy()], axis=1))
    counts["partner_edges"] = int(len(p))
    groups = [(name, tables[name].rename(columns={"grp": "g"})[["target", "g"]]) for name in MEMBERSHIPS if name in tables]
    groups.append(("disease", tables["disease"].rename(columns={"disease": "g"})[["target", "g"]]))
    for name, m in groups:
        m = m[m.target.isin(ti)].drop_duplicates()
        m = m[m.groupby("g").target.transform("size") >= 2]
        same = m.g.isin(di) if name == "disease" else pd.Series(False, index=m.index)
        ids = {g: n + i for i, g in enumerate(sorted(set(m.g[~same])))}
        node = np.where(same, m.g.map(di), m.g.map(ids)).astype(np.int64) if len(m) else np.zeros(0, np.int64)
        edges.append(np.stack([m.target.map(ti).to_numpy(), node], axis=1) if len(m) else np.zeros((0, 2), np.int64))
        counts[f"{name}_edges"], counts[f"{name}_nodes"] = int(len(m)), len(ids)
        if name == "disease":
            counts["disease_edges_onto_a_trait_node"] = int(same.sum())
        n += len(ids)
    return np.concatenate(edges).astype(np.int64), n, counts


def load_context(release: str, which: str) -> dict[str, pd.DataFrame]:
    C = paths.context(release)
    out = {"partners": pd.read_parquet(C / "partners.parquet", columns=["a", "b"]), "disease": pd.read_parquet(C / f"{which}.parquet", columns=["target", "disease"])}
    out.update({name: pd.read_parquet(C / f"{name}.parquet") for name in MEMBERSHIPS})
    return out


def adjacency(n: int, edges: np.ndarray, dev: str):
    """The symmetric, degree-normalised adjacency of an undirected graph, with no self loops."""
    import torch
    e = edges[edges[:, 0] != edges[:, 1]]
    e = np.unique(np.concatenate([e, e[:, ::-1]]), axis=0)
    deg = np.bincount(e[:, 0], minlength=n).astype(np.float64)
    w = 1.0 / np.sqrt(deg[e[:, 0]] * deg[e[:, 1]])
    return torch.sparse_coo_tensor(torch.as_tensor(e.T.copy()), torch.as_tensor(w, dtype=torch.float32), (n, n), check_invariants=False).coalesce().to(dev)


# ---- the model ---------------------------------------------------------------------------------------------

def build(n_nodes: int, T: int, D: int, setting: dict, adj, seed: int, dev: str):
    import torch

    class Net(torch.nn.Module):
        def __init__(self) -> None:
            super().__init__()
            g = torch.Generator().manual_seed(seed)
            self.vec = torch.nn.Parameter(torch.randn(n_nodes, setting["dim"], generator=g) * 0.1)
            self.bias = torch.nn.Parameter(torch.zeros(T + D)) if adj is None else None      # M1 only

        def vectors(self):
            if adj is None:
                return self.vec
            layers = [self.vec]
            for _ in range(setting["layers"]):
                layers.append(torch.sparse.mm(adj, layers[-1]))
            return torch.stack(layers).mean(dim=0)

        def pair(self, v, t, d):
            s = (v[t] * v[T + d]).sum(dim=1)
            return s if self.bias is None else s + self.bias[t] + self.bias[T + d]

        def matrix(self) -> np.ndarray:
            with torch.no_grad():
                v = self.vectors()
                s = v[:T] @ v[T:T + D].T
                if self.bias is not None:
                    s = s + self.bias[:T, None] + self.bias[None, T:T + D]
                return s.float().cpu().numpy()

    return Net().to(dev)


def draw_negatives(A, d, rng, T: int):
    """For each trait in `d`, a target with no link to it, uniformly."""
    import torch
    neg = torch.randint(0, T, d.shape, generator=rng, device=d.device)
    for _ in range(20):
        hit = A[neg, d]
        if not bool(hit.any()):
            break
        neg[hit] = torch.randint(0, T, (int(hit.sum()),), generator=rng, device=d.device)
    return neg


def fit(A: np.ndarray, setting: dict, seed: int, dev: str, passes: int, edges: np.ndarray | None = None, n_nodes: int | None = None,
        validate=None, log=None) -> tuple[np.ndarray, list[dict]]:
    """Fit on the links in `A` for `passes` passes. With `validate`, a function of the score matrix, it is
    called every EVERY passes, fitting stops PATIENCE passes after its best value, and the matrix
    returned is the best one's. Returns (score matrix, history)."""
    import torch
    T, D = A.shape
    torch.manual_seed(seed)
    adj = None if edges is None else adjacency(n_nodes, edges, dev)
    net = build(n_nodes or T + D, T, D, setting, adj, seed, dev)
    opt = torch.optim.Adam(net.parameters(), lr=setting["lr"])
    At = torch.as_tensor(A, device=dev)
    pos_t, pos_d = (torch.as_tensor(x, device=dev) for x in np.nonzero(A))
    rng = torch.Generator(device=dev).manual_seed(seed)
    history, best, best_S = [], (-1.0, 0), None
    for p in range(1, passes + 1):
        order = torch.randperm(len(pos_t), generator=rng, device=dev)
        total = 0.0
        for k in range(0, len(order), BATCH):
            b = order[k:k + BATCH]
            t, d = pos_t[b], pos_d[b]
            neg = draw_negatives(At, d, rng, T)
            v = net.vectors()
            gap = net.pair(v, t, d) - net.pair(v, neg, d)
            reg = (net.vec[t].pow(2).sum() + net.vec[T + d].pow(2).sum() + net.vec[neg].pow(2).sum()) / len(b)
            loss = torch.nn.functional.softplus(-gap).mean() + setting["l2"] * reg
            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(loss.detach()) * len(b)
        if validate is not None and (p % EVERY == 0 or p == passes):
            S = net.matrix()
            v = float(validate(S))
            history.append({"passes": p, "loss": total / len(order), "validation": v})
            if log:
                log(f"      pass {p:3d}  loss {total / len(order):.4f}  validation {v:.4f}")
            if v > best[0]:
                best, best_S = (v, p), S
            if p - best[1] >= PATIENCE:
                break
    return (net.matrix() if best_S is None else best_S), history


# ---- the stages --------------------------------------------------------------------------------------------

def graph_for(release: str, model: str, M: dict) -> tuple[np.ndarray | None, int | None, dict]:
    """The edges and node count a model is fitted on for the links in M["train"], and what was joined."""
    A = M["train"]
    if not MODELS[model]["graph"]:
        return None, None, {}
    e = link_edges(A)
    if MODELS[model]["context"] is None:
        return e, A.shape[0] + A.shape[1], {"link_edges": int(len(e))}
    ce, n, counts = context_edges(load_context(release, MODELS[model]["context"]), M["targets"], M["traits"])
    return np.concatenate([e, ce]), n, {"link_edges": int(len(e)), **counts, "nodes": n}


def inner_matrices(cutoff: int) -> tuple[dict, dict]:
    tr = split.load(cutoff)["train"]
    inner = split.split(tr, INNER)
    return baselines.matrices({"train": inner["train"], "test": inner["test"]}), {k: {"links": len(v), "sha256": split.digest(v)} for k, v in inner.items()}


def primary(S: np.ndarray, A: np.ndarray, Y: np.ndarray, seeds: tuple[int, ...] = (0,)) -> float:
    return float(np.nanmean(score.by_axis(S, A, Y, axis=1, seeds=seeds)["ap"]))


def tune(release: str, cutoff: int, model: str, only: int | None = None) -> dict:
    import torch
    M, parts = inner_matrices(cutoff)
    A, Y = M["train"], M["test"]
    dev = device()
    edges, n_nodes, joined = graph_for(release, model, M)
    base = {k: primary(v, A, Y) for k, v in baselines.scorers(A).items() if k in ("B0", "B1", "B3")}
    print(f"  {model} on {dev}: fit {int(A.sum()):,} links, validate {int(Y.sum()):,}; baselines on validation {base}; joined {joined}", flush=True)
    rows = []
    for i, s in enumerate(settings(model)[:only]):
        t0 = time.time()
        _, hist = fit(A, s, 0, dev, MAX_PASSES, edges, n_nodes, validate=lambda S: primary(S, A, Y), log=None)
        best = max(hist, key=lambda h: h["validation"])
        rows.append({**s, "validation": best["validation"], "passes": best["passes"], "passes_run": hist[-1]["passes"], "seconds": round(time.time() - t0, 1)})
        print(f"    {i + 1:2d}/{len(settings(model))}  {s}  validation {best['validation']:.4f} at {best['passes']} passes  ({rows[-1]['seconds']} s)", flush=True)
    top = max(rows, key=lambda r: r["validation"])          # the first of equals, in grid order
    out = {"model": model, "release": release, "cutoff": cutoff, "inner_cutoff": INNER, "inner_split": parts, "baselines_on_validation": base, "joined": joined,
           "grid": {k: list(v) for k, v in GRID.items()}, "tried": rows, "complete": only is None,
           "chosen": {k: top[k] for k in ("dim", "layers", "lr", "l2", "passes")}, "validation": top["validation"],
           "torch": torch.__version__, "device": dev}
    d = out_dir(cutoff, model)
    d.mkdir(parents=True, exist_ok=True)
    (d / "tune.json").write_text(json.dumps(out, indent=1) + "\n")
    if only is None:
        chosen_path(cutoff, model).write_text(json.dumps(out, indent=1) + "\n")
    return out


def final(release: str, cutoff: int, model: str) -> dict:
    import torch
    cp, d = chosen_path(cutoff, model), out_dir(cutoff, model)
    if not cp.is_file():
        raise SystemExit(f"no committed setting for {model}; run:  python -m pharos_gap_kg.model tune {release} --cutoff {cutoff} --model {model}  and commit {cp.name}")
    if (d / "final.json").is_file():
        raise SystemExit(f"REFUSED: {model} has been scored on the test links ({d / 'final.json'}); the protocol scores a model once")
    ch = json.loads(cp.read_text())
    if not ch.get("complete"):
        raise SystemExit(f"REFUSED: {cp.name} is from a partial run of the grid")
    M = baselines.matrices(split.load(cutoff))
    A, Y = M["train"], M["test"]
    dev = device()
    edges, n_nodes, joined = graph_for(release, model, M)
    s = ch["chosen"]
    mats, by_seed = [], {}
    for seed in SEEDS:
        S, _ = fit(A, s, seed, dev, s["passes"], edges, n_nodes)
        mats.append(S)
        by_seed[seed] = primary(S, A, Y, score.SEEDS)
        print(f"    seed {seed}: primary {by_seed[seed]:.4f}", flush=True)
    S = np.mean(mats, axis=0).astype(np.float32)
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "scores.npy", S)
    out = {"model": model, "release": release, "cutoff": cutoff, "setting": s, "joined": joined, "primary_by_seed": by_seed,
           "primary_of_the_mean_matrix": primary(S, A, Y, score.SEEDS), "scores_sha256": graph.file_sha(d / "scores.npy"),
           "chosen_sha256": graph.file_sha(cp), "torch": torch.__version__, "device": dev}
    (d / "final.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def scored(release: str, cutoff: int, model: str) -> dict:
    """Every measure and control for the model's matrix, beside the baselines computed the same way."""
    from pharos_gap_kg import controls, ontology, position, registry
    d = out_dir(cutoff, model)
    fin = json.loads((d / "final.json").read_text())
    if graph.file_sha(d / "scores.npy") != fin["scores_sha256"]:
        raise SystemExit(f"REFUSED: {d / 'scores.npy'} is not the matrix final.json describes")
    M = baselines.matrices(split.load(cutoff))
    A, Y = M["train"], M["test"]
    S = {"model": np.load(d / "scores.npy"), **{k: v for k, v in baselines.scorers(A).items() if k in ("B1", "B3")}}
    terms = ontology.read(registry.require("efo_obo"))
    R = controls.lineage_matrix(M["traits"], terms)
    tg = position.target_genes(release, position.genes(registry.require("ensembl_gtf_113")))
    pairs = position.neighbour_pairs(tg, {t: i for i, t in enumerate(M["targets"])})
    linked = (((A.astype(np.float32) @ R.astype(np.float32).T) > 0) | (position.neighbour_count(pairs, A) > 0)) & ~A
    new = Y & ~linked
    per, res = {}, {}
    for k, s in S.items():
        per[k] = {"full": score.by_axis(s, A, Y, axis=1)["ap"], "low": score.by_axis(s, A, new, axis=1, exclude=Y & linked)["ap"],
                  "neutral": score.by_axis(s, A, new, axis=1, exclude=linked)["ap"], "locus": controls.collapsed(s, A, Y, pairs)["ap"]}
        res[k] = {v: score.interval(a) for v, a in per[k].items()}
    stored = json.loads((paths.DOCS_DATA / f"controls_gene_trait_{cutoff}.json").read_text())["results"]
    for k, names in (("B1", "B1"), ("B3", "B3")):
        was = [stored[names][v][0] for v in ("primary", "new_links_lower_bound", "new_links_neutral", "locus_collapsed")]
        now = [res[k][v][0] for v in ("full", "low", "neutral", "locus")]
        if not np.allclose(was, now, atol=1e-9):
            raise SystemExit(f"REFUSED: {k} recomputed here is {now}, and the committed controls say {was}; the model would be compared with a different baseline")
    paired = {}
    for other in ("B3", "B1"):
        paired[f"model-{other}"] = {}
        for v in ("full", "low", "neutral", "locus"):
            m, lo, hi = score.interval(per["model"][v] - per[other][v])
            paired[f"model-{other}"][v] = {"difference": m, "low": lo, "high": hi, "half_width": (hi - lo) / 2}
    years = sorted(int(y) for y in np.unique(M["year"][Y]) if (M["year"][Y] == y).sum() >= 100)
    full = baselines.evaluate(S["model"], M, years)
    out = {"model": model, "cutoff": cutoff, "setting": fin["setting"], "primary_by_seed": fin["primary_by_seed"], "results": res, "paired": paired,
           "secondary": {k: full[k] for k in ("trait_auc", "trait_recall_at_100", "target_ap", "target_auc", "overall_ap_by_seed")},
           "by_year": {str(y): v for y, v in full["by_year"].items()}, "scores_sha256": fin["scores_sha256"]}
    (paths.DOCS_DATA / f"model_{model}_{cutoff}_scored.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def _selftest() -> int:
    try:
        import torch
    except ImportError:
        print("  FAIL  torch is not installed: pip install -e '.[model]'")
        return 1
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    ck("the grid is every combination, and the model without a graph has no steps", len(settings("lightgcn")) == 24 and len(settings("mf")) == 8 and {s["layers"] for s in settings("mf")} == {0})
    A = np.zeros((4, 3), bool)
    A[[0, 1, 2, 3, 0], [0, 0, 1, 2, 1]] = True
    e = link_edges(A)
    ck("a link is an edge from its target node to its trait node, numbered after the targets", sorted(map(tuple, e.tolist())) == [(0, 4), (0, 5), (1, 4), (2, 5), (3, 6)])
    adj = adjacency(7, e, "cpu").to_dense().numpy()
    ck("the adjacency is symmetric and each edge is weighted by its two degrees", np.allclose(adj, adj.T) and abs(adj[0, 4] - 1 / np.sqrt(2 * 2)) < 1e-6 and abs(adj[3, 6] - 1.0) < 1e-6 and adj.diagonal().sum() == 0)
    tables = {"partners": pd.DataFrame({"a": ["t0", "t1", "zz"], "b": ["t1", "t2", "t0"]}),
              "disease": pd.DataFrame({"target": ["t0", "t1", "t2", "t3", "t0"], "disease": ["d1", "d1", "x9", "x9", "solo"]}),
              "pathways": pd.DataFrame({"target": ["t0", "t1", "t2", "zz"], "grp": ["p1", "p1", "p2", "p2"]})}
    ce, n, counts = context_edges(tables, ["t0", "t1", "t2", "t3"], ["d0", "d1", "d2"])
    ck("partners join two training targets, and a protein outside them is left out", counts["partner_edges"] == 2 and (0, 1) in map(tuple, ce.tolist()))
    ck("a group needs two training targets to become a node", counts["pathways_nodes"] == 1 and counts["pathways_edges"] == 2)
    ck("a disease with a trait's identifier is that trait's node, and another shared disease is a new node",
       counts["disease_edges_onto_a_trait_node"] == 2 and counts["disease_nodes"] == 1 and (0, 4 + 1) in map(tuple, ce.tolist()) and n == 4 + 3 + 1 + 1)
    At = torch.as_tensor(np.array([[True, False], [True, False], [False, False]]))
    neg = draw_negatives(At, torch.zeros(500, dtype=torch.long), torch.Generator().manual_seed(0), 3)
    ck("a negative is never a target linked to the trait", set(neg.tolist()) == {2})

    # three groups of 20 targets and 4 traits each; a target links to nearly all of its own group's traits,
    # so a held-out link is far likelier inside a group than a random pair is
    rng = np.random.default_rng(0)
    own = (np.arange(60)[:, None] // 20) == (np.arange(12)[None, :] // 4)
    full = own & (rng.random((60, 12)) < 0.9)
    held = full & (rng.random((60, 12)) < 0.3)
    A = full & ~held
    chance = primary(np.zeros(A.shape, np.float32), A, held)
    val = lambda S: primary(S, A, held)
    s = {"dim": 16, "layers": 2, "lr": 0.01, "l2": 1e-6}
    S_g, hist = fit(A, s, 0, "cpu", 60, link_edges(A), 72, validate=val)
    S_m, _ = fit(A, {**s, "layers": 0}, 0, "cpu", 60, validate=val)
    ck(f"both models rank held-out links of a planted structure well above chance ({chance:.2f} -> graph {val(S_g):.2f}, tables {val(S_m):.2f})", val(S_g) > 3 * chance and val(S_m) > 2 * chance)
    ck("validation is recorded every ten passes and the best matrix is the one returned", [h["passes"] for h in hist][:3] == [10, 20, 30] and abs(val(S_g) - max(h["validation"] for h in hist)) < 1e-9)
    S_again, _ = fit(A, s, 0, "cpu", 60, link_edges(A), 72, validate=val)
    ck("the same seed gives the same matrix", np.allclose(S_g, S_again, atol=1e-5))
    S_fixed, h2 = fit(A, s, 1, "cpu", 20, link_edges(A), 72)
    ck("without validation the model is fitted for exactly the passes asked and nothing is recorded", h2 == [] and S_fixed.shape == A.shape)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=("tune", "final", "score"))
    ap.add_argument("release")
    ap.add_argument("--cutoff", type=int, required=True)
    ap.add_argument("--model", choices=list(MODELS), required=True)
    ap.add_argument("--only", type=int, help="tune: try the first N settings only (a dry run; chooses nothing)")
    a = ap.parse_args(argv)
    if a.stage == "tune":
        out = tune(a.release, a.cutoff, a.model, a.only)
        print(f"  best: {out['chosen']} validation {out['validation']:.4f} (neighbours baseline {out['baselines_on_validation']['B3']:.4f})" + ("" if out["complete"] else "  [partial grid: nothing chosen]"))
    elif a.stage == "final":
        out = final(a.release, a.cutoff, a.model)
        print(f"  {a.model}: primary {out['primary_of_the_mean_matrix']:.4f} for the mean of five seeds")
    else:
        out = scored(a.release, a.cutoff, a.model)
        for other, views in out["paired"].items():
            print("  " + other + ": " + "  ".join(f"{v} {x['difference']:+.4f} [{x['low']:+.4f}, {x['high']:+.4f}]" for v, x in views.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
