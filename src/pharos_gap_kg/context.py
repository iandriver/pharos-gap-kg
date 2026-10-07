"""Do Pharos's undated disease associations already hold the links first reported after the cutoff?

The `disease` table gives target and disease pairs from twelve sources, as of 2022, with no dates.
Used as context for a model it could carry the held-out links in under another source's name. This
counts, per source, how many candidates it asserts and how many of those are later links, and scores
each source as if it were a predictor. It also reports the primary measure by kind of trait.

Run:  python -m pharos_gap_kg.context pharos319 --cutoff 2018
      python -m pharos_gap_kg.context --selftest
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import baselines, ontology, paths, registry, score, split

DISEASE_LIKE = ("MONDO:", "Orphanet:", "HP:")


def kinds(traits: list[str], terms: dict) -> tuple[list[str | None], np.ndarray]:
    """Each trait's live term, and whether that term is disease-like."""
    canon = [ontology.canonical(t, terms) for t in traits]
    return canon, np.array([bool(c) and c.startswith(DISEASE_LIKE) for c in canon])


def source_matrices(rows: pd.DataFrame, targets: list[str], canon: list[str | None]) -> dict[str, np.ndarray]:
    """source -> matrix[target, trait]: the source has a row for the target under the trait's live term.

    `rows` has columns source, target, mondoid."""
    ti = {t: i for i, t in enumerate(targets)}
    by_term: dict[str, list[int]] = {}
    for j, c in enumerate(canon):
        if c:
            by_term.setdefault(c, []).append(j)
    rows = rows[rows.target.isin(ti) & rows.mondoid.isin(by_term)]
    out = {}
    for s, g in rows.groupby("source"):
        C = np.zeros((len(targets), len(canon)), bool)
        r = g.target.map(ti).to_numpy()
        for i, term in zip(r, g.mondoid):
            C[i, by_term[term]] = True
        out[s] = C
    return out


def audit(C: np.ndarray, A: np.ndarray, Y: np.ndarray, cols: np.ndarray) -> dict:
    cand = ~A[:, cols]
    c, y = C[:, cols], Y[:, cols]
    n_in, n_pos = int((c & cand).sum()), int((c & y).sum())
    r = score.by_axis(C.astype(np.float32)[:, cols], A[:, cols], y, axis=1)
    return {"candidates_in_source": n_in, "positives_in_source": n_pos,
            "precision": (n_pos / n_in) if n_in else None, "recall": n_pos / int(y.sum()),
            "training_links_in_source": int((c & A[:, cols]).sum()), "primary": score.interval(r["ap"]), "per_trait": r["ap"]}


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    terms = ontology.load("""[Term]
id: MONDO:1
name: a disease

[Term]
id: efo:EFO_0000002
name: obsolete_a disease
is_obsolete: true
replaced_by: http://purl.obolibrary.org/obo/MONDO_1

[Term]
id: efo:EFO_0000003
name: a measurement
""".splitlines())
    traits = ["EFO:0000002", "EFO:0000003", "MONDO:1"]
    canon, dis = kinds(traits, terms)
    ck("a retired id is disease-like when its replacement is, and a measurement is not", canon == ["MONDO:1", "EFO:0000003", "MONDO:1"] and dis.tolist() == [True, False, True])
    rows = pd.DataFrame({"source": ["S", "S", "T", "S"], "target": ["A", "B", "A", "Z"], "mondoid": ["MONDO:1", "MONDO:9", "MONDO:1", "MONDO:1"]})
    C = source_matrices(rows, ["A", "B", "C"], canon)
    ck("a row marks every trait that stands for its term, for that target only", C["S"].tolist() == [[True, False, True], [False, False, False], [False, False, False]])
    ck("each source has its own matrix, and unknown targets and terms are dropped", set(C) == {"S", "T"} and C["T"].sum() == 2)
    A = np.array([[1, 0, 0], [0, 0, 1], [0, 0, 0]], bool)
    Y = np.array([[0, 0, 1], [1, 0, 0], [0, 0, 0]], bool)
    a = audit(C["S"], A, Y, np.array([0, 2]))
    # candidates on traits 0 and 2: (B,0),(C,0),(A,2),(C,2). Source asserts (A,0) [a training link] and (A,2) [a positive].
    ck("candidates, positives and training links in a source are counted apart", (a["candidates_in_source"], a["positives_in_source"], a["training_links_in_source"]) == (1, 1, 1))
    ck("precision is positives over candidates in the source; recall is over all positives", a["precision"] == 1.0 and a["recall"] == 0.5)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("release")
    ap.add_argument("--cutoff", type=int, required=True)
    a = ap.parse_args(argv)
    M = baselines.matrices(split.load(a.cutoff))
    A, Y = M["train"], M["test"]
    terms = ontology.read(registry.require("efo_obo"))
    canon, dis = kinds(M["traits"], terms)
    scored = Y.any(axis=0)
    S = baselines.scorers(A)
    res = {"by_kind": {}}
    for name, mask in (("disease_like", dis), ("measurement_or_other", ~dis)):
        cols = np.flatnonzero(mask & scored)
        res["by_kind"][name] = {"traits": int(len(cols)), "positives": int(Y[:, cols].sum()), "candidates": int((~A[:, cols]).sum())}
        per = {}
        for k in ("B0", "B1", "B3"):
            r = score.by_axis(S[k][:, cols], A[:, cols], Y[:, cols], axis=1)
            per[k] = r["ap"]
            res["by_kind"][name][k] = score.interval(r["ap"])
        m, lo, hi = score.interval(per["B3"] - per["B1"])
        res["by_kind"][name]["B3-B1"] = {"difference": m, "low": lo, "high": hi}
        b = res["by_kind"][name]
        print(f"{name}: {b['traits']} traits, {b['positives']:,} positives; chance {b['B0'][0]:.4f}, degree {b['B1'][0]:.4f}, "
              f"neighbours {b['B3'][0]:.4f}; neighbours minus degree {m:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    d = paths.tables(a.release)
    rows = pq.read_table(d / "disease.parquet", columns=["dtype", "protein_id", "mondoid"]).to_pandas().dropna(subset=["mondoid"])
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot"]).to_pandas().rename(columns={"id": "protein_id"})
    rows = rows.merge(prot, on="protein_id")
    rows["target"] = "UniProtKB:" + rows.uniprot
    rows = rows.rename(columns={"dtype": "source"})[["source", "target", "mondoid"]].drop_duplicates()
    mondo = np.array([bool(c) and c.startswith("MONDO:") for c in canon])
    cols = np.flatnonzero(mondo & scored)
    C = source_matrices(rows, M["targets"], canon)
    C["any source"] = np.any(list(C.values()), axis=0)
    base = {k: score.by_axis(S[k][:, cols], A[:, cols], Y[:, cols], axis=1)["ap"] for k in ("B0", "B1", "B3")}
    res["audit"] = {"traits": int(len(cols)), "positives": int(Y[:, cols].sum()), "candidates": int((~A[:, cols]).sum()),
                    "base_rate": float(Y[:, cols].sum() / (~A[:, cols]).sum()),
                    "baselines": {k: score.interval(v) for k, v in base.items()}, "sources": {}}
    au = res["audit"]
    print(f"\naudit over {au['traits']} scored traits whose live term is a MONDO id: {au['positives']:,} positives among {au['candidates']:,} candidates "
          f"(base rate {au['base_rate']:.4f}); degree {au['baselines']['B1'][0]:.4f}, neighbours {au['baselines']['B3'][0]:.4f}")
    print(f"  {'source':34s} {'candidates in it':>16s} {'positives':>10s} {'precision':>10s} {'recall':>7s} {'training links':>15s} {'as a scorer':>12s}")
    for s in sorted(C, key=lambda k: -int((C[k][:, cols] & Y[:, cols]).sum())):
        r = audit(C[s], A, Y, cols)
        per = r.pop("per_trait")
        m, lo, hi = score.interval(per - base["B1"])
        r["minus_degree"] = {"difference": m, "low": lo, "high": hi}
        au["sources"][s] = r
        prec = f"{r['precision']:.3f}" if r["precision"] is not None else "-"
        print(f"  {s:34s} {r['candidates_in_source']:>16,} {r['positives_in_source']:>10,} {prec:>10s} {r['recall']:>7.3f} "
              f"{r['training_links_in_source']:>15,} {r['primary'][0]:>12.4f}")
    frozen = json.loads(split.frozen_path(a.cutoff).read_text())
    out = paths.DOCS_DATA / f"context_gene_trait_{a.cutoff}.json"
    out.write_text(json.dumps({"cutoff": a.cutoff, "release": a.release, "split": {k: frozen[k]["sha256"] for k in split.PARTS}, **res}, indent=1) + "\n")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
