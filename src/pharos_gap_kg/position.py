"""Where each target's gene sits, which targets are neighbours, and which positives share a locus.

Coordinates come from an Ensembl GTF (gene lines on primary chromosomes). A Pharos target is tied to
its Ensembl gene ids by the snapshot's own cross-references. Two targets are neighbours when genes of
theirs lie on one chromosome with at most `WINDOW` bases between them; overlapping genes are 0 apart.

Run:  python -m pharos_gap_kg.position --selftest
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import paths

WINDOW = 500_000
PRIMARY = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
GENE_ID = re.compile(r'gene_id "([^"]+)"')


def genes(path: Path) -> pd.DataFrame:
    """ensg, chr, start, end for every gene on a primary chromosome."""
    rows = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 9)
            if c[2] == "gene" and c[0] in PRIMARY:
                rows.append((GENE_ID.search(c[8]).group(1), c[0], int(c[3]), int(c[4])))
    return pd.DataFrame(rows, columns=["ensg", "chr", "start", "end"])


def target_genes(release: str, g: pd.DataFrame) -> pd.DataFrame:
    """target (UniProtKB CURIE), ensg, chr, start, end: one row per gene a target is tied to."""
    d = paths.tables(release)
    prot = pq.read_table(d / "protein.parquet", columns=["id", "uniprot"]).to_pandas().rename(columns={"id": "protein_id"})
    x = pq.read_table(d / "xref.parquet", columns=["xtype", "protein_id", "value"]).to_pandas()
    x = x[x.xtype == "ENSG"][["protein_id", "value"]].rename(columns={"value": "ensg"})
    t = pq.read_table(d / "tiga.parquet", columns=["protein_id", "ensg"]).to_pandas()
    m = pd.concat([x, t]).dropna()
    m["ensg"] = m.ensg.str.split(".").str[0]
    m = m.drop_duplicates().merge(g, on="ensg").merge(prot, on="protein_id")
    m["target"] = "UniProtKB:" + m.uniprot
    return m[["target", "ensg", "chr", "start", "end"]].drop_duplicates()


def neighbour_pairs(tg: pd.DataFrame, index: dict[str, int], window: int = WINDOW) -> np.ndarray:
    """Pairs (i, j), i != j, of target indices that are neighbours. Both directions are listed."""
    tg = tg[tg.target.isin(index)]
    pairs: set[tuple[int, int]] = set()
    for _, grp in tg.groupby("chr", sort=False):
        grp = grp.sort_values("start")
        idx = grp.target.map(index).to_numpy()
        start, end = grp.start.to_numpy(), grp.end.to_numpy()
        for a in range(len(grp)):
            reach = end[a] + window
            b = a + 1
            while b < len(grp) and start[b] <= reach:
                if idx[a] != idx[b]:
                    pairs.add((idx[a], idx[b]))
                    pairs.add((idx[b], idx[a]))
                b += 1
    return np.array(sorted(pairs), dtype=np.int64).reshape(-1, 2)


def neighbour_table(tg: pd.DataFrame, targets: list[str], window: int = WINDOW) -> pd.DataFrame:
    """a, b, gap: every ordered pair of neighbouring targets and the fewest bases between genes of theirs."""
    keep = set(targets)
    tg = tg[tg.target.isin(keep)]
    best: dict[tuple[str, str], int] = {}
    for _, grp in tg.groupby("chr", sort=False):
        grp = grp.sort_values("start")
        name, start, end = grp.target.to_numpy(), grp.start.to_numpy(), grp.end.to_numpy()
        for a in range(len(grp)):
            reach = end[a] + window
            b = a + 1
            while b < len(grp) and start[b] <= reach:
                if name[a] != name[b]:
                    gap = max(0, int(start[b] - end[a]))
                    for k in ((name[a], name[b]), (name[b], name[a])):
                        if gap < best.get(k, window + 1):
                            best[k] = gap
                b += 1
    return pd.DataFrame([(a, b, g) for (a, b), g in sorted(best.items())], columns=["a", "b", "gap"])


def neighbour_count(pairs: np.ndarray, A: np.ndarray) -> np.ndarray:
    """For each (target, trait): how many of the target's neighbours have a training link to the trait."""
    out = np.zeros(A.shape, np.float32)
    if len(pairs):
        order = np.argsort(pairs[:, 0], kind="stable")
        rows, cols = pairs[order, 0], pairs[order, 1]
        cut = np.flatnonzero(np.diff(rows)) + 1
        for i, js in zip(rows[np.r_[0, cut]], np.split(cols, cut)):
            out[i] = A[js].sum(axis=0)
    return out


def components(pairs: np.ndarray, members: np.ndarray) -> np.ndarray:
    """A locus label for each of `members` (target indices): connected components of the neighbour
    relation restricted to them. A member with no neighbour among them is a locus of one."""
    label = {int(m): int(m) for m in members}

    def find(x: int) -> int:
        while label[x] != x:
            label[x] = label[label[x]]
            x = label[x]
        return x

    inside = np.isin(pairs[:, 0], members) & np.isin(pairs[:, 1], members) if len(pairs) else np.zeros(0, bool)
    for i, j in pairs[inside]:
        a, b = find(int(i)), find(int(j))
        if a != b:
            label[max(a, b)] = min(a, b)
    return np.array([find(int(m)) for m in members])


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    tg = pd.DataFrame({"target": ["A", "B", "C", "D", "E", "E", "F"], "ensg": list("abcdefg"),
                       "chr": ["1", "1", "1", "2", "1", "2", "1"],
                       "start": [1_000, 400_000, 1_200_000, 1_000, 1_650_000, 900_000_000, 5_000_000],
                       "end": [2_000, 410_000, 1_210_000, 2_000, 1_660_000, 900_001_000, 5_001_000]})
    index = {t: i for i, t in enumerate("ABCDEFG")}
    P = neighbour_pairs(tg, index)
    have = {(int(a), int(b)) for a, b in P}
    ck("genes within the window are neighbours, in both directions", {(0, 1), (1, 0)} <= have)
    ck("a gap is measured from one gene's end to the next one's start", (1, 2) not in have)          # 1,200,000 - 410,000 > window
    ck("a target with two genes is a neighbour through either", (2, 4) in have and (4, 2) in have)    # 1,650,000 - 1,210,000 <= window
    ck("different chromosomes are never neighbours", (0, 3) not in have and (3, 4) not in have)
    ck("a target in no pair is absent, and no target neighbours itself", 5 not in {a for a, _ in have} and all(a != b for a, b in have))
    A = np.zeros((7, 2), bool)
    A[1, 0] = A[2, 1] = True
    N = neighbour_count(P, A)
    ck("a neighbour's training link is counted for the target, not the target's own", N[0, 0] == 1 and N[1, 0] == 0 and N[4, 1] == 1 and N.sum() == 2)
    ck("members chained by neighbours share a locus; the rest stand alone",
       len(set(components(P, np.array([0, 1, 5])))) == 2 and len(set(components(P, np.array([0, 2])))) == 2
       and len(set(components(P, np.array([2, 4, 3])))) == 2)
    ck("no pairs gives every member its own locus", len(set(components(np.zeros((0, 2), np.int64), np.array([3, 4])))) == 2)
    NT = neighbour_table(tg, list("ABCDEFG"))
    gap = {(r.a, r.b): r.gap for r in NT.itertuples()}
    ck("the table lists the same pairs as the index form, with the bases between them",
       {(index[a], index[b]) for a, b in gap} == have and gap[("A", "B")] == 398_000 and gap[("B", "A")] == 398_000)
    ck("overlapping genes are 0 apart", neighbour_table(pd.DataFrame({"target": ["P", "Q"], "ensg": ["p", "q"], "chr": ["1", "1"], "start": [10, 50], "end": [100, 60]}), ["P", "Q"]).gap.tolist() == [0, 0])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else 0)
