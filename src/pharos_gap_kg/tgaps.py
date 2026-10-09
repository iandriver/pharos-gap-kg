"""Target-level gaps: what a target lacks for the next development level, and whether that predicted what came.

Pharos grades every target by explicit evidence. So for a target at one level, what it lacks for the
next is a matter of record, and so is how far it is from having it. This module reads that from a
database release (docs/PROTOCOL.md, "Target-level gaps"):

  THE DARK GAP    a Tdark target meets at least two of three criteria (PubMed score under 5, at most 3
                  GeneRIFs, at most 50 antibodies) and has no experimental GO leaf annotation. It leaves
                  by clearing one or two of them. The gap says how many, how far each is, and which is
                  nearest: that is what it proposes.
  THE LIGAND GAP  a target with no qualifying compound lacks one at its family's potency cutoff. The
                  release holds only compounds that pass, so this gap has no distance.

and then asks Pharos 4.0.0, two to four years on, what happened: which Tdark targets were promoted,
through which criterion, and which targets gained a ligand.

Run:  python -m pharos_gap_kg.tgaps pharos319
      python -m pharos_gap_kg.tgaps --selftest
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from pharos_gap_kg import graph, levels, paths, prune2, registry, score

CLEARS = {"pubmed_score": 5.0, "generifs": 4.0, "antibodies": 51.0}      # the least value that is no longer dark
CRITERIA = list(CLEARS)
RECENT = (2019, 2020, 2021)                                              # the last three complete years of the score in 3.19
CUTOFF = {"Kinase": "30 nM", "GPCR": "100 nM", "oGPCR": "100 nM", "NR": "100 nM", "IC": "10 uM"}
CUTOFF_OTHER = "1 uM"
AGREE_AT_LEAST = 0.995
N_BOOT = 2000


# ---- the evidence behind a level, as one release holds it ---------------------------------------------------

def evidence(release: str) -> pd.DataFrame:
    """One row per target: uniprot, sym, tdl, fam, pubmed_score, generifs, antibodies, go (an experimental
    leaf annotation), recent (score added in RECENT), has_rows (a ligand activity row), partners_with_ligand."""
    d = paths.tables(release)
    rd = lambda t, cols: pq.read_table(d / f"{t}.parquet", columns=cols).to_pandas()
    x = levels.old_levels(release).merge(rd("t2tc", ["target_id", "protein_id"]), on="target_id")
    ti = rd("tdl_info", ["itype", "protein_id", "number_value", "integer_value"])
    pick = lambda itype, col: ti[ti.itype == itype].dropna(subset=["protein_id"]).set_index("protein_id")[col]
    x["pubmed_score"] = x.protein_id.map(pick("JensenLab PubMed Score", "number_value"))
    x["antibodies"] = x.protein_id.map(pick("Ab Count", "integer_value"))
    if x.pubmed_score.isna().any() or x.antibodies.isna().any():
        raise SystemExit(f"REFUSED {release}: {int(x.pubmed_score.isna().sum())} targets have no PubMed score and {int(x.antibodies.isna().sum())} no antibody count")
    x["generifs"] = x.protein_id.map(rd("generif", ["protein_id"]).protein_id.value_counts()).fillna(0).astype(int)
    x["go"] = x.protein_id.isin(set(ti[ti.itype == "Experimental MF/BP Leaf Term GOA"].protein_id))
    pm = rd("pmscore", ["protein_id", "year", "score"])
    x["recent"] = x.protein_id.map(pm[pm.year.isin(RECENT)].groupby("protein_id").score.sum()).fillna(0.0)
    x["has_rows"] = x.target_id.isin(set(rd("ncats_ligand_activity", ["target_id"]).target_id))
    holders = set("UniProtKB:" + x.uniprot[x.tdl.isin(levels.HAS_LIGAND)])
    p = pd.read_parquet(paths.context(release) / "partners.parquet", columns=["a", "b"])
    both = pd.concat([p.rename(columns={"a": "t", "b": "o"}), p.rename(columns={"b": "t", "a": "o"})])
    x["partners_with_ligand"] = ("UniProtKB:" + x.uniprot).map(both[both.o.isin(holders)].groupby("t").o.nunique()).fillna(0).astype(int)
    return x.drop(columns=["protein_id"])


def dark(x: pd.DataFrame) -> pd.DataFrame:
    """Which of the three dark criteria each row meets."""
    return pd.DataFrame({"pubmed_score": x.pubmed_score < CLEARS["pubmed_score"], "generifs": x.generifs < CLEARS["generifs"],
                         "antibodies": x.antibodies < CLEARS["antibodies"]}, index=x.index)


def regrade(x: pd.DataFrame, go: pd.Series) -> pd.Series:
    """The level the rule gives a target that holds no qualifying compound."""
    return pd.Series(np.where((dark(x).sum(axis=1) >= 2) & ~go, "Tdark", "Tbio"), index=x.index)


def rule_check(x: pd.DataFrame) -> dict:
    low = x[x.tdl.isin(["Tdark", "Tbio"])]
    mine = regrade(low, low.go)
    out = {"targets": int(len(low)), "agree": float((mine == low.tdl).mean()), "differ": int((mine != low.tdl).sum()),
           "rule_dark_release_bio": int(((mine == "Tdark") & (low.tdl == "Tbio")).sum()), "rule_bio_release_dark": int(((mine == "Tbio") & (low.tdl == "Tdark")).sum())}
    if out["agree"] < AGREE_AT_LEAST:
        raise SystemExit(f"REFUSED: the rule gives the release's own grade for {out['agree']:.4f} of its Tdark and Tbio targets; the gaps would be read from a rule the release does not follow")
    return out


# ---- the dark gap ------------------------------------------------------------------------------------------

def need(criterion: str, value: float) -> str:
    """What clearing one criterion takes, in the units a person would act on."""
    if criterion == "pubmed_score":
        return f"PubMed score up by {CLEARS[criterion] - value:.1f} (now {value:.1f}, clears at 5)"
    n = int(round(CLEARS[criterion] - value))
    if criterion == "generifs":
        return f"{n} more GeneRIF{'s' if n != 1 else ''} (now {int(value)}, clears at 4)"
    return f"{n} more antibod{'ies' if n != 1 else 'y'} (now {int(value)}, clears at 51)"


def dark_gaps(x: pd.DataFrame) -> pd.DataFrame:
    """For each Tdark target whose evidence agrees with its grade: flips, cost, route, proposal, and its
    shortfall on each criterion (NaN where the criterion is already clear)."""
    t = x[x.tdl == "Tdark"]
    met = dark(t)
    t, met = t[met.sum(axis=1) >= 2], met[met.sum(axis=1) >= 2]
    short = pd.DataFrame({c: (CLEARS[c] - t[c]) / CLEARS[c] for c in CRITERIA}, index=t.index).where(met)
    flips = (met.sum(axis=1) - 1).astype(int)
    order = np.argsort(short.fillna(np.inf).to_numpy(), axis=1, kind="stable")
    vals = np.take_along_axis(short.fillna(np.inf).to_numpy(), order, axis=1)
    out = t.copy()
    out["flips"] = flips
    out["cost"] = [float(v[:k].sum()) for v, k in zip(vals, flips)]
    out["route"] = ["|".join(CRITERIA[j] for j in o[:k]) for o, k in zip(order, flips)]
    out["proposal"] = ["; and ".join(need(c, r[c]) for c in route.split("|")) for route, (_, r) in zip(out.route, out.iterrows())]
    for c in CRITERIA:
        out[f"short_{c}"] = short[c]
    return out


def bands(g: pd.DataFrame) -> pd.Series:
    """`1 flip, nearest third` ... `2 flips, farthest third`: thirds of route cost within each number of flips."""
    name = {0: "nearest third", 1: "middle third", 2: "farthest third"}
    out = pd.Series("", index=g.index)
    for k, sub in g.groupby("flips"):
        third = np.minimum((sub.cost.rank(method="first") - 1) * 3 // len(sub), 2).astype(int)
        out.loc[sub.index] = [f"{k} flip{'s' if k != 1 else ''}, {name[t]}" for t in third]
    return out


# ---- what Pharos 4.0.0 shows -------------------------------------------------------------------------------

def later() -> pd.DataFrame:
    """uniprot, tdl_new and the evidence counts behind it, as numbers."""
    cols = {"uniprot_id": "uniprot", "tdl": "tdl_new", "tdl_pm_score": "pubmed_score", "tdl_generif_count": "generifs", "tdl_antibody_count": "antibodies",
            "tdl_go_term_count": "go_terms", "tdl_ligand_count": "ligands", "tdl_drug_count": "drugs"}
    d = pd.read_csv(registry.require("pharos400_canonicals"), dtype=str, keep_default_na=False, usecols=list(cols)).rename(columns=cols)
    for c in list(cols.values())[2:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def later_check(n: pd.DataFrame) -> dict:
    """Does 4.0.0's level follow from its own evidence columns, by the rule in its release notes?"""
    low = n[n.tdl_new.isin(["Tdark", "Tbio"]) & n[CRITERIA].notna().all(axis=1)]
    mine = regrade(low, low.go_terms.fillna(0) > 0)
    return {"targets": int(len(low)), "agree": float((mine == low.tdl_new).mean()), "rule_dark_release_bio": int(((mine == "Tdark") & (low.tdl_new == "Tbio")).sum()),
            "rule_bio_release_dark": int(((mine == "Tbio") & (low.tdl_new == "Tdark")).sum()),
            "without_all_three_counts": int((n.tdl_new.isin(["Tdark", "Tbio"]) & ~n[CRITERIA].notna().all(axis=1)).sum())}


# ---- measures ----------------------------------------------------------------------------------------------

def measure(value: np.ndarray, label: np.ndarray, seed: int = 0, n_boot: int = N_BOOT) -> dict:
    """ROC area and average precision of ranking by `value`, each with a 95% interval over targets."""
    value, label = np.asarray(value, float), np.asarray(label, bool)
    a = prune2.auc_interval(value, label, seed=seed, n_boot=n_boot)
    ap = float(np.mean([score.ranked(value, label, np.random.default_rng(s))[0] for s in score.SEEDS]))
    rng = np.random.default_rng(seed)
    b = []
    for _ in range(n_boot):
        i = rng.integers(0, len(value), len(value))
        if label[i].any() and not label[i].all():
            b.append(score.ranked(value[i], label[i], rng)[0])
    lo, hi = np.quantile(b, [0.025, 0.975])
    return {"auc": a["auc"], "auc_low": a["low"], "auc_high": a["high"], "ap": ap, "ap_low": float(lo), "ap_high": float(hi), "base_rate": float(label.mean())}


def paired_auc(a: np.ndarray, b: np.ndarray, label: np.ndarray, seed: int = 0, n_boot: int = N_BOOT) -> dict:
    """ROC area of `a` minus that of `b` on the same targets, with a 95% interval."""
    label = np.asarray(label, bool)
    ca, cb = (np.unique(np.asarray(v, float), return_inverse=True)[1] for v in (a, b))
    rng = np.random.default_rng(seed)
    n = len(label)
    d = []
    for _ in range(n_boot):
        w = rng.multinomial(n, np.full(n, 1 / n)).astype(float)
        d.append(prune2.wauc(ca, label, w) - prune2.wauc(cb, label, w))
    lo, hi = np.nanquantile(d, [0.025, 0.975])
    return {"difference": prune2.wauc(ca, label) - prune2.wauc(cb, label), "low": float(lo), "high": float(hi)}


def proposal_check(g: pd.DataFrame, seed: int = 0, n_boot: int = N_BOOT) -> dict:
    """For promoted targets that needed one flip: was the criterion the gap named among those that
    cleared? `g` has route, the 3.19 dark flags `was_<c>` and the 4.0.0 values `<c>_new`."""
    one = g[(g.flips == 1) & g.promoted]
    cleared = pd.DataFrame({c: one[f"was_{c}"] & (one[f"{c}_new"] >= CLEARS[c]) for c in CRITERIA})
    through = cleared.any(axis=1)
    c, r = cleared[through], one.route[through]
    right = np.array([bool(row[name]) for name, (_, row) in zip(r, c.iterrows())], float)
    chance = c.sum(axis=1).to_numpy() / 2.0
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(right), (n_boot, len(right))) if len(right) else np.zeros((0, 0), int)
    b = right[idx].mean(axis=1) - chance[idx].mean(axis=1) if len(right) else np.array([np.nan])
    sh = prune2.share(right.astype(bool))
    return {"promoted_needing_one_flip": int(len(one)), "left_through_the_counts": int(through.sum()), "left_another_way": int((~through).sum()),
            "both_criteria_cleared": int((c.sum(axis=1) == 2).sum()), "named_criterion_cleared": sh, "by_naming_at_random": float(chance.mean()) if len(chance) else float("nan"),
            "above_random": {"difference": float(right.mean() - chance.mean()) if len(right) else float("nan"), "low": float(np.nanquantile(b, 0.025)), "high": float(np.nanquantile(b, 0.975))},
            "named": {k: int(v) for k, v in r.value_counts().items()}, "cleared": {k: int(v) for k, v in c.sum().items()}}


# ---- the back-test -----------------------------------------------------------------------------------------

DARK_SCORES = {"pubmed score": lambda g: g.pubmed_score, "GeneRIF count": lambda g: g.generifs, "antibody count": lambda g: g.antibodies,
               "recent attention": lambda g: g.recent, "flips needed": lambda g: -g.flips, "route cost": lambda g: -g.cost}
LIGAND_SCORES = {"pubmed score": lambda g: g.pubmed_score, "GeneRIF count": lambda g: g.generifs, "antibody count": lambda g: g.antibodies,
                 "recent attention": lambda g: g.recent, "partners with a ligand": lambda g: g.partners_with_ligand}
SINGLE = ("pubmed score", "GeneRIF count", "antibody count")


def run(release: str) -> dict:
    x = evidence(release)
    out = {"release": release, "later": "pharos400_canonicals", "rule_on_this_release": rule_check(x)}
    n = later()
    out["rule_on_the_later_release"] = later_check(n)
    g = dark_gaps(x)
    out["dark_gap"] = {"tdark_targets": int((x.tdl == "Tdark").sum()), "with_a_gap": int(len(g)), "grade_and_evidence_disagree": int((x.tdl == "Tdark").sum() - len(g)),
                       "flips": {str(k): int(v) for k, v in g.flips.value_counts().sort_index().items()},
                       "route": {k: int(v) for k, v in g.route.value_counts().items()}}
    met = dark(g)
    j = g.assign(**{f"was_{c}": met[c] for c in CRITERIA}, band=bands(g)).merge(n.rename(columns={c: f"{c}_new" for c in CRITERIA}), on="uniprot")
    j["promoted"] = j.tdl_new != "Tdark"
    y = j.promoted.to_numpy()
    d = out["dark_gap"]
    d.update({"in_both_releases": int(len(j)), "promoted": int(y.sum()), "to": {k: int(v) for k, v in j.tdl_new[j.promoted].value_counts().items()}})
    d["scores"] = {k: measure(f(j).to_numpy(), y) for k, f in DARK_SCORES.items()}
    best = max(SINGLE, key=lambda k: d["scores"][k]["auc"])
    d["best_single_count"] = best
    d["route_cost_against_best_single_count"] = paired_auc(DARK_SCORES["route cost"](j), DARK_SCORES[best](j), y)
    d["route_cost_against_each"] = {k: paired_auc(DARK_SCORES["route cost"](j), DARK_SCORES[k](j), y) for k in DARK_SCORES if k != "route cost"}
    d["bands"] = {b: prune2.share(s.promoted.to_numpy()) for b, s in sorted(j.groupby("band"))}
    d["by_flips"] = {str(k): prune2.share(s.promoted.to_numpy()) for k, s in j.groupby("flips")}
    d["proposal"] = proposal_check(j)
    stay = j[~j.promoted]
    d["counts_between_releases_for_targets_still_dark"] = {
        c: {"targets": int(stay[f"{c}_new"].notna().sum()), "fell": float((stay[f"{c}_new"] < stay[c] - 1e-9).mean()), "rose": float((stay[f"{c}_new"] > stay[c] + 1e-9).mean()),
            "rank_correlation": float(stay[[c, f"{c}_new"]].corr(method="spearman").iloc[0, 1])} for c in CRITERIA}
    d["promoted_how"] = {"gained_a_go_term": int((j.promoted & (j.go_terms.fillna(0) > 0)).sum()), "gained_a_ligand_or_drug": int((j.promoted & ((j.ligands.fillna(0) > 0) | (j.drugs.fillna(0) > 0))).sum()),
                         **{f"cleared_{c}": int((j.promoted & j[f"was_{c}"] & (j[f"{c}_new"] >= CLEARS[c])).sum()) for c in CRITERIA}}

    cand = x[~x.tdl.isin(levels.HAS_LIGAND) & ~x.has_rows].merge(n[["uniprot", "tdl_new"]], on="uniprot")
    gained = cand.tdl_new.isin(levels.HAS_LIGAND).to_numpy()
    stored = json.loads((paths.DOCS_DATA / f"levels_{release}_pharos400.json").read_text())["ligand_label"]
    if (len(cand), int(gained.sum())) != (stored["candidates"], stored["gained"]):
        raise SystemExit(f"REFUSED: {len(cand):,} candidates and {int(gained.sum()):,} gains here, and the committed label says {stored['candidates']:,} and {stored['gained']:,}")
    out["ligand_gap"] = {"candidates": int(len(cand)), "gained": int(gained.sum()), "scores": {k: measure(f(cand).to_numpy(), gained) for k, f in LIGAND_SCORES.items()},
                         "cutoff_by_family": {k: int(v) for k, v in cand.fam.fillna("none").map(lambda f: CUTOFF.get(f, CUTOFF_OTHER)).value_counts().items()}}
    lb = max(SINGLE, key=lambda k: out["ligand_gap"]["scores"][k]["auc"])
    out["ligand_gap"]["best_single_count"] = lb
    out["ligand_gap"]["partners_against_best_single_count"] = paired_auc(LIGAND_SCORES["partners with a ligand"](cand), LIGAND_SCORES[lb](cand), gained)

    C = paths.context(release)
    table = g.assign(band=bands(g))[["uniprot", "sym", "fam", "pubmed_score", "generifs", "antibodies", "recent", "flips", "cost", "route", "proposal", "band"]]
    table.to_parquet(C / "target_gaps_dark.parquet", index=False)
    lig = x[~x.tdl.isin(levels.HAS_LIGAND) & ~x.has_rows][["uniprot", "sym", "tdl", "fam", "pubmed_score", "partners_with_ligand"]]
    lig = lig.assign(cutoff=lig.fam.map(lambda f: CUTOFF.get(f, CUTOFF_OTHER)))
    lig.to_parquet(C / "target_gaps_ligand.parquet", index=False)
    out["files"] = {f: graph.file_sha(C / f) for f in ("target_gaps_dark.parquet", "target_gaps_ligand.parquet")}
    out["inputs"] = {"tables_sha256": json.loads((paths.tables(release) / "_manifest.json").read_text())["sha256"],
                     "later_sha256": json.loads(registry.SOURCES["pharos400_canonicals"].path.with_name(registry.SOURCES["pharos400_canonicals"].filename + ".lineage.json").read_text())["sha256"]}
    (paths.DOCS_DATA / f"target_gaps_{release}.json").write_text(prune2.dump(out))
    return out


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    x = pd.DataFrame({"uniprot": list("abcdefg"), "sym": list("ABCDEFG"), "fam": [None] * 7,
                      "tdl": ["Tdark", "Tdark", "Tdark", "Tbio", "Tbio", "Tdark", "Tchem"],
                      "pubmed_score": [4.0, 1.0, 4.9, 80.0, 1.0, 9.0, 1.0], "generifs": [3, 0, 9, 40, 0, 2, 0], "antibodies": [200, 10, 50, 300, 0, 400, 0],
                      "go": [False, False, False, False, True, False, False], "recent": [0.0] * 7})
    m = dark(x)
    ck("a criterion is dark under 5 for the score, at 3 or fewer GeneRIFs, at 50 or fewer antibodies",
       m.pubmed_score.tolist() == [True, True, True, False, True, False, True] and m.generifs.tolist() == [True, True, False, False, True, True, True] and m.antibodies.tolist()[:3] == [False, True, True])
    ck("the rule grades a target dark on two criteria and no GO annotation, and bio otherwise", regrade(x, x.go).tolist() == ["Tdark", "Tdark", "Tdark", "Tbio", "Tbio", "Tbio", "Tdark"])
    try:
        rule_check(x)
        ck("a rule the release does not follow is refused", False)
    except SystemExit:
        ck("a rule the release does not follow is refused", True)      # f is graded dark and the rule says bio: 5 of 6 agree
    g = dark_gaps(x).set_index("sym")
    ck("a gap is given to dark targets whose evidence agrees with the grade", list(g.index) == ["A", "B", "C"])
    ck("flips needed are the dark criteria met less one", g.flips.to_dict() == {"A": 1, "B": 2, "C": 1})
    ck("the route is the nearest criterion, and its cost that criterion's shortfall over its threshold",
       g.route["A"] == "pubmed_score" and abs(g.cost["A"] - 1 / 5) < 1e-12 and g.route["C"] == "antibodies" and abs(g.cost["C"] - 1 / 51) < 1e-12)
    ck("two flips take the two nearest criteria and add their shortfalls", g.route["B"] == "pubmed_score|antibodies" and abs(g.cost["B"] - (4 / 5 + 41 / 51)) < 1e-12)
    ck("a criterion already clear has no shortfall", np.isnan(g.short_antibodies["A"]) and abs(g.short_generifs["A"] - 0.25) < 1e-12)
    ck("the proposal says what clearing the route takes", g.proposal["A"] == "PubMed score up by 1.0 (now 4.0, clears at 5)" and g.proposal["C"] == "1 more antibody (now 50, clears at 51)"
       and need("generifs", 0) == "4 more GeneRIFs (now 0, clears at 4)")
    big = pd.DataFrame({"flips": [1] * 6 + [2] * 3, "cost": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0, 1.5, 2.0]})
    ck("bands are thirds of cost within each number of flips", bands(big).tolist() == ["1 flip, nearest third"] * 2 + ["1 flip, middle third"] * 2 + ["1 flip, farthest third"] * 2
       + ["2 flips, nearest third", "2 flips, middle third", "2 flips, farthest third"])
    n = pd.DataFrame({"tdl_new": ["Tdark", "Tbio", "Tbio", "Tdark"], "pubmed_score": [1.0, 9.0, 1.0, np.nan], "generifs": [0, 9, 0, 0], "antibodies": [0, 90, 0, 0], "go_terms": [0, 0, 2, 0]})
    lc = later_check(n)
    ck("the later release is checked against its own counts, and a target without them is set aside", lc["targets"] == 3 and lc["agree"] == 1.0 and lc["without_all_three_counts"] == 1)
    v, y = np.r_[np.arange(50.0), np.arange(50.0) + 25], np.r_[np.zeros(50, bool), np.ones(50, bool)]
    ms = measure(v, y, n_boot=200)
    ck("a score that separates is above a half on both measures, with the interval above the base rate", ms["auc"] > 0.8 and ms["auc_low"] > 0.5 and ms["ap"] > ms["base_rate"] and ms["ap_low"] > 0.5)
    pa = paired_auc(v, -v, y, n_boot=200)
    ck("a score is compared with another on the same targets", pa["difference"] > 0.6 and pa["low"] > 0)
    pr = pd.DataFrame({"flips": [1, 1, 1, 1, 2], "promoted": [True, True, True, False, True], "route": ["pubmed_score", "generifs", "pubmed_score", "generifs", "generifs|antibodies"],
                       "was_pubmed_score": [True, True, True, True, True], "was_generifs": [True, True, False, True, True], "was_antibodies": [False, False, True, False, True],
                       "pubmed_score_new": [6.0, 7.0, 2.0, 1.0, 9.0], "generifs_new": [1, 9, 0, 0, 9], "antibodies_new": [0, 0, 30, 0, 99]})
    pc = proposal_check(pr, n_boot=200)
    # three promoted on one flip: the first cleared the score (named: right), the second cleared both (named GeneRIFs: right), the third cleared neither
    ck("the proposal is checked on targets that needed one flip and left through the counts", pc["promoted_needing_one_flip"] == 3 and pc["left_through_the_counts"] == 2 and pc["left_another_way"] == 1)
    ck("naming is right when the named criterion cleared, and chance is the share of the two that cleared",
       pc["named_criterion_cleared"]["share"] == 1.0 and pc["both_criteria_cleared"] == 1 and abs(pc["by_naming_at_random"] - 0.75) < 1e-12)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) != 1:
        raise SystemExit(__doc__)
    out = run(argv[0])
    r, l = out["rule_on_this_release"], out["rule_on_the_later_release"]
    print(f"  rule against the release's own grade: {r['agree']:.4f} of {r['targets']:,}; against 4.0.0's: {l['agree']:.4f} of {l['targets']:,}")
    d = out["dark_gap"]
    print(f"  dark gap: {d['with_a_gap']:,} targets, flips {d['flips']}, {d['in_both_releases']:,} in both releases, {d['promoted']:,} promoted ({d['promoted'] / d['in_both_releases']:.1%}) to {d['to']}")
    for k, m in d["scores"].items():
        print(f"    {k:18s} AUC {m['auc']:.3f} [{m['auc_low']:.3f}, {m['auc_high']:.3f}]  AP {m['ap']:.3f} [{m['ap_low']:.3f}, {m['ap_high']:.3f}]")
    p = d["route_cost_against_best_single_count"]
    print(f"    route cost minus {d['best_single_count']}: {p['difference']:+.3f} [{p['low']:+.3f}, {p['high']:+.3f}]")
    for b, s in d["bands"].items():
        print(f"    {b:26s} promoted {s['share']:.1%} [{s['low']:.1%}, {s['high']:.1%}] of {s['n']:,}")
    q = d["proposal"]
    print(f"    proposal: {q['named_criterion_cleared']['share']:.1%} of {q['left_through_the_counts']:,} named the criterion that cleared; at random {q['by_naming_at_random']:.1%}; "
          f"above random {q['above_random']['difference']:+.3f} [{q['above_random']['low']:+.3f}, {q['above_random']['high']:+.3f}]; left another way {q['left_another_way']:,}")
    print(f"    promoted how: {d['promoted_how']}")
    print(f"    counts between releases, targets still dark: {d['counts_between_releases_for_targets_still_dark']}")
    g = out["ligand_gap"]
    print(f"  ligand gap: {g['candidates']:,} candidates, {g['gained']:,} gained ({g['gained'] / g['candidates']:.1%})")
    for k, m in g["scores"].items():
        print(f"    {k:22s} AUC {m['auc']:.3f} [{m['auc_low']:.3f}, {m['auc_high']:.3f}]  AP {m['ap']:.3f} [{m['ap_low']:.3f}, {m['ap_high']:.3f}]")
    p = g["partners_against_best_single_count"]
    print(f"    partners with a ligand minus {g['best_single_count']}: {p['difference']:+.3f} [{p['low']:+.3f}, {p['high']:+.3f}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
