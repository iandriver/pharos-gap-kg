"""The self-test of `ask`: a six-node graph written to a temporary bundle, and every question put to it."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pandas as pd

from pharos_gap_kg import ask, graph


def tiny(d: Path, gaps: bool = True, scorer: bool = True, no_link: bool = True) -> None:
    T = lambda i, sym, lvl="Tbio": {"id": f"UniProtKB:P{i}", "kind": "target", "category": "biolink:Protein", "name": sym, "full_name": f"protein {sym}",
                                    "level": lvl, "family": "", "xref": "", "live_term": "", "disease_like": False,
                                    "chr": "6" if sym == "BBB" else "1", "start": 30_000_000 if sym == "BBB" else 1000 * i, "end": 30_000_500 if sym == "BBB" else 1000 * i + 500}
    D = lambda i, name, live=None, dis=True: {"id": i, "kind": "trait", "category": "biolink:Disease", "name": name, "full_name": "", "level": "",
                                              "family": "", "xref": "", "live_term": live or i, "disease_like": dis, "chr": "", "start": 0, "end": 0}
    nodes = pd.DataFrame([T(1, "AAA"), T(2, "BBB"), T(3, "CCC"), T(4, "DDD", "Tdark"), T(5, "AAB"), T(6, "AAA"),
                          D("MONDO:1", "heart disease"), D("MONDO:2", "arrhythmia"), D("MONDO:3", "lung disease"),
                          D("EFO:9", "heart disease", live="MONDO:1"), D("EFO:7", "body height", dis=False)])
    L = [("P1", "MONDO:1", 2010, 2021, 3), ("P1", "EFO:7", 2012, 2012, 1), ("P2", "MONDO:1", 2015, 2015, 1), ("P2", "EFO:7", 2011, 2011, 1),
         ("P2", "MONDO:3", 2020, 2020, 1), ("P3", "MONDO:2", 2016, 2016, 1), ("P3", "EFO:7", 2013, 2013, 1), ("P1", "MONDO:2", 2020, 2020, 1),
         ("P3", "EFO:9", 2014, 2014, 1), ("P6", "MONDO:1", 2010, 2010, 1), ("P6", "EFO:7", 2012, 2012, 1)]
    links = pd.DataFrame([{"edge_id": f"x:e{i}", "target": f"UniProtKB:{t}", "trait": tr, "first_year": a, "last_year": b, "n_studies": n, "n_papers": n}
                          for i, (t, tr, a, b, n) in enumerate(L)])
    ev = [("P1", "MONDO:1", "S1", 11, 2010), ("P1", "MONDO:1", "S2", 12, 2014), ("P1", "MONDO:1", "S3", 13, 2021)]
    ev += [(t, tr, f"S{i + 10}", 100 + i, a) for i, (t, tr, a, b, n) in enumerate(L) if not (t == "P1" and tr == "MONDO:1")]
    evidence = pd.DataFrame([{"target": f"UniProtKB:{t}", "trait": tr, "study": s, "pmid": p, "year": y} for t, tr, s, p, y in ev])
    lin = pd.DataFrame([("MONDO:1", "MONDO:2", "narrower"), ("MONDO:2", "MONDO:1", "broader"), ("MONDO:1", "EFO:9", "same term"),
                        ("EFO:9", "MONDO:1", "same term"), ("EFO:9", "MONDO:2", "narrower"), ("MONDO:2", "EFO:9", "broader")], columns=["a", "b", "relation"])
    nb = pd.DataFrame([("UniProtKB:P3", "UniProtKB:P4", 1200), ("UniProtKB:P4", "UniProtKB:P3", 1200)], columns=["a", "b", "gap"])
    for name, df in (("nodes", nodes), ("links", links), ("evidence", evidence), ("lineage", lin), ("neighbours", nb)):
        df.to_parquet(d / f"{name}.parquet", index=False)
    cal = {g: {"candidates": 10.0, "positives": 2.0, "share": 0.2} for g in ("all", "linked", "not_linked", "disease_like_not_linked")}
    bt = {"headline": {"cutoff": 2018, "data_ends": 2021, "chance": 0.02, "degree": 0.07, "neighbours": 0.1},
          "calibration": {"cutoff": 2018, "data_ends": 2021, "base_rate": 0.01,
                          "by_trait": {f"{a}-{b}": cal for a, b in ask.TRAIT_BINS}, "by_target": {f"{a}-{b}": cal for a, b in ask.TARGET_BINS}}}
    if scorer:
        bt["scorer"] = {"name": "combine[B1+B3+B5]:F1", "weights": {"B1": 0.03125, "B3": 1.0, "B5": 0.25}, "weights_fitted_on_links_to": 2021, "served_by_target": "combiner"}
        bt["prospective"] = {"cutoff": 2021, "first_year": 2022, "data_ends": 2025, "traits_scored": 40, "test_links": 900, "combiner": 0.11, "neighbours": 0.10, "degree": 0.08,
                             "combiner_minus_neighbours_all": [0.01, 0.002, 0.018], "by_target": {"combiner": 0.07, "neighbours": 0.06}}
        bt["calibration"].update({"cutoff": 2021, "first_year": 2022, "data_ends": 2025, "prospective": True, "scorer": bt["scorer"]["name"],
                                  "pooled_over": {"traits": 80, "targets": 600}})
        if no_link:
            side = lambda n, who: {who: n, "position": 0.14, "tie_breaker_alone": 0.12, "chance": 0.01, "position_minus_tie_breaker": [0.02, 0.005, 0.04]}
            one = {"all": {"candidates": 100.0, "positives": 4.0, "share": 0.04}}
            bt["no_link"] = {"headline": {"cutoff": 2021, "first_year": 2022, "data_ends": 2025, "later_links": 300, "later_links_listed": 60, "targets_with_no_link": 50,
                                          "of_them_listed_for_some_trait": 30, "by_trait": side(20, "traits"), "by_target": side(25, "targets"), "earlier_periods": {}},
                             "calibration": {"bins": [f"{a}-{b}" for a, b in ask.NO_LINK_BINS], "by_trait": {f"{a}-{b}": one for a, b in ask.NO_LINK_BINS},
                                             "by_target": {f"{a}-{b}": one for a, b in ask.NO_LINK_BINS}, "first_year": 2022, "data_ends": 2025,
                                             "base_rate_among_listed": 0.03, "base_rate_among_all": 0.001, "listed_pairs": 2000}}
    (d / "backtest.json").write_text(json.dumps(bt))
    files = ["nodes.parquet", "links.parquet", "evidence.parquet", "lineage.parquet", "neighbours.parquet", "backtest.json"]
    if gaps:
        # DDD is Tdark, one criterion short; CCC is Tbio with no compound
        tg = pd.DataFrame([{"id": "UniProtKB:P4", "level": "Tdark", "family": "Kinase", "cutoff": "30 nM", "partners_with_ligand": 2, "lacks_ligand": True,
                            "pubmed_score": 4.2, "generifs": 3, "antibodies": 90.0, "flips": 1, "cost": 0.16, "route": "pubmed_score",
                            "proposal": "PubMed score up by 0.8 (now 4.2, clears at 5)", "band": "1 flip, nearest third", "dark_gap_id": "x:g1", "ligand_gap_id": "x:g2"},
                           {"id": "UniProtKB:P3", "level": "Tbio", "family": None, "cutoff": "1 uM", "partners_with_ligand": 0, "lacks_ligand": True,
                            "pubmed_score": None, "generifs": None, "antibodies": None, "flips": None, "cost": None, "route": None, "proposal": None, "band": None,
                            "dark_gap_id": None, "ligand_gap_id": None}])
        tg["flips"] = tg.flips.astype("Int64")
        tg.to_parquet(d / "target_gaps.parquet", index=False)
        sh = lambda p, n: {"n": n, "share": p, "low": p - 0.05, "high": p + 0.05}
        rates = {"earlier": "Tiny 1", "later": "Tiny 2", "dark": {"targets": 1, "in_both_releases": 1, "promoted": 1, "route_cost_auc": 0.86,
                                                                    "best_single_count": "pubmed score", "best_single_count_auc": 0.79,
                                                                    "bands": {"1 flip, nearest third": {"promoted": sh(0.5, 40), "fell_in_the_mirror": sh(0.06, 30)}}},
                 "ligand": {"candidates": 200, "gained": 20}}
        (d / "target_gaps.json").write_text(json.dumps(rates))
        files += ["target_gaps.parquet", "target_gaps.json"]
    man = {"release": "tiny", "source_file": "tiny.sql.gz", "source_url": "https://example.org/tiny.sql.gz", "source_published": "2024-01-01",
           "source_sha256": "0" * 64, "source_tables": ["tiga"], "ontology_sha256": "1" * 64, "coordinates_sha256": "2" * 64, "window": 500000,
           "first_year": 2010, "latest_year": 2021, "targets": 6, "targets_with_a_link": 4, "traits": 5, "links": len(links), "studies": 11,
           "papers": 11, "lineage_pairs": 3, "neighbour_pairs": 1, "conformance": "10 of 10 checks passed", "dated_graphs": [],
           "files": {f: graph.file_sha(d / f) for f in files}}
    (d / "manifest.json").write_text(json.dumps(man))


def run() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        tiny(d)
        G = ask.Graph(d)

        # names
        ck("an accession and a CURIE each resolve to the one target", all(G.resolve(x, "target").get("i") == 0 for x in ("P1", "UniProtKB:P1")))
        r = G.resolve("aaa", "target")
        ck("a symbol shared by two Pharos entries resolves to the one with more links, and says so", r["ok"] and r["i"] == 0 and "2 entries" in r["note"])
        ck("a trait resolves by id, by id written with an underscore, and by exact name", G.resolve("MONDO:2")["i"] == 1 and G.resolve("MONDO_2")["i"] == 1 and G.resolve("Arrhythmia")["i"] == 1)
        r = G.resolve("heart disease")
        ck("one term under two ids resolves to the id with more links, and says so", r["ok"] and G.d.id[r["i"]] == "MONDO:1" and "EFO:9" in r["note"])
        r = G.resolve("disease")
        ck("a word that fits several traits is ambiguous and lists them", not r["ok"] and r["problem"] == "ambiguous" and len(r["candidates"]) >= 3)
        r = G.resolve("ZZZ9")
        ck("a name that matches nothing is not in the graph", not r["ok"] and r["problem"] == "not in the graph" and r["candidates"] == [])
        ck("two letters are too few to search symbols by", G.find("AA")["found"] is False)
        f = G.find("disease")
        ck("find says when no name matches exactly, and to ask rather than pick", f["exact_match"] is False and f["ask_which_is_meant"] is True and "Do not pick" in f["answer"])
        f = G.find("arrhythmia")
        ck("find says when a name matches exactly", f["exact_match"] is True and f["ask_which_is_meant"] is False and f["matches"][0]["id"] == "MONDO:2")
        ck("a three-letter prefix finds the symbols that start with it", {m["symbol"] for m in G.find("AAA", "target")["matches"]} == {"AAA"} and {m["symbol"] for m in G.find("aab", "target")["matches"]} == {"AAB"})

        # known
        k = G.known("AAA")
        ck("known counts a target's links and splits diseases from measurements", k["links"] == 3 and "2 to diseases" in k["answer"] and k["top"][0]["id"] == "MONDO:1")
        k = G.known("AAA", as_of=2015)
        ck("as of a year, later links are left out and counted apart", k["links"] == 2 and k["first_reported_after_as_of"] == 1 and k["as_of"] == 2015)
        k = G.known("DDD")
        ck("a target with no link is in the graph and says absence is not a negative", k["found"] and k["links"] == 0 and "not evidence" in k["answer"])
        k = G.known("MONDO:1")
        ck("known on a trait counts every link and lists each gene symbol once", k["links"] == 3 and [x["symbol"] for x in k["top"]] == ["AAA", "BBB"]
           and k["top"][0]["pharos_entries_with_this_symbol"] == 2 and k["more"] == 0)

        # why
        w = G.why("AAA", "MONDO:1", as_of=2015)
        ck("why returns the studies dated by the year, and counts the later ones", w["found"] and [s["study"] for s in w["studies"]] == ["S1", "S2"] and w["studies_after_as_of"] == 1 and w["edge"] == "x:e0")
        ck("why names the file the link was read from and its digest", w["source"]["file"] == "tiny.sql.gz" and len(w["source"]["sha256"]) == 64)
        w = G.why("AAA", "arrhythmia", as_of=2018)
        ck("before its first study a link is absent, with the year it was reported", not w["found"] and w["reported_later"] == 2020)
        ck("the absence names the broader trait the target already had", [x["id"] for x in w["same_trait_under_another_wording"]] == ["MONDO:1"])
        w = G.why("DDD", "arrhythmia")
        ck("an absent link names the neighbouring gene that has it", not w["found"] and w["reported_later"] is None and w["neighbouring_genes_with_this_link"][0]["symbol"] == "CCC")
        w = G.why("BBB", "arrhythmia", as_of=2018)
        ck("an absent link to a gene in the MHC region says a high rank there means little", "MHC" in w.get("caution", "") and w["candidates_for_this_trait"] == 5)
        ck("an unknown argument stops the question", G.why("QQQ", "arrhythmia")["found"] is False and G.why("QQQ", "arrhythmia")["problem"] == "not in the graph")

        # gaps
        g = G.gaps("DDD")
        ck("a dark target is told how many criteria must clear, the nearest route, and how often that closed",
           g["found"] and [x["kind"] for x in g["gaps"]] == ["dark", "ligand"] and g["gaps"][0]["criteria_to_clear"] == 1
           and "PubMed score up by 0.8" in g["answer"] and "50% (45% to 55%)" in g["answer"] and "6% of 30 targets fell the other way" in g["answer"])
        ck("its ligand gap names the cutoff for its family and says the distance is not recorded", "30 nM" in g["gaps"][1]["lacks"] and "not recorded" in g["gaps"][1]["distance"]
           and "20 (10%)" in g["gaps"][1]["how_often_this_closed"])
        ck("each gap cites the graph edge that states it, where one was exported", [x.get("graph_edge") for x in g["gaps"]] == ["x:g1", "x:g2"])
        g = G.gaps("CCC")
        ck("a target with biology and no compound has only the ligand gap, and no edge is cited where none was exported",
           [x["kind"] for x in g["gaps"]] == ["ligand"] and "1 uM" in g["answer"] and "graph_edge" not in g["gaps"][0])
        ck("every gap answer says what a gap is and is not", "not the same as learning its function" in g["reading"])
        ck("a target graded without a compound that the ligand table gives one is told the two disagree", G.gaps("BBB")["gaps"] == [] and "disagree" in G.gaps("BBB")["answer"])
        ck("an unknown name stops the question", G.gaps("QQQ")["found"] is False and G.gaps("QQQ")["problem"] == "not in the graph")
        ls = G.gaps(top=5)
        ck("with no target, the dark targets nearest to leaving are listed with their routes", [t["symbol"] for t in ls["targets"]] == ["DDD"] and "86%" in ls["how_this_ranking_did"])
        ck("a family narrows the list, and a family with none says which exist", G.gaps(family="kinase")["targets"][0]["symbol"] == "DDD" and "Kinase" in G.gaps(family="GPCR")["answer"])
        with tempfile.TemporaryDirectory() as tmp2:
            tiny(Path(tmp2), gaps=False)
            ck("a bundle built without the gaps says how to add them", ask.Graph(Path(tmp2)).gaps("DDD")["problem"] == "not built")

        # connect
        c = G.connect("arrhythmia", "AAA", as_of=2018)
        ck("connect takes its two things in either order", c["a"]["symbol"] == "AAA" and c["direct_link"] is None and c["through_trait_wording"][0]["id"] == "MONDO:1")
        ck("connect lists the target's traits that share targets with the trait", [x["id"] for x in c["through_shared_traits"]] == ["EFO:7"] and "share targets with this one: 1." in c["answer"])
        c = G.connect("AAA", "BBB")
        ck("two targets are tied by the traits they share", c["shared_traits"] == 2 and c["traits_a"] == 3 and c["traits_b"] == 3)
        ck("two neighbouring targets are told their distance", G.connect("CCC", "DDD")["bases_apart_if_neighbours"] == 1200)
        c = G.connect("MONDO:1", "MONDO:2")
        ck("two traits are tied by shared targets and by the ontology", c["shared_targets"] == 1 and "narrower" in c["ontology_relation"])

        # predict
        p = G.predict("arrhythmia", as_of=2018)
        syms = [x["symbol"] for x in p["candidates"]]
        ck("predict ranks the targets without the link and leaves out the one that has it", syms == ["AAA", "BBB"] and p["candidates"][0]["rank"] == 1)
        ck("predict says how many of its candidates were later reported", "Reported between 2019 and 2021: 1 of these 2" in p["answer"])
        ck("a gene held as two Pharos entries is one candidate, and says it stands for two", p["candidates"][0]["pharos_entries_with_this_symbol"] == 2)
        ck("a candidate in the MHC region carries the caution, and the answer counts them", "MHC" in p["candidates"][1]["caution"] and p["candidates"][0]["caution"] == ""
           and "best read as one lead: 1 of them" in p["answer"])
        ck("why says when the gene is in the MHC region", "MHC" in G.why("BBB", "MONDO:1")["caveat"] and "MHC" not in G.why("CCC", "MONDO:2")["caveat"])
        ck("a candidate reported after the year says when", {x["symbol"]: x["reported_after_as_of"] for x in p["candidates"]}["AAA"] == 2020)
        ck("a candidate an earlier link explains says how", any("heart disease" in e for x in p["candidates"] if x["symbol"] == "AAA" for e in x["explained_by"]))
        ck("each candidate quotes how often its rank has held on later links, over every trait ranked, and that it is a lower bound",
           all("links first reported from 2022 to 2025" in x["how_often_this_held"] and "counted over every trait ranked" in x["how_often_this_held"]
               and "not known to be absent" in x["how_often_this_held"] for x in p["candidates"]))
        import math
        w = G.weights
        ck("a candidate's score is the weighted sum of its three parts, and the parts are shown",
           all(abs(x["score"] - (w["B3"] * math.log1p(x["score_parts"]["shared_trait_neighbours"]) + w["B1"] * math.log1p(x["score_parts"]["traits_already_linked"])
                                 + w["B5"] * math.log1p(x["score_parts"]["genes_nearby_with_this_link"]))) < 2e-3 for x in p["candidates"]))
        ck("the method names the three parts, the test on later links, and that an earlier year is an illustration",
           "three parts" in p["method"] and "first reported from 2022 to 2025" in p["method"] and "0.110" in p["method"] and "illustration and not a test" in p["method"]
           and "illustration" not in G.predict("lung disease")["method"])
        ck("a target with no link is not in the main list, though a linked gene sits beside it, and the answer says how to rank such targets",
           "DDD" not in [x["symbol"] for x in G.predict("arrhythmia")["candidates"]] and G.why("DDD", "arrhythmia")["rank_among_candidates"] is None
           and G.predict("arrhythmia")["targets_with_no_link_beside_a_linked_gene"] == 1 and "unlinked=true ranks those" in G.predict("arrhythmia")["answer"])

        # targets with no link, ranked by position alone
        u = G.predict("arrhythmia", unlinked=True)
        ck("asked for, a trait's targets with no link are ranked by the linked genes beside them",
           [(x["rank"], x["symbol"], x["genes_nearby_with_this_link"], x["links_held_by_genes_nearby"]) for x in u["candidates"]] == [(1, "DDD", 1, 3)]
           and u["candidates"][0]["nearest"] == [{"symbol": "CCC", "bases_apart": 1200}] and "out of 1" in u["answer"])
        ck("that list says nothing is fitted, how the ranking held, and what share of listed pairs were later linked",
           "Nothing is fitted" in u["method"] and "0.140 (mean average precision over 20 traits) against 0.120" in u["method"] and "60 of those 300 links" in u["method"]
           and "4.0% of targets with no link listed at ranks 1-5 for a trait" in u["candidates"][0]["how_often_this_held"] and "Of all listed pairs 3.0% were" in u["candidates"][0]["how_often_this_held"]
           and "same association signal" in u["caveat"])
        ck("as of a year before the neighbour's link, the target is not listed", G.predict("arrhythmia", as_of=2015, unlinked=True)["candidates"] == []
           and "No target without a link" in G.predict("arrhythmia", as_of=2015, unlinked=True)["answer"])
        d4 = G.predict("DDD")
        ck("a target with no link is given the traits its neighbouring genes are linked to, commonest trait first, wordings of one trait as one lead",
           [(x["rank"], x["name"], x["genes_nearby_with_this_link"], x["targets_linked_to_this_trait"]) for x in d4["candidates"]] == [(1, "body height", 1, 4), (2, "arrhythmia", 1, 2)]
           and d4["candidates"][1]["other_wordings_of_this_trait"] == ["heart disease"] and d4["answer"].startswith("DDD has no link as of 2021. The 2 highest-ranked traits")
           and "4.0% of traits listed at ranks 1-5 for a target with no link" in d4["candidates"][0]["how_often_this_held"])
        wd = G.why("DDD", "arrhythmia")
        ck("why gives a target with no link its place among such targets, and the share that held",
           wd["rank_among_targets_with_no_link"] == 1 and wd["targets_with_no_link_listed_for_this_trait"] == 1 and "ranks 1 of the 1 targets with no link" in wd["answer"]
           and "listed at ranks 1-5 for a trait" in wd["calibration"])
        ck("status says how targets with no link are ranked", G.status()["ranking"]["targets_with_no_link"]["ranked_by"] == "position alone, nothing fitted")
        ck("status names the scorer and the test it was put to", G.status()["ranking"]["scorer"]["name"] == "combine[B1+B3+B5]:F1" and G.status()["ranking"]["tested_on_later_links"]["first_year"] == 2022)
        q = G.predict("arrhythmia", as_of=2018, unexplained_only=True)
        ck("unexplained-only drops every candidate an earlier link explains", q["candidates"] == [])
        pt = G.predict("BBB", as_of=2018)["candidates"]
        ck("predict for a target ranks the traits it does not have, with wordings of one trait as one lead",
           [x["id"] for x in pt] == ["MONDO:2"] and pt[0]["other_wordings_of_this_trait"] == ["heart disease"])
        ck("with no link and no linked gene nearby there is nothing to rank from",
           G.predict("AAB")["candidates"] == [] and "no gene within 500,000 bases of it has one, so there is nothing to rank from" in G.predict("AAB")["answer"])
        ck("ranked by target, the method says that reading was checked after the test", "not declared before it" in G.predict("BBB", as_of=2018)["method"] and "score_parts" in pt[0])

        # a bundle that names no scorer is ranked by neighbours alone, as before 11 Oct 2026
        with tempfile.TemporaryDirectory() as tmp0:
            tiny(Path(tmp0), scorer=False)
            G0 = ask.Graph(Path(tmp0))
            p0 = G0.predict("arrhythmia", as_of=2018)
            ck("with no record of how position rankings held, a target with no link is not ranked",
               G0.predict("DDD")["candidates"] == [] and G0.predict("DDD")["answer"].startswith("DDD has no link as of 2021, so there is nothing to rank from")
               and "does not record" in G0.predict("arrhythmia", unlinked=True)["answer"] and "rank_among_targets_with_no_link" not in G0.why("DDD", "arrhythmia"))
            ck("and a share from the earlier backtest says it is counted over the traits that had a later link",
               "counted over the traits that had a later link" in p0["candidates"][0]["how_often_this_held"])
            ck("with no scorer recorded the ranking is shared-trait neighbours, and says so",
               G0.weights is None and "No model is fitted" in p0["method"] and "score_parts" not in p0["candidates"][0] and p0["candidates"][0]["score"] == 0.5
               and "In the 2018 backtest" in p0["candidates"][0]["how_often_this_held"] and "ranking" not in G0.status()
               and "by shared-trait neighbours BBB ranks" in G0.why("BBB", "arrhythmia", as_of=2018)["answer"])

        # leads
        near = {1: {2}, 2: {1, 3}, 3: {2}, 7: set(), 8: set(), 9: set()}
        L = ask.group_leads([(1, 1), (2, 7), (3, 2), (4, 3), (5, 8), (6, 9)], lambda i: near[i], lambda i: "s8" if i in (8, 9) else f"s{i}", top=3)
        ck("neighbours chain into one lead, at the rank of its best member", [(x["rank"], x["first"], sorted(x["members"])) for x in L] == [(1, 1, [1, 2, 3]), (2, 7, [7]), (5, 8, [8, 9])])
        ck("a lead lists the other labels in it once, and items sharing a label are one lead", L[0]["also"] == ["s2", "s3"] and L[2]["also"] == [])
        ck("grouping stops when another lead would exceed the number asked for", len(ask.group_leads([(1, 1), (2, 7), (3, 8)], lambda i: set(), str, top=2)) == 2)

        # whats_new, status
        n = G.whats_new("AAA", 2015)
        ck("whats_new counts links first reported after the year, by year", n["new_links"] == 1 and n["by_year"] == {2020: 1})
        ck("status says what the graph is and what it leaves out", "11 GWAS links" in G.status()["answer"] and len(G.status()["left_out"]) == 3)
        try:
            G.known("AAA", as_of=1990)
            ck("a year before any study is refused", False)
        except ValueError:
            ck("a year before any study is refused", True)

        # the bundle is checked
        (d / "links.parquet").write_bytes((d / "links.parquet").read_bytes() + b" ")
        try:
            ask.Graph(d)
            ck("a bundle file that is not the one its manifest describes is refused", False)
        except SystemExit:
            ck("a bundle file that is not the one its manifest describes is refused", True)
    return 0 if ok else 1
