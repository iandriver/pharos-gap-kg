"""The self-test of `ask`: a six-node graph written to a temporary bundle, and every question put to it."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pandas as pd

from pharos_gap_kg import ask, graph


def tiny(d: Path) -> None:
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
    (d / "backtest.json").write_text(json.dumps(bt))
    files = ["nodes.parquet", "links.parquet", "evidence.parquet", "lineage.parquet", "neighbours.parquet", "backtest.json"]
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
        ck("an unknown argument stops the question", G.why("QQQ", "arrhythmia")["found"] is False and G.why("QQQ", "arrhythmia")["problem"] == "not in the graph")

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
        ck("each candidate quotes how often its rank has held, and that it is a lower bound", all("backtest" in x["how_often_this_held"] and "not known to be absent" in x["how_often_this_held"] for x in p["candidates"]))
        q = G.predict("arrhythmia", as_of=2018, unexplained_only=True)
        ck("unexplained-only drops every candidate an earlier link explains", q["candidates"] == [])
        ck("predict for a target ranks the traits it does not have", [x["id"] for x in G.predict("BBB", as_of=2018)["candidates"]] == ["MONDO:2", "EFO:9"])
        ck("with nothing linked there is nothing to rank from", G.predict("DDD")["candidates"] == [] and "nothing to rank" in G.predict("DDD")["answer"])

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
