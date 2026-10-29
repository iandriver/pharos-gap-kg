"""Proxy cases in the study-size table: the findings read against the record of what the rule changed."""
import json
import re

import pytest

from pharos_gap_kg import paths, studysize

REC = paths.DOCS_DATA / "proxy_floors_pharos319.json"
if not REC.is_file() or paths.ROOT is None:
    pytest.skip("no record of the proxy rule", allow_module_level=True)
D = json.loads(REC.read_text())
F = re.sub(r"\s+", " ", (paths.ROOT / "findings" / "DGAPS.md").read_text().split("## Added 29 Oct 2026: proxy cases")[1])
S = json.loads((paths.DOCS_DATA / "study_sizes_pharos319.json").read_text())


def has(*wanted):
    for w in wanted:
        assert w in F, f"findings/DGAPS.md (29 Oct) no longer says: {w}"


def test_the_studies_with_proxies():
    px = D["catalog_studies_with_proxies"]
    mine = sorted(k for k, v in px.items() if v["in_the_release"])
    assert mine == sorted(S["proxies"]["release_studies_with_proxies"]) and len(mine) == D["release_studies_with_proxies"]
    has(f"Of the {S['catalog_studies']:,} studies in the Catalog file, {len(px)} have a line with proxies, {len(mine)} of them in the release",
        ", ".join(mine[:-1]) + " and " + mine[-1], f"The {D['catalog_studies_with_no_proxy_words_read_as_before']:,} studies whose line holds no proxy words are read as before")
    g = px["GCST006035"]
    assert not g["in_the_release"] and g["proxy_cases"] == 967
    has(f"reads {g['proxy_cases']} of these cases as proxy cases")


def test_what_changed():
    t = {r["trait"]: r for r in D["traits"]}
    assert set(t) == {"MONDO:0004975", "EFO:0009268"} and D["traits_changed"] == 2
    has(f"The floor changes for {D['traits_changed']} of the {D['traits_with_a_sized_study']:,} traits")
    ad = t["MONDO:0004975"]["after"]
    assert (ad["study"], ad["smallest_odds_ratio_it_could_find"], ad["smallest_odds_ratio_with_proxies_as_cases"]) == ("GCST007320", 1.077, 1.046)
    assert ad["smallest_odds_ratio_with_proxies_as_cases"] == studysize.smallest_effect(24087, 383378, 0, 47793)["smallest_odds_ratio_with_proxies_as_cases"]
    has(f"floor of {ad['smallest_odds_ratio_it_could_find']}", f"with its {ad['proxy_cases']:,} proxy cases counted as cases, {ad['smallest_odds_ratio_with_proxies_as_cases']}")
    fh = t["EFO:0009268"]
    b, a, e = fh["before"], fh["after"], fh["earlier_largest_study_under_the_rule"]
    assert fh["kind"] == "the largest study changed" and e["design"] == "proxies only"
    has(f"Its largest study was {b['study']}, read before as {b['individuals']:,} individuals", f"{a['study']} ({a['cases']:,} cases, {a['controls']:,} controls)", f"floor of {b['smallest_effect_it_could_find_in_standard_deviations_per_allele']} standard deviations",
        f"{e['proxy_cases']:,} people with a parental history and {e['proxy_controls']:,} without", f"now {a['study']}", f"with a floor of {a['smallest_odds_ratio_it_could_find']}",
        f"{e['study']}'s floor is {e['smallest_odds_ratio_with_proxies_as_cases']:.3f}")
    assert e["proxy_cases"] + e["proxy_controls"] == b["individuals"]
    for s in ("GCST005922", "GCST90104087"):
        assert all(r["before"]["study"] != s for r in D["traits"])
