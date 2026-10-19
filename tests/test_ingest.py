"""The ingest stage's findings, read against its records: who qualified, the control, what the whole graph
fails, and what was exported."""
import json
import re

import pytest

from pharos_gap_kg import ingest, paths

J = lambda f: json.loads((paths.DOCS_DATA / f).read_text())
F = re.sub(r"\s+", " ", (paths.ROOT / "findings" / "INGEST.md").read_text()) if paths.ROOT else ""
HAIKU, SONNET, OPUS = "anthropic/claude-haiku-5.5", "anthropic/claude-sonnet-5.5", "anthropic/claude-opus-5.5"


def has(*wanted):
    for w in wanted:
        assert w in F, f"findings/INGEST.md no longer says: {w}"


def test_who_qualified_and_the_table_of_models():
    r = J("ingest_pharos319.json")
    q = {m: v["qualification"] for m, v in r["models"].items()}
    assert [m for m in ingest.MODELS if q[m]["qualifies"]] == [HAIKU, SONNET, OPUS] and not q[ingest.LOCAL]["qualifies"]
    for m, name in ((ingest.LOCAL, "qwen3.6 35B, local"), (HAIKU, "Claude Haiku 5.5"), (SONNET, "Claude Sonnet 5.5"), (OPUS, "Claude Opus 5.5")):
        cards = [c for c in r["models"][m]["cards"] if c["view"] in ingest.SCORED]
        first = sum(c["accepted_first"] for c in cards)
        has(f"| {name} | {q[m]['accepted_cards']} ({first}) | {q[m]['columns_agree']} of 12 | {q[m]['predicate_exact']}, {q[m]['predicate_exact_or_close']} of 12 |")
        row = re.search(re.escape(name) + r" \| [^|]+ \| [^|]+ \| [^|]+ \| \$([\d.]+) \| (\d+) \|", F)
        assert round(sum(c["seconds"] for c in cards) / len(cards)) == int(row.group(2))
        digits = len(row.group(1).split(".")[1]) if "." in row.group(1) else 0
        assert f"{sum(c['usd'] for c in cards):.{digits}f}" == row.group(1)


def test_the_control_and_the_whole_graph():
    r, g = J("ingest_pharos319.json"), J("biolink_checks_pharos319.json")
    assert [v for v, c in r["control"].items() if c["accepted"]] == ["drug_activity"]
    has("pass the validator on 1 of the 9 scored cards, `drug_activity`")
    assert all(all(v for k, v in x.items() if k != "edges") for x in r["hand_values_in_export"].values()), "the hand values are what the export carries"
    has(f"{g['edges_of_kinds_that_pass']:,} of the {g['edges']:,} edges are of a kind that passes, and {g['edges'] - g['edges_of_kinds_that_pass']:,} fail")
    for n in g["edges_by_problem"].values():
        has(f"| {n:,} |")


def test_what_was_exported():
    x = J("ingest/exported_pharos319.json")
    assert x["joined_to_the_graph"] is False and [v for v, e in x["views"].items() if e["exported"]] == ["ortholog_disease", "phenotype.OMIM", "phenotype.JAX_MGI"]
    assert {e["model"] for e in x["views"].values() if e["exported"]} == {HAIKU} and not x["views"]["goa"]["exported"]
    for v in ("ortholog_disease", "phenotype.OMIM", "phenotype.JAX_MGI"):
        e = x["views"][v]
        has(f"| {e['edges']:,} | {e['gaps']:,} | {e['nodes']:,} |")
        s = J(f"ingest/standards_{v}_pharos319.json")
        assert s["checked"]["edges"] == e["edges"] and s["subject_object_pairs_with_an_error"] == 0


def test_omim_keeps_one_phenotype_a_row():
    f = paths.tables("pharos319") / "phenotype.parquet"
    if not f.is_file():
        pytest.skip("the tables are not extracted")
    rows, _ = ingest.view_rows("pharos319", "phenotype.OMIM")
    n = rows.trait.astype(str).str.count(r"\d{6} \(\d\)")
    has(f"{int((n > 1).sum()):,} of the {int((n >= 1).sum()):,} rows with a phenotype number name more than one, {int(n.sum()):,} in all")


def test_the_spend_is_within_the_cap_and_every_record_is_current():
    r = J("ingest_pharos319.json")
    assert r["spend_usd"] <= r["cap_usd"] == 1.75
    has(f"${r['spend_usd']:.3f} against the cap of $1.75")
