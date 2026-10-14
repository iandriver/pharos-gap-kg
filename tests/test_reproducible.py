"""What a rebuild from a clean clone depends on, checked without the data where possible.

The PubMed and NCBI answers the graph was built from are committed under `recorded/`. They must be the
files their manifest names, and they must be the answers the committed record of the lift cut was
computed from. A sample drawn again that is not the committed one must be refused, never recorded."""
import json

import pandas as pd
import pytest

from pharos_gap_kg import paths, recorded

R = paths.RECORDED / "pharos319"


def test_the_recorded_answers_are_the_files_their_manifest_names():
    man = json.loads((R / "manifest.json").read_text())["files"]
    assert set(man) == set(recorded.places("pharos319"))
    assert all(recorded.sha(R / n) == h for n, h in man.items())


def test_the_recorded_answers_are_the_ones_the_lift_cut_was_computed_from():
    rec = json.loads((paths.DOCS_DATA / "pruning2_applied_pharos319.json").read_text())
    applied = pd.read_parquet(R / "applied.parquet")
    assert json.loads((R / "applied.json").read_text()) == rec["pubmed"]
    by = applied[applied.lit_withheld.astype(bool)].groupby("source").size().to_dict()
    assert by == {s: v["withheld_by_the_cut"] for s, v in rec["sources"].items() if v.get("withheld_by_the_cut")}
    sample = json.loads((paths.DOCS_DATA / "pruning2_sample_pharos319.json").read_text())["files"]["sample.parquet"]
    assert json.loads((R / "search.json").read_text())["sample_sha256"] == sample


def test_a_sample_drawn_again_that_differs_is_refused_and_not_recorded(tmp_path, monkeypatch):
    from pharos_gap_kg import prune2
    if not (paths.context("pharos319") / "disease_links.parquet").is_file():
        pytest.skip("the context tables are not built")
    rec = paths.DOCS_DATA / "pruning2_sample_pharos319.json"
    j = json.loads(rec.read_text())
    j["files"]["sample.parquet"] = "0" * 64
    (tmp_path / rec.name).write_text(json.dumps(j))
    monkeypatch.setattr(paths, "DOCS_DATA", tmp_path)
    with pytest.raises(SystemExit, match="REFUSED: the sample just drawn"):
        prune2.freeze("pharos319")
    assert json.loads((tmp_path / rec.name).read_text())["files"]["sample.parquet"] == "0" * 64
