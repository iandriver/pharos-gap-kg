"""Each module's own self-test also runs under pytest, so one command covers both."""
from pharos_gap_kg import applycut, ask_selftest, baselines, bundle, calibrate, context, controls, ctxgraph, ctxtables, dated, enrich, evidence, litcheck, planted, prune2, usability, graph, levels, ontology, paths, position, pubyear, score, snapshots, split, sqldump


def test_paths_selftest():
    assert paths._selftest() == 0


def test_sqldump_selftest():
    assert sqldump._selftest() == 0


def test_snapshots_selftest():
    assert snapshots._selftest() == 0


def test_levels_selftest():
    assert levels._selftest() == 0


def test_pubyear_selftest():
    assert pubyear._selftest() == 0


def test_dated_selftest():
    assert dated._selftest() == 0


def test_split_selftest():
    assert split._selftest() == 0


def test_graph_selftest():
    assert graph._selftest() == 0


def test_score_selftest():
    assert score._selftest() == 0


def test_baselines_selftest():
    assert baselines._selftest() == 0


def test_ontology_selftest():
    assert ontology._selftest() == 0


def test_position_selftest():
    assert position._selftest() == 0


def test_controls_selftest():
    assert controls._selftest() == 0


def test_context_selftest():
    assert context._selftest() == 0


def test_bundle_selftest():
    assert bundle._selftest() == 0


def test_calibrate_selftest():
    assert calibrate._selftest() == 0


def test_ask_selftest():
    assert ask_selftest.run() == 0


def test_usability_selftest():
    assert usability._selftest() == 0


def test_ctxtables_selftest():
    assert ctxtables._selftest() == 0


def test_enrich_selftest():
    assert enrich._selftest() == 0


def test_planted_selftest():
    assert planted._selftest() == 0


def test_ctxgraph_selftest():
    assert ctxgraph._selftest() == 0


def test_evidence_selftest():
    assert evidence._selftest() == 0


def test_litcheck_selftest():
    assert litcheck._selftest() == 0


def test_prune2_selftest():
    assert prune2._selftest() == 0


def test_applycut_selftest():
    assert applycut._selftest() == 0


def test_fulltext_selftest():
    from pharos_gap_kg import fulltext
    assert fulltext._selftest() == 0


def test_model_selftest():
    import pytest
    pytest.importorskip("torch")
    from pharos_gap_kg import model
    assert model._selftest() == 0
