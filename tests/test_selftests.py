"""Each module's own self-test also runs under pytest, so one command covers both."""
from pharos_gap_kg import ask_selftest, baselines, bundle, calibrate, context, controls, dated, usability, graph, levels, ontology, paths, position, pubyear, score, snapshots, split, sqldump


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
