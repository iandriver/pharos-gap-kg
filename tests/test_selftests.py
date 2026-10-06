"""Each module's own self-test also runs under pytest, so one command covers both."""
from pharos_gaps import baselines, dated, graph, levels, paths, pubyear, score, snapshots, split, sqldump


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
