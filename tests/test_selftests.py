"""Each module's own self-test also runs under pytest, so one command covers both."""
from pharos_gaps import paths, sqldump


def test_paths_selftest():
    assert paths._selftest() == 0


def test_sqldump_selftest():
    assert sqldump._selftest() == 0
