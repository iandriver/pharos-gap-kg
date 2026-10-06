PY ?= $(firstword $(wildcard $(CURDIR)/.venv/bin/python) python3)

.PHONY: check paths selftest test
check: paths selftest test

paths:
	$(PY) -m pharos_gaps.paths --selftest

selftest:
	$(PY) -m pharos_gaps.sqldump --selftest
	$(PY) -m pharos_gaps.snapshots --selftest

test:
	$(PY) -m pytest -q
