PY ?= $(firstword $(wildcard $(CURDIR)/.venv/bin/python) python3)

.PHONY: check paths selftest test
check: paths selftest test

paths:
	$(PY) -m pharos_gaps.paths --selftest

selftest:
	$(PY) -m pharos_gaps.sqldump --selftest
	$(PY) -m pharos_gaps.snapshots --selftest
	$(PY) -m pharos_gaps.levels --selftest
	$(PY) -m pharos_gaps.pubyear --selftest
	$(PY) -m pharos_gaps.dated --selftest
	$(PY) -m pharos_gaps.split --selftest
	$(PY) -m pharos_gaps.graph --selftest
	$(PY) -m pharos_gaps.score --selftest
	$(PY) -m pharos_gaps.baselines --selftest
	$(PY) -m pharos_gaps.ontology --selftest
	$(PY) -m pharos_gaps.position --selftest
	$(PY) -m pharos_gaps.controls --selftest

test:
	$(PY) -m pytest -q
