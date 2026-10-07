PY ?= $(firstword $(wildcard $(CURDIR)/.venv/bin/python) python3)

.PHONY: check paths selftest test
check: paths selftest test

paths:
	$(PY) -m pharos_gap_kg.paths --selftest

selftest:
	$(PY) -m pharos_gap_kg.sqldump --selftest
	$(PY) -m pharos_gap_kg.snapshots --selftest
	$(PY) -m pharos_gap_kg.levels --selftest
	$(PY) -m pharos_gap_kg.pubyear --selftest
	$(PY) -m pharos_gap_kg.dated --selftest
	$(PY) -m pharos_gap_kg.split --selftest
	$(PY) -m pharos_gap_kg.graph --selftest
	$(PY) -m pharos_gap_kg.score --selftest
	$(PY) -m pharos_gap_kg.baselines --selftest
	$(PY) -m pharos_gap_kg.ontology --selftest
	$(PY) -m pharos_gap_kg.position --selftest
	$(PY) -m pharos_gap_kg.controls --selftest

test:
	$(PY) -m pytest -q
