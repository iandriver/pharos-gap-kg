PY ?= $(firstword $(wildcard $(CURDIR)/.venv/bin/python) python3)

.PHONY: check paths selftest test model
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
	$(PY) -m pharos_gap_kg.context --selftest
	$(PY) -m pharos_gap_kg.bundle --selftest
	$(PY) -m pharos_gap_kg.calibrate --selftest
	$(PY) -m pharos_gap_kg.ask --selftest
	$(PY) -m pharos_gap_kg.usability --selftest
	$(PY) -m pharos_gap_kg.ctxtables --selftest
	$(PY) -m pharos_gap_kg.enrich --selftest
	$(PY) -m pharos_gap_kg.planted --selftest
	$(PY) -m pharos_gap_kg.ctxgraph --selftest
	$(PY) -m pharos_gap_kg.evidence --selftest
	$(PY) -m pharos_gap_kg.litcheck --selftest
	$(PY) -m pharos_gap_kg.prune2 --selftest
	$(PY) -m pharos_gap_kg.applycut --selftest
	$(PY) -m pharos_gap_kg.fulltext --selftest
	$(PY) -m pharos_gap_kg.tgaps --selftest
	$(PY) -m pharos_gap_kg.tgapgraph --selftest
	$(PY) -m pharos_gap_kg.combine --selftest
	$(PY) -m pharos_gap_kg.rehearse --selftest
	$(PY) -m pharos_gap_kg.coldstart --selftest
	$(PY) -m pharos_gap_kg.gwascat --selftest

test:
	$(PY) -m pytest -q

# the fitted models need the [model] extra (torch), so their self-test is not part of `check`
model:
	$(PY) -m pharos_gap_kg.model --selftest
