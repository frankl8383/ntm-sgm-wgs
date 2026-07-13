PYTHON ?= python

.PHONY: validate test

validate:
	$(PYTHON) -m compileall -q scripts tests
	bash -n scripts/run_recovery.sh
	$(PYTHON) -m pytest -q

test:
	$(PYTHON) -m pytest -q
