PYTHON ?= python

.PHONY: validate test figures figure-check

validate:
	$(PYTHON) -m compileall -q scripts tests figures
	bash -n scripts/run_recovery.sh
	$(PYTHON) -m pytest -q

test:
	$(PYTHON) -m pytest -q

figures:
	cd figures && $(PYTHON) run_all.py

figure-check:
	cd figures && $(PYTHON) qa_check.py
