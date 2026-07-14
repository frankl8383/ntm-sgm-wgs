PYTHON ?= python

.PHONY: validate test figures figure-check figure-check-fallback

validate:
	$(PYTHON) -m compileall -q scripts tests figures
	bash -n scripts/run_recovery.sh
	$(PYTHON) -m pytest -q

test:
	$(PYTHON) -m pytest -q

figures:
	cd figures && $(PYTHON) run_all.py
	$(MAKE) figure-check-fallback PYTHON="$(PYTHON)"

figure-check:
	cd figures && $(PYTHON) qa_check.py

figure-check-fallback:
	cd figures && NTM_FIGURE_FONT="DejaVu Sans" NTM_FIGURE_OUTPUT_DIR="/tmp/ntm_sgm_wgs_figure_fallback" $(PYTHON) qa_check.py
