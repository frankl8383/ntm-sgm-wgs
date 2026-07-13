#!/usr/bin/env python3
"""Regenerate and QA-gate all seven figures from packaged source data.

Usage:
    python run_all.py            # regenerate all figures (PDF + PNG) in place
    QA_CROPS=1 python run_all.py # also emit per-panel crops for Figure 2

The standalone QA runner executes every plotting script, checks each rendered
figure with figstyle.qa_report(), and exits nonzero on a collision or script
error. Plot scripts still write their PDF and PNG outputs in place.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    result = subprocess.run(
        [sys.executable, os.path.join(HERE, "qa_check.py")],
        cwd=HERE,
    )
    if result.returncode == 0:
        print("ALL FIGURES REGENERATED AND QA-CLEAN")
    else:
        print("FIGURE REGENERATION OR LAYOUT QA FAILED", file=sys.stderr)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
