#!/usr/bin/env python3
"""Standalone collision / layout QA for the figure scripts.

The core detector is figstyle.qa_report(fig), the render-then-verify gate every
figure script calls before saving. It inspects a rendered Matplotlib Figure and
flags four classes of layout defect that a human reviewer typically catches only
by eye:

  OOB          text drawn past the figure canvas bounds (clipped / cut off)
  CROSS_PANEL  text from one axes overlapping a *foreign* axes box (>15 px^2)
  TXT_TXT      two text objects overlapping within the same axes
  TXT_ON_LINE  text crossing a thin reference line (axvline/axhline). Note:
               a zero-thickness reference line has window width/height == 0, for
               which Matplotlib's Bbox.overlaps() returns False, so the detector
               pads the thin dimension +/-1.5 px before testing.

qa_report(fig) returns a list of problem dicts; an empty list == clean.

Automated QA is necessary but NOT sufficient: the flagship figure needed three
rounds of human-caught fixes on top of a passing automated gate. Always also do
a per-panel perceptual crop review (set QA_CROPS=1 when running plot_figure2.py
to emit per-panel crops, and read the rendered PNG at full size).

Usage:
    python qa_check.py                 # rebuild every figure and QA-gate it
    python qa_check.py plot_figure3.py # QA-gate one script's figure
"""
import os, sys, importlib.util, runpy

HERE = os.path.dirname(os.path.abspath(__file__))

def _load_figstyle():
    spec = importlib.util.spec_from_file_location("figstyle", os.path.join(HERE, "figstyle.py"))
    fs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fs)
    return fs

def check_script(script, fs):
    """Run a plot script and QA every figure it leaves open."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.close("all")
    runpy.run_path(os.path.join(HERE, script), run_name="__qa__")
    figure_numbers = plt.get_fignums()
    if not figure_numbers:
        raise RuntimeError(f"{script} left no rendered figure for layout QA")
    problems = {}
    for num in figure_numbers:
        fig = plt.figure(num)
        probs = fs.qa_report(fig)
        if probs:
            problems[f"{script}#fig{num}"] = probs
    return problems

def main():
    fs = _load_figstyle()
    scripts = sys.argv[1:] or [
        "plot_figure1.py", "plot_figure2.py", "plot_figure3.py", "plot_figure4.py",
        "plot_figure5.py", "plot_supplementary_s1.py", "plot_supplementary_s2.py",
    ]
    all_clean = True
    for s in scripts:
        try:
            probs = check_script(s, fs)
        except Exception as e:
            print(f"ERROR {s}: {e}"); all_clean = False; continue
        if probs:
            all_clean = False
            print(f"FAIL {s}: {sum(len(v) for v in probs.values())} problem(s)")
            for k, v in probs.items():
                for p in v:
                    print(f"     {k}: {p}")
        else:
            print(f"OK   {s}: 0 layout problems")
    print("\n" + ("ALL FIGURES QA-CLEAN" if all_clean else "QA PROBLEMS FOUND"))
    sys.exit(0 if all_clean else 1)

if __name__ == "__main__":
    main()
