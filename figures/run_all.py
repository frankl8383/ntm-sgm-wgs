#!/usr/bin/env python3
"""Regenerate all seven redesigned figures from packaged source data.

Usage:
    python run_all.py            # regenerate all figures (PDF + PNG) in place
    QA_CROPS=1 python run_all.py # also emit per-panel crops for Figure 2

Each plot_*.py is self-contained: it chdir's to its own directory, adds it to
sys.path, imports figstyle.py, reads its inputs from source_data/ and
phylogeny/, runs the figstyle collision QA gate, and writes Figure_N.pdf/.png.
"""
import subprocess, sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [
    "plot_figure1.py", "plot_figure2.py", "plot_figure3.py", "plot_figure4.py",
    "plot_figure5.py", "plot_supplementary_s1.py", "plot_supplementary_s2.py",
]
def main():
    fails = []
    for s in SCRIPTS:
        t0 = time.time()
        r = subprocess.run([sys.executable, os.path.join(HERE, s)],
                           capture_output=True, text=True, cwd=HERE)
        dt = time.time() - t0
        ok = r.returncode == 0
        print(f"{'OK ' if ok else 'FAIL'} {s:28s} {dt:5.1f}s")
        if not ok:
            fails.append(s)
            print(r.stderr[-800:])
    print("\n" + ("ALL FIGURES REGENERATED" if not fails else f"FAILED: {fails}"))
    sys.exit(1 if fails else 0)
if __name__ == "__main__":
    main()
