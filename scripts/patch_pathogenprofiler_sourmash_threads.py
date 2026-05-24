#!/usr/bin/env python3
"""Patch pathogenprofiler sourmash taxonomy compatibility in local envs.

NTM-Profiler 0.8.1 calls pathogenprofiler taxonomy classification with a
``threads=`` keyword. Some pathogenprofiler releases accept this keyword for
Sylph but not for Sourmash, causing species-only FASTA profiling to fail before
any biological result is produced. This patch keeps Sourmash behavior unchanged
and only lets the method ignore unknown keyword arguments.
"""

from __future__ import annotations

import importlib
from pathlib import Path


def main() -> int:
    taxonomy = importlib.import_module("pathogenprofiler.taxonomy")
    path = Path(taxonomy.__file__)
    text = path.read_text()

    old = (
        "def classify(self, ref_db: str, "
        "intersect_bp: int=500000,f_match_threshold: float=0.1) "
        "-> List[TaxonomicHit]:"
    )
    new = (
        "def classify(self, ref_db: str, "
        "intersect_bp: int=500000,f_match_threshold: float=0.1, **kwargs) "
        "-> List[TaxonomicHit]:"
    )

    if new in text:
        print(f"already_patched\t{path}")
        return 0
    if old not in text:
        raise SystemExit(f"expected SourmashSig.classify signature not found in {path}")

    backup = path.with_suffix(path.suffix + ".bak_ntm_sgm")
    if not backup.exists():
        backup.write_text(text)
    path.write_text(text.replace(old, new, 1))
    print(f"patched\t{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
