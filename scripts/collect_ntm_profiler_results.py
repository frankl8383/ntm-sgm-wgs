#!/usr/bin/env python3
"""Collect NTM-Profiler JSON species calls into a TSV table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def collect_one(path: Path) -> dict[str, object]:
    data = json.loads(path.read_text())
    sample_id = data.get("id") or path.name.split(".")[0]
    taxa = data.get("taxa") or []
    top = taxa[0] if taxa else {}
    db = (data.get("pipeline") or {}).get("species_db_version") or {}
    return {
        "sample_id": sample_id,
        "ntm_profiler_species": top.get("species", "NA"),
        "ntm_profiler_ani": top.get("ani", "NA"),
        "ntm_profiler_relative_abundance": top.get("relative_abundance", "NA"),
        "ntm_profiler_accession": top.get("accession", "NA"),
        "ntm_profiler_ncbi_organism_name": top.get("ncbi_organism_name", "NA"),
        "ntm_profiler_prediction_method": top.get("prediction_method", "NA"),
        "ntm_profiler_qc_fail_taxa_count": len(data.get("qc_fail_taxa") or []),
        "ntm_profiler_db_name": db.get("name", "NA"),
        "ntm_profiler_db_commit": db.get("commit", "NA"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob("*/*.results.json"))
    rows = [collect_one(path) for path in paths]
    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(
            columns=[
                "sample_id",
                "ntm_profiler_species",
                "ntm_profiler_ani",
                "ntm_profiler_relative_abundance",
                "ntm_profiler_accession",
                "ntm_profiler_ncbi_organism_name",
                "ntm_profiler_prediction_method",
                "ntm_profiler_qc_fail_taxa_count",
                "ntm_profiler_db_name",
                "ntm_profiler_db_commit",
            ]
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
