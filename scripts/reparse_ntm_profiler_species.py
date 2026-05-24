#!/usr/bin/env python3
"""Reparse NTM-Profiler species JSON with explicit accession-name warnings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def has_name_conflict(predicted: str, organism: str) -> bool:
    predicted = (predicted or "NA").strip()
    organism = (organism or "NA").strip()
    if predicted in {"", "NA"} or organism in {"", "NA"}:
        return False
    return predicted.lower() not in organism.lower()


def collect_one(path: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    data = json.loads(path.read_text())
    sample_id = data.get("id") or path.name.split(".")[0]
    pipeline = data.get("pipeline") or {}
    species_db = pipeline.get("species_db_version") or {}
    taxa = data.get("taxa") or []
    qc_fail_taxa = data.get("qc_fail_taxa") or []
    qc = data.get("qc") or {}
    top = taxa[0] if taxa else {}
    predicted = top.get("species", "NA")
    organism = top.get("ncbi_organism_name", "NA")
    summary = {
        "sample_id": sample_id,
        "ntm_profiler_predicted_species": predicted,
        "ntm_profiler_top_ani": top.get("ani", "NA"),
        "ntm_profiler_top_relative_abundance": top.get("relative_abundance", "NA"),
        "ntm_profiler_closest_accession": top.get("accession", "NA"),
        "ntm_profiler_closest_accession_organism_name": organism,
        "ntm_profiler_prediction_method": top.get("prediction_method", "NA"),
        "ntm_profiler_species_vs_accession_name_conflict": has_name_conflict(predicted, organism),
        "ntm_profiler_notes": ";".join(top.get("notes") or []) if top else "NA",
        "ntm_profiler_qc_fail_taxa_count": len(qc_fail_taxa),
        "ntm_profiler_num_sequences": qc.get("num_sequences", "NA"),
        "ntm_profiler_num_bases": qc.get("num_bases", "NA"),
        "ntm_profiler_n50": qc.get("n50", "NA"),
        "ntm_profiler_db_name": species_db.get("name", "NA"),
        "ntm_profiler_db_commit": species_db.get("commit", "NA"),
        "ntm_profiler_db_date": species_db.get("date", "NA"),
    }
    long_rows: list[dict[str, object]] = []
    for rank, taxon in enumerate(taxa, start=1):
        pred = taxon.get("species", "NA")
        org = taxon.get("ncbi_organism_name", "NA")
        long_rows.append(
            {
                "sample_id": sample_id,
                "rank": rank,
                "predicted_species": pred,
                "ani": taxon.get("ani", "NA"),
                "relative_abundance": taxon.get("relative_abundance", "NA"),
                "abundance": taxon.get("abundance", "NA"),
                "coverage": taxon.get("coverage", "NA"),
                "closest_accession": taxon.get("accession", "NA"),
                "closest_accession_organism_name": org,
                "prediction_method": taxon.get("prediction_method", "NA"),
                "species_vs_accession_name_conflict": has_name_conflict(pred, org),
                "notes": ";".join(taxon.get("notes") or []),
            }
        )
    return summary, long_rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    parser.add_argument("--long-output", required=True, type=Path)
    args = parser.parse_args()

    summaries = []
    long_rows = []
    for path in sorted(args.input_dir.glob("*/*.results.json")):
        summary, rows = collect_one(path)
        summaries.append(summary)
        long_rows.extend(rows)

    summary_df = pd.DataFrame(summaries)
    long_df = pd.DataFrame(long_rows)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(args.summary_output, sep="\t", index=False)
    long_df.to_csv(args.long_output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
