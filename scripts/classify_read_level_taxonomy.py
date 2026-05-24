#!/usr/bin/env python3
"""Create a practical read-level taxonomy judgement table from Kraken/Bracken summary."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml


OUTPUT_FIELDS = [
    "sample_id",
    "top_genus",
    "top_species",
    "mycobacterium_fraction",
    "top_species_fraction",
    "second_species",
    "second_species_fraction",
    "read_level_taxonomy_status",
    "reason",
    "recommended_next_step",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, help="read_qc_taxonomy_summary.tsv")
    parser.add_argument("--thresholds", required=True, help="config/qc_thresholds.yaml")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def is_mycobacterium_name(name: str) -> bool:
    return name.lower().startswith("mycobacterium")


def classify(row: dict[str, str], min_mycobacterium: float, min_top_species: float, mixed_warning: float) -> tuple[str, str, str]:
    top_genus = row["top_genus"]
    top_species = row["top_species"]
    second_species = row["second_species"]
    myco = as_float(row["mycobacterium_fraction"])
    top_frac = as_float(row["top_species_fraction"])
    second_frac = as_float(row["second_species_fraction"])

    if myco < min_mycobacterium:
        if top_genus == "Mycobacterium" or is_mycobacterium_name(top_species):
            return (
                "possible_low_purity_mycobacterium",
                f"Mycobacterium fraction {myco:.3f} is below {min_mycobacterium:.2f}, but the top genus/species is Mycobacterium.",
                "Keep for assembly QC, but treat as possible mixed or low-purity NTM until CheckM2/GUNC and ANI confirm.",
            )
        if myco >= 0.10:
            return (
                "dominant_non_mycobacterium_with_mycobacterium_signal",
                f"Top species is {top_species}; Mycobacterium fraction is {myco:.3f}.",
                "Do not include in SGM main analysis yet; verify culture purity and inspect assembly/binning evidence.",
            )
        return (
            "dominant_non_mycobacterium",
            f"Top species is {top_species}; Mycobacterium fraction is only {myco:.3f}.",
            "Treat as likely non-NTM or strong contamination at read level; verify sample identity before downstream NTM analysis.",
        )

    if second_frac >= mixed_warning or top_frac < min_top_species:
        if second_species != "NA" and not is_mycobacterium_name(second_species):
            return (
                "mycobacterium_dominant_with_non_mycobacterium_signal",
                f"Mycobacterium fraction {myco:.3f}, but second species {second_species} is {second_frac:.3f}.",
                "Inspect assembly QC and contamination tools carefully before inclusion.",
            )
        return (
            "mycobacterium_dominant_species_level_ambiguous",
            f"Mycobacterium fraction {myco:.3f}, but top/second species fractions are {top_frac:.3f}/{second_frac:.3f}.",
            "Keep as NTM candidate; resolve species with assembly-level ANI, NTM-Profiler, and marker genes.",
        )

    return (
        "mycobacterium_dominant_single_species_support",
        f"Mycobacterium fraction {myco:.3f}; top species fraction {top_frac:.3f}; second species fraction {second_frac:.3f}.",
        "Proceed to assembly QC and multi-evidence species calling.",
    )


def main() -> int:
    args = parse_args()
    with Path(args.thresholds).open("r", encoding="utf-8") as handle:
        thresholds = yaml.safe_load(handle) or {}
    tax_thresholds = thresholds.get("taxonomy_reads", {})
    min_mycobacterium = float(tax_thresholds.get("min_mycobacterium_fraction_for_ntm_candidate", 0.70))
    min_top_species = float(tax_thresholds.get("min_top_species_fraction_for_single_species", 0.60))
    mixed_warning = float(tax_thresholds.get("mixed_species_warning_fraction", 0.10))

    with Path(args.summary).open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    output_rows = []
    for row in rows:
        status, reason, next_step = classify(row, min_mycobacterium, min_top_species, mixed_warning)
        output_rows.append(
            {
                "sample_id": row["sample_id"],
                "top_genus": row["top_genus"],
                "top_species": row["top_species"],
                "mycobacterium_fraction": row["mycobacterium_fraction"],
                "top_species_fraction": row["top_species_fraction"],
                "second_species": row["second_species"],
                "second_species_fraction": row["second_species_fraction"],
                "read_level_taxonomy_status": status,
                "reason": reason,
                "recommended_next_step": next_step,
            }
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {len(output_rows)} rows: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
