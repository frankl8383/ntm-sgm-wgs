#!/usr/bin/env python3
"""Integrate read-level taxonomy and assembly QC into an initial genome-level judgement."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "sample_id",
    "top_species",
    "read_level_taxonomy_status",
    "assembly_size",
    "gc_content",
    "num_contigs",
    "n50",
    "checkm2_completeness",
    "checkm2_contamination",
    "gunc_pass_or_fail",
    "gunc_taxonomic_level",
    "gunc_contamination_portion",
    "gunc_clade_separation_score",
    "gunc_reference_representation_score",
    "assembly_qc_status",
    "initial_genome_level_status",
    "reason",
    "recommended_next_step",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read-judgement", required=True)
    parser.add_argument("--assembly-qc", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["sample_id"]: row for row in reader if row.get("sample_id")}


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def as_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def classify(
    read_status: str,
    assembly_status: str,
    size: int,
    gc: float,
    contigs: int,
    n50: int,
    checkm2_completeness: float,
    checkm2_contamination: float,
    gunc_status: str,
) -> tuple[str, str, str]:
    read_non_mycobacterium = read_status.startswith("dominant_non_mycobacterium")
    read_low_purity = read_status == "possible_low_purity_mycobacterium"
    read_mycobacterium = read_status.startswith("mycobacterium_dominant")
    assembly_severe = size > 8_000_000 or size < 4_000_000 or contigs > 2000 or n50 < 20_000
    checkm2_low_quality = 0 < checkm2_completeness < 90
    checkm2_high_contamination = checkm2_contamination > 10
    checkm2_warning_contamination = 5 < checkm2_contamination <= 10
    gunc_failed = gunc_status.lower() in {"false", "fail", "failed"}
    genome_qc_severe = assembly_severe or checkm2_low_quality or checkm2_high_contamination or gunc_failed
    genome_qc_warning = checkm2_warning_contamination
    assembly_warning = assembly_status != "pass_basic_assembly_thresholds"

    if read_non_mycobacterium:
        if genome_qc_severe:
            return (
                "likely_non_ntm_or_mixed_contaminated",
                "Read-level taxonomy is dominated by non-Mycobacterium and genome-level QC is abnormal.",
                "Exclude from SGM main analysis unless sample identity is rescued by CheckM2/GUNC/ANI.",
            )
        return (
            "likely_non_ntm_read_level",
            "Read-level taxonomy is dominated by non-Mycobacterium.",
            "Do not include in SGM main analysis; retain for contamination/misidentification reporting.",
        )

    if read_low_purity:
        return (
            "possible_mixed_or_low_purity_ntm",
            "Read-level Mycobacterium signal is below the NTM-candidate threshold and/or genome-level QC is suspicious.",
            "Run CheckM2/GUNC and ANI, then decide whether to exclude or resequence/re-isolate.",
        )

    if read_mycobacterium and genome_qc_severe:
        return (
            "mycobacterium_reads_but_assembly_suspicious",
            "Reads are Mycobacterium-dominant, but assembly/CheckM2/GUNC evidence suggests mixture, contamination, or low-quality assembly.",
            "Prioritize CheckM2/GUNC and inspect ANI against type/reference genomes before downstream inclusion.",
        )

    if read_mycobacterium and (assembly_warning or genome_qc_warning):
        return (
            "mycobacterium_candidate_with_assembly_warning",
            "Reads are Mycobacterium-dominant but genome-level QC has warnings.",
            "Keep as NTM candidate; use CheckM2/GUNC, ANI, and marker genes for final inclusion.",
        )

    if read_mycobacterium:
        return (
            "mycobacterium_candidate_basic_qc_pass",
            "Reads are Mycobacterium-dominant and basic assembly QC passes configured thresholds.",
            "Proceed to ANI, NTM-Profiler, marker genes, and CheckM2/GUNC.",
        )

    return (
        "unresolved_initial_status",
        "Read and assembly evidence did not match a predefined initial category.",
        "Inspect manually before downstream inclusion.",
    )


def main() -> int:
    args = parse_args()
    reads = read_tsv(Path(args.read_judgement))
    assemblies = read_tsv(Path(args.assembly_qc))
    output_rows = []
    for sample_id, assembly in assemblies.items():
        read = reads.get(sample_id, {})
        size = as_int(assembly.get("assembly_size", "0"))
        gc = as_float(assembly.get("gc_content", "0"))
        contigs = as_int(assembly.get("num_contigs", "0"))
        n50 = as_int(assembly.get("n50", "0"))
        checkm2_completeness = as_float(assembly.get("checkm2_completeness", "0"))
        checkm2_contamination = as_float(assembly.get("checkm2_contamination", "0"))
        gunc_status = assembly.get("gunc_pass_or_fail", "NA")
        read_status = read.get("read_level_taxonomy_status", "NA")
        assembly_status = assembly.get("assembly_qc_status", "NA")
        status, reason, next_step = classify(
            read_status,
            assembly_status,
            size,
            gc,
            contigs,
            n50,
            checkm2_completeness,
            checkm2_contamination,
            gunc_status,
        )
        output_rows.append(
            {
                "sample_id": sample_id,
                "top_species": read.get("top_species", "NA"),
                "read_level_taxonomy_status": read_status,
                "assembly_size": assembly.get("assembly_size", "NA"),
                "gc_content": assembly.get("gc_content", "NA"),
                "num_contigs": assembly.get("num_contigs", "NA"),
                "n50": assembly.get("n50", "NA"),
                "checkm2_completeness": assembly.get("checkm2_completeness", "NA"),
                "checkm2_contamination": assembly.get("checkm2_contamination", "NA"),
                "gunc_pass_or_fail": gunc_status,
                "gunc_taxonomic_level": assembly.get("gunc_taxonomic_level", "NA"),
                "gunc_contamination_portion": assembly.get("gunc_contamination_portion", "NA"),
                "gunc_clade_separation_score": assembly.get("gunc_clade_separation_score", "NA"),
                "gunc_reference_representation_score": assembly.get(
                    "gunc_reference_representation_score", "NA"
                ),
                "assembly_qc_status": assembly_status,
                "initial_genome_level_status": status,
                "reason": reason,
                "recommended_next_step": next_step,
            }
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {len(output_rows)} rows: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
