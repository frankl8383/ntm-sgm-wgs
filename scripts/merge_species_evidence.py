#!/usr/bin/env python3
"""Merge read, assembly, FastANI, NTM-Profiler, and marker evidence."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd


MAC_SPECIES = {
    "Mycobacterium avium",
    "Mycobacterium intracellulare",
    "Mycobacterium paraintracellulare",
    "Mycobacterium colombiense",
}


def marker_majority(row: pd.Series) -> str:
    species = []
    for col in row.index:
        if col.endswith("_organism_name") and pd.notna(row[col]) and str(row[col]) != "NA":
            species.append(str(row[col]))
    if not species:
        return "NA"
    return Counter(species).most_common(1)[0][0]


def marker_species_set(row: pd.Series) -> str:
    species = sorted(
        {
            str(row[col])
            for col in row.index
            if col.endswith("_organism_name") and pd.notna(row[col]) and str(row[col]) != "NA"
        }
    )
    return ";".join(species) if species else "NA"


def same_complex(*species: str) -> bool:
    called = [s for s in species if s and s != "NA"]
    return bool(called) and all(s in MAC_SPECIES for s in called)


def growth_type(species: str) -> str:
    return "SGM" if species in MAC_SPECIES else "unknown"


def confidence(row: pd.Series) -> str:
    ani = float(row.get("fastani_best_hit_ani", 0) or 0)
    af = float(row.get("fastani_best_hit_alignment_fraction", 0) or 0)
    if ani >= 96 and af >= 0.5 and same_complex(
        row.get("fastani_best_hit_species", "NA"),
        row.get("ntm_profiler_species", "NA"),
        row.get("marker_majority_species", "NA"),
    ):
        return "High"
    if ani >= 95 and af >= 0.5:
        return "Moderate"
    return "Low"


def conflict_flag(row: pd.Series) -> str:
    final = row.get("fastani_best_hit_species", "NA")
    ntm = row.get("ntm_profiler_species", "NA")
    markers = row.get("marker_species_set", "NA")
    if final == ntm and final in markers:
        return "none_detected"
    if same_complex(final, ntm) and all(s in MAC_SPECIES for s in markers.split(";") if s and s != "NA"):
        return "MAC_species_granularity_conflict"
    return "major_species_conflict"


def downstream_reason(row: pd.Series) -> str:
    status = row.get("initial_genome_level_status", "NA")
    assembly_status = row.get("assembly_qc_status", "NA")
    gunc_status = str(row.get("gunc_pass_or_fail", "NA"))
    if gunc_status.lower() in {"false", "fail", "failed"}:
        return f"assembly_suspicious_GUNC_fail: {assembly_status}"
    if status == "mycobacterium_candidate_with_assembly_warning":
        return f"assembly_warning: {assembly_status}"
    if status == "mycobacterium_candidate_basic_qc_pass":
        return "NA"
    return f"not_priority_for_species_matrix: {status}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--read-assembly", required=True, type=Path)
    parser.add_argument("--ntm-profiler", required=True, type=Path)
    parser.add_argument("--fastani-best", required=True, type=Path)
    parser.add_argument("--marker-wide", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    read_assembly = pd.read_csv(args.read_assembly, sep="\t")
    ntm = pd.read_csv(args.ntm_profiler, sep="\t")
    fastani = pd.read_csv(args.fastani_best, sep="\t")
    marker = pd.read_csv(args.marker_wide, sep="\t")
    marker["marker_majority_species"] = marker.apply(marker_majority, axis=1)
    marker["marker_species_set"] = marker.apply(marker_species_set, axis=1)

    df = (
        read_assembly.merge(ntm, on="sample_id", how="inner")
        .merge(fastani, on="sample_id", how="inner")
        .merge(marker, on="sample_id", how="left")
    )
    df["final_wgs_species_call"] = df["fastani_best_hit_species"]
    df["final_species_complex_call"] = df["final_wgs_species_call"].map(
        lambda s: "Mycobacterium avium complex" if s in MAC_SPECIES else "NA"
    )
    df["final_growth_type_call"] = df["final_wgs_species_call"].map(growth_type)
    df["confidence_level"] = df.apply(confidence, axis=1)
    df["evidence_conflict_flag"] = df.apply(conflict_flag, axis=1)
    df["exclude_from_downstream_reason"] = df.apply(downstream_reason, axis=1)
    df["gunc_status"] = df.get("gunc_pass_or_fail", "NA")

    keep = [
        "sample_id",
        "top_species",
        "read_level_taxonomy_status",
        "assembly_size",
        "gc_content",
        "num_contigs",
        "n50",
        "checkm2_completeness",
        "checkm2_contamination",
        "gunc_status",
        "gunc_contamination_portion",
        "gunc_clade_separation_score",
        "gunc_reference_representation_score",
        "assembly_qc_status",
        "initial_genome_level_status",
        "ntm_profiler_species",
        "ntm_profiler_ani",
        "ntm_profiler_accession",
        "ntm_profiler_ncbi_organism_name",
        "fastani_best_hit_species",
        "fastani_best_hit_accession",
        "fastani_best_hit_strain",
        "fastani_best_hit_ani",
        "fastani_best_hit_alignment_fraction",
        "marker_majority_species",
        "marker_species_set",
        "final_wgs_species_call",
        "final_species_complex_call",
        "final_growth_type_call",
        "confidence_level",
        "evidence_conflict_flag",
        "exclude_from_downstream_reason",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df[keep].to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
