#!/usr/bin/env python3
"""Build a conservative species evidence and conflict-resolution table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


MAC_TERMS = [
    ("paraintracellulare", "Mycobacterium paraintracellulare"),
    ("chimaera", "Mycobacterium chimaera"),
    ("yongonense", "Mycobacterium yongonense"),
    ("colombiense", "Mycobacterium colombiense"),
    ("marseillense", "Mycobacterium marseillense"),
    ("timonense", "Mycobacterium timonense"),
    ("arosiense", "Mycobacterium arosiense"),
    ("bouchedurhonense", "Mycobacterium bouchedurhonense"),
    ("mantenii", "Mycobacterium mantenii"),
    ("intracellulare", "Mycobacterium intracellulare"),
    ("avium", "Mycobacterium avium"),
]


def norm_species(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    lower = text.lower()
    for term, species in MAC_TERMS:
        if term in lower:
            return species
    if "gordonae" in lower:
        return "Mycobacterium gordonae"
    if "kansasii" in lower:
        return "Mycobacterium kansasii"
    if "simiae" in lower:
        return "Mycobacterium simiae"
    if "malmoense" in lower:
        return "Mycobacterium malmoense"
    if "marinum" in lower:
        return "Mycobacterium marinum"
    if "ulcerans" in lower:
        return "Mycobacterium ulcerans"
    if "mycobacter" in lower:
        return "Mycobacterium_unresolved"
    return "NA"


def parse_marker_species_set(row: pd.Series) -> str:
    species = []
    for value in row.values:
        if pd.isna(value) or str(value).strip() == "":
            continue
        species.append(norm_species(str(value).split("|", 1)[0]))
    species = sorted(set(s for s in species if s != "NA"))
    return ";".join(species) if species else "NA"


def mac_confidence(row: pd.Series) -> str:
    if row.get("gunc_status") is False or str(row.get("gunc_status")).lower() == "false":
        return "MAC_supported_but_genome_QC_warning"
    evidence = {
        row.get("public_best_norm"),
        row.get("type_panel_best_norm"),
        row.get("ntm_predicted_norm"),
        row.get("ntm_accession_norm"),
    }
    marker_set = str(row.get("independent_marker_top1_species_set", ""))
    evidence |= set(marker_set.split(";")) if marker_set and marker_set != "NA" else set()
    non_na = {e for e in evidence if e and e != "NA"}
    if non_na and all(e.startswith("Mycobacterium") and e not in {"Mycobacterium_unresolved"} for e in non_na):
        return "High_MAC_or_SGM_confidence"
    return "Review_needed"


def species_confidence(row: pd.Series) -> str:
    public = row.get("public_best_norm")
    type_best = row.get("type_panel_best_norm")
    ntm_pred = row.get("ntm_predicted_norm")
    ntm_acc = row.get("ntm_accession_norm")
    markers = set(str(row.get("independent_marker_top1_species_set", "NA")).split(";"))
    markers.discard("NA")
    core = {public, type_best, ntm_pred, ntm_acc} | markers
    core = {x for x in core if x and x != "NA"}
    if len(core) == 1:
        return "High_species_confidence"
    if public == type_best and (not markers or public in markers):
        return "High_with_minor_tool_or_marker_ambiguity"
    if public == ntm_pred and type_best != public:
        return "Moderate_public_NTM_support_type_panel_conflict"
    if public == type_best and ntm_pred != public:
        return "Moderate_ANI_support_NTM_profiler_granularity_conflict"
    return "Species_level_conflict_requires_review"


def review_action(row: pd.Series) -> str:
    conf = row["species_level_confidence_revised"]
    if conf == "High_species_confidence":
        return "retain_current_call"
    if "type_panel_conflict" in conf:
        return "review_type_strain_vs_public_context; do not revise automatically"
    if "NTM_profiler_granularity" in conf:
        return "retain_ANI_call_with_MAC_granularity_note"
    return "manual_review_before_final_species_claim"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-final", required=True, type=Path)
    parser.add_argument("--type-best", required=True, type=Path)
    parser.add_argument("--ntm-reparsed", required=True, type=Path)
    parser.add_argument("--marker-wide", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    current = pd.read_csv(args.current_final, sep="\t")
    type_best = pd.read_csv(args.type_best, sep="\t")
    ntm = pd.read_csv(args.ntm_reparsed, sep="\t")
    markers = pd.read_csv(args.marker_wide, sep="\t")

    type_best = type_best.rename(
        columns={
            "fastani_best_hit_species": "type_panel_best_species",
            "fastani_best_hit_accession": "type_panel_best_accession",
            "fastani_best_hit_strain": "type_panel_best_strain",
            "fastani_best_hit_ani": "type_panel_best_ani",
            "fastani_best_hit_alignment_fraction": "type_panel_best_alignment_fraction",
        }
    )
    marker_cols = [c for c in markers.columns if c.endswith("_top1_hit")]
    markers["independent_marker_top1_species_set"] = markers[marker_cols].apply(parse_marker_species_set, axis=1)

    df = (
        current.merge(
            type_best[
                [
                    "sample_id",
                    "type_panel_best_species",
                    "type_panel_best_accession",
                    "type_panel_best_strain",
                    "type_panel_best_ani",
                    "type_panel_best_alignment_fraction",
                ]
            ],
            on="sample_id",
            how="left",
        )
        .merge(ntm, on="sample_id", how="left")
        .merge(markers, on="sample_id", how="left")
    )

    df["public_best_norm"] = df["fastani_best_hit_species"].map(norm_species)
    df["type_panel_best_norm"] = df["type_panel_best_species"].map(norm_species)
    df["ntm_predicted_norm"] = df["ntm_profiler_predicted_species"].map(norm_species)
    df["ntm_accession_norm"] = df["ntm_profiler_closest_accession_organism_name"].map(norm_species)
    df["current_final_norm"] = df["final_wgs_species_call"].map(norm_species)
    df["MAC_complex_confidence"] = df.apply(mac_confidence, axis=1)
    df["species_level_confidence_revised"] = df.apply(species_confidence, axis=1)
    df["downstream_inclusion_status_revised"] = df["exclude_from_downstream_reason"].apply(
        lambda v: "included_high_confidence_downstream"
        if pd.isna(v) or str(v).strip() in {"", "NA", "nan"}
        else f"excluded_or_warning:{v}"
    )
    df["recommended_review_action"] = df.apply(review_action, axis=1)

    columns = [
        "sample_id",
        "current_final_norm",
        "public_best_norm",
        "fastani_best_hit_accession",
        "fastani_best_hit_strain",
        "fastani_best_hit_ani",
        "fastani_best_hit_alignment_fraction",
        "type_panel_best_norm",
        "type_panel_best_accession",
        "type_panel_best_strain",
        "type_panel_best_ani",
        "type_panel_best_alignment_fraction",
        "ntm_predicted_norm",
        "ntm_accession_norm",
        "ntm_profiler_species_vs_accession_name_conflict",
        "independent_marker_top1_species_set",
        *marker_cols,
        "gunc_status",
        "checkm2_contamination",
        "num_contigs",
        "MAC_complex_confidence",
        "species_level_confidence_revised",
        "downstream_inclusion_status_revised",
        "recommended_review_action",
        "evidence_conflict_flag",
        "exclude_from_downstream_reason",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df[[c for c in columns if c in df.columns]].to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
