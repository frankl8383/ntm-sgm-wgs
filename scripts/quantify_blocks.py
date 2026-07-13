#!/usr/bin/env python3
"""Freeze manuscript-facing MI accessory evidence and figure source tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


CORE_BLOCKS = [
    "TMI_aromatic_catabolism_associated_block",
    "TMI_methyltransferase_hydrolase_cupin_block",
    "MP_MIP_nitrogen_redox_associated_block",
    "MP_MIP_oxidoreductase_associated_block",
]
SECONDARY_SIGNAL = "MP_MIP_ArsN1_family_B_signal"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotated-candidates", required=True)
    parser.add_argument("--full-family-validation", required=True)
    parser.add_argument("--block-members", required=True)
    parser.add_argument("--block-validation", required=True)
    parser.add_argument("--block-per-genome", required=True)
    parser.add_argument("--panel-with-project", required=True)
    parser.add_argument("--secondary-signal", default="")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    annotations = pd.read_csv(args.annotated_candidates, sep="\t", dtype=str).fillna("")
    family_validation = pd.read_csv(
        args.full_family_validation, sep="\t", dtype=str
    ).fillna("")
    block_members = pd.read_csv(args.block_members, sep="\t", dtype=str).fillna("")
    block_validation = pd.read_csv(
        args.block_validation, sep="\t", dtype=str
    ).fillna("")
    expected_signals = set(CORE_BLOCKS) | {SECONDARY_SIGNAL}
    observed_signals = set(block_validation.block_id)
    if observed_signals != expected_signals:
        raise SystemExit(
            "Expected the frozen five-signal multiplicity family; "
            f"observed={sorted(observed_signals)}"
        )
    if set(block_validation.multiplicity_family_size) != {"5"}:
        raise SystemExit("Block statistics were not corrected across five tracked signals")
    per_genome = pd.read_csv(args.block_per_genome, sep="\t", dtype=str).fillna("")
    panel = pd.read_csv(args.panel_with_project, sep="\t", dtype=str).fillna("")

    member_map = (
        block_members.groupby("cluster_id").block_id.apply(lambda x: ";".join(sorted(x)))
    )
    final = annotations.merge(
        family_validation,
        on="primary_cluster_id",
        how="left",
        validate="one_to_one",
    )
    final["curated_block_membership"] = final.primary_cluster_id.map(member_map).fillna("")
    supported = final.full_atlas_direction_supported.str.lower().eq("true")
    core_member = final.curated_block_membership.map(
        lambda value: any(block in CORE_BLOCKS for block in value.split(";") if block)
    )
    secondary = final.curated_block_membership.str.contains(SECONDARY_SIGNAL, regex=False)
    final["manuscript_evidence_tier"] = "not_retained_after_full_atlas_sensitivity"
    final.loc[supported, "manuscript_evidence_tier"] = "full_atlas_supported_family"
    final.loc[supported & core_member, "manuscript_evidence_tier"] = (
        "core_multigene_block_member"
    )
    final.loc[secondary, "manuscript_evidence_tier"] = (
        "secondary_stress_associated_signal"
    )
    preferred = [
        "primary_cluster_id",
        "association_direction",
        "conservative_product",
        "annotation_status",
        "curated_block_membership",
        "manuscript_evidence_tier",
        "full_atlas_direction_supported",
        "full_public_tmi_present",
        "full_public_tmi_total",
        "full_public_tmi_prevalence",
        "full_public_mp_present",
        "full_public_mp_total",
        "full_public_mp_prevalence",
        "full_local_tmi_present",
        "full_local_tmi_total",
        "full_local_mp_present",
        "full_local_mp_total",
        "full_public_fisher_fdr",
        "primary_fisher_fdr",
    ]
    remaining = [column for column in final.columns if column not in preferred]
    final = final[preferred + remaining]
    final.to_csv(
        output_dir / "stable77_accessory_families_final_evidence.tsv",
        sep="\t",
        index=False,
    )
    final[supported].to_csv(
        output_dir / "full_atlas_supported73_accessory_families.tsv",
        sep="\t",
        index=False,
    )
    final[supported & core_member].to_csv(
        output_dir / "core_multigene_block_family_members.tsv",
        sep="\t",
        index=False,
    )

    core_block_validation = block_validation[
        block_validation.block_id.isin(CORE_BLOCKS)
    ].copy()
    core_block_validation.to_csv(
        output_dir / "core_multigene_blocks_final_evidence.tsv",
        sep="\t",
        index=False,
    )
    if args.secondary_signal:
        secondary_signal = pd.read_csv(
            args.secondary_signal, sep="\t", dtype=str
        ).fillna("")
    else:
        secondary_signal = block_validation[
            block_validation.block_id.eq(SECONDARY_SIGNAL)
        ].copy()
    secondary_signal.to_csv(
        output_dir / "secondary_arsn1_signal_evidence.tsv", sep="\t", index=False
    )

    primary_genomes = panel[
        panel.source.isin(["public", "local"])
        & panel.accessory_lineage.isin(["MI_TMI_lineage", "MI_MP_MIP_lineage"])
    ].copy()
    block_matrix = per_genome[per_genome.block_id.isin(CORE_BLOCKS)].copy()
    block_matrix["syntenic_block_present"] = (
        block_matrix.syntenic_block_present.str.lower().eq("true").astype(int)
    )
    block_matrix = block_matrix.pivot(
        index="tree_id", columns="block_id", values="syntenic_block_present"
    ).reset_index()
    figure_source = primary_genomes[
        [
            "tree_id",
            "source",
            "sample_id",
            "reporting_label",
            "fine_label",
            "accessory_lineage",
            "cohort_source",
            "bioproject",
            "high_qc_public_sensitivity_include",
        ]
    ].merge(block_matrix, on="tree_id", how="left", validate="one_to_one")
    figure_source[CORE_BLOCKS] = figure_source[CORE_BLOCKS].fillna(0).astype(int)
    figure_source.to_csv(
        output_dir / "figure_source_mi_core_tree_accessory_blocks.tsv",
        sep="\t",
        index=False,
    )

    tier_counts = final.manuscript_evidence_tier.value_counts().to_dict()
    report = [
        "# Frozen MI accessory-genome evidence",
        "",
        f"Threshold-stable families entering full-atlas validation: {final.shape[0]}",
        f"Families retained in the full public MI atlas: {int(supported.sum())}",
        f"Families in four retained multigene blocks: {int((supported & core_member).sum())}",
        f"Other full-atlas-supported families: {tier_counts.get('full_atlas_supported_family', 0)}",
        "Multiplicity family: four candidate multigene intervals plus the exploratory ArsN1 family B signal (five tests).",
        "Secondary ArsN1 family B signal: retained separately as a single-family stress-associated observation; it does not meet the multigene syntenic-interval definition and is not a clinical resistance determinant.",
        "",
        "The four retained blocks are supported by balanced public discovery, local directional concordance, full-atlas sensitivity, stricter assembly-quality sensitivity, mixed-BioProject stratification and leave-one-BioProject-out analysis. Functional names are conservative homology annotations and do not establish phenotype, adaptation, horizontal transfer or causality.",
    ]
    (output_dir / "mi_accessory_evidence_freeze_report.md").write_text(
        "\n".join(report) + "\n"
    )
    print(
        f"families={final.shape[0]} supported={int(supported.sum())} "
        f"core_block_members={int((supported & core_member).sum())}"
    )


if __name__ == "__main__":
    main()
