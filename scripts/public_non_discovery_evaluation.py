#!/usr/bin/env python3
"""Evaluate five frozen lineage signals in public genomes excluded from discovery."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact


CORE_BLOCKS = (
    "TMI_aromatic_catabolism_associated_block",
    "TMI_methyltransferase_hydrolase_cupin_block",
    "MP_MIP_nitrogen_redox_associated_block",
    "MP_MIP_oxidoreductase_associated_block",
)
SECONDARY_SIGNAL = "MP_MIP_ArsN1_family_B_signal"
TRACKED_SIGNALS = CORE_BLOCKS + (SECONDARY_SIGNAL,)
LINEAGES = ("MI_TMI_lineage", "MI_MP_MIP_lineage")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def bh_adjust(p_values: list[float]) -> list[float]:
    ordered = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * len(p_values)
    running = 1.0
    for rank_index in range(len(ordered) - 1, -1, -1):
        original_index, p_value = ordered[rank_index]
        rank = rank_index + 1
        running = min(running, p_value * len(p_values) / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    base = project / "analysis_global_mac_upgrade" / "results"
    discovery_manifest = base / "17_mi_accessory_upgrade" / "panel" / "mi_accessory_panel_manifest.tsv"
    full_manifest = base / "17_mi_accessory_upgrade" / "full141_panel" / "mi_accessory_panel_manifest.tsv"
    block_calls = (
        base
        / "17_mi_accessory_upgrade"
        / "full141_block_validation"
        / "curated_blocks_full_atlas_per_genome.tsv"
    )
    output = base / "24_public_non_discovery_evaluation"

    discovery_ids = {
        row["tree_id"]
        for row in read_tsv(discovery_manifest)
        if row["source"] == "public" and row["accessory_lineage"] in LINEAGES
    }
    full_rows = {
        row["tree_id"]: row
        for row in read_tsv(full_manifest)
        if row["source"] == "public" and row["accessory_lineage"] in LINEAGES
    }
    evaluation_ids = set(full_rows) - discovery_ids
    if len(discovery_ids) != 60 or len(evaluation_ids) != 47:
        raise ValueError(
            f"Unexpected public split: discovery={len(discovery_ids)}, evaluation={len(evaluation_ids)}"
        )

    membership_rows: list[dict[str, object]] = []
    for tree_id in sorted(full_rows):
        row = full_rows[tree_id]
        membership_rows.append(
            {
                "tree_id": tree_id,
                "accessory_lineage": row["accessory_lineage"],
                "cohort_role": "discovery" if tree_id in discovery_ids else "non_discovery_evaluation",
            }
        )
    write_tsv(
        output / "public_discovery_evaluation_membership.tsv",
        membership_rows,
        ["tree_id", "accessory_lineage", "cohort_role"],
    )

    calls: dict[tuple[str, str], bool] = {}
    evaluation_call_rows: list[dict[str, object]] = []
    for row in read_tsv(block_calls):
        if row["tree_id"] not in evaluation_ids or row["block_id"] not in TRACKED_SIGNALS:
            continue
        present = row["syntenic_block_present"].lower() == "true"
        calls[(row["block_id"], row["tree_id"])] = present
        evaluation_call_rows.append(
            {
                "block_id": row["block_id"],
                "signal_class": (
                    "single_family_exploratory_signal"
                    if row["block_id"] == SECONDARY_SIGNAL
                    else "candidate_multigene_syntenic_interval"
                ),
                "association_direction": row["association_direction"],
                "tree_id": row["tree_id"],
                "accessory_lineage": full_rows[row["tree_id"]]["accessory_lineage"],
                "syntenic_block_present": present,
            }
        )
    write_tsv(
        output / "non_discovery_public_signal_calls.tsv",
        evaluation_call_rows,
        [
            "block_id",
            "signal_class",
            "association_direction",
            "tree_id",
            "accessory_lineage",
            "syntenic_block_present",
        ],
    )

    stats: list[dict[str, object]] = []
    p_values: list[float] = []
    for block_id in TRACKED_SIGNALS:
        expected = "MI_TMI_lineage" if block_id.startswith("TMI_") else "MI_MP_MIP_lineage"
        opposite = "MI_MP_MIP_lineage" if expected == "MI_TMI_lineage" else "MI_TMI_lineage"
        grouped: dict[str, list[bool]] = defaultdict(list)
        for tree_id in evaluation_ids:
            grouped[full_rows[tree_id]["accessory_lineage"]].append(calls[(block_id, tree_id)])
        expected_present = sum(grouped[expected])
        expected_absent = len(grouped[expected]) - expected_present
        opposite_present = sum(grouped[opposite])
        opposite_absent = len(grouped[opposite]) - opposite_present
        odds_ratio, p_value = fisher_exact(
            [[expected_present, expected_absent], [opposite_present, opposite_absent]],
            alternative="two-sided",
        )
        p_values.append(float(p_value))
        stats.append(
            {
                "block_id": block_id,
                "signal_class": (
                    "single_family_exploratory_signal"
                    if block_id == SECONDARY_SIGNAL
                    else "candidate_multigene_syntenic_interval"
                ),
                "multiplicity_family_size": 5,
                "expected_lineage": expected,
                "tmi_present": sum(grouped["MI_TMI_lineage"]),
                "tmi_total": len(grouped["MI_TMI_lineage"]),
                "mp_mip_present": sum(grouped["MI_MP_MIP_lineage"]),
                "mp_mip_total": len(grouped["MI_MP_MIP_lineage"]),
                "expected_present": expected_present,
                "expected_total": len(grouped[expected]),
                "opposite_present": opposite_present,
                "opposite_total": len(grouped[opposite]),
                "prevalence_difference": expected_present / len(grouped[expected])
                - opposite_present / len(grouped[opposite]),
                "odds_ratio": float(odds_ratio),
                "odds_ratio_display": "Inf" if math.isinf(odds_ratio) else f"{odds_ratio:.3g}",
                "fisher_p": float(p_value),
            }
        )
    for row, q_value in zip(stats, bh_adjust(p_values), strict=True):
        row["bh_fdr_five_signals"] = q_value
    fields = [
        "block_id",
        "signal_class",
        "multiplicity_family_size",
        "expected_lineage",
        "tmi_present",
        "tmi_total",
        "mp_mip_present",
        "mp_mip_total",
        "expected_present",
        "expected_total",
        "opposite_present",
        "opposite_total",
        "prevalence_difference",
        "odds_ratio",
        "odds_ratio_display",
        "fisher_p",
        "bh_fdr_five_signals",
    ]
    write_tsv(
        output / "non_discovery_public_signal_statistics.tsv",
        stats,
        fields,
    )
    core_stats = [row for row in stats if row["block_id"] in CORE_BLOCKS]
    secondary_stats = [row for row in stats if row["block_id"] == SECONDARY_SIGNAL]
    core_calls = [row for row in evaluation_call_rows if row["block_id"] in CORE_BLOCKS]
    write_tsv(
        output / "non_discovery_public_block_statistics.tsv", core_stats, fields
    )
    write_tsv(
        output / "non_discovery_public_arsn1_statistics.tsv", secondary_stats, fields
    )
    write_tsv(
        output / "non_discovery_public_block_calls.tsv",
        core_calls,
        [
            "block_id",
            "signal_class",
            "association_direction",
            "tree_id",
            "accessory_lineage",
            "syntenic_block_present",
        ],
    )

    summary = output / "non_discovery_public_block_evaluation.md"
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("# Frozen non-discovery public evaluation\n\n")
        handle.write(
            "The original public discovery panel contained 32 TMI and 28 MP-MIP genomes. "
            "The later public atlas contributed 32 additional TMI and 15 additional MP-MIP genomes "
            "that were not used for family discovery or local interval curation. Five frozen signals "
            "were evaluated together: four candidate multigene intervals and one exploratory ArsN1 "
            "family signal.\n\n"
        )
        for row in core_stats:
            handle.write(
                f"- {row['block_id']}: {row['expected_present']}/{row['expected_total']} in the "
                f"expected lineage versus {row['opposite_present']}/{row['opposite_total']} in the "
                f"opposite lineage; five-signal FDR={row['bh_fdr_five_signals']:.6g}.\n"
            )
        row = secondary_stats[0]
        handle.write(
            f"- Exploratory ArsN1 signal: {row['expected_present']}/{row['expected_total']} in the "
            f"expected lineage versus {row['opposite_present']}/{row['opposite_total']} in the "
            f"opposite lineage; five-signal FDR={row['bh_fdr_five_signals']:.6g}.\n"
        )
        handle.write(
            "\nThis is a frozen public non-discovery evaluation, not a population-prevalence estimate "
            "or a clinical validation cohort.\n"
        )

    print(output / "non_discovery_public_block_statistics.tsv")


if __name__ == "__main__":
    main()
