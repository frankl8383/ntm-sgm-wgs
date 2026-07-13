#!/usr/bin/env python3
"""Summarize assembly and accessory metrics for direct and recovered genomes."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from statistics import median


CORE_SIGNALS = {
    "TMI_aromatic_catabolism_associated_block",
    "TMI_methyltransferase_hydrolase_cupin_block",
    "MP_MIP_nitrogen_redox_associated_block",
    "MP_MIP_oxidoreductase_associated_block",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "analysis_global_mac_upgrade/results/32_submission_freeze/"
            "recovered_direct_bias"
        ),
    )
    return parser.parse_args()


def read_tsv(path: Path, fieldnames: list[str] | None = None) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", fieldnames=fieldnames))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def fasta_stats(path: Path) -> dict[str, float | int]:
    lengths: list[int] = []
    gc = 0
    total = 0
    current = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
                continue
            sequence = line.strip().upper()
            current += len(sequence)
            total += len(sequence)
            gc += sequence.count("G") + sequence.count("C")
    if current:
        lengths.append(current)
    target = total / 2
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            n50 = length
            break
    return {
        "assembly_size_bp": total,
        "contigs": len(lengths),
        "n50_bp": n50,
        "gc_percent": 100 * gc / total,
    }


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    output = args.output_dir if args.output_dir.is_absolute() else project / args.output_dir
    results = project / "analysis_global_mac_upgrade/results"

    cohort = read_tsv(results / "14_rescue_full_batch/cohort/updated_interpretable_21_manifest.tsv")
    residual = {
        row["sample_id"]: row
        for row in read_tsv(
            results
            / "19_review_resolution/cohort21_residual_screen/"
            "cohort21_selected_assembly_residual_burden.tsv"
        )
    }
    panel = read_tsv(
        results / "17_mi_accessory_upgrade/full141_panel/mi_accessory_panel_manifest.tsv"
    )
    local_panel = {row["sample_id"]: row for row in panel if row["source"] == "local"}
    protein_counts = {
        row["tree_id"]: int(row["predicted_proteins"])
        for row in read_tsv(
            results
            / "17_mi_accessory_upgrade/full141_prodigal/"
            "uniform_prodigal_genome_counts.tsv"
        )
    }
    sequence_meta = {
        row["sequence_id"]: row["tree_id"]
        for row in read_tsv(
            results
            / "17_mi_accessory_upgrade/full141_prodigal/"
            "uniform_prodigal_sequence_metadata.tsv"
        )
    }
    search_rows = read_tsv(
        results
        / "17_mi_accessory_upgrade/full141_candidate_search/stable77_vs_full141.tsv",
        fieldnames=["query", "target", "pident", "alnlen", "qcov", "tcov", "evalue", "bits"],
    )
    family_sets: dict[str, set[str]] = {}
    for row in search_rows:
        tree_id = sequence_meta.get(row["target"])
        if tree_id:
            family_sets.setdefault(tree_id, set()).add(row["query"])

    block_rows = read_tsv(
        results
        / "32_submission_freeze/accessory_fdr5/"
        "curated_blocks_full_atlas_per_genome.tsv"
    )
    blocks_by_tree: dict[str, int] = {}
    for row in block_rows:
        if row["block_id"] in CORE_SIGNALS and truth(row["syntenic_block_present"]):
            blocks_by_tree[row["tree_id"]] = blocks_by_tree.get(row["tree_id"], 0) + 1
    route_concordance = {
        row["sample_id"]: row
        for row in read_tsv(
            results
            / "19_review_resolution/recovered_accessory_concordance/"
            "stable77_route_concordance_summary.tsv"
        )
    }

    rows: list[dict[str, object]] = []
    for cohort_row in cohort:
        sample = cohort_row["sample_id"]
        assembly = Path(cohort_row["selected_assembly_path"])
        stats = fasta_stats(assembly)
        residual_row = residual[sample]
        panel_row = local_panel.get(sample)
        tree_id = panel_row["tree_id"] if panel_row else ""
        route = route_concordance.get(sample, {})
        rows.append(
            {
                "sample_id": sample,
                "cohort_group": "recovered" if cohort_row["cohort_source"].startswith("rescued") else "direct",
                "cohort_source": cohort_row["cohort_source"],
                "analysis_lineage": panel_row["accessory_lineage"] if panel_row else "outside_MI_accessory_panel",
                **stats,
                "checkm2_completeness_percent": float(cohort_row["checkm2_completeness"]),
                "checkm2_contamination_percent": float(cohort_row["checkm2_contamination"]),
                "callable_fraction": (
                    int(residual_row["callable_positions_depth_ge_20"])
                    / int(residual_row["total_reference_positions"])
                ),
                "median_effective_depth": int(residual_row["median_effective_depth"]),
                "residual_20_80_sites_per_mbp": float(
                    residual_row["mixed_sites_20_80_per_mbp_callable"]
                ),
                "predicted_proteins": protein_counts.get(tree_id, ""),
                "stable77_families_detected": len(family_sets.get(tree_id, set())) if tree_id else "",
                "candidate_multigene_intervals_present": blocks_by_tree.get(tree_id, 0) if tree_id else "",
                "strict_meta_family_call_concordance": route.get("family_call_concordance", ""),
                "strict_meta_presence_jaccard": route.get("presence_jaccard", ""),
            }
        )

    write_tsv(output / "local21_recovered_direct_metrics.tsv", rows)
    mi_rows = [row for row in rows if row["analysis_lineage"] != "outside_MI_accessory_panel"]
    write_tsv(output / "mi17_recovered_direct_accessory_metrics.tsv", mi_rows)

    assembly_metrics = [
        "assembly_size_bp",
        "contigs",
        "n50_bp",
        "checkm2_completeness_percent",
        "checkm2_contamination_percent",
        "callable_fraction",
        "median_effective_depth",
        "residual_20_80_sites_per_mbp",
    ]
    accessory_metrics = [
        "predicted_proteins",
        "stable77_families_detected",
        "candidate_multigene_intervals_present",
    ]
    summaries: list[dict[str, object]] = []
    for scope, scope_rows, metrics in (
        ("local21", rows, assembly_metrics),
        ("MI17_accessory_panel", mi_rows, assembly_metrics + accessory_metrics),
    ):
        for metric in metrics:
            summary: dict[str, object] = {"scope": scope, "metric": metric}
            for group in ("direct", "recovered"):
                values = [
                    float(row[metric])
                    for row in scope_rows
                    if row["cohort_group"] == group and row[metric] != ""
                ]
                summary[f"{group}_n"] = len(values)
                summary[f"{group}_median"] = median(values) if values else ""
                summary[f"{group}_minimum"] = min(values) if values else ""
                summary[f"{group}_maximum"] = max(values) if values else ""
            summaries.append(summary)
    write_tsv(output / "recovered_direct_group_summary.tsv", summaries)
    report = [
        "# Recovered-versus-direct bias audit",
        "",
        "Metrics are descriptive. Direct and recovered groups differ in ascertainment and lineage composition, so no causal or exchangeable-group test is reported.",
        "",
        f"Local genomes: {len(rows)} ({sum(row['cohort_group'] == 'direct' for row in rows)} direct; {sum(row['cohort_group'] == 'recovered' for row in rows)} recovered).",
        f"M. intracellulare-complex accessory panel: {len(mi_rows)} genomes.",
    ]
    (output / "recovered_direct_bias_audit.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(f"local21={len(rows)} mi_panel={len(mi_rows)}")


if __name__ == "__main__":
    main()
