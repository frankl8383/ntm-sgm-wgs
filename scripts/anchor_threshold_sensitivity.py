#!/usr/bin/env python3
"""Evaluate external principal-lineage membership across fixed ANI floors."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


THRESHOLDS = (98.0, 98.5, 98.6, 99.0)
BLOCK_LINEAGE = {
    "TMI_aromatic_catabolism_associated_block": "TMI",
    "TMI_methyltransferase_hydrolase_cupin_block": "TMI",
    "MP_MIP_nitrogen_redox_associated_block": "MP_MIP",
    "MP_MIP_oxidoreductase_associated_block": "MP_MIP",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "analysis_global_mac_upgrade/results/32_submission_freeze/"
            "anchor_threshold_sensitivity"
        ),
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    output = args.output_dir if args.output_dir.is_absolute() else project / args.output_dir
    external = (
        project
        / "analysis_global_mac_upgrade/results/27_postreview_validation/"
        "external_PRJNA983112"
    )
    lineage = read_tsv(external / "external_lineage_ani.tsv")
    calls = read_tsv(external / "external_block_calls.tsv")
    calls_by_sample_block = {
        (row["sample_id"], row["block_id"]): truth(row["syntenic_block_present"])
        for row in calls
    }

    sample_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    for threshold in THRESHOLDS:
        qualified: list[tuple[dict[str, str], str]] = []
        for row in lineage:
            anchor_max = float(row["principal_lineage_anchor_max_ani"])
            is_qualified = anchor_max >= threshold
            assigned = ""
            if is_qualified:
                assigned = (
                    "TMI"
                    if float(row["tmi_ani"]) >= float(row["mp_mip_ani"])
                    else "MP_MIP"
                )
                qualified.append((row, assigned))
            sample_rows.append(
                {
                    "threshold_ani_percent": threshold,
                    "sample_id": row["sample_id"],
                    "principal_anchor_max_ani": anchor_max,
                    "qualified": is_qualified,
                    "assigned_lineage": assigned or "outside_principal_comparison",
                    "published_expected_lineage": row["frozen_expected_lineage"],
                    "expected_lineage_concordant": (
                        is_qualified and assigned == row["frozen_expected_lineage"]
                    ),
                    "external_qc_pass": truth(row["external_qc_pass"]),
                }
            )

        tmi = [row for row, assigned in qualified if assigned == "TMI"]
        mp = [row for row, assigned in qualified if assigned == "MP_MIP"]
        summary_rows.append(
            {
                "threshold_ani_percent": threshold,
                "qualified_total": len(qualified),
                "qualified_tmi": len(tmi),
                "qualified_mp_mip": len(mp),
                "qualified_qc_pass": sum(truth(row["external_qc_pass"]) for row, _ in qualified),
                "expected_lineage_concordant": sum(
                    assigned == row["frozen_expected_lineage"]
                    for row, assigned in qualified
                ),
                "all_four_block_directions_retained": True,
            }
        )
        samples_by_lineage = {
            "TMI": [row["sample_id"] for row, assigned in qualified if assigned == "TMI"],
            "MP_MIP": [row["sample_id"] for row, assigned in qualified if assigned == "MP_MIP"],
        }
        for block, expected_lineage in BLOCK_LINEAGE.items():
            other_lineage = "MP_MIP" if expected_lineage == "TMI" else "TMI"
            expected_samples = samples_by_lineage[expected_lineage]
            other_samples = samples_by_lineage[other_lineage]
            expected_present = sum(
                calls_by_sample_block[(sample, block)] for sample in expected_samples
            )
            other_present = sum(
                calls_by_sample_block[(sample, block)] for sample in other_samples
            )
            direction_retained = (
                expected_present / len(expected_samples) > other_present / len(other_samples)
                if expected_samples and other_samples
                else False
            )
            block_rows.append(
                {
                    "threshold_ani_percent": threshold,
                    "block_id": block,
                    "expected_lineage": expected_lineage,
                    "expected_lineage_present": expected_present,
                    "expected_lineage_total": len(expected_samples),
                    "other_lineage_present": other_present,
                    "other_lineage_total": len(other_samples),
                    "direction_retained": direction_retained,
                }
            )
            if not direction_retained:
                summary_rows[-1]["all_four_block_directions_retained"] = False

    write_tsv(output / "external_anchor_threshold_sample_calls.tsv", sample_rows)
    write_tsv(output / "external_anchor_threshold_summary.tsv", summary_rows)
    write_tsv(output / "external_anchor_threshold_block_direction.tsv", block_rows)
    report = [
        "# External anchor-threshold sensitivity",
        "",
        "The 98.60% ANI value is an operational inclusion floor, not a species boundary.",
        "",
        "| ANI floor (%) | Qualified | TMI | MP-MIP | QC pass | Expected-lineage concordant | Four block directions retained |",
        "|---:|---:|---:|---:|---:|---:|---|",
    ]
    report.extend(
        "| {threshold_ani_percent:.1f} | {qualified_total} | {qualified_tmi} | "
        "{qualified_mp_mip} | {qualified_qc_pass} | {expected_lineage_concordant} | "
        "{all_four_block_directions_retained} |".format(**row)
        for row in summary_rows
    )
    (output / "external_anchor_threshold_sensitivity.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(f"thresholds={len(THRESHOLDS)} samples={len(lineage)}")


if __name__ == "__main__":
    main()
