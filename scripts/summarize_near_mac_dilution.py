#!/usr/bin/env python3
"""Summarize staged near-MAC dilution results without retuning the primary rule."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


PRIMARY_20_80_THRESHOLD = 161.8
CALIBRATION_PAIRS = {"pairA", "pairB"}
EXPECTED_CONDITIONS_PER_PAIR = 36


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--allow-partial-calibration", action="store_true")
    return parser.parse_args()


def load_condition_rows(root: Path) -> pd.DataFrame:
    pattern = re.compile(
        r"^(pair[A-D])__(.+)_major__(.+)_minor__p([0-9]+)__s([0-9]+)$"
    )
    rows: list[dict[str, object]] = []
    for path in root.glob("conditions/*/*.minor_allele_burden.tsv"):
        match = pattern.match(path.parent.name)
        if not match:
            raise ValueError(f"Unexpected condition name: {path.parent.name}")
        result = pd.read_csv(path, sep="\t").iloc[0]
        callable_positions = int(result["callable_positions_depth_ge_20"])
        mixed_10_90 = int(result["mixed_sites_maf_0.10_0.90"])
        rows.append(
            {
                "pair_id": match.group(1),
                "major": match.group(2),
                "minor": match.group(3),
                "minor_percent": int(match.group(4)),
                "seed": int(match.group(5)),
                "callable_positions_depth_ge_20": callable_positions,
                "median_effective_depth": int(result["median_effective_depth"]),
                "mixed_sites_10_90_per_mbp_callable": (
                    mixed_10_90 * 1_000_000 / callable_positions if callable_positions else 0.0
                ),
                "mixed_sites_20_80_per_mbp_callable": float(
                    result["mixed_sites_20_80_per_mbp_callable"]
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["pair_id", "major", "minor", "minor_percent", "seed"]
    )


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    root = project / "results" / "near_mac_dilution" / "stage1_clean_reference"
    rows = load_condition_rows(root)
    if rows.empty:
        raise ValueError("No completed near-MAC dilution conditions were found")
    calibration_counts = rows.loc[rows["pair_id"].isin(CALIBRATION_PAIRS)].groupby("pair_id").size()
    calibration_complete = all(
        calibration_counts.get(pair_id, 0) == EXPECTED_CONDITIONS_PER_PAIR
        for pair_id in CALIBRATION_PAIRS
    )
    if not calibration_complete and not args.allow_partial_calibration:
        raise ValueError(
            "Both calibration pairs must be complete before freezing the secondary threshold"
        )
    available_calibration = set(calibration_counts[calibration_counts.gt(0)].index)
    pure_controls = rows.loc[
        rows["pair_id"].isin(available_calibration) & rows["minor_percent"].eq(0)
    ]
    secondary_threshold = 5 * pure_controls["mixed_sites_10_90_per_mbp_callable"].max()
    threshold_status = (
        "frozen_from_pairA_and_pairB_pure_controls"
        if calibration_complete
        else "provisional_partial_calibration_not_for_validation_claims"
    )
    rows["primary_20_80_positive"] = rows["mixed_sites_20_80_per_mbp_callable"].gt(
        PRIMARY_20_80_THRESHOLD
    )
    rows["secondary_10_90_positive"] = rows["mixed_sites_10_90_per_mbp_callable"].gt(
        secondary_threshold
    )
    rows["combined_positive"] = rows["primary_20_80_positive"] | rows[
        "secondary_10_90_positive"
    ]
    rows["primary_20_80_threshold"] = PRIMARY_20_80_THRESHOLD
    rows["secondary_10_90_threshold"] = secondary_threshold
    rows["secondary_threshold_status"] = threshold_status

    outdir = root / "summary_analysis"
    outdir.mkdir(parents=True, exist_ok=True)
    rows.to_csv(outdir / "nearmac_dilution_condition_metrics.tsv", sep="\t", index=False)
    aggregate = (
        rows.groupby(["pair_id", "major", "minor", "minor_percent"])
        .agg(
            replicates=("seed", "size"),
            median_depth=("median_effective_depth", "median"),
            burden_10_90_mean=("mixed_sites_10_90_per_mbp_callable", "mean"),
            burden_10_90_min=("mixed_sites_10_90_per_mbp_callable", "min"),
            burden_10_90_max=("mixed_sites_10_90_per_mbp_callable", "max"),
            burden_20_80_mean=("mixed_sites_20_80_per_mbp_callable", "mean"),
            burden_20_80_min=("mixed_sites_20_80_per_mbp_callable", "min"),
            burden_20_80_max=("mixed_sites_20_80_per_mbp_callable", "max"),
            primary_detection_rate=("primary_20_80_positive", "mean"),
            secondary_detection_rate=("secondary_10_90_positive", "mean"),
            combined_detection_rate=("combined_positive", "mean"),
        )
        .reset_index()
    )
    aggregate.to_csv(outdir / "nearmac_dilution_aggregate.tsv", sep="\t", index=False)

    report = [
        "# Near-MAC dilution summary",
        "",
        f"- Completed conditions: {len(rows)}.",
        f"- Primary frozen 20-80% threshold: {PRIMARY_20_80_THRESHOLD:.1f} sites/Mb.",
        f"- Secondary 10-90% threshold: {secondary_threshold:.2f} sites/Mb.",
        f"- Secondary threshold status: {threshold_status}.",
        "",
    ]
    for row in aggregate.itertuples(index=False):
        report.append(
            f"- {row.pair_id} {row.major} major/{row.minor} minor, {row.minor_percent}%: "
            f"n={row.replicates}, 10-90 burden mean {row.burden_10_90_mean:.1f}, "
            f"20-80 burden mean {row.burden_20_80_mean:.1f}, combined detection {row.combined_detection_rate:.0%}."
        )
    report.extend(
        [
            "",
            "Stage 1 maps mixtures to the clean major-source assembly. It calibrates the residual-mixture statistic but does not replace complete two-route reconstruction of selected conditions.",
        ]
    )
    (outdir / "nearmac_dilution_report.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print("\n".join(report))


if __name__ == "__main__":
    main()
