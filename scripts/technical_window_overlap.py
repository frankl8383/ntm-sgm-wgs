#!/usr/bin/env python3
"""Summarize dependence among the seeded near-MAC technical read windows."""

from __future__ import annotations

import argparse
import csv
from itertools import combinations
from pathlib import Path
from statistics import mean, median


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    return parser.parse_args()


def read_one(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one row in {path}, found {len(rows)}")
    return rows[0]


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    ordered = sorted(intervals)
    total = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    window_root = (
        project
        / "analysis_global_mac_upgrade/results/20_nearmac_dilution/"
        "stage1_clean_reference/seeded_windows"
    )
    output = (
        project
        / "analysis_global_mac_upgrade/results/32_submission_freeze/"
        "technical_window_dependence"
    )

    by_sample: dict[str, list[dict[str, object]]] = {}
    for path in sorted(window_root.glob("*/seed_*/window.tsv")):
        row = read_one(path)
        sample = row["sample_id"]
        start = int(row["window_start_zero_based"])
        length = int(row["window_pairs"])
        by_sample.setdefault(sample, []).append(
            {
                "sample_id": sample,
                "seed": int(row["seed"]),
                "total_source_pairs": int(row["total_source_pairs"]),
                "window_start_zero_based": start,
                "window_pairs": length,
                "window_end_exclusive": start + length,
            }
        )

    if len(by_sample) != 8:
        raise ValueError(f"Expected eight source genomes, found {len(by_sample)}")

    pairwise_rows: list[dict[str, object]] = []
    source_rows: list[dict[str, object]] = []
    for sample, windows in sorted(by_sample.items()):
        windows = sorted(windows, key=lambda row: int(row["seed"]))
        if len(windows) != 3:
            raise ValueError(f"Expected three windows for {sample}, found {len(windows)}")
        intervals = [
            (int(row["window_start_zero_based"]), int(row["window_end_exclusive"]))
            for row in windows
        ]
        union_pairs = interval_union_length(intervals)
        nominal_pairs = sum(int(row["window_pairs"]) for row in windows)
        source_rows.append(
            {
                "sample_id": sample,
                "source_pairs": windows[0]["total_source_pairs"],
                "technical_windows": len(windows),
                "nominal_window_pairs": nominal_pairs,
                "unique_union_pairs": union_pairs,
                "unique_fraction_of_nominal_windows": union_pairs / nominal_pairs,
            }
        )
        for left, right in combinations(windows, 2):
            overlap = max(
                0,
                min(int(left["window_end_exclusive"]), int(right["window_end_exclusive"]))
                - max(int(left["window_start_zero_based"]), int(right["window_start_zero_based"])),
            )
            denominator = min(int(left["window_pairs"]), int(right["window_pairs"]))
            pairwise_rows.append(
                {
                    "sample_id": sample,
                    "seed_a": left["seed"],
                    "seed_b": right["seed"],
                    "start_a_zero_based": left["window_start_zero_based"],
                    "start_b_zero_based": right["window_start_zero_based"],
                    "window_pairs": denominator,
                    "overlap_pairs": overlap,
                    "overlap_fraction_of_window": overlap / denominator,
                }
            )

    overlaps = [float(row["overlap_fraction_of_window"]) for row in pairwise_rows]
    union_fractions = [
        float(row["unique_fraction_of_nominal_windows"]) for row in source_rows
    ]
    summary_rows = [
        {
            "source_genomes": len(source_rows),
            "windows_per_source": 3,
            "pairwise_window_comparisons": len(pairwise_rows),
            "mean_pairwise_overlap_fraction": mean(overlaps),
            "median_pairwise_overlap_fraction": median(overlaps),
            "minimum_pairwise_overlap_fraction": min(overlaps),
            "maximum_pairwise_overlap_fraction": max(overlaps),
            "mean_unique_fraction_of_nominal_windows": mean(union_fractions),
            "minimum_unique_fraction_of_nominal_windows": min(union_fractions),
            "maximum_unique_fraction_of_nominal_windows": max(union_fractions),
            "interpretation": "Dependent deterministic technical subsamples; not independent replicates",
        }
    ]

    write_tsv(output / "technical_window_pairwise_overlap.tsv", pairwise_rows)
    write_tsv(output / "technical_window_source_summary.tsv", source_rows)
    write_tsv(output / "technical_window_overlap_summary.tsv", summary_rows)

    summary = summary_rows[0]
    report = [
        "# Technical-window dependence",
        "",
        f"Eight source genomes contributed three deterministic two-million-pair windows each. The {len(pairwise_rows)} within-source pairwise comparisons had a mean overlap of {100 * float(summary['mean_pairwise_overlap_fraction']):.1f}% and a maximum overlap of {100 * float(summary['maximum_pairwise_overlap_fraction']):.1f}%.",
        "",
        "These windows are dependent technical subsamples and are not independent biological or computational replicates.",
    ]
    (output / "technical_window_dependence.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(output)


if __name__ == "__main__":
    main()
