#!/usr/bin/env python3
"""Summarize substitution-only intermediate-frequency alleles from mpileup."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


def base_counts(ref: str, pile: str) -> Counter[str]:
    counts: Counter[str] = Counter()
    index = 0
    while index < len(pile):
        char = pile[index]
        if char == "^":
            index += 2
            continue
        if char == "$":
            index += 1
            continue
        if char in "+-":
            index += 1
            start = index
            while index < len(pile) and pile[index].isdigit():
                index += 1
            length = int(pile[start:index]) if index > start else 0
            index += length
            continue
        if char in ".,":
            counts[ref.upper()] += 1
        elif char.upper() in {"A", "C", "G", "T"}:
            counts[char.upper()] += 1
        index += 1
    return counts


def histogram_median(histogram: Counter[int], total: int) -> int:
    midpoint = (total + 1) // 2
    running = 0
    for depth in sorted(histogram):
        running += histogram[depth]
        if running >= midpoint:
            return depth
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--route", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-depth", type=int, default=20)
    parser.add_argument("--min-minor-count", type=int, default=5)
    args = parser.parse_args()

    callable_positions = 0
    mixed_20_80 = 0
    mixed_10_90 = 0
    total_positions = 0
    depth_hist: Counter[int] = Counter()
    for line in sys.stdin:
        fields = line.rstrip().split("\t")
        if len(fields) < 5:
            continue
        total_positions += 1
        counts = base_counts(fields[2], fields[4])
        effective_depth = sum(counts.values())
        if effective_depth < args.min_depth:
            continue
        callable_positions += 1
        depth_hist[effective_depth] += 1
        ordered = sorted(counts.values(), reverse=True)
        if len(ordered) < 2 or ordered[1] < args.min_minor_count:
            continue
        minor_fraction = ordered[1] / effective_depth
        if 0.10 <= minor_fraction <= 0.90:
            mixed_10_90 += 1
        if 0.20 <= minor_fraction <= 0.80:
            mixed_20_80 += 1

    row = {
        "sample_id": args.sample,
        "route": args.route,
        "total_reference_positions": total_positions,
        "callable_positions_depth_ge_20": callable_positions,
        "median_effective_depth": histogram_median(depth_hist, callable_positions),
        "mixed_sites_maf_0.10_0.90": mixed_10_90,
        "mixed_sites_maf_0.20_0.80": mixed_20_80,
        "mixed_sites_20_80_per_mbp_callable": (mixed_20_80 * 1_000_000 / callable_positions) if callable_positions else 0,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
