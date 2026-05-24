#!/usr/bin/env python3
"""Collect per-sample fastp JSON metrics into a read-QC summary table."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import yaml


REQUIRED_SAMPLE_COLUMNS = ("sample_id",)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, help="Input TSV samplesheet.")
    parser.add_argument("--fastp-dir", required=True, help="Directory with {sample}.fastp.json files.")
    parser.add_argument("--thresholds", required=True, help="QC thresholds YAML.")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument(
        "--assumed-genome-size-bp",
        type=float,
        default=6_000_000,
        help="Genome size used for rough depth estimates when no assembly exists yet.",
    )
    return parser.parse_args()


def read_samplesheet(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [col for col in REQUIRED_SAMPLE_COLUMNS if col not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"Samplesheet is missing required columns: {', '.join(missing)}")
        return [row["sample_id"] for row in reader if row.get("sample_id")]


def load_thresholds(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return data


def pct(value: float | int | None) -> str:
    if value is None:
        return "NA"
    return f"{float(value) * 100:.4f}"


def number(value: float | int | None, digits: int = 4) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.{digits}f}"


def main() -> None:
    args = parse_args()
    samplesheet = Path(args.samplesheet)
    fastp_dir = Path(args.fastp_dir)
    thresholds = load_thresholds(Path(args.thresholds))

    read_thresholds = thresholds.get("read_qc", {})
    assembly_thresholds = thresholds.get("assembly_qc", {})
    min_reads = int(read_thresholds.get("min_total_reads", 500_000))
    min_depth = float(read_thresholds.get("min_estimated_depth", 20))
    preferred_depth = float(read_thresholds.get("preferred_estimated_depth", 30))
    min_q30_rate = float(read_thresholds.get("min_q30_rate", 0.75))
    gc_min = float(assembly_thresholds.get("mycobacterium_gc_min", 60))
    gc_max = float(assembly_thresholds.get("mycobacterium_gc_max", 72))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "sample_id",
        "raw_reads",
        "clean_reads",
        "read_pass_rate",
        "raw_bases",
        "clean_bases",
        "estimated_depth_assuming_6mb",
        "q20_rate",
        "q30_rate",
        "gc_content",
        "duplication_rate",
        "adapter_trimmed_reads",
        "low_quality_reads",
        "too_many_N_reads",
        "too_short_reads",
        "gc_outside_mycobacterium_range",
        "read_qc_status",
        "warnings",
    ]

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for sample_id in read_samplesheet(samplesheet):
            fastp_json = fastp_dir / f"{sample_id}.fastp.json"
            warnings: list[str] = []
            if not fastp_json.exists():
                writer.writerow(
                    {
                        "sample_id": sample_id,
                        "read_qc_status": "FAIL",
                        "warnings": f"missing_fastp_json:{fastp_json}",
                    }
                )
                continue

            data = json.loads(fastp_json.read_text(encoding="utf-8"))
            before = data.get("summary", {}).get("before_filtering", {})
            after = data.get("summary", {}).get("after_filtering", {})
            filtering = data.get("filtering_result", {})
            duplication = data.get("duplication", {})
            adapter = data.get("adapter_cutting", {})

            raw_reads = int(before.get("total_reads", 0) or 0)
            clean_reads = int(after.get("total_reads", 0) or 0)
            raw_bases = int(before.get("total_bases", 0) or 0)
            clean_bases = int(after.get("total_bases", 0) or 0)
            q20_rate = float(after.get("q20_rate", 0) or 0)
            q30_rate = float(after.get("q30_rate", 0) or 0)
            gc_content = float(after.get("gc_content", 0) or 0)
            read_pass_rate = clean_reads / raw_reads if raw_reads else 0.0
            estimated_depth = clean_bases / float(args.assumed_genome_size_bp)
            gc_percent = gc_content * 100
            gc_warning = gc_percent < gc_min or gc_percent > gc_max

            if clean_reads < min_reads:
                warnings.append("clean_reads_below_minimum")
            if q30_rate < min_q30_rate:
                warnings.append("q30_rate_below_minimum")
            if estimated_depth < min_depth:
                warnings.append("estimated_depth_below_minimum")
            elif estimated_depth < preferred_depth:
                warnings.append("estimated_depth_below_preferred")
            if gc_warning:
                warnings.append("gc_outside_mycobacterium_range")

            critical_warnings = {
                "clean_reads_below_minimum",
                "q30_rate_below_minimum",
                "estimated_depth_below_minimum",
            }
            if critical_warnings.intersection(warnings):
                status = "WARN"
            elif warnings == ["gc_outside_mycobacterium_range"]:
                status = "PASS_WITH_GC_WARNING"
            elif warnings:
                status = "PASS_WITH_WARNING"
            else:
                status = "PASS"

            writer.writerow(
                {
                    "sample_id": sample_id,
                    "raw_reads": raw_reads,
                    "clean_reads": clean_reads,
                    "read_pass_rate": number(read_pass_rate),
                    "raw_bases": raw_bases,
                    "clean_bases": clean_bases,
                    "estimated_depth_assuming_6mb": number(estimated_depth, 2),
                    "q20_rate": pct(q20_rate),
                    "q30_rate": pct(q30_rate),
                    "gc_content": pct(gc_content),
                    "duplication_rate": pct(duplication.get("rate")),
                    "adapter_trimmed_reads": adapter.get("adapter_trimmed_reads", "NA"),
                    "low_quality_reads": filtering.get("low_quality_reads", "NA"),
                    "too_many_N_reads": filtering.get("too_many_N_reads", "NA"),
                    "too_short_reads": filtering.get("too_short_reads", "NA"),
                    "gc_outside_mycobacterium_range": str(gc_warning).upper(),
                    "read_qc_status": status,
                    "warnings": ";".join(warnings) if warnings else "none",
                }
            )


if __name__ == "__main__":
    main()
