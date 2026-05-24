#!/usr/bin/env python3
"""Collect assembly FASTA, CheckM2, and optional GUNC metrics into one QC table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import yaml


FIELDS = [
    "sample_id",
    "assembly_size",
    "num_contigs",
    "n50",
    "gc_content",
    "largest_contig",
    "quast_status",
    "checkm2_completeness",
    "checkm2_contamination",
    "gunc_pass_or_fail",
    "gunc_taxonomic_level",
    "gunc_contamination_portion",
    "gunc_clade_separation_score",
    "gunc_reference_representation_score",
    "assembly_qc_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--assembly-dir", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--quast-report")
    parser.add_argument("--checkm2-report")
    parser.add_argument("--gunc-report")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_samples(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row["sample_id"] for row in reader if row.get("sample_id")]


def read_checkm2(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["Name"]: row for row in reader if row.get("Name")}


def read_gunc(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            sample = row.get("genome") or row.get("Name") or row.get("sample_id")
            if not sample:
                continue
            sample = Path(sample).stem
            pass_value = row.get("pass.GUNC") or row.get("pass_gunc") or row.get("pass")
            rows[sample] = {
                "gunc_pass_or_fail": pass_value if pass_value not in (None, "") else "NA",
                "gunc_taxonomic_level": row.get("taxonomic_level", "NA") or "NA",
                "gunc_contamination_portion": row.get("contamination_portion", "NA") or "NA",
                "gunc_clade_separation_score": row.get("clade_separation_score", "NA") or "NA",
                "gunc_reference_representation_score": row.get("reference_representation_score", "NA") or "NA",
            }
        return rows


def fasta_lengths(path: Path) -> tuple[list[int], int, int]:
    lengths: list[int] = []
    gc = 0
    current = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
                continue
            seq = line.upper()
            current += len(seq)
            gc += seq.count("G") + seq.count("C")
    if current:
        lengths.append(current)
    total = sum(lengths)
    return lengths, total, gc


def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    half = sum(lengths) / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0


def status_for_metrics(metrics: dict[str, Any], thresholds: dict[str, Any]) -> tuple[str, str]:
    assembly_qc = thresholds.get("assembly_qc", {})
    max_contigs = int(assembly_qc.get("max_contigs_preferred", 500))
    min_n50 = int(assembly_qc.get("min_n50_preferred", 20000))
    gc_min = float(assembly_qc.get("mycobacterium_gc_min", 60))
    gc_max = float(assembly_qc.get("mycobacterium_gc_max", 72))
    min_completeness = float(assembly_qc.get("min_completeness", 90))
    max_contamination_preferred = float(assembly_qc.get("max_contamination_preferred", 5))
    max_contamination_warning = float(assembly_qc.get("max_contamination_warning", 10))

    warnings: list[str] = []
    if metrics["num_contigs"] > max_contigs:
        warnings.append(f"contigs>{max_contigs}")
    if metrics["n50"] < min_n50:
        warnings.append(f"n50<{min_n50}")
    if not (gc_min <= metrics["gc_content"] <= gc_max):
        warnings.append(f"gc_outside_{gc_min:g}_{gc_max:g}")
    checkm2_completeness = metrics.get("checkm2_completeness")
    checkm2_contamination = metrics.get("checkm2_contamination")
    if isinstance(checkm2_completeness, float) and checkm2_completeness < min_completeness:
        warnings.append(f"checkm2_completeness<{min_completeness:g}")
    if isinstance(checkm2_contamination, float):
        if checkm2_contamination > max_contamination_warning:
            warnings.append(f"checkm2_contamination>{max_contamination_warning:g}")
        elif checkm2_contamination > max_contamination_preferred:
            warnings.append(f"checkm2_contamination>{max_contamination_preferred:g}")
    if str(metrics.get("gunc_pass_or_fail", "NA")).lower() in {"false", "fail", "failed"}:
        warnings.append("gunc_fail")

    if warnings:
        return "WARN", ";".join(warnings)
    return "PASS", "pass_basic_assembly_thresholds"


def format_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def main() -> int:
    args = parse_args()
    samples = read_samples(Path(args.samplesheet))
    assembly_dir = Path(args.assembly_dir)
    with Path(args.thresholds).open("r", encoding="utf-8") as handle:
        thresholds = yaml.safe_load(handle) or {}
    checkm2 = read_checkm2(Path(args.checkm2_report) if args.checkm2_report else None)
    gunc = read_gunc(Path(args.gunc_report) if args.gunc_report else None)
    quast_report = Path(args.quast_report) if args.quast_report else None
    quast_status = "QUAST_5.3.0_RUN" if quast_report and quast_report.exists() else "NOT_RUN_BASIC_FASTA_STATS_ONLY"

    rows: list[dict[str, Any]] = []
    for sample in samples:
        fasta = assembly_dir / sample / f"{sample}.assembly.fasta"
        if not fasta.exists():
            raise FileNotFoundError(f"Missing assembly FASTA: {fasta}")
        lengths, total, gc = fasta_lengths(fasta)
        gc_content = (gc / total * 100) if total else 0.0
        metrics: dict[str, Any] = {
            "sample_id": sample,
            "assembly_size": total,
            "num_contigs": len(lengths),
            "n50": n50(lengths),
            "gc_content": gc_content,
            "largest_contig": max(lengths) if lengths else 0,
        }
        checkm2_row = checkm2.get(sample, {})
        if checkm2_row:
            metrics["checkm2_completeness"] = float(checkm2_row.get("Completeness", "nan"))
            metrics["checkm2_contamination"] = float(checkm2_row.get("Contamination", "nan"))
        else:
            metrics["checkm2_completeness"] = "NA"
            metrics["checkm2_contamination"] = "NA"
        gunc_row = gunc.get(sample, {})
        metrics["gunc_pass_or_fail"] = gunc_row.get("gunc_pass_or_fail", "NA")
        metrics["gunc_taxonomic_level"] = gunc_row.get("gunc_taxonomic_level", "NA")
        metrics["gunc_contamination_portion"] = gunc_row.get("gunc_contamination_portion", "NA")
        metrics["gunc_clade_separation_score"] = gunc_row.get("gunc_clade_separation_score", "NA")
        metrics["gunc_reference_representation_score"] = gunc_row.get(
            "gunc_reference_representation_score", "NA"
        )
        _, assembly_qc_status = status_for_metrics(metrics, thresholds)
        metrics.update(
            {
                "quast_status": quast_status,
                "assembly_qc_status": assembly_qc_status,
            }
        )
        rows.append(metrics)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field, "NA")) for field in FIELDS})
    print(f"Wrote {len(rows)} rows: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
