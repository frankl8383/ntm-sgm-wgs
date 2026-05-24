#!/usr/bin/env python3
"""Validate the NTM WGS samplesheet before running workflow modules."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


REQUIRED_COLUMNS = [
    "sample_id",
    "fastq_r1",
    "fastq_r2",
    "initial_species_label",
    "initial_growth_type",
    "collection_date",
    "isolation_source",
    "patient_id_anonymized",
    "sequencing_platform",
    "read_length",
    "notes",
]

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, help="Input samplesheet TSV.")
    parser.add_argument("--output", required=True, help="Output validation report TSV.")
    parser.add_argument("--project-root", default=".", help="Project root for resolving relative FASTQ paths.")
    return parser.parse_args()


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return project_root / path


def is_fastq(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in FASTQ_SUFFIXES)


def read_token(path: Path) -> str | None:
    name = path.name.lower()
    r1_tokens = ("_r1", ".r1", "_1.", ".1.", "clean.1", "read1")
    r2_tokens = ("_r2", ".r2", "_2.", ".2.", "clean.2", "read2")
    if any(token in name for token in r1_tokens):
        return "R1"
    if any(token in name for token in r2_tokens):
        return "R2"
    return None


def add_report(reports: list[dict[str, str]], sample_id: str, check: str, status: str, severity: str, message: str, value: str = "") -> None:
    reports.append(
        {
            "sample_id": sample_id,
            "check": check,
            "status": status,
            "severity": severity,
            "message": message,
            "value": value,
        }
    )


def validate_rows(rows: list[dict[str, str]], fieldnames: list[str] | None, project_root: Path) -> tuple[list[dict[str, str]], int]:
    reports: list[dict[str, str]] = []
    critical_errors = 0

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in (fieldnames or [])]
    if missing_columns:
        add_report(
            reports,
            "GLOBAL",
            "required_columns",
            "FAIL",
            "ERROR",
            "Missing required columns.",
            ",".join(missing_columns),
        )
        return reports, 1
    add_report(reports, "GLOBAL", "required_columns", "PASS", "INFO", "All required columns are present.")

    if not rows:
        add_report(reports, "GLOBAL", "non_empty_samplesheet", "FAIL", "ERROR", "Samplesheet has no sample rows.")
        return reports, 1
    add_report(reports, "GLOBAL", "non_empty_samplesheet", "PASS", "INFO", f"Found {len(rows)} sample rows.")

    seen: dict[str, int] = {}
    for row_number, row in enumerate(rows, start=2):
        sample_id = (row.get("sample_id") or "").strip()
        display_id = sample_id or f"ROW_{row_number}"

        if not sample_id:
            add_report(reports, display_id, "sample_id_present", "FAIL", "ERROR", "sample_id is empty.", str(row_number))
            critical_errors += 1
            continue

        seen[sample_id] = seen.get(sample_id, 0) + 1
        if SAMPLE_ID_RE.match(sample_id):
            add_report(reports, sample_id, "sample_id_format", "PASS", "INFO", "sample_id format is valid.", sample_id)
        else:
            add_report(
                reports,
                sample_id,
                "sample_id_format",
                "FAIL",
                "ERROR",
                "sample_id may contain only letters, numbers, dot, underscore, and hyphen.",
                sample_id,
            )
            critical_errors += 1

        r1_value = (row.get("fastq_r1") or "").strip()
        r2_value = (row.get("fastq_r2") or "").strip()
        r1_path = resolve_path(project_root, r1_value)
        r2_path = resolve_path(project_root, r2_value)

        for label, value, path in (("fastq_r1", r1_value, r1_path), ("fastq_r2", r2_value, r2_path)):
            if not value:
                add_report(reports, sample_id, f"{label}_present", "FAIL", "ERROR", f"{label} is empty.")
                critical_errors += 1
                continue
            if not is_fastq(path):
                add_report(reports, sample_id, f"{label}_suffix", "FAIL", "ERROR", f"{label} does not look like FASTQ.", value)
                critical_errors += 1
            else:
                add_report(reports, sample_id, f"{label}_suffix", "PASS", "INFO", f"{label} suffix is FASTQ-like.", value)
            if not path.exists():
                add_report(reports, sample_id, f"{label}_exists", "FAIL", "ERROR", f"{label} path does not exist.", str(path))
                critical_errors += 1
            elif not path.is_file():
                add_report(reports, sample_id, f"{label}_exists", "FAIL", "ERROR", f"{label} exists but is not a file.", str(path))
                critical_errors += 1
            else:
                add_report(reports, sample_id, f"{label}_exists", "PASS", "INFO", f"{label} exists.", str(path))

        r1_token = read_token(r1_path)
        r2_token = read_token(r2_path)
        sample_in_names = sample_id.lower() in r1_path.name.lower() and sample_id.lower() in r2_path.name.lower()
        if r1_token == "R1" and r2_token == "R2" and sample_in_names:
            add_report(reports, sample_id, "paired_reads", "PASS", "INFO", "R1/R2 filenames are paired and include sample_id.")
        else:
            add_report(
                reports,
                sample_id,
                "paired_reads",
                "FAIL",
                "ERROR",
                "Could not confirm paired R1/R2 filenames containing sample_id.",
                f"{r1_path.name};{r2_path.name}",
            )
            critical_errors += 1

    duplicates = sorted(sample_id for sample_id, count in seen.items() if count > 1)
    if duplicates:
        add_report(reports, "GLOBAL", "sample_id_unique", "FAIL", "ERROR", "Duplicate sample_id values found.", ",".join(duplicates))
        critical_errors += len(duplicates)
    else:
        add_report(reports, "GLOBAL", "sample_id_unique", "PASS", "INFO", "All sample_id values are unique.")

    return reports, critical_errors


def main() -> int:
    args = parse_args()
    samplesheet = Path(args.samplesheet)
    project_root = Path(args.project_root).resolve()
    output = Path(args.output)

    if not samplesheet.is_absolute():
        samplesheet = project_root / samplesheet
    if not output.is_absolute():
        output = project_root / output

    if not samplesheet.exists():
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["sample_id", "check", "status", "severity", "message", "value"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(
                {
                    "sample_id": "GLOBAL",
                    "check": "samplesheet_exists",
                    "status": "FAIL",
                    "severity": "ERROR",
                    "message": "Samplesheet does not exist.",
                    "value": str(samplesheet),
                }
            )
        print(f"Samplesheet does not exist: {samplesheet}", file=sys.stderr)
        return 1

    with samplesheet.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    reports, critical_errors = validate_rows(rows, fieldnames, project_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "check", "status", "severity", "message", "value"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(reports)

    if critical_errors:
        print(f"Samplesheet validation failed with {critical_errors} critical error(s). See {output}", file=sys.stderr)
        return 1
    print(f"Samplesheet validation passed for {len(rows)} sample(s). Report: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
