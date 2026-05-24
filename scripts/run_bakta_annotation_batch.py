#!/usr/bin/env python3
"""Run Bakta annotation for downstream-retained local MAC/SGM assemblies."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def read_samples(path: Path) -> list[str]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "sample_id" not in reader.fieldnames:
            raise SystemExit(f"Missing sample_id column in {path}")
        samples: list[str] = []
        for row in reader:
            status = row.get("downstream_inclusion_status_revised", "included_high_confidence_downstream")
            if status and status != "included_high_confidence_downstream":
                continue
            samples.append(row["sample_id"])
    if not samples:
        raise SystemExit(f"No downstream-retained samples found in {path}")
    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, default=Path("results/tables/table1_high_confidence_isolates.tsv"))
    parser.add_argument("--assemblies-dir", type=Path, default=Path("results/assemblies"))
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/annotation/bakta"))
    parser.add_argument("--logs-dir", type=Path, default=Path("results/logs/bakta"))
    parser.add_argument("--summary", type=Path, default=Path("results/tables/bakta_annotation_run_summary.tsv"))
    parser.add_argument("--bakta", default="bakta")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    samples = read_samples(args.samples)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.logs_dir.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for sample_id in samples:
        assembly = args.assemblies_dir / sample_id / f"{sample_id}.assembly.fasta"
        out_dir = args.output_dir / sample_id
        log_path = args.logs_dir / f"{sample_id}.bakta.log"
        json_path = out_dir / f"{sample_id}.json"
        gff_path = out_dir / f"{sample_id}.gff3"
        tsv_path = out_dir / f"{sample_id}.tsv"

        if not assembly.exists():
            rows.append(
                {
                    "sample_id": sample_id,
                    "assembly": str(assembly),
                    "status": "missing_assembly",
                    "returncode": "",
                    "output_dir": str(out_dir),
                    "gff3": str(gff_path),
                    "tsv": str(tsv_path),
                    "json": str(json_path),
                    "log": str(log_path),
                }
            )
            continue

        if json_path.exists() and gff_path.exists() and tsv_path.exists() and not args.force:
            rows.append(
                {
                    "sample_id": sample_id,
                    "assembly": str(assembly),
                    "status": "already_done",
                    "returncode": 0,
                    "output_dir": str(out_dir),
                    "gff3": str(gff_path),
                    "tsv": str(tsv_path),
                    "json": str(json_path),
                    "log": str(log_path),
                }
            )
            continue

        cmd = [
            args.bakta,
            "--db",
            str(args.db),
            "--output",
            str(out_dir),
            "--prefix",
            sample_id,
            "--threads",
            str(args.threads),
            "--genus",
            "Mycobacterium",
            str(assembly),
        ]
        if args.force:
            cmd.insert(1, "--force")

        with log_path.open("w") as log_handle:
            proc = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, text=True)

        rows.append(
            {
                "sample_id": sample_id,
                "assembly": str(assembly),
                "status": "ok" if proc.returncode == 0 else "failed",
                "returncode": proc.returncode,
                "output_dir": str(out_dir),
                "gff3": str(gff_path),
                "tsv": str(tsv_path),
                "json": str(json_path),
                "log": str(log_path),
            }
        )

    fieldnames = ["sample_id", "assembly", "status", "returncode", "output_dir", "gff3", "tsv", "json", "log"]
    with args.summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    failed = [row for row in rows if row["status"] in {"failed", "missing_assembly"}]
    if failed:
        raise SystemExit(f"Bakta failed or missing assembly for {len(failed)} sample(s); see {args.summary}")


if __name__ == "__main__":
    main()
