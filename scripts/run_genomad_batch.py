#!/usr/bin/env python3
"""Run geNomad end-to-end for downstream-retained local assemblies."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def load_samples(path: Path) -> list[str]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    samples = [
        row["sample_id"]
        for row in rows
        if row.get("downstream_inclusion_status_revised", "included_high_confidence_downstream")
        == "included_high_confidence_downstream"
    ]
    if not samples:
        raise SystemExit(f"No downstream-retained samples in {path}")
    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, default=Path("results/tables/table1_high_confidence_isolates.tsv"))
    parser.add_argument("--assemblies-dir", type=Path, default=Path("results/assemblies"))
    parser.add_argument("--database", type=Path, default=Path("data/databases/genomad_latest/genomad_db"))
    parser.add_argument("--outdir", type=Path, default=Path("results/mobilome/genomad"))
    parser.add_argument("--logs-dir", type=Path, default=Path("results/logs/genomad"))
    parser.add_argument("--summary", type=Path, default=Path("results/tables/genomad_run_summary.tsv"))
    parser.add_argument("--genomad", default="genomad")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--splits", type=int, default=8)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    samples = load_samples(args.samples)
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.logs_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for sample_id in samples:
        assembly = args.assemblies_dir / sample_id / f"{sample_id}.assembly.fasta"
        sample_out = args.outdir / sample_id
        prefix = f"{sample_id}.assembly"
        plasmid_summary = sample_out / f"{prefix}_summary" / f"{prefix}_plasmid_summary.tsv"
        virus_summary = sample_out / f"{prefix}_summary" / f"{prefix}_virus_summary.tsv"
        provirus_summary = sample_out / f"{prefix}_find_proviruses" / f"{prefix}_provirus.tsv"
        log_path = args.logs_dir / f"{sample_id}.genomad.log"

        if not assembly.exists():
            rows.append(
                {
                    "sample_id": sample_id,
                    "status": "missing_assembly",
                    "returncode": "",
                    "assembly": str(assembly),
                    "output_dir": str(sample_out),
                    "plasmid_summary": str(plasmid_summary),
                    "virus_summary": str(virus_summary),
                    "provirus_summary": str(provirus_summary),
                    "log": str(log_path),
                }
            )
            continue

        if plasmid_summary.exists() and virus_summary.exists() and provirus_summary.exists() and not args.force:
            rows.append(
                {
                    "sample_id": sample_id,
                    "status": "already_done",
                    "returncode": 0,
                    "assembly": str(assembly),
                    "output_dir": str(sample_out),
                    "plasmid_summary": str(plasmid_summary),
                    "virus_summary": str(virus_summary),
                    "provirus_summary": str(provirus_summary),
                    "log": str(log_path),
                }
            )
            continue

        cmd = [
            args.genomad,
            "end-to-end",
            "--cleanup",
            "--conservative",
            "--splits",
            str(args.splits),
            "--threads",
            str(args.threads),
            str(assembly),
            str(sample_out),
            str(args.database),
        ]
        with log_path.open("w") as log:
            proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)

        rows.append(
            {
                "sample_id": sample_id,
                "status": "ok" if proc.returncode == 0 else "failed",
                "returncode": proc.returncode,
                "assembly": str(assembly),
                "output_dir": str(sample_out),
                "plasmid_summary": str(plasmid_summary),
                "virus_summary": str(virus_summary),
                "provirus_summary": str(provirus_summary),
                "log": str(log_path),
            }
        )

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "status",
        "returncode",
        "assembly",
        "output_dir",
        "plasmid_summary",
        "virus_summary",
        "provirus_summary",
        "log",
    ]
    with args.summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    failed = [row for row in rows if row["status"] in {"failed", "missing_assembly"}]
    if failed:
        raise SystemExit(f"geNomad failed for {len(failed)} sample(s); see {args.summary}")


if __name__ == "__main__":
    main()
