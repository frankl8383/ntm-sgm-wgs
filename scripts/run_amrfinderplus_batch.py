#!/usr/bin/env python3
"""Run AMRFinderPlus for selected local NTM genomes."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-table", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--logdir", required=True, type=Path)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--amrfinder-bin", default="amrfinder")
    return parser.parse_args()


def run_one(args: argparse.Namespace, row: pd.Series) -> None:
    sample_id = str(row["sample_id"])
    protein_fasta = Path(row["protein_fasta"])
    if not protein_fasta.exists():
        raise FileNotFoundError(protein_fasta)

    out_tsv = args.outdir / f"{sample_id}.amrfinderplus.tsv"
    out_prot = args.outdir / f"{sample_id}.amrfinderplus.proteins.faa"
    log_file = args.logdir / f"{sample_id}.amrfinderplus.log"
    cmd = [
        args.amrfinder_bin,
        "--protein",
        str(protein_fasta),
        "--plus",
        "--name",
        sample_id,
        "--threads",
        str(args.threads),
        "--output",
        str(out_tsv),
        "--protein_output",
        str(out_prot),
        "--log",
        str(log_file),
    ]
    if args.database:
        cmd.extend(["--database", str(args.database)])
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.logdir.mkdir(parents=True, exist_ok=True)

    samples = pd.read_csv(args.sample_table, sep="\t")
    if samples.empty:
        raise SystemExit(f"No samples in {args.sample_table}")
    for row in samples.itertuples(index=False):
        run_one(args, pd.Series(row._asdict()))


if __name__ == "__main__":
    main()
