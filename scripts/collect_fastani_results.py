#!/usr/bin/env python3
"""Summarize FastANI all-vs-reference results."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def sample_id_from_path(path: str) -> str:
    p = Path(path)
    name = p.name
    for suffix in [".assembly.fasta", ".fasta", ".fa", ".fna"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return p.stem


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastani", required=True, type=Path)
    parser.add_argument("--reference-metadata", required=True, type=Path)
    parser.add_argument("--all-output", required=True, type=Path)
    parser.add_argument("--best-output", required=True, type=Path)
    args = parser.parse_args()

    cols = ["query_fasta", "reference_fasta", "ani", "matched_fragments", "total_fragments"]
    df = pd.read_csv(args.fastani, sep="\t", names=cols)
    if df.empty:
        raise SystemExit(f"FastANI result is empty: {args.fastani}")

    df["sample_id"] = df["query_fasta"].map(sample_id_from_path)
    df["alignment_fraction"] = df["matched_fragments"] / df["total_fragments"]

    meta = pd.read_csv(args.reference_metadata, sep="\t")
    df = df.merge(meta, on="reference_fasta", how="left")
    ordered_cols = [
        "sample_id",
        "query_fasta",
        "accession",
        "organism_name",
        "strain",
        "ani",
        "alignment_fraction",
        "matched_fragments",
        "total_fragments",
        "assembly_level",
        "source_database",
        "source_dataset",
        "type_material",
        "reference_fasta",
    ]
    df = df[ordered_cols].sort_values(["sample_id", "ani", "alignment_fraction"], ascending=[True, False, False])

    best = df.groupby("sample_id", as_index=False, sort=False).head(1).copy()
    best = best.rename(
        columns={
            "organism_name": "fastani_best_hit_species",
            "accession": "fastani_best_hit_accession",
            "strain": "fastani_best_hit_strain",
            "ani": "fastani_best_hit_ani",
            "alignment_fraction": "fastani_best_hit_alignment_fraction",
        }
    )

    args.all_output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.all_output, sep="\t", index=False)
    best.to_csv(args.best_output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
