#!/usr/bin/env python3
"""Select high-confidence local genomes for AMR/virulence/mobilome screening."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conflict-table", required=True, type=Path)
    parser.add_argument("--gene-call-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.conflict_table, sep="\t")
    required = {
        "sample_id",
        "current_final_norm",
        "downstream_inclusion_status_revised",
        "gunc_status",
        "checkm2_contamination",
    }
    missing = required.difference(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {args.conflict_table}: {sorted(missing)}")

    included = df[df["downstream_inclusion_status_revised"].eq("included_high_confidence_downstream")].copy()
    included["protein_fasta"] = included["sample_id"].map(
        lambda sample_id: args.gene_call_dir / f"{sample_id}.genecalls.faa"
    )
    included["protein_fasta_exists"] = included["protein_fasta"].map(lambda p: p.exists())
    included["amr_module_species"] = included["current_final_norm"]
    included["amr_module_note"] = (
        "High-confidence downstream genome; screen with AMRFinderPlus --plus. "
        "Interpret NTM resistance only after species-level review and AST context."
    )

    out_cols = [
        "sample_id",
        "amr_module_species",
        "protein_fasta",
        "protein_fasta_exists",
        "gunc_status",
        "checkm2_contamination",
        "species_level_confidence_revised",
        "evidence_conflict_flag",
        "recommended_review_action",
        "amr_module_note",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    included[out_cols].to_csv(args.output, sep="\t", index=False)

    missing_fastas = included.loc[~included["protein_fasta_exists"], "sample_id"].tolist()
    if missing_fastas:
        raise SystemExit(f"Missing protein FASTA files for AMR module: {', '.join(missing_fastas)}")


if __name__ == "__main__":
    main()
