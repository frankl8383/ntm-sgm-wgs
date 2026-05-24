#!/usr/bin/env python3
"""Parse all-vs-reference FastANI output into best-hit tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastani", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--best-out", required=True, type=Path)
    parser.add_argument("--all-out", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cols = ["query_fasta", "reference_fasta", "ani", "matched_fragments", "total_fragments"]
    hits = pd.read_csv(args.fastani, sep="\t", names=cols)
    hits["alignment_fraction"] = hits["matched_fragments"] / hits["total_fragments"]
    meta = pd.read_csv(args.metadata, sep="\t", keep_default_na=False)
    local = meta[meta["sample_type"].eq("local")][["sample_id", "context_fasta", "species"]].rename(
        columns={"sample_id": "sample_id", "context_fasta": "query_fasta", "species": "local_species"}
    )
    ref = meta[meta["sample_type"].eq("public")][
        [
            "sample_id",
            "accession",
            "species",
            "strain",
            "assembly_level",
            "source_database",
            "source_dataset",
            "type_material",
            "context_fasta",
        ]
    ].rename(
        columns={
            "sample_id": "reference_sample_id",
            "context_fasta": "reference_fasta",
            "species": "fastani_best_hit_species",
            "accession": "fastani_best_hit_accession",
            "strain": "fastani_best_hit_strain",
        }
    )
    annotated = hits.merge(local, on="query_fasta", how="left").merge(ref, on="reference_fasta", how="left")
    annotated = annotated.sort_values(
        ["sample_id", "ani", "alignment_fraction", "matched_fragments"],
        ascending=[True, False, False, False],
    )
    args.all_out.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(args.all_out, sep="\t", index=False)

    best = annotated.dropna(subset=["sample_id"]).groupby("sample_id", as_index=False).head(1).copy()
    best = best.rename(
        columns={
            "ani": "fastani_best_hit_ani",
            "alignment_fraction": "fastani_best_hit_alignment_fraction",
        }
    )
    best_cols = [
        "sample_id",
        "query_fasta",
        "local_species",
        "fastani_best_hit_accession",
        "fastani_best_hit_species",
        "fastani_best_hit_strain",
        "fastani_best_hit_ani",
        "fastani_best_hit_alignment_fraction",
        "matched_fragments",
        "total_fragments",
        "assembly_level",
        "source_database",
        "source_dataset",
        "type_material",
        "reference_fasta",
    ]
    best[[c for c in best_cols if c in best.columns]].to_csv(args.best_out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
