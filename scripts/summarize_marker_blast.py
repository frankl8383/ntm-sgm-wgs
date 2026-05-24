#!/usr/bin/env python3
"""Summarize marker gene BLAST hits for local assemblies."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast-dir", required=True, type=Path)
    parser.add_argument("--marker-metadata", required=True, type=Path)
    parser.add_argument("--output-long", required=True, type=Path)
    parser.add_argument("--output-wide", required=True, type=Path)
    args = parser.parse_args()

    cols = [
        "query_contig",
        "marker_sequence_id",
        "pident",
        "alignment_length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "query_length",
        "subject_length",
    ]
    frames = []
    for path in sorted(args.blast_dir.glob("*.blast.tsv")):
        sample_id = path.name.replace(".blast.tsv", "")
        if path.stat().st_size == 0:
            continue
        frame = pd.read_csv(path, sep="\t", names=cols)
        frame["sample_id"] = sample_id
        frames.append(frame)

    if frames:
        hits = pd.concat(frames, ignore_index=True)
    else:
        hits = pd.DataFrame(columns=cols + ["sample_id"])

    marker_meta = pd.read_csv(args.marker_metadata, sep="\t")
    hits = hits.merge(marker_meta, on="marker_sequence_id", how="left")
    hits["subject_coverage"] = hits["alignment_length"] / hits["subject_length"]
    hits = hits[hits["subject_coverage"] >= 0.8].copy()
    hits = hits.sort_values(
        ["sample_id", "marker", "pident", "bitscore", "subject_coverage"],
        ascending=[True, True, False, False, False],
    )
    best = hits.groupby(["sample_id", "marker"], as_index=False, sort=False).head(1).copy()
    best["marker_best_hit"] = (
        best["organism_name"].fillna("NA")
        + " | "
        + best["accession"].fillna("NA")
        + " | "
        + best["strain"].fillna("NA")
    )
    args.output_long.parent.mkdir(parents=True, exist_ok=True)
    best.to_csv(args.output_long, sep="\t", index=False)

    wide = best.pivot_table(
        index="sample_id",
        columns="marker",
        values=["organism_name", "accession", "pident", "subject_coverage"],
        aggfunc="first",
    )
    wide.columns = [f"{marker}_{field}" for field, marker in wide.columns]
    wide = wide.reset_index()
    wide.to_csv(args.output_wide, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
