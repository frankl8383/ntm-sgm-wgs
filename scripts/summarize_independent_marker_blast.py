#!/usr/bin/env python3
"""Summarize independent type-panel marker BLAST hits for local assemblies."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]


def sample_id_from_path(path: Path) -> str:
    name = path.name
    for suffix in [".independent_marker_blast.tsv", ".blast.tsv", ".tsv"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast-dir", required=True, type=Path)
    parser.add_argument("--marker-metadata", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--output-top-hits", required=True, type=Path)
    parser.add_argument("--output-wide", required=True, type=Path)
    args = parser.parse_args()

    marker_meta = pd.read_csv(args.marker_metadata, sep="\t")
    marker_meta = marker_meta.rename(
        columns={"organism_name": "marker_reference_species", "length": "marker_reference_length"}
    )

    rows: list[pd.DataFrame] = []
    for path in sorted(args.blast_dir.glob("*.tsv")):
        if path.stat().st_size == 0:
            continue
        df = pd.read_csv(path, sep="\t", names=BLAST_COLUMNS)
        if df.empty:
            continue
        df["sample_id"] = sample_id_from_path(path)
        rows.append(df)

    if rows:
        hits = pd.concat(rows, ignore_index=True)
        hits = hits.merge(marker_meta, left_on="sseqid", right_on="marker_sequence_id", how="left")
        hits["subject_coverage"] = hits["length"].astype(float) / hits["slen"].astype(float)
        hits["marker"] = hits["marker"].fillna("unknown_marker")
        hits["marker_reference_species"] = hits["marker_reference_species"].fillna("NA")
        hits["accession"] = hits["accession"].fillna("NA")
        hits["strain"] = hits["strain"].fillna("NA")
        hits = hits.sort_values(
            ["sample_id", "marker", "bitscore", "pident", "subject_coverage", "evalue"],
            ascending=[True, True, False, False, False, True],
        )
        top = hits.groupby(["sample_id", "marker"], as_index=False, sort=False).head(args.top_n).copy()
        top["rank_within_marker"] = top.groupby(["sample_id", "marker"]).cumcount() + 1
        top["marker_hit_label"] = (
            top["marker_reference_species"].astype(str)
            + "|"
            + top["accession"].astype(str)
            + "|"
            + top["strain"].astype(str)
        )
    else:
        top = pd.DataFrame(
            columns=[
                "sample_id",
                "marker",
                "rank_within_marker",
                "marker_hit_label",
                "marker_reference_species",
                "accession",
                "strain",
                "pident",
                "bitscore",
                "subject_coverage",
            ]
        )

    keep_cols = [
        "sample_id",
        "marker",
        "rank_within_marker",
        "marker_hit_label",
        "marker_reference_species",
        "accession",
        "strain",
        "pident",
        "length",
        "slen",
        "subject_coverage",
        "evalue",
        "bitscore",
        "qseqid",
        "sseqid",
    ]
    top = top[[c for c in keep_cols if c in top.columns]]

    best = top[top["rank_within_marker"].eq(1)].copy()
    wide = best.pivot(index="sample_id", columns="marker", values="marker_hit_label")
    wide.columns = [f"{col}_top1_hit" for col in wide.columns]
    wide = wide.reset_index()

    args.output_top_hits.parent.mkdir(parents=True, exist_ok=True)
    top.to_csv(args.output_top_hits, sep="\t", index=False)
    wide.to_csv(args.output_wide, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
