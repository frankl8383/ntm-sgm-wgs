#!/usr/bin/env python3
"""Summarize AMRFinderPlus AMR/STRESS/VIRULENCE results across samples."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-table", required=True, type=Path)
    parser.add_argument("--amrfinder-dir", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--sample-summary", required=True, type=Path)
    return parser.parse_args()


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip().lower().replace(" ", "_").replace("/", "_") for c in df.columns]
    return df


def load_hits(sample_table: pd.DataFrame, amrfinder_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for sample_id in sample_table["sample_id"].astype(str):
        path = amrfinder_dir / f"{sample_id}.amrfinderplus.tsv"
        if not path.exists() or path.stat().st_size == 0:
            continue
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        df = normalize_columns(df)
        if "name" in df.columns:
            df = df.rename(columns={"name": "sample_id"})
        else:
            df.insert(0, "sample_id", sample_id)
        df["sample_id"] = df["sample_id"].fillna(sample_id).astype(str)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def pick_first_available(row: pd.Series, columns: list[str]) -> str:
    for col in columns:
        value = row.get(col, "")
        if pd.notna(value) and str(value).strip():
            return str(value)
    return "unclassified_element"


def parse_node_metrics(protein_id: str) -> tuple[float | None, float | None]:
    """Parse Shovill/SPAdes-style NODE length and coverage from a Prodigal protein id."""
    text = str(protein_id)
    length_match = re.search(r"_length_([0-9.]+)", text)
    cov_match = re.search(r"_cov_([0-9.]+)", text)
    length = float(length_match.group(1)) if length_match else None
    coverage = float(cov_match.group(1)) if cov_match else None
    return length, coverage


def main() -> None:
    args = parse_args()
    sample_table = pd.read_csv(args.sample_table, sep="\t")
    hits = load_hits(sample_table, args.amrfinder_dir)

    for out in [args.summary, args.matrix, args.sample_summary]:
        out.parent.mkdir(parents=True, exist_ok=True)

    if hits.empty:
        empty_cols = [
            "sample_id",
            "amr_module_species",
            "element_class",
            "element_name",
            "method",
            "scope",
            "sequence_name",
            "percent_identity",
            "coverage_of_reference_sequence",
        ]
        pd.DataFrame(columns=empty_cols).to_csv(args.summary, sep="\t", index=False)
        pd.DataFrame({"sample_id": sample_table["sample_id"]}).to_csv(args.matrix, sep="\t", index=False)
        sample_table.assign(
            amrfinderplus_total_hits=0,
            amrfinderplus_amr_hits=0,
            amrfinderplus_stress_hits=0,
            amrfinderplus_virulence_hits=0,
        ).to_csv(args.sample_summary, sep="\t", index=False)
        return

    hits = hits.merge(
        sample_table[["sample_id", "amr_module_species"]],
        on="sample_id",
        how="left",
    )
    hits["element_class"] = hits.apply(
        lambda row: pick_first_available(row, ["element_type", "type", "subtype", "scope"]),
        axis=1,
    )
    hits["element_class"] = hits["element_class"].str.upper().replace({"PLUS": "PLUS"})
    hits["element_name"] = hits.apply(
        lambda row: pick_first_available(row, ["gene_symbol", "element_symbol", "allele", "hmm_id", "protein_identifier"]),
        axis=1,
    )
    hits["element_name"] = hits["element_name"].str.replace(r"\s+", "_", regex=True)
    hits["element_label"] = hits["element_class"] + ":" + hits["element_name"]
    metrics = hits.get("protein_id", pd.Series([""] * len(hits))).map(parse_node_metrics)
    hits["contig_length_from_id"] = [value[0] for value in metrics]
    hits["contig_coverage_from_id"] = [value[1] for value in metrics]
    method_upper = hits.get("method", pd.Series([""] * len(hits))).fillna("").astype(str).str.upper()
    hits["hit_warning"] = "none"
    hits.loc[method_upper.str.contains("PARTIAL", na=False), "hit_warning"] = "partial_hit"
    hits.loc[hits["contig_coverage_from_id"].fillna(9999) < 5, "hit_warning"] = (
        hits.loc[hits["contig_coverage_from_id"].fillna(9999) < 5, "hit_warning"]
        .replace("none", "")
        .map(lambda value: f"{value};low_contig_coverage".strip(";"))
    )
    hits.loc[hits["contig_length_from_id"].fillna(999999) < 1000, "hit_warning"] = (
        hits.loc[hits["contig_length_from_id"].fillna(999999) < 1000, "hit_warning"]
        .replace("none", "")
        .map(lambda value: f"{value};short_contig".strip(";"))
    )
    hits["hit_warning"] = hits["hit_warning"].replace("", "none")

    preferred_cols = [
        "sample_id",
        "amr_module_species",
        "element_class",
        "element_name",
        "element_label",
        "method",
        "scope",
        "sequence_name",
        "protein_identifier",
        "gene_symbol",
        "element_symbol",
        "class",
        "subclass",
        "percent_identity",
        "coverage_of_reference_sequence",
        "alignment_length",
        "accession_of_closest_sequence",
        "name_of_closest_sequence",
        "hmm_id",
        "hmm_description",
        "contig_length_from_id",
        "contig_coverage_from_id",
        "hit_warning",
    ]
    cols = [col for col in preferred_cols if col in hits.columns]
    hits[cols].to_csv(args.summary, sep="\t", index=False)

    matrix = (
        hits.assign(present=1)
        .pivot_table(index="sample_id", columns="element_label", values="present", aggfunc="max", fill_value=0)
        .reset_index()
    )
    matrix = sample_table[["sample_id", "amr_module_species"]].merge(matrix, on="sample_id", how="left").fillna(0)
    element_cols = [c for c in matrix.columns if c not in {"sample_id", "amr_module_species"}]
    matrix[element_cols] = matrix[element_cols].astype(int)
    matrix.to_csv(args.matrix, sep="\t", index=False)

    sample_summary = sample_table.copy()
    counts = hits.groupby(["sample_id", "element_class"]).size().unstack(fill_value=0)
    sample_summary = sample_summary.merge(counts, left_on="sample_id", right_index=True, how="left").fillna(0)
    sample_summary["amrfinderplus_total_hits"] = hits.groupby("sample_id").size().reindex(sample_summary["sample_id"]).fillna(0).astype(int).values
    for cls in ["AMR", "STRESS", "VIRULENCE"]:
        source_col = cls if cls in sample_summary.columns else None
        target_col = f"amrfinderplus_{cls.lower()}_hits"
        sample_summary[target_col] = sample_summary[source_col].astype(int) if source_col else 0
    sample_summary.to_csv(args.sample_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
