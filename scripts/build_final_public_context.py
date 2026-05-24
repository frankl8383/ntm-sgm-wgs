#!/usr/bin/env python3
"""Build the final MAC public-context metadata with expanded M. paraintracellulare."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-metadata", required=True, type=Path)
    parser.add_argument("--expanded-paraintra-metadata", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    base = pd.read_csv(args.base_metadata, sep="\t", keep_default_na=False)
    expanded = pd.read_csv(args.expanded_paraintra_metadata, sep="\t", keep_default_na=False)

    # Retain all local genomes from the base set, retain non-paraintracellulare public genomes
    # from the original context, and replace the old 3-public paraintracellulare background
    # with the curated expanded 11-public panel.
    base_keep = base[
        (base["sample_type"].eq("local"))
        | ~(
            base["sample_type"].eq("public")
            & base["species"].eq("Mycobacterium paraintracellulare")
        )
    ].copy()
    expanded_public = expanded[expanded["sample_type"].eq("public")].copy()

    combined = pd.concat([base_keep, expanded_public], ignore_index=True)
    combined = combined.drop_duplicates(subset=["sample_id"], keep="first")
    species_order = {
        "Mycobacterium avium": 0,
        "Mycobacterium colombiense": 1,
        "Mycobacterium intracellulare": 2,
        "Mycobacterium paraintracellulare": 3,
    }
    type_order = {"local": 0, "public": 1}
    combined["_species_order"] = combined["species"].map(species_order).fillna(99).astype(int)
    combined["_type_order"] = combined["sample_type"].map(type_order).fillna(9).astype(int)
    combined = combined.sort_values(["_species_order", "_type_order", "sample_id"]).drop(
        columns=["_species_order", "_type_order"]
    )

    missing = [p for p in combined["context_fasta"].astype(str) if not Path(p).exists()]
    if missing:
        raise FileNotFoundError("Missing context FASTA files:\n" + "\n".join(missing[:20]))

    combined.to_csv(args.outdir / "mac_public_context_metadata_final.tsv", sep="\t", index=False)
    (args.outdir / "mac_public_context_fasta_list_final.txt").write_text(
        "\n".join(combined["context_fasta"].astype(str)) + "\n", encoding="utf-8"
    )

    local = combined[combined["sample_type"].eq("local")]
    public = combined[combined["sample_type"].eq("public")]
    public_ref_list = public[["sample_id", "context_fasta"]].copy()
    public_ref_list.to_csv(args.outdir / "mac_public_context_public_reference_list_final.tsv", sep="\t", index=False)
    (args.outdir / "mac_public_context_public_fastani_reference_list_final.txt").write_text(
        "\n".join(public_ref_list["context_fasta"].astype(str)) + "\n", encoding="utf-8"
    )
    (args.outdir / "mac_public_context_local_query_list_final.txt").write_text(
        "\n".join(local["context_fasta"].astype(str)) + "\n", encoding="utf-8"
    )

    summary = (
        combined.groupby(["species", "sample_type"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .rename(columns={"local": "n_local", "public": "n_public"})
    )
    for col in ["n_local", "n_public"]:
        if col not in summary.columns:
            summary[col] = 0
    summary["n_total"] = summary["n_local"] + summary["n_public"]
    summary.to_csv(args.outdir / "mac_public_context_summary_final.tsv", sep="\t", index=False)

    overall = pd.DataFrame(
        [
            {"metric": "local_genomes", "value": int(local.shape[0])},
            {"metric": "public_genomes", "value": int(public.shape[0])},
            {"metric": "combined_context_genomes", "value": int(combined.shape[0])},
        ]
    )
    overall.to_csv(args.outdir / "mac_public_context_counts_final.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
