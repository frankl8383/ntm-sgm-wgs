#!/usr/bin/env python3
"""Plot local priority-isolate FastANI values against the type-strain panel."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SPECIES_ORDER = [
    "Mycobacterium avium",
    "Mycobacterium intracellulare",
    "Mycobacterium paraintracellulare",
    "Mycobacterium colombiense",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ani-all-hits", required=True, type=Path)
    parser.add_argument("--best-hits", required=True, type=Path)
    parser.add_argument("--conflict-table", required=True, type=Path)
    parser.add_argument("--output-png", required=True, type=Path)
    parser.add_argument("--output-pdf", required=True, type=Path)
    return parser.parse_args()


def short_taxon(name: str) -> str:
    if not name:
        return "NA"
    return (
        name.replace("Mycobacterium ", "M. ")
        .replace("Mycolicibacterium ", "Mycolicibacterium ")
        .replace(" subsp. ", " subsp. ")
    )


def ref_label(row: pd.Series) -> str:
    organism = short_taxon(str(row["organism_name"]))
    strain = str(row.get("strain", "")).strip()
    accession = str(row["accession"])
    if strain and strain.lower() != "nan":
        return f"{organism} | {strain} | {accession}"
    return f"{organism} | {accession}"


def main() -> int:
    args = parse_args()
    all_hits = pd.read_csv(args.ani_all_hits, sep="\t", keep_default_na=False)
    best_hits = pd.read_csv(args.best_hits, sep="\t", keep_default_na=False)
    conflicts = pd.read_csv(args.conflict_table, sep="\t", keep_default_na=False)

    best_accessions = best_hits["fastani_best_hit_accession"].drop_duplicates().tolist()
    extra_close = (
        all_hits[all_hits["ani"].astype(float) >= 98.0]
        .sort_values(["organism_name", "accession"])["accession"]
        .drop_duplicates()
        .tolist()
    )
    selected_accessions = list(dict.fromkeys(best_accessions + extra_close))

    selected = all_hits[all_hits["accession"].isin(selected_accessions)].copy()
    ref_meta = (
        selected[["accession", "organism_name", "strain"]]
        .drop_duplicates("accession")
        .assign(label=lambda df: df.apply(ref_label, axis=1))
    )
    selected = selected.merge(ref_meta[["accession", "label"]], on="accession", how="left")

    sample_order = (
        conflicts.assign(
            species_rank=lambda df: df["current_final_norm"].map(
                {
                    "avium": 0,
                    "intracellulare": 1,
                    "paraintracellulare": 2,
                    "colombiense": 3,
                    "chimaera": 4,
                    "yongonense": 5,
                    "unresolved": 6,
                }
            ).fillna(99)
        )
        .sort_values(["species_rank", "sample_id"])["sample_id"]
        .tolist()
    )
    ref_order = (
        ref_meta.assign(
            taxon_key=lambda df: df["organism_name"].str.replace("Mycobacterium ", "", regex=False),
            best_rank=lambda df: df["accession"].map({acc: i for i, acc in enumerate(best_accessions)}).fillna(999),
        )
        .sort_values(["taxon_key", "best_rank", "accession"])["label"]
        .tolist()
    )

    matrix = (
        selected.pivot_table(index="sample_id", columns="label", values="ani", aggfunc="max")
        .reindex(index=sample_order, columns=ref_order)
        .astype(float)
    )

    values = matrix.to_numpy()
    plotted_values = values.copy()
    plotted_values[plotted_values < 95.0] = np.nan
    masked = np.ma.masked_invalid(plotted_values)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#F2F2F2")

    fig_width = max(12, len(ref_order) * 0.32)
    fig_height = max(5.5, len(sample_order) * 0.32)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(masked, aspect="auto", interpolation="nearest", cmap=cmap, vmin=95.0, vmax=100.0)

    ax.set_xticks(np.arange(len(ref_order)))
    ax.set_xticklabels(ref_order, rotation=90, ha="center", fontsize=5.7)
    ax.set_yticks(np.arange(len(sample_order)))
    ax.set_yticklabels(sample_order, fontsize=8)
    ax.set_title(
        "Priority local isolates vs comprehensive MAC/SGM type-strain panel\n"
        "Cells below 95% ANI are blanked; values >=98% are printed.",
        loc="left",
        fontsize=11,
    )
    ax.set_xlabel("Type-material / representative reference genomes")
    ax.set_ylabel("Local priority isolate")

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            if np.isfinite(value) and value >= 98.0:
                ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=4.8, color="white")

    cbar = fig.colorbar(image, ax=ax, fraction=0.018, pad=0.012)
    cbar.set_label("FastANI (%)")
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_png, dpi=300, bbox_inches="tight")
    fig.savefig(args.output_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output_png} and {args.output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
