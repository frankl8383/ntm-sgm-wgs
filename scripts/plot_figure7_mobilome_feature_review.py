#!/usr/bin/env python3
"""Plot conservative mobilome feature-review heatmap."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch


SPECIES_COLORS = {
    "Mycobacterium avium": "#4C78A8",
    "Mycobacterium intracellulare": "#F58518",
    "Mycobacterium paraintracellulare": "#54A24B",
    "Mycobacterium colombiense": "#B279A2",
}


def short_clade(text: str) -> str:
    for species in SPECIES_COLORS:
        if species in str(text):
            return species.replace("Mycobacterium ", "M. ")
    return "MAC boundary"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=Path("results/tables/mobilome_feature_summary.tsv"))
    parser.add_argument("--genomad-summary", type=Path, default=Path("results/tables/genomad_mobilome_sample_summary.tsv"))
    parser.add_argument("--png", type=Path, default=Path("results/figures/Figure7_mobilome_feature_review.png"))
    parser.add_argument("--pdf", type=Path, default=Path("results/figures/Figure7_mobilome_feature_review.pdf"))
    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep="\t")
    if df.empty:
        raise SystemExit(f"No rows in {args.summary}")
    genomad = pd.read_csv(args.genomad_summary, sep="\t") if args.genomad_summary.exists() else pd.DataFrame()

    marker_cols = [
        "integrase_marker_count",
        "transposase_marker_count",
        "recombinase_marker_count",
        "relaxase_mobilization_marker_count",
        "conjugation_t4ss_marker_count",
        "phage_structural_marker_count",
        "plasmid_associated_marker_count",
    ]
    labels = [
        "Integrase",
        "Transposase/IS",
        "Recombinase",
        "Relaxase/mob",
        "Conjugation/T4SS",
        "Phage structural",
        "Plasmid-associated",
    ]
    for col in marker_cols:
        if col not in df.columns:
            df[col] = 0

    df["clade_short"] = df["public_context_clade"].map(short_clade)
    clade_order = {
        "M. avium": 0,
        "M. intracellulare": 1,
        "M. paraintracellulare": 2,
        "M. colombiense": 3,
        "MAC boundary": 4,
    }
    df = df.sort_values(["clade_short", "sample_id"], key=lambda s: s.map(clade_order).fillna(9) if s.name == "clade_short" else s)
    if not genomad.empty and "sample_id" in genomad.columns:
        keep_cols = [
            "sample_id",
            "n_conservative_plasmids",
            "n_conservative_viruses",
            "n_find_provirus_regions",
        ]
        for col in keep_cols:
            if col not in genomad.columns:
                genomad[col] = 0
        df = df.merge(genomad[keep_cols], on="sample_id", how="left")
    for col in ["n_conservative_plasmids", "n_conservative_viruses", "n_find_provirus_regions"]:
        if col not in df.columns:
            df[col] = 0
        df[col] = df[col].fillna(0).astype(float)

    matrix = df[marker_cols].fillna(0).astype(float).to_numpy()
    display = np.log1p(matrix)
    cmap = LinearSegmentedColormap.from_list("mobilome_blues", ["#F7FBFF", "#9ECAE1", "#08519C"])

    fig = plt.figure(figsize=(18.2, 7.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[3.2, 8.0, 3.0, 2.2], wspace=0.18)
    ax_meta = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_genomad = fig.add_subplot(gs[0, 2])
    ax_bar = fig.add_subplot(gs[0, 3])

    y = np.arange(len(df))
    clade_colors = []
    for clade in df["clade_short"]:
        full = clade.replace("M. ", "Mycobacterium ")
        clade_colors.append(SPECIES_COLORS.get(full, "#8C8C8C"))

    ax_meta.set_xlim(0, 1.6)
    ax_meta.set_ylim(-0.5, len(df) - 0.5)
    ax_meta.invert_yaxis()
    for i, (_, row) in enumerate(df.iterrows()):
        ax_meta.add_patch(plt.Rectangle((0.05, i - 0.36), 0.12, 0.72, facecolor=clade_colors[i], edgecolor="none"))
        ax_meta.text(0.22, i, row["sample_id"], va="center", ha="left", fontsize=9, fontweight="bold")
        ax_meta.text(0.62, i, row["clade_short"], va="center", ha="left", fontsize=7)
    ax_meta.set_axis_off()

    im = ax.imshow(display, aspect="auto", cmap=cmap, vmin=0)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_title("Bakta annotation keyword screen", fontsize=11, fontweight="bold", pad=10)
    ax.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(df), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = int(matrix[i, j])
            if val > 0:
                ax.text(j, i, str(val), ha="center", va="center", fontsize=7, color="black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.012)
    cbar.set_label("log1p count", fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    genomad_cols = [
        "n_conservative_plasmids",
        "n_conservative_viruses",
        "n_find_provirus_regions",
    ]
    genomad_labels = ["Plasmid", "Virus", "Provirus"]
    genomad_matrix = df[genomad_cols].fillna(0).astype(float).to_numpy()
    genomad_display = np.log1p(genomad_matrix)
    genomad_cmap = LinearSegmentedColormap.from_list("genomad_gold", ["#FFFDF2", "#F6D36B", "#9A5B00"])
    im2 = ax_genomad.imshow(genomad_display, aspect="auto", cmap=genomad_cmap, vmin=0)
    ax_genomad.set_yticks(y)
    ax_genomad.set_yticklabels([])
    ax_genomad.set_xticks(np.arange(len(genomad_labels)))
    ax_genomad.set_xticklabels(genomad_labels, rotation=35, ha="right", fontsize=8)
    ax_genomad.set_title("geNomad conservative calls", fontsize=10, fontweight="bold", pad=10)
    ax_genomad.set_xticks(np.arange(-0.5, len(genomad_labels), 1), minor=True)
    ax_genomad.set_yticks(np.arange(-0.5, len(df), 1), minor=True)
    ax_genomad.grid(which="minor", color="white", linewidth=1.0)
    ax_genomad.tick_params(which="minor", bottom=False, left=False)
    for i in range(genomad_matrix.shape[0]):
        for j in range(genomad_matrix.shape[1]):
            val = int(genomad_matrix[i, j])
            if val > 0:
                ax_genomad.text(j, i, str(val), ha="center", va="center", fontsize=7, color="black")
    cbar2 = fig.colorbar(im2, ax=ax_genomad, fraction=0.052, pad=0.018)
    cbar2.set_label("log1p count", fontsize=8)
    cbar2.ax.tick_params(labelsize=8)

    totals = df["total_mobilome_marker_count"].fillna(0).astype(float)
    ax_bar.barh(y, totals, color="#4C78A8", edgecolor="white", linewidth=0.8)
    ax_bar.set_ylim(-0.5, len(df) - 0.5)
    ax_bar.invert_yaxis()
    ax_bar.set_yticks([])
    ax_bar.set_xlabel("Total marker\ncount", fontsize=8)
    ax_bar.set_ylabel("")
    ax_bar.set_title("Per-isolate burden", fontsize=10, fontweight="bold")
    ax_bar.spines[["top", "right", "left"]].set_visible(False)
    ax_bar.tick_params(axis="x", labelsize=8)

    legend_handles = [Patch(facecolor=color, label=species.replace("Mycobacterium ", "M. ")) for species, color in SPECIES_COLORS.items()]
    fig.legend(handles=legend_handles, loc="lower center", ncol=4, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle(
        "Figure 7. Conservative mobilome-associated feature review",
        y=0.985,
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.935,
        "Bakta/geNomad short-read draft-assembly screen; counts are not evidence for complete plasmids, horizontal transfer, or transmission.",
        ha="center",
        fontsize=9,
        color="#444444",
    )
    args.png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.png, dpi=300, bbox_inches="tight")
    fig.savefig(args.pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
