#!/usr/bin/env python3
"""Plot a compact AMRFinderPlus element heatmap for local high-confidence isolates."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap


SPECIES_COLORS = {
    "Mycobacterium avium": "#0072B2",
    "Mycobacterium intracellulare": "#E69F00",
    "Mycobacterium paraintracellulare": "#009E73",
    "Mycobacterium colombiense": "#CC79A7",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--max-elements", type=int, default=45)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    matrix = pd.read_csv(args.matrix, sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")
    element_cols = [c for c in matrix.columns if c not in {"sample_id", "amr_module_species"}]

    args.png.parent.mkdir(parents=True, exist_ok=True)
    args.pdf.parent.mkdir(parents=True, exist_ok=True)

    if not element_cols:
        fig, ax = plt.subplots(figsize=(7.5, 3.0))
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No AMRFinderPlus AMR/STRESS/VIRULENCE hits detected",
            ha="center",
            va="center",
            fontsize=11,
        )
        fig.savefig(args.png, dpi=300, bbox_inches="tight")
        fig.savefig(args.pdf, bbox_inches="tight")
        plt.close(fig)
        return

    prevalence = matrix[element_cols].sum(axis=0).sort_values(ascending=False)
    selected = prevalence.head(args.max_elements).index.tolist()
    plot_df = matrix[["sample_id", "amr_module_species", *selected]].copy()
    plot_df = plot_df.sort_values(["amr_module_species", "sample_id"])

    data = plot_df[selected].to_numpy(dtype=float)
    fig_height = max(4.6, 0.32 * len(plot_df) + 1.5)
    fig_width = max(6.2, 0.8 * len(selected) + 3.2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    cmap = ListedColormap(["#f2f2f2", "#222222"])
    ax.imshow(data, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)

    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["sample_id"], fontsize=8)
    ax.set_xticks(range(len(selected)))
    warning_labels: set[str] = set()
    if "element_label" in summary.columns and "hit_warning" in summary.columns:
        warning_labels = set(summary.loc[summary["hit_warning"].fillna("none").ne("none"), "element_label"].astype(str))
    short_labels = []
    for label in selected:
        short = label.replace("VIRULENCE:", "V:").replace("STRESS:", "S:").replace("AMR:", "A:")
        if label in warning_labels:
            short = f"{short}*"
        short_labels.append(short)
    ax.set_xticklabels(short_labels, fontsize=8, rotation=45, ha="right")
    ax.set_xlabel("AMRFinderPlus element", fontsize=9)
    ax.set_ylabel("")
    ax.set_title("AMRFinderPlus --plus screening of high-confidence local MAC/SGM isolates", loc="left", fontsize=10)
    ax.set_xticks([x - 0.5 for x in range(1, len(selected))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(plot_df))], minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    for y, species in enumerate(plot_df["amr_module_species"]):
        ax.scatter(
            [-0.75],
            [y],
            s=45,
            marker="s",
            color=SPECIES_COLORS.get(species, "#777777"),
            clip_on=False,
            edgecolor="white",
            linewidth=0.3,
        )
    ax.set_xlim(-1.1, len(selected) - 0.5)

    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=color, label=species.replace("Mycobacterium ", "M. "))
        for species, color in SPECIES_COLORS.items()
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.text(
        0,
        -0.28,
        "Black cells indicate element presence. Asterisked elements have low-confidence hit warnings. These are genomic features, not clinical resistance predictions without AST/species-specific review.",
        transform=ax.transAxes,
        fontsize=8,
        va="top",
    )
    fig.tight_layout()
    fig.savefig(args.png, dpi=300, bbox_inches="tight")
    fig.savefig(args.pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
