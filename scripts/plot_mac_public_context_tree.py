#!/usr/bin/env python3
"""Plot a MAC public-context tree with local isolates highlighted."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo


SPECIES_COLORS = {
    "Mycobacterium avium": "#4C78A8",
    "Mycobacterium intracellulare": "#F58518",
    "Mycobacterium paraintracellulare": "#54A24B",
    "Mycobacterium colombiense": "#B279A2",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-png", required=True, type=Path)
    parser.add_argument("--output-pdf", required=True, type=Path)
    parser.add_argument("--title", default="MAC public-context tree with local high-confidence SGM isolates")
    parser.add_argument("--xlabel", default="Branch length")
    return parser.parse_args()


def leaf_key(name: str) -> str:
    return Path(name).name.replace(".fna", "")


def main() -> int:
    args = parse_args()
    tree = Phylo.read(args.tree, "newick")
    metadata = pd.read_csv(args.metadata, sep="\t", keep_default_na=False)
    metadata["tree_label"] = metadata["context_fasta"].map(leaf_key)
    meta = metadata.set_index("tree_label").to_dict("index")

    terminals = tree.get_terminals()
    n_leaves = len(terminals)
    fig_height = max(9, n_leaves * 0.07)
    fig, ax = plt.subplots(figsize=(13, fig_height))
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda clade: "")
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel("")
    ax.set_title(args.title)

    y_positions = {terminal: idx + 1 for idx, terminal in enumerate(terminals)}
    max_x = max((tree.distance(terminal) for terminal in terminals), default=1.0)
    label_x = max_x * 1.02 if max_x > 0 else 0.02

    for terminal in terminals:
        label = leaf_key(terminal.name)
        info = meta.get(label, {})
        sample = info.get("sample_id", label)
        species = info.get("species", "NA")
        sample_type = info.get("sample_type", "public")
        color = SPECIES_COLORS.get(species, "#777777")
        y = y_positions[terminal]
        x = tree.distance(terminal)
        marker = "*" if sample_type == "local" else "o"
        size = 65 if sample_type == "local" else 18
        ax.scatter([x], [y], c=[color], s=size, marker=marker, edgecolors="black", linewidths=0.35, zorder=5)
        if sample_type == "local":
            ax.text(label_x, y, sample, va="center", fontsize=7.5, fontweight="bold", color="black")

    handles = []
    for species, color in SPECIES_COLORS.items():
        handles.append(
            plt.Line2D([0], [0], marker="o", color="w", label=species, markerfacecolor=color, markersize=7)
        )
    handles.append(
        plt.Line2D([0], [0], marker="*", color="w", label="local high-confidence SGM", markerfacecolor="white", markeredgecolor="black", markersize=10)
    )
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    ax.set_xlim(left=0, right=label_x * 1.18)
    fig.tight_layout()
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_png, dpi=300)
    fig.savefig(args.output_pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
