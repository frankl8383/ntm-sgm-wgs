#!/usr/bin/env python3
"""Create public-context Figures 4 and 5 for MAC analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


SPECIES_COLORS = {
    "Mycobacterium avium": "#0072B2",
    "Mycobacterium intracellulare": "#E69F00",
    "Mycobacterium paraintracellulare": "#009E73",
    "Mycobacterium colombiense": "#CC79A7",
}

SPECIES_SHORT = {
    "Mycobacterium avium": "M. avium",
    "Mycobacterium intracellulare": "M. intracellulare",
    "Mycobacterium paraintracellulare": "M. paraintracellulare",
    "Mycobacterium colombiense": "M. colombiense",
}

GROUP_TITLES = {
    "M_avium": "M. avium",
    "M_colombiense": "M. colombiense",
    "M_intracellulare": "M. intracellulare",
    "M_paraintracellulare": "M. paraintracellulare",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--fastani-best", required=True, type=Path)
    parser.add_argument("--mashtree", required=True, type=Path)
    parser.add_argument("--core-summary", required=True, type=Path)
    parser.add_argument("--core-dir", required=True, type=Path)
    parser.add_argument("--figure4-png", required=True, type=Path)
    parser.add_argument("--figure4-pdf", required=True, type=Path)
    parser.add_argument("--figure5-png", required=True, type=Path)
    parser.add_argument("--figure5-pdf", required=True, type=Path)
    parser.add_argument("--expanded-paraintra-tree", type=Path)
    parser.add_argument("--expanded-paraintra-metadata", type=Path)
    parser.add_argument("--expanded-paraintra-sites", type=int)
    return parser.parse_args()


def tree_label(path_or_name: str) -> str:
    return Path(path_or_name).name.replace(".fna", "")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="bottom",
        ha="right",
    )


def add_species_legend(fig: plt.Figure, loc: tuple[float, float]) -> None:
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=color, label=SPECIES_SHORT[species], markersize=7)
        for species, color in SPECIES_COLORS.items()
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            marker="*",
            color="w",
            markerfacecolor="white",
            markeredgecolor="black",
            label="Local high-confidence SGM",
            markersize=10,
        )
    )
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=loc, frameon=False, fontsize=9)


def draw_tree_panel(
    ax: plt.Axes,
    tree_path: Path,
    metadata: pd.DataFrame,
    title: str,
    xlabel: str,
    local_fontsize: float = 6.0,
    show_xaxis: bool = True,
) -> None:
    tree = Phylo.read(tree_path, "newick")
    meta = metadata.copy()
    meta["tree_label"] = meta["context_fasta"].map(tree_label)
    meta_by_label = meta.set_index("tree_label").to_dict("index")

    terminals = tree.get_terminals()
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda clade: "")
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    if show_xaxis:
        ax.set_xlabel(xlabel, fontsize=8)
    else:
        ax.set_xlabel("")
        ax.set_xticklabels([])
    ax.tick_params(axis="both", labelsize=7, length=2)
    ax.grid(False)

    y_positions = {terminal: idx + 1 for idx, terminal in enumerate(terminals)}
    max_x = max((tree.distance(terminal) for terminal in terminals), default=1.0)
    label_offset = max(max_x * 0.018, 0.001)
    label_right = max_x + (label_offset * 10)
    local_points: list[tuple[float, float, str]] = []
    for terminal in terminals:
        label = tree_label(terminal.name)
        info = meta_by_label.get(label, {})
        species = info.get("species", "NA")
        sample_type = info.get("sample_type", "public")
        sample_id = info.get("sample_id", label)
        color = SPECIES_COLORS.get(species, "#777777")
        x = tree.distance(terminal)
        y = y_positions[terminal]
        if sample_type == "local":
            ax.scatter([x], [y], c=["white"], edgecolors="black", s=58, marker="*", linewidths=0.8, zorder=5)
            local_points.append((x, y, sample_id))
        else:
            ax.scatter([x], [y], c=[color], edgecolors="white", s=14, marker="o", linewidths=0.25, zorder=4)
    dense_many_local = len(terminals) > 30 and len(local_points) > 3
    if dense_many_local:
        local_text = "Local isolates: " + ", ".join(sample_id for _, _, sample_id in sorted(local_points, key=lambda item: item[2]))
        ax.text(
            0.98,
            0.95,
            local_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=6.2,
            bbox=dict(facecolor="white", edgecolor="#999999", linewidth=0.4, alpha=0.85),
        )
    else:
        min_label_sep = 0.9
        previous_y: float | None = None
        for x, y, sample_id in sorted(local_points, key=lambda item: item[1]):
            text_y = y
            if previous_y is not None and text_y - previous_y < min_label_sep:
                text_y = previous_y + min_label_sep
            previous_y = text_y
            ax.plot([x, x + label_offset * 0.8], [y, text_y], color="#333333", lw=0.35, zorder=5)
            ax.text(x + label_offset, text_y, sample_id, va="center", ha="left", fontsize=local_fontsize, fontweight="bold")
    ax.set_xlim(left=0, right=label_right)


def figure4(args: argparse.Namespace, metadata: pd.DataFrame, fastani: pd.DataFrame) -> None:
    local_meta = metadata[metadata["sample_type"] == "local"].set_index("sample_id")
    fastani = fastani.copy()
    fastani["local_species"] = fastani["sample_id"].map(local_meta["species"])
    fastani["sample_order_key"] = fastani["local_species"].map(lambda s: list(SPECIES_COLORS).index(s))
    fastani = fastani.sort_values(["sample_order_key", "sample_id"]).reset_index(drop=True)

    fig = plt.figure(figsize=(14.5, 9.5), constrained_layout=False)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.05, 1.45], height_ratios=[1.35, 0.9], wspace=0.25, hspace=0.35)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[:, 1])

    y = np.arange(len(fastani))
    colors = fastani["local_species"].map(SPECIES_COLORS)
    sizes = 70 + 260 * (fastani["fastani_best_hit_alignment_fraction"].astype(float).clip(0.80, 1.0) - 0.80) / 0.20
    ax_a.scatter(fastani["fastani_best_hit_ani"], y, c=colors, s=sizes, edgecolors="black", linewidths=0.4)
    ax_a.axvline(95, color="#666666", linestyle="--", linewidth=0.8)
    ax_a.axvline(96, color="#999999", linestyle=":", linewidth=0.8)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(fastani["sample_id"], fontsize=8)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Best public FastANI (%)")
    ax_a.set_title("Nearest public genome per local isolate", fontsize=10, loc="left")
    xmin = 94.8
    ax_a.set_xlim(xmin, 100.15)
    ax_a.grid(axis="x", color="#D0D0D0", linewidth=0.5)
    for _, row in fastani.iterrows():
        label = f"{row['fastani_best_hit_accession']} ({row['fastani_best_hit_strain']})"
        ax_a.text(100.22, row.name, label, va="center", ha="left", fontsize=6.6, clip_on=False)
    ax_a.set_xlim(xmin, 100.7)
    panel_label(ax_a, "A")

    comp = metadata.groupby(["species", "sample_type"]).size().unstack(fill_value=0).reindex(SPECIES_COLORS.keys())
    x = np.arange(len(comp))
    width = 0.38
    ax_b.bar(
        x - width / 2,
        comp.get("local", 0),
        width=width,
        color="#333333",
        edgecolor="#333333",
        label="Local",
    )
    ax_b.bar(
        x + width / 2,
        comp.get("public", 0),
        width=width,
        color="#BDBDBD",
        edgecolor="#555555",
        linewidth=0.4,
        label="Public",
    )
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([SPECIES_SHORT[s] for s in comp.index], rotation=28, ha="right", fontsize=8)
    ax_b.set_ylabel("Genomes")
    ax_b.set_title("Public-context composition", fontsize=10, loc="left")
    ax_b.legend(frameon=False, fontsize=8)
    ax_b.grid(axis="y", color="#D0D0D0", linewidth=0.5)
    panel_label(ax_b, "B")

    draw_tree_panel(
        ax_c,
        args.mashtree,
        metadata,
        "MAC public-context MashTree overview",
        "Mash distance",
        local_fontsize=6.5,
        show_xaxis=True,
    )
    panel_label(ax_c, "C")
    add_species_legend(fig, (0.81, 0.98))
    fig.suptitle("Figure 4. Public-genome contextualization of high-confidence local MAC/SGM isolates", fontsize=13, y=0.995)
    args.figure4_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure4_png, dpi=300, bbox_inches="tight")
    fig.savefig(args.figure4_pdf, bbox_inches="tight")
    plt.close(fig)


def figure5(args: argparse.Namespace, metadata: pd.DataFrame, core_summary: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(14.5, 10.5), constrained_layout=False)
    gs = GridSpec(2, 2, figure=fig, wspace=0.24, hspace=0.32)
    groups = ["M_avium", "M_colombiense", "M_intracellulare", "M_paraintracellulare"]
    labels = ["A", "B", "C", "D"]
    for idx, group in enumerate(groups):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        title = GROUP_TITLES[group]
        row = core_summary[core_summary["species_group"] == group].iloc[0]
        subtitle = (
            f"{title}: {int(row['n_local_genomes'])} local + {int(row['n_public_genomes'])} public; "
            f"{int(row['alignment_sites']):,} alignment sites"
        )
        tree_path = args.core_dir / group / f"{group}.iqtree_gtr.treefile"
        panel_metadata = metadata
        xlabel = "ML branch length (SKA alignment; GTR)"
        if group == "M_paraintracellulare":
            if args.expanded_paraintra_tree and args.expanded_paraintra_metadata:
                tree_path = args.expanded_paraintra_tree
                panel_metadata = pd.read_csv(args.expanded_paraintra_metadata, sep="\t", keep_default_na=False)
                sites = args.expanded_paraintra_sites or int(row["alignment_sites"])
                subtitle = f"{title}: 4 local + 11 curated public; {sites:,} SNP-only sites"
                xlabel = "ML branch length (SNP-only SKA alignment; GTR+ASC)"
            else:
                asc_tree_path = args.core_dir / group / f"{group}.iqtree_gtr_asc.treefile"
                if asc_tree_path.exists():
                    tree_path = asc_tree_path
                    xlabel = "ML branch length (SNP-only SKA alignment; GTR+ASC)"
        if group == "M_paraintracellulare" and not (args.expanded_paraintra_tree and args.expanded_paraintra_metadata):
            asc_tree_path = args.core_dir / group / f"{group}.iqtree_gtr_asc.treefile"
            if asc_tree_path.exists():
                tree_path = asc_tree_path
                xlabel = "ML branch length (SNP-only SKA alignment; GTR+ASC)"

        draw_tree_panel(
            ax,
            tree_path,
            panel_metadata,
            subtitle,
            xlabel,
            local_fontsize=6.2,
            show_xaxis=True,
        )
        panel_label(ax, labels[idx])
    add_species_legend(fig, (0.77, 0.99))
    fig.suptitle("Figure 5. Species-specific core-SKA/IQ-TREE phylogenies", fontsize=13, y=0.995)
    args.figure5_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure5_png, dpi=300, bbox_inches="tight")
    fig.savefig(args.figure5_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    metadata = pd.read_csv(args.metadata, sep="\t", keep_default_na=False)
    fastani = pd.read_csv(args.fastani_best, sep="\t", keep_default_na=False)
    core_summary = pd.read_csv(args.core_summary, sep="\t", keep_default_na=False)
    figure4(args, metadata, fastani)
    figure5(args, metadata, core_summary)
    print(f"Wrote {args.figure4_png} and {args.figure4_pdf}")
    print(f"Wrote {args.figure5_png} and {args.figure5_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
