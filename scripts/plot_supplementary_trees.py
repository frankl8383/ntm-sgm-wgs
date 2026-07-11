#!/usr/bin/env python3
"""Plot complete public-context trees for the upgraded supplement."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo


MM = 1 / 25.4
COLORS = {
    "branch": "#B7B7B7",
    "ink": "#242424",
    "muted": "#777777",
    "tmi": "#3E6FA3",
    "mp": "#B85C55",
    "chimaera": "#8A78A8",
    "yongonense": "#D19A32",
    "marseillense": "#66A99A",
    "avium": "#3E6FA3",
    "timonense": "#D19A32",
    "colombiense": "#9A6D9E",
    "other": "#A9A9A9",
    "anchor": "#202020",
}


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.7,
            "axes.titlesize": 7.6,
            "axes.labelsize": 6.8,
            "xtick.labelsize": 6.1,
            "ytick.labelsize": 6.1,
            "axes.linewidth": 0.7,
            "legend.fontsize": 6.1,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.06,
        1.02,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def clade_coordinates(tree: Phylo.BaseTree.Tree) -> tuple[dict, dict]:
    x = tree.depths()
    if not max(x.values()):
        x = tree.depths(unit_branch_lengths=True)
    terminals = tree.get_terminals()
    y = {terminal: index for index, terminal in enumerate(reversed(terminals))}

    def assign(clade):
        for child in clade.clades:
            assign(child)
        if clade not in y:
            y[clade] = float(np.mean([y[child] for child in clade.clades]))

    assign(tree.root)
    return x, y


def branch_plot(ax: plt.Axes, tree: Phylo.BaseTree.Tree) -> tuple[dict, dict]:
    x, y = clade_coordinates(tree)
    parent_by_child = {
        child: parent
        for parent in tree.find_clades(order="preorder")
        for child in parent.clades
    }
    for clade in tree.find_clades(order="preorder"):
        if clade.clades:
            child_y = [y[child] for child in clade.clades]
            ax.plot(
                [x[clade], x[clade]],
                [min(child_y), max(child_y)],
                color=COLORS["branch"],
                linewidth=0.45,
                zorder=1,
            )
        if clade in parent_by_child:
            parent = parent_by_child[clade]
            ax.plot(
                [x[parent], x[clade]],
                [y[clade], y[clade]],
                color=COLORS["branch"],
                linewidth=0.45,
                zorder=1,
            )
    return x, y


def mi_color(lineage: str) -> str:
    if lineage == "MI_TMI_lineage":
        return COLORS["tmi"]
    if lineage == "MI_MP_MIP_lineage":
        return COLORS["mp"]
    if "chimaera" in lineage:
        return COLORS["chimaera"]
    if "yongonense" in lineage:
        return COLORS["yongonense"]
    if "marseillense" in lineage:
        return COLORS["marseillense"]
    return COLORS["other"]


def context_color(panel: str, row: pd.Series) -> str:
    if row.source == "anchor":
        return COLORS["anchor"]
    label = row.reporting_label.lower()
    if panel == "avium":
        if "timonense" in label:
            return COLORS["timonense"]
        if "avium" in label:
            return COLORS["avium"]
    else:
        if "colombiense" in label:
            return COLORS["colombiense"]
    return COLORS["other"]


def draw_context_tree(
    ax: plt.Axes,
    tree_path: Path,
    panel: pd.DataFrame,
    title: str,
    subtitle: str,
    panel_kind: str,
    label_all_local: bool = True,
    exclude_lineages: set[str] | None = None,
    exclude_tree_ids: set[str] | None = None,
) -> pd.DataFrame:
    tree = Phylo.read(tree_path, "newick")
    if exclude_tree_ids:
        for tip in list(tree.get_terminals()):
            if tip.name in exclude_tree_ids:
                tree.prune(tip)
        panel = panel[~panel.tree_id.isin(exclude_tree_ids)].copy()
    if exclude_lineages:
        panel_index_for_pruning = panel.set_index("tree_id")
        for tip in list(tree.get_terminals()):
            lineage = panel_index_for_pruning.loc[tip.name].get(
                "accessory_lineage", ""
            )
            if lineage in exclude_lineages:
                tree.prune(tip)
        panel = panel[~panel.accessory_lineage.isin(exclude_lineages)].copy()
    anchor_rows = panel[panel.source.eq("anchor")]
    if not anchor_rows.empty:
        anchor_name = anchor_rows.tree_id.iat[0]
        anchor = next((tip for tip in tree.get_terminals() if tip.name == anchor_name), None)
        if anchor is not None:
            tree.root_with_outgroup(anchor)
    tree.ladderize(reverse=True)
    x, y = branch_plot(ax, tree)
    panel_index = panel.set_index("tree_id")
    max_x = max(x.values())
    source_rows = []
    local_positions = []

    for tip in tree.get_terminals():
        row = panel_index.loc[tip.name]
        if panel_kind == "mi":
            color = mi_color(row.accessory_lineage)
            lineage = row.accessory_lineage
        else:
            color = context_color(panel_kind, row)
            lineage = row.analysis_lineage
        marker = "*" if row.source == "local" else "D" if row.source == "anchor" else "o"
        size = 20 if row.source == "local" else 10 if row.source == "anchor" else 6
        ax.scatter(
            x[tip],
            y[tip],
            marker=marker,
            s=size,
            facecolor=color,
            edgecolor="white" if row.source == "local" else "none",
            linewidth=0.35,
            zorder=3,
            clip_on=False,
        )
        if row.source == "local":
            local_positions.append((x[tip], y[tip], row.sample_id, color, row.cohort_source))
        source_rows.append(
            {
                "tree_id": tip.name,
                "source": row.source,
                "sample_id": row.sample_id,
                "reporting_label": row.reporting_label,
                "fine_label": row.fine_label,
                "analysis_lineage": lineage,
                "tree_x": x[tip],
                "tree_y": y[tip],
            }
        )

    if label_all_local:
        offset = max(max_x * 0.012, 0.001)
        for tip_x, tip_y, sample_id, color, cohort_source in local_positions:
            label = sample_id
            if cohort_source == "rescued_from_mixed_or_failed_original_assembly":
                label = f"{sample_id} (R)"
            ax.text(
                tip_x + offset,
                tip_y,
                label,
                fontsize=5.2 if panel_kind == "mi" else 5.8,
                fontweight="bold",
                color=color,
                ha="left",
                va="center",
                clip_on=False,
            )

    ax.set_ylim(-1, len(tree.get_terminals()))
    ax.set_xlim(-max_x * 0.01, max_x * (1.12 if panel_kind == "mi" else 1.18))
    ax.set_yticks([])
    ax.set_xlabel("Substitutions per complete-core SNP site")
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.tick_params(axis="x", length=2.5)
    ax.set_title(title, loc="left", pad=8)
    ax.text(
        0,
        1.005,
        subtitle,
        transform=ax.transAxes,
        fontsize=5.8,
        color=COLORS["muted"],
        ha="left",
        va="bottom",
    )
    return pd.DataFrame(source_rows)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mi-tree", required=True)
    parser.add_argument("--mi-panel", required=True)
    parser.add_argument("--local-labels", required=True)
    parser.add_argument("--avium-tree", required=True)
    parser.add_argument("--avium-panel", required=True)
    parser.add_argument("--colombiense-tree", required=True)
    parser.add_argument("--colombiense-panel", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    configure()
    output_dir = Path(args.output_dir)
    source_dir = output_dir / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)

    local_labels = pd.read_csv(args.local_labels, sep="\t", dtype=str).fillna("")
    cohort_source = local_labels.set_index("sample_id").cohort_source.to_dict()

    def add_cohort_source(panel: pd.DataFrame) -> pd.DataFrame:
        panel = panel.copy()
        if "cohort_source" not in panel:
            panel["cohort_source"] = panel.sample_id.map(cohort_source).fillna(panel.source)
        return panel

    mi_panel = add_cohort_source(pd.read_csv(args.mi_panel, sep="\t", dtype=str).fillna(""))
    fig, ax = plt.subplots(figsize=(183 * MM, 220 * MM))
    mi_source = draw_context_tree(
        ax,
        Path(args.mi_tree),
        mi_panel,
        "Current public context of the local M. intracellulare complex",
        "90 genomes displayed from the 94-genome tree; four distant M. marseillense records omitted for readability",
        "mi",
        exclude_lineages={"MI_marseillense_lineage"},
    )
    panel_label(ax, "a")
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["tmi"], markeredgecolor="none", label="TMI"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["mp"], markeredgecolor="none", label="MP-MIP"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["chimaera"], markeredgecolor="none", label="Chimaera-adjacent"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["yongonense"], markeredgecolor="none", label="Yongonense boundary"),
        mpl.lines.Line2D([], [], marker="*", linestyle="", markerfacecolor="#555555", markeredgecolor="none", label="Local genome; R, recovered"),
    ]
    ax.legend(handles=handles, loc="upper left", ncol=2, handletextpad=0.3, columnspacing=0.8)
    save_figure(fig, output_dir / "Supplementary_Figure_S1_complete_MI_context")
    mi_source.to_csv(source_dir / "Supplementary_Figure_S1_source.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 150 * MM), gridspec_kw={"wspace": 0.28})
    avium_panel = add_cohort_source(pd.read_csv(args.avium_panel, sep="\t", dtype=str).fillna(""))
    avium_source = draw_context_tree(
        axes[0],
        Path(args.avium_tree),
        avium_panel,
        "M. avium / M. timonense public context",
        "63/64 genomes shown; 8,921 sites; one distant public tip omitted",
        "avium",
        exclude_tree_ids={"PUB_GCF_010723675.1"},
    )
    panel_label(axes[0], "a")
    colombiense_panel = add_cohort_source(pd.read_csv(args.colombiense_panel, sep="\t", dtype=str).fillna(""))
    colombiense_source = draw_context_tree(
        axes[1],
        Path(args.colombiense_tree),
        colombiense_panel,
        "M. colombiense public context",
        "29 genomes; 851 complete-core variable sites",
        "colombiense",
    )
    panel_label(axes[1], "b")
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["avium"], markeredgecolor="none", label="M. avium label"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["timonense"], markeredgecolor="none", label="M. timonense label"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["colombiense"], markeredgecolor="none", label="M. colombiense label"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["other"], markeredgecolor="none", label="Other public label"),
        mpl.lines.Line2D([], [], marker="*", linestyle="", markerfacecolor="#555555", markeredgecolor="none", label="Local genome; R, recovered"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.01), frameon=False)
    fig.subplots_adjust(bottom=0.12, top=0.94, left=0.055, right=0.985)
    save_figure(fig, output_dir / "Supplementary_Figure_S2_non_MI_contexts")
    pd.concat(
        [avium_source.assign(panel="avium_timonense"), colombiense_source.assign(panel="colombiense")],
        ignore_index=True,
    ).to_csv(source_dir / "Supplementary_Figure_S2_source.tsv", sep="\t", index=False)

    print(output_dir)


if __name__ == "__main__":
    main()
