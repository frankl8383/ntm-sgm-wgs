#!/usr/bin/env python3
"""Plot Figure 2: current MAC atlas and local ANI context."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle


MM = 1 / 25.4
GROUP_ORDER = [
    "M_intracellulare_complex",
    "M_avium_timonense_boundary",
    "M_colombiense",
]
GROUP_LABELS = {
    "M_intracellulare_complex": "M. intracellulare complex",
    "M_avium_timonense_boundary": "M. avium public context",
    "M_colombiense": "M. colombiense",
}
COLORS = {
    "ink": "#242424",
    "muted": "#707070",
    "grid": "#D7D7D7",
    "mi": "#3E6FA3",
    "avium": "#D29A32",
    "colombiense": "#8D6A9F",
    "direct": "#555555",
    "rescued": "#2F7F7A",
    "neutral": "#E5E5E5",
}


def group_color(group: str) -> str:
    return {
        GROUP_ORDER[0]: COLORS["mi"],
        GROUP_ORDER[1]: COLORS["avium"],
        GROUP_ORDER[2]: COLORS["colombiense"],
    }.get(group, COLORS["muted"])


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.titlesize": 7.5,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.2,
            "ytick.labelsize": 6.2,
            "axes.linewidth": 0.7,
            "legend.fontsize": 6.2,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.08) -> None:
    ax.text(
        x,
        1.035,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def plot_atlas_flow(
    ax: plt.Axes, atlas: pd.DataFrame, raw_count: int, qc_pass_count: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    panel_label(ax, "a", x=-0.035)
    ax.set_title("Current NCBI MAC atlas used for public context", loc="left", pad=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        (0.00, "Current records", raw_count, "#ECECEC"),
        (0.20, "Assembly-QC\npass", qc_pass_count, "#DDE6EE"),
        (0.43, "QC-first paired/\nBioSample dedup.", len(atlas), "#C7D9E8"),
    ]
    for x, title, count, color in boxes:
        ax.add_patch(Rectangle((x, 0.30), 0.16, 0.47, facecolor=color, edgecolor="none"))
        ax.text(x + 0.08, 0.60, title, ha="center", va="center", fontsize=6.7)
        ax.text(x + 0.08, 0.40, f"n = {count}", ha="center", va="center", fontsize=8, fontweight="bold")
    for start, end in [((0.16, 0.54), (0.195, 0.54)), ((0.36, 0.54), (0.425, 0.54))]:
        ax.add_patch(
            FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=8, color=COLORS["muted"], linewidth=0.8)
        )

    composition = atlas.reporting_species.value_counts().rename_axis("reporting_species").reset_index(name="n")
    top = composition.iloc[:4].copy()
    other = int(composition.iloc[4:].n.sum())
    if other:
        top = pd.concat(
            [top, pd.DataFrame([{"reporting_species": "Other MAC labels", "n": other}])],
            ignore_index=True,
        )
    total = top.n.sum()
    left = 0.66
    width = 0.32
    species_colors = ["#4F7FAE", "#8FAECD", "#A87991", "#7AB0A1", "#C9C9C9"]
    cursor = left
    for row, color in zip(top.itertuples(index=False), species_colors, strict=False):
        segment = width * row.n / total
        ax.add_patch(Rectangle((cursor, 0.47), segment, 0.19, facecolor=color, edgecolor="white", linewidth=0.5))
        if segment > 0.035:
            ax.text(cursor + segment / 2, 0.565, str(row.n), ha="center", va="center", fontsize=6.2, color="white" if color != "#C9C9C9" else COLORS["ink"])
        cursor += segment
    ax.text(left, 0.75, f"Current reporting labels (n = {len(atlas)})", fontsize=6.7, ha="left")
    legend_y = 0.30
    for index, (row, color) in enumerate(zip(top.itertuples(index=False), species_colors, strict=False)):
        col = index % 2
        rr = index // 2
        x = left + col * 0.17
        y = legend_y - rr * 0.105
        ax.scatter(x, y, s=14, color=color)
        label = str(row.reporting_species).replace("Mycobacterium ", "M. ")
        ax.text(x + 0.012, y, f"{label} ({row.n})", va="center", ha="left", fontsize=5.4)

    flow = pd.DataFrame(
        [
            ("current_records", raw_count),
            ("assembly_qc_pass", qc_pass_count),
            ("qc_first_pair_and_biosample_deduplicated", len(atlas)),
        ],
        columns=["stage", "n"],
    )
    return flow, composition


def plot_local_composition(ax: plt.Axes, labels: pd.DataFrame) -> pd.DataFrame:
    panel_label(ax, "b")
    ax.set_title("Updated local cohort", loc="left", pad=8)
    table = (
        labels.groupby(["broad_analysis_panel", "cohort_source"]).size().unstack(fill_value=0).reindex(GROUP_ORDER).fillna(0)
    )
    direct_column = "directly_retained_original"
    rescued_column = "rescued_from_mixed_or_failed_original_assembly"
    direct = table.get(direct_column, pd.Series(0, index=table.index)).astype(int)
    rescued = table.get(rescued_column, pd.Series(0, index=table.index)).astype(int)
    y = np.arange(len(table))[::-1]
    ax.barh(y, direct, color="#A9A9A9", height=0.55, label="Direct")
    ax.barh(y, rescued, left=direct, color=COLORS["rescued"], height=0.55, label="Recovered")
    for index, group in enumerate(table.index):
        total = int(direct[group] + rescued[group])
        ax.text(total + 0.25, y[index], str(total), va="center", fontsize=6.3, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([GROUP_LABELS[group] for group in table.index])
    ax.set_xlim(0, 19)
    ax.set_xlabel("Local genomes")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower right", ncol=1, handlelength=1.2)
    source = table.reset_index()
    source["direct"] = source.pop(direct_column) if direct_column in source else 0
    source["recovered"] = source.pop(rescued_column) if rescued_column in source else 0
    return source


def plot_best_ani(ax: plt.Axes, merged: pd.DataFrame) -> pd.DataFrame:
    panel_label(ax, "c", x=-0.17)
    ax.set_title("Nearest public genome and type anchor", loc="left", pad=8)
    data = merged.copy()
    data["group_order"] = data.broad_analysis_panel.map({group: index for index, group in enumerate(GROUP_ORDER)})
    data = data.sort_values(["group_order", "provisional_genomic_lineage", "best_public_ani"], ascending=[True, True, False]).reset_index(drop=True)
    y = np.arange(len(data))[::-1]
    for index, row in data.iterrows():
        color = group_color(row.broad_analysis_panel)
        public = float(row.best_public_ani)
        type_ani = float(row.best_type_ani)
        ax.plot([public, type_ani], [y[index], y[index]], color="#C5C5C5", linewidth=0.75)
        ax.scatter(public, y[index], s=19, color=color, marker="o", zorder=3)
        ax.scatter(type_ani, y[index], s=20, facecolor="white", edgecolor=color, marker="D", linewidth=0.8, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(data.sample_id)
    for tick, group in zip(ax.get_yticklabels(), data.broad_analysis_panel, strict=True):
        tick.set_color(group_color(group))
        tick.set_fontweight("bold")
    ax.set_xlim(97.4, 100.1)
    ax.set_ylim(-0.8, len(data) + 2.2)
    ax.set_xticks([97.5, 98, 98.5, 99, 99.5, 100])
    ax.set_xlabel("ANI (%)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", color="#666666", label="Nearest current public genome"),
        mpl.lines.Line2D([], [], marker="D", linestyle="", markerfacecolor="white", markeredgecolor="#666666", label="Best direct type anchor"),
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        borderaxespad=0.35,
        handletextpad=0.3,
        labelspacing=0.35,
    )
    return data


def plot_margins(ax: plt.Axes, merged: pd.DataFrame) -> pd.DataFrame:
    panel_label(ax, "d", x=-0.16)
    ax.set_title("Separation from the nearest distinct label", loc="left", pad=8)
    data = merged.copy()
    x = data.public_distinct_species_ani_margin.astype(float)
    y = data.type_distinct_broad_group_ani_margin.astype(float)
    for group in GROUP_ORDER:
        subset = data[data.broad_analysis_panel.eq(group)]
        ax.scatter(
            subset.public_distinct_species_ani_margin.astype(float),
            subset.type_distinct_broad_group_ani_margin.astype(float),
            s=24,
            color=group_color(group),
            edgecolor="white",
            linewidth=0.4,
            label=GROUP_LABELS[group],
        )
    label_offsets = {
        "Ma49": (3, 6),
        "Mi8": (3, 5),
        "Mi22": (3, -9),
        "Mi23": (3, 5),
        "Ma9": (3, 5),
        "Mi2": (-18, -9),
        "Mix16B": (3, -9),
        "Mi24": (3, 3),
    }
    for row in data.itertuples(index=False):
        public_margin = float(row.public_distinct_species_ani_margin)
        type_margin = float(row.type_distinct_broad_group_ani_margin)
        if public_margin < 0.035 or type_margin < 0.12 or row.sample_id in {"Mi23", "Mi24", "Mix16B"}:
            ax.annotate(
                row.sample_id,
                (public_margin, type_margin),
                xytext=label_offsets.get(row.sample_id, (3, 3)),
                textcoords="offset points",
                fontsize=5.4,
            )
    ax.axvline(0.15, color="#A5A5A5", linestyle="--", linewidth=0.7)
    ax.axhline(0.15, color="#A5A5A5", linestyle="--", linewidth=0.7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.0025, 13)
    ax.set_ylim(0.05, 13)
    ax.set_xlabel("Public-label ANI margin (percentage points)")
    ax.set_ylabel("Type-group ANI margin (percentage points)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", fontsize=5.4, handletextpad=0.3, borderaxespad=0.35)
    ax.text(0.0032, 0.17, "narrow public-label margin", fontsize=5.1, color=COLORS["muted"], rotation=90, va="bottom")
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas-metadata", required=True)
    parser.add_argument("--local-labels", required=True)
    parser.add_argument("--ani-summary", required=True)
    parser.add_argument("--raw-count", type=int, default=1283)
    parser.add_argument("--qc-pass-count", type=int, default=957)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    configure()
    atlas = pd.read_csv(args.atlas_metadata, sep="\t", dtype=str).fillna("")
    labels = pd.read_csv(args.local_labels, sep="\t", dtype=str).fillna("")
    ani = pd.read_csv(args.ani_summary, sep="\t", dtype=str).fillna("")
    merged = labels.merge(ani, on=["sample_id", "cohort_source"], how="left", validate="one_to_one", suffixes=("", "_ani"))

    output_dir = Path(args.output_dir)
    source_dir = output_dir / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(183 * MM, 145 * MM))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[0.48, 1.0],
        left=0.065,
        right=0.985,
        top=0.96,
        bottom=0.10,
        hspace=0.36,
    )
    ax_flow = fig.add_subplot(outer[0])
    bottom = outer[1].subgridspec(1, 3, width_ratios=[0.85, 1.05, 1.18], wspace=0.52)
    ax_comp = fig.add_subplot(bottom[0, 0])
    ax_ani = fig.add_subplot(bottom[0, 1])
    ax_margin = fig.add_subplot(bottom[0, 2])

    flow_source, atlas_composition = plot_atlas_flow(
        ax_flow, atlas, args.raw_count, args.qc_pass_count
    )
    local_source = plot_local_composition(ax_comp, labels)
    ani_source = plot_best_ani(ax_ani, merged)
    margin_source = plot_margins(ax_margin, merged)

    flow_source.to_csv(source_dir / "Figure2a_atlas_selection_flow.tsv", sep="\t", index=False)
    atlas_composition.to_csv(source_dir / "Figure2a_atlas_reporting_label_composition.tsv", sep="\t", index=False)
    local_source.to_csv(source_dir / "Figure2b_local_group_composition.tsv", sep="\t", index=False)
    ani_source.to_csv(source_dir / "Figure2c_local_public_type_ani.tsv", sep="\t", index=False)
    margin_source.to_csv(source_dir / "Figure2d_local_ani_margins.tsv", sep="\t", index=False)

    stem = output_dir / "Figure2_current_MAC_atlas_and_local_ANI_context"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(stem)


if __name__ == "__main__":
    main()
