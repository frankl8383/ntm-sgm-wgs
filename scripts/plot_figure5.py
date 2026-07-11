#!/usr/bin/env python3
"""Plot structural and sampling-sensitivity checks for the four frozen blocks."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


BLOCK_ORDER = [
    "TMI_aromatic_catabolism_associated_block",
    "TMI_methyltransferase_hydrolase_cupin_block",
    "MP_MIP_nitrogen_redox_associated_block",
    "MP_MIP_oxidoreductase_associated_block",
]

BLOCK_LABELS = {
    "TMI_aromatic_catabolism_associated_block": "TMI aromatic\n26 families",
    "TMI_methyltransferase_hydrolase_cupin_block": "TMI methyl/hydrolase\n5 families",
    "MP_MIP_nitrogen_redox_associated_block": "MP-MIP N/redox\n3 families",
    "MP_MIP_oxidoreductase_associated_block": "MP-MIP oxidoreductase\n5 families",
}

TMI = "#4C78A8"
MP = "#C65D57"
INK = "#2F2F2F"
MID = "#777777"
LIGHT = "#E6E6E6"
PALE = "#F4F4F4"
ACCENT = "#D7A928"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    return parser.parse_args()


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.titlesize": 7.5,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.3,
            "ytick.labelsize": 6.3,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        fontweight="bold",
        color=INK,
    )


def load_prevalence(project: Path) -> pd.DataFrame:
    base = project / "results"
    non_discovery = pd.read_csv(
        base
        / "24_public_non_discovery_evaluation"
        / "non_discovery_public_block_statistics.tsv",
        sep="\t",
    )
    complete = pd.read_csv(
        base / "22_complete_genome_block_validation" / "complete_chromosome_block_prevalence.tsv",
        sep="\t",
    )
    project_one = pd.read_csv(
        base / "23_clone_project_sensitivity" / "one_per_bioproject_lineage_block_statistics.tsv",
        sep="\t",
    )
    clone = pd.read_csv(
        base / "23_clone_project_sensitivity" / "strict_d20_kmer002_block_statistics.tsv",
        sep="\t",
    )

    rows: list[dict[str, object]] = []
    specifications = [
        (
            "Non-discovery",
            non_discovery,
            "tmi_present",
            "tmi_total",
            "mp_mip_present",
            "mp_mip_total",
        ),
        (
            "Complete",
            complete,
            "complete_chromosome_tmi_present",
            "complete_chromosome_tmi_total",
            "complete_chromosome_mp_mip_present",
            "complete_chromosome_mp_mip_total",
        ),
        (
            "One/project",
            project_one,
            "tmi_present",
            "tmi_total",
            "mp_mip_present",
            "mp_mip_total",
        ),
        (
            "Near-clone",
            clone,
            "tmi_present",
            "tmi_total",
            "mp_mip_present",
            "mp_mip_total",
        ),
    ]
    for dataset, frame, tp, tt, mp, mt in specifications:
        indexed = frame.set_index("block_id")
        for block in BLOCK_ORDER:
            row = indexed.loc[block]
            for lineage, present_key, total_key in (
                ("TMI", tp, tt),
                ("MP-MIP", mp, mt),
            ):
                present = int(row[present_key])
                total = int(row[total_key])
                rows.append(
                    {
                        "dataset": dataset,
                        "block_id": block,
                        "lineage": lineage,
                        "present": present,
                        "total": total,
                        "prevalence": present / total,
                    }
                )
    result = pd.DataFrame(rows)
    expected = len(BLOCK_ORDER) * 4 * 2
    if len(result) != expected:
        raise ValueError(f"Expected {expected} prevalence rows, found {len(result)}")
    return result


def plot_prevalence_matrix(ax: plt.Axes, data: pd.DataFrame) -> None:
    datasets = ["Non-discovery", "Complete", "One/project", "Near-clone"]
    lineages = ["TMI", "MP-MIP"]
    ax.set_xlim(-0.5, 7.5)
    ax.set_ylim(3.5, -1.20)
    ax.set_xticks([])
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels([BLOCK_LABELS[item] for item in BLOCK_ORDER], ha="right")
    ax.tick_params(axis="y", length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for dataset_index, dataset in enumerate(datasets):
        x0 = dataset_index * 2
        ax.text(
            x0 + 0.5,
            -1.00,
            dataset,
            ha="center",
            va="bottom",
            fontsize=6.5,
            fontweight="bold",
            color=INK,
            clip_on=False,
        )
        for lineage_index, lineage in enumerate(lineages):
            x = x0 + lineage_index
            ax.text(
                x,
                -0.60,
                "T" if lineage == "TMI" else "M",
                ha="center",
                va="bottom",
                fontsize=6,
                color=TMI if lineage == "TMI" else MP,
                fontweight="bold",
                clip_on=False,
            )
            for block_index, block in enumerate(BLOCK_ORDER):
                row = data.loc[
                    data["dataset"].eq(dataset)
                    & data["lineage"].eq(lineage)
                    & data["block_id"].eq(block)
                ].iloc[0]
                prevalence = float(row["prevalence"])
                color = TMI if lineage == "TMI" else MP
                face = mpl.colors.to_rgba(color, 0.10 + 0.90 * prevalence)
                rect = Rectangle(
                    (x - 0.43, block_index - 0.40),
                    0.86,
                    0.80,
                    facecolor=face,
                    edgecolor="white",
                    linewidth=0.8,
                )
                ax.add_patch(rect)
                text_color = "white" if prevalence >= 0.70 else INK
                ax.text(
                    x,
                    block_index,
                    f"{int(row['present'])}/{int(row['total'])}",
                    ha="center",
                    va="center",
                    fontsize=5.7,
                    color=text_color,
                    fontweight="bold" if prevalence >= 0.70 else "normal",
                )
        if dataset_index < len(datasets) - 1:
            ax.axvline(x0 + 1.5, color="#D6D6D6", linewidth=0.6)
    ax.set_title("Frozen public evaluation and sensitivity checks", loc="left", pad=3, fontweight="bold")
    ax.text(
        1.0,
        -0.17,
        "T, TMI; M, MP-MIP. Cell labels show present/total.",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.7,
        color=MID,
    )


def draw_locus_row(
    ax: plt.Axes,
    y: float,
    label: str,
    block_color: str,
    block_present: bool,
    count_label: str,
    partial: bool = False,
) -> None:
    ax.text(0.01, y, label, ha="left", va="center", fontsize=6.4, color=INK)
    left_x, left_w = 0.28, 0.17
    block_x, block_w = 0.47, 0.23
    right_x, right_w = 0.72, 0.17
    for x, width in ((left_x, left_w), (right_x, right_w)):
        ax.add_patch(
            FancyBboxPatch(
                (x, y - 0.055),
                width,
                0.11,
                boxstyle="round,pad=0.004,rounding_size=0.008",
                linewidth=0.6,
                edgecolor=MID,
                facecolor="#CFCFCF",
            )
        )
    if block_present:
        ax.add_patch(
            FancyBboxPatch(
                (block_x, y - 0.075),
                block_w,
                0.15,
                boxstyle="round,pad=0.004,rounding_size=0.008",
                linewidth=0.7,
                edgecolor=block_color,
                facecolor=mpl.colors.to_rgba(block_color, 0.82),
            )
        )
    else:
        ax.add_patch(
            Rectangle(
                (block_x, y - 0.06),
                block_w,
                0.12,
                linewidth=0.8,
                edgecolor=MID,
                facecolor=PALE,
                linestyle=(0, (2, 1.5)),
            )
        )
        ax.text(block_x + block_w / 2, y, "absent" if not partial else "partial", ha="center", va="center", fontsize=5.5, color=MID)
    ax.text(0.92, y, count_label, ha="right", va="center", fontsize=6.1, color=INK)


def plot_structural_context(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Conserved flanks distinguish structural absence", loc="left", pad=3, fontweight="bold")
    ax.text(0.01, 0.86, "MP-MIP oxidoreductase block", fontsize=6.7, fontweight="bold", color=MP)
    draw_locus_row(ax, 0.72, "MP-MIP", MP, True, "10/10 full")
    draw_locus_row(ax, 0.54, "TMI", MP, False, "26/26 flanks")
    ax.text(0.01, 0.37, "TMI methyl/hydrolase block", fontsize=6.7, fontweight="bold", color=TMI)
    draw_locus_row(ax, 0.23, "TMI", TMI, True, "22/26 full")
    draw_locus_row(ax, 0.07, "MP-MIP", TMI, False, "8/10 flanks")
    ax.text(0.585, 0.95, "lineage-associated block", ha="center", va="center", fontsize=5.5, color=MID)
    ax.annotate("", xy=(0.585, 0.88), xytext=(0.585, 0.93), arrowprops={"arrowstyle": "-", "color": MID, "lw": 0.6})


def plot_mi31_exception(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Mi31 separates block content from one-contig synteny", loc="left", pad=3, fontweight="bold")
    labels = ["Strict assembly", "Selected meta assembly", "Original reads"]
    ys = [0.72, 0.47, 0.22]
    for label, y in zip(labels, ys):
        ax.text(0.01, y, label, ha="left", va="center", fontsize=6.3, color=INK)
    ax.add_patch(Rectangle((0.34, 0.66), 0.48, 0.12, facecolor=mpl.colors.to_rgba(TMI, 0.82), edgecolor=TMI, linewidth=0.7))
    ax.text(0.58, 0.72, "26/26 families, one contig", ha="center", va="center", fontsize=5.7, color="white", fontweight="bold")
    ax.add_patch(Rectangle((0.34, 0.41), 0.30, 0.12, facecolor=mpl.colors.to_rgba(TMI, 0.82), edgecolor=TMI, linewidth=0.7))
    ax.add_patch(Rectangle((0.68, 0.41), 0.14, 0.12, facecolor=mpl.colors.to_rgba(TMI, 0.82), edgecolor=TMI, linewidth=0.7))
    ax.plot([0.645, 0.675], [0.47, 0.47], color=ACCENT, linewidth=1.2, linestyle=(0, (2, 1.4)))
    ax.text(0.49, 0.47, "26/26 families", ha="center", va="center", fontsize=5.5, color="white", fontweight="bold")
    ax.text(0.86, 0.47, "split", ha="left", va="center", fontsize=5.7, color=ACCENT, fontweight="bold")
    ax.add_patch(Rectangle((0.34, 0.16), 0.48, 0.12, facecolor=mpl.colors.to_rgba(TMI, 0.25), edgecolor=TMI, linewidth=0.7))
    ax.text(0.58, 0.22, "100% breadth at depth >=5", ha="center", va="center", fontsize=5.7, color=INK)
    ax.text(0.99, 0.02, "Reads support content, not long-range gene order.", ha="right", va="bottom", fontsize=5.5, color=MID)


def plot_minimum_changes(ax: plt.Axes, transitions: pd.DataFrame) -> None:
    frame = transitions.loc[transitions["scenario"].eq("frozen_one_contig_synteny")].set_index("block_id")
    ax.set_xlim(-0.2, 2.6)
    ax.set_ylim(3.6, -0.6)
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels([BLOCK_LABELS[item] for item in BLOCK_ORDER])
    ax.set_xticks([0, 1, 2])
    ax.set_xlabel("Minimum number of changes")
    ax.grid(axis="x", color="#E5E5E5", linewidth=0.6, zorder=0)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for index, block in enumerate(BLOCK_ORDER):
        row = frame.loc[block]
        gains = int(row["inferred_gains"])
        losses = int(row["inferred_losses"])
        ax.barh(index - 0.12, gains, height=0.20, color="#6F9E73", edgecolor="none", label="gain" if index == 0 else None, zorder=2)
        ax.barh(index + 0.12, losses, height=0.20, color="#B78978", edgecolor="none", label="loss" if index == 0 else None, zorder=2)
        ax.text(gains + 0.08, index - 0.12, str(gains), va="center", ha="left", fontsize=5.7, color=INK)
        ax.text(losses + 0.08, index + 0.12, str(losses), va="center", ha="left", fontsize=5.7, color=INK)
    ax.legend(loc="lower right", ncol=2, handlelength=1.2, columnspacing=0.9, fontsize=5.8)
    ax.set_title("Minimum changes on the selected tree", loc="left", pad=3, fontweight="bold")
    ax.text(1.0, -0.20, "Descriptive state history; mechanism is not inferred.", transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=MID)


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    results = project / "results"
    outdir = results / "figures"
    source = outdir / "source_data"
    outdir.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)

    prevalence = load_prevalence(project)
    transitions = pd.read_csv(
        results / "24_block_ancestral_reconstruction" / "block_minimum_change_summary.tsv",
        sep="\t",
    )
    flanks = pd.read_csv(
        results / "25_block_flanking_context" / "block_flanking_context_summary.tsv",
        sep="\t",
    )
    route = pd.read_csv(
        results / "19_review_resolution" / "recovered_accessory_concordance" / "four_block_calls_by_route.tsv",
        sep="\t",
    )
    reads = pd.read_csv(
        results / "19_review_resolution" / "recovered_block_read_support" / "four_block_read_support.tsv",
        sep="\t",
    )
    mi31_block = "TMI_aromatic_catabolism_associated_block"
    mi31_route = route.loc[route["sample_id"].eq("Mi31") & route["block_id"].eq(mi31_block)].copy()
    mi31_reads = reads.loc[reads["sample_id"].eq("Mi31") & reads["block_id"].eq(mi31_block)].copy()
    if set(mi31_route["route"]) != {"strict", "meta"} or len(mi31_reads) != 1:
        raise ValueError("Mi31 route/read evidence is incomplete")

    prevalence.to_csv(source / "Figure5a_prevalence_sensitivity.tsv", sep="\t", index=False)
    flanks.to_csv(source / "Figure5b_flanking_context.tsv", sep="\t", index=False)
    pd.concat(
        [
            mi31_route.assign(evidence_type="assembly_route"),
            mi31_reads.assign(evidence_type="original_reads"),
        ],
        ignore_index=True,
        sort=False,
    ).to_csv(source / "Figure5c_Mi31_route_read_evidence.tsv", sep="\t", index=False)
    transitions.to_csv(source / "Figure5d_minimum_change_summary.tsv", sep="\t", index=False)

    set_style()
    width_in = 183 / 25.4
    height_in = 157 / 25.4
    fig = plt.figure(figsize=(width_in, height_in), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        left=0.16,
        right=0.985,
        bottom=0.10,
        top=0.965,
        width_ratios=[1.12, 1.0],
        height_ratios=[1.05, 0.90],
        wspace=0.36,
        hspace=0.40,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    plot_prevalence_matrix(ax_a, prevalence)
    plot_structural_context(ax_b)
    plot_mi31_exception(ax_c)
    plot_minimum_changes(ax_d, transitions)
    for ax, label in zip((ax_a, ax_b, ax_c, ax_d), "abcd"):
        panel_label(ax, label)

    stem = outdir / "Figure5_structural_block_validation"
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=600,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)
    print(f"figure\t{stem}")
    print(f"source_data\t{source}")


if __name__ == "__main__":
    main()
