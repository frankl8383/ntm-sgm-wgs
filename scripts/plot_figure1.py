#!/usr/bin/env python3
"""Plot upgraded Figure 1: cohort recovery and mixture safeguards."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle


MM = 1 / 25.4
COLORS = {
    "ink": "#242424",
    "muted": "#6A6A6A",
    "grid": "#D9D9D9",
    "neutral": "#E8E8E8",
    "direct": "#2F7F7A",
    "candidate": "#D6A43A",
    "recovered": "#2F7F7A",
    "failed": "#B95A55",
    "final": "#315F8C",
    "tmi": "#3E6FA3",
}


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.titlesize": 7.5,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "legend.fontsize": 6.5,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.055,
        1.035,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def draw_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    count: int,
    color: str,
    text_color: str = "#242424",
) -> None:
    ax.add_patch(
        Rectangle(
            (x, y),
            width,
            height,
            facecolor=color,
            edgecolor="none",
            linewidth=0,
        )
    )
    ax.text(
        x + width / 2,
        y + height * 0.62,
        title,
        ha="center",
        va="center",
        fontsize=6.7,
        color=text_color,
        linespacing=1.05,
    )
    ax.text(
        x + width / 2,
        y + height * 0.22,
        f"n = {count}",
        ha="center",
        va="center",
        fontsize=8,
        fontweight="bold",
        color=text_color,
    )


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=0.8,
            color=COLORS["muted"],
            shrinkA=1,
            shrinkB=1,
        )
    )


def plot_flow(ax: plt.Axes) -> pd.DataFrame:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "a")
    ax.text(0, 1.035, "Cohort reconstruction after targeted recovery", va="bottom")

    draw_box(ax, 0.00, 0.40, 0.14, 0.34, "Presumed\nSGM-NTM", 38, COLORS["neutral"])
    draw_box(ax, 0.24, 0.64, 0.15, 0.28, "Directly\ninterpretable", 13, "#C8E0DD")
    draw_box(ax, 0.24, 0.31, 0.15, 0.28, "Recovery\ncandidates", 11, "#F1DFB8")
    draw_box(ax, 0.24, -0.02, 0.15, 0.28, "Not advanced\nto recovery", 14, "#EEEEEE")
    draw_box(ax, 0.51, 0.46, 0.15, 0.28, "Recovered by\ntwo routes", 8, "#C8E0DD")
    draw_box(ax, 0.51, 0.10, 0.15, 0.28, "Residual\nwithin-MAC mixture", 3, "#E8C8C5")
    draw_box(ax, 0.82, 0.46, 0.18, 0.28, "Interpretable\nMAC/SGM genomes", 21, "#BED2E4")

    arrow(ax, (0.14, 0.57), (0.23, 0.78))
    arrow(ax, (0.14, 0.57), (0.23, 0.45))
    arrow(ax, (0.14, 0.57), (0.23, 0.12))
    arrow(ax, (0.39, 0.45), (0.50, 0.60))
    arrow(ax, (0.39, 0.45), (0.50, 0.24))
    arrow(ax, (0.39, 0.78), (0.81, 0.64))
    arrow(ax, (0.66, 0.60), (0.81, 0.60))
    ax.text(
        0.585,
        0.82,
        "13 direct + 8 recovered",
        ha="center",
        color=COLORS["muted"],
        fontsize=6.3,
    )
    ax.text(
        0.585,
        0.035,
        "Excluded before lineage and accessory analyses",
        ha="center",
        color=COLORS["muted"],
        fontsize=6.1,
    )
    return pd.DataFrame(
        [
            ("initial_presumed_sgm_ntm", 38),
            ("directly_interpretable", 13),
            ("recovery_candidates", 11),
            ("not_advanced_to_recovery", 14),
            ("recovered_two_routes", 8),
            ("residual_within_mac_mixture", 3),
            ("final_interpretable", 21),
        ],
        columns=["stage", "n"],
    )


def plot_benchmark(ax: plt.Axes, benchmark: pd.DataFrame) -> pd.DataFrame:
    panel_label(ax, "b")
    ax.set_title("Synthetic-mixture recruitment", loc="left", pad=7)
    strict = benchmark[
        benchmark.route.eq("strict_paired_bait_recruitment_spades")
    ].copy()
    strict["sensitivity_percent"] = strict.target_pair_sensitivity.astype(float) * 100
    strict["precision_percent"] = strict.recruited_pair_precision.astype(float) * 100
    cross = strict[strict.benchmark.str.startswith("cross")].copy()
    near = strict[strict.benchmark.eq("nearMAC50")].copy()
    ax.scatter(
        cross.sensitivity_percent,
        cross.precision_percent,
        s=28,
        facecolor=COLORS["direct"],
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
        label="Cross-genus mixture",
    )
    label_offsets = {25: (7, -14), 65: (7, -2), 95: (7, 10)}
    for row in cross.itertuples(index=False):
        target = int(round(float(row.target_fraction) * 100))
        ax.annotate(
            f"{target}% target",
            (row.sensitivity_percent, row.precision_percent),
            xytext=label_offsets[target],
            textcoords="offset points",
            fontsize=6.1,
            ha="left",
            arrowprops={"arrowstyle": "-", "color": "#9A9A9A", "lw": 0.45},
        )
    ax.scatter(
        near.sensitivity_percent,
        near.precision_percent,
        s=32,
        marker="D",
        facecolor=COLORS["failed"],
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
        label="50:50 within-MAC control",
    )
    if not near.empty:
        row = near.iloc[0]
        ax.annotate(
            "near-neighbour\ncontrol",
            (row.sensitivity_percent, row.precision_percent),
            xytext=(5, 0),
            textcoords="offset points",
            fontsize=6.1,
            ha="left",
            va="center",
        )
    ax.set_xlim(84.5, 92.5)
    ax.set_ylim(40, 103)
    ax.set_xlabel("Target-pair sensitivity (%)")
    ax.set_ylabel("Recruited-pair precision (%)")
    ax.set_yticks([40, 60, 80, 100])
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(
        0.98,
        0.95,
        "cross-genus mixtures",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.1,
        color=COLORS["direct"],
    )
    return strict[
        [
            "benchmark",
            "route",
            "target_fraction",
            "sensitivity_percent",
            "precision_percent",
            "background_pair_false_positive_rate",
            "source_ani",
            "source_fragment_recovery",
            "benchmark_interpretation",
        ]
    ]


def plot_mixture_burden(
    ax: plt.Axes, outcomes: pd.DataFrame, mixture_threshold: float
) -> pd.DataFrame:
    panel_label(ax, "c")
    ax.set_title("Residual within-MAC mixture after recovery", loc="left", pad=7)
    outcomes = outcomes.copy()
    outcomes["recovered"] = outcomes.final_rescue_decision.str.startswith(
        "rescued_interpretable"
    )
    outcomes = outcomes.sort_values(
        ["recovered", "meta_mixed_sites_per_mbp"], ascending=[False, True]
    ).reset_index(drop=True)
    x = np.arange(len(outcomes))
    strict_values = outcomes.strict_mixed_sites_per_mbp.astype(float).to_numpy()
    meta_values = outcomes.meta_mixed_sites_per_mbp.astype(float).to_numpy()
    for index, row in outcomes.iterrows():
        color = COLORS["recovered"] if row.recovered else COLORS["failed"]
        ax.plot(
            [index, index],
            [strict_values[index], meta_values[index]],
            color="#B8B8B8",
            linewidth=0.8,
            zorder=1,
        )
        ax.scatter(
            index,
            strict_values[index],
            s=19,
            marker="o",
            facecolor=color,
            edgecolor="white",
            linewidth=0.45,
            zorder=3,
        )
        ax.scatter(
            index,
            meta_values[index],
            s=22,
            marker="s",
            facecolor=color,
            edgecolor="white",
            linewidth=0.45,
            zorder=3,
        )
    threshold = float(mixture_threshold)
    outcomes["benchmark_derived_mixture_threshold"] = threshold
    ax.axhline(threshold, color=COLORS["candidate"], linestyle="--", linewidth=1)
    ax.axvline(7.5, color="#D0D0D0", linewidth=0.7)
    ax.text(
        len(outcomes) - 0.45,
        threshold * 1.12,
        f"fixed threshold {threshold:.1f}",
        color="#8A651E",
        fontsize=6.1,
        ha="right",
        va="bottom",
    )
    ax.set_yscale("log")
    ax.set_ylim(10, 6000)
    ax.set_ylabel("Intermediate-frequency sites per Mb")
    ax.set_xticks(x)
    ax.set_xticklabels(outcomes.sample_id, rotation=55, ha="right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", length=0)
    route_handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor="#777777", markeredgecolor="none", label="Strict route"),
        mpl.lines.Line2D([], [], marker="s", linestyle="", markerfacecolor="#777777", markeredgecolor="none", label="MetaSPAdes bin"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["recovered"], markeredgecolor="none", label="Recovered"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["failed"], markeredgecolor="none", label="Residual mixture"),
    ]
    ax.legend(
        handles=route_handles,
        loc="upper left",
        ncol=2,
        handletextpad=0.3,
        columnspacing=0.7,
        borderaxespad=0.35,
        labelspacing=0.35,
    )
    return outcomes[
        [
            "sample_id",
            "strict_mixed_sites_per_mbp",
            "meta_mixed_sites_per_mbp",
            "benchmark_derived_mixture_threshold",
            "reporting_name_or_anchor",
            "final_rescue_decision",
            "recovered",
        ]
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-summary", required=True)
    parser.add_argument("--rescue-outcomes", required=True)
    parser.add_argument("--mixture-threshold", type=float, default=161.826223)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    configure()
    output_dir = Path(args.output_dir)
    source_dir = output_dir / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    benchmark = pd.read_csv(args.benchmark_summary, sep="\t")
    outcomes = pd.read_csv(args.rescue_outcomes, sep="\t")

    fig = plt.figure(figsize=(183 * MM, 139 * MM))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.92, 1.28],
        width_ratios=[0.82, 1.18],
        left=0.065,
        right=0.985,
        top=0.965,
        bottom=0.12,
        hspace=0.52,
        wspace=0.34,
    )
    ax_flow = fig.add_subplot(grid[0, :])
    ax_benchmark = fig.add_subplot(grid[1, 0])
    ax_burden = fig.add_subplot(grid[1, 1])
    flow_source = plot_flow(ax_flow)
    benchmark_source = plot_benchmark(ax_benchmark, benchmark)
    burden_source = plot_mixture_burden(
        ax_burden, outcomes, args.mixture_threshold
    )

    flow_source.to_csv(source_dir / "Figure1a_cohort_flow.tsv", sep="\t", index=False)
    benchmark_source.to_csv(
        source_dir / "Figure1b_synthetic_mixture_benchmark.tsv", sep="\t", index=False
    )
    burden_source.to_csv(
        source_dir / "Figure1c_residual_mixture_burden.tsv", sep="\t", index=False
    )

    stem = output_dir / "Figure1_recovery_and_mixture_safeguards"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(stem)


if __name__ == "__main__":
    main()
