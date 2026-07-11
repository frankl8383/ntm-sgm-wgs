#!/usr/bin/env python3
"""Plot the reciprocal near-MAC dilution benchmark and route confirmation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PRIMARY_THRESHOLD = 161.8
AUXILIARY_THRESHOLD = 952.42
FRACTIONS = [0, 5, 10, 20, 30, 50]
PAIR_ORDER = ["pairA", "pairB", "pairC", "pairD"]
PAIR_SAMPLES = {
    "pairA": ("Mi1", "Ma20"),
    "pairB": ("Mi2", "Mi22"),
    "pairC": ("Mi23", "Mi8"),
    "pairD": ("Mi25", "Mi4"),
}
PAIR_LABELS = {
    "pairA": "A  Mi1 / Ma20",
    "pairB": "B  Mi2 / Mi22",
    "pairC": "C  Mi23 / Mi8",
    "pairD": "D  Mi25 / Mi4",
}
PAIR_COLORS = {
    "pairA": "#3B6FB6",
    "pairB": "#29926D",
    "pairC": "#D79A24",
    "pairD": "#A95A8B",
}
INK = "#292929"
MID = "#6D6D6D"
GRID = "#D9D9D9"
PALE = "#F2F2F2"


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
            "xtick.labelsize": 6.2,
            "ytick.labelsize": 6.2,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        fontweight="bold",
        color=INK,
        clip_on=False,
    )


def sample_from_path(value: str) -> str:
    return Path(value).name.split(".")[0]


def load_pair_ani(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["query", "reference", "ani", "matched", "total"],
    )
    raw["query_sample"] = raw["query"].map(sample_from_path)
    raw["reference_sample"] = raw["reference"].map(sample_from_path)
    rows: list[dict[str, object]] = []
    for pair_id in PAIR_ORDER:
        left, right = PAIR_SAMPLES[pair_id]
        selected = raw.loc[
            raw["query_sample"].isin([left, right])
            & raw["reference_sample"].isin([left, right])
            & raw["query_sample"].ne(raw["reference_sample"])
        ]
        if len(selected) != 2:
            raise ValueError(f"Expected two reciprocal ANI rows for {pair_id}, found {len(selected)}")
        rows.append(
            {
                "pair_id": pair_id,
                "pair_label": PAIR_LABELS[pair_id],
                "design_role": "calibration" if pair_id in {"pairA", "pairB"} else "untouched validation",
                "left_sample": left,
                "right_sample": right,
                "ani_mean": selected["ani"].mean(),
                "ani_min": selected["ani"].min(),
                "ani_max": selected["ani"].max(),
            }
        )
    return pd.DataFrame(rows)


def load_inputs(project: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    root = project / "results" / "near_mac_dilution"
    aggregate = pd.read_csv(
        root / "stage1_clean_reference" / "summary_analysis" / "nearmac_dilution_aggregate.tsv",
        sep="\t",
    )
    stage2_path = root / "stage2_two_route" / "selected_two_route_summary.tsv"
    if not stage2_path.exists():
        raise FileNotFoundError(
            f"Stage 2 summary is not yet available: {stage2_path}. "
            "Complete 72_run_nearmac_selected_two_route.py first."
        )
    stage2 = pd.read_csv(stage2_path, sep="\t")
    pair_ani = load_pair_ani(root / "source_pair_all_vs_all_fastani.tsv")
    return aggregate, stage2, pair_ani


def direction_rank(frame: pd.DataFrame) -> pd.Series:
    ranks: list[int] = []
    for row in frame.itertuples(index=False):
        left, right = PAIR_SAMPLES[row.pair_id]
        ranks.append(0 if (row.major, row.minor) == (left, right) else 1)
    return pd.Series(ranks, index=frame.index)


def plot_pair_ani(ax: plt.Axes, data: pd.DataFrame) -> None:
    y = np.arange(len(data))
    for index, row in data.iterrows():
        color = PAIR_COLORS[row["pair_id"]]
        ax.plot(
            [row["ani_min"], row["ani_max"]],
            [index, index],
            color=color,
            linewidth=1.4,
            solid_capstyle="round",
        )
        ax.scatter(row["ani_mean"], index, s=24, color=color, edgecolor="white", linewidth=0.6, zorder=3)
        ax.text(
            row["ani_max"] + 0.015,
            index,
            f"{row['ani_mean']:.2f}%",
            ha="left",
            va="center",
            fontsize=6,
            color=INK,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(data["pair_label"])
    ax.invert_yaxis()
    ax.set_xlim(98.82, 99.62)
    ax.set_xlabel("Reciprocal FastANI (%)")
    ax.set_title("Challenge-pair similarity", loc="left", fontweight="bold", pad=3)
    ax.grid(axis="x", color=GRID, linewidth=0.5)
    ax.tick_params(axis="y", length=0)
    ax.text(
        0.0,
        -0.24,
        "A/B calibration; C/D untouched validation",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.7,
        color=MID,
    )


def plot_dilution_curves(
    ax: plt.Axes,
    data: pd.DataFrame,
    mean_column: str,
    min_column: str,
    max_column: str,
    threshold: float,
    title: str,
) -> None:
    ranked = data.copy()
    ranked["direction_rank"] = direction_rank(ranked)
    for (pair_id, major, minor), group in ranked.groupby(["pair_id", "major", "minor"], sort=False):
        group = group.sort_values("minor_percent")
        rank = int(group["direction_rank"].iloc[0])
        color = PAIR_COLORS[pair_id]
        linestyle = "-" if rank == 0 else (0, (3, 1.5))
        ax.plot(
            group["minor_percent"],
            group[mean_column],
            color=color,
            linewidth=1.25,
            linestyle=linestyle,
            marker="o",
            markersize=2.7,
            markeredgewidth=0,
        )
        ax.fill_between(
            group["minor_percent"],
            group[min_column],
            group[max_column],
            color=color,
            alpha=0.10,
            linewidth=0,
        )
    ax.axhline(threshold, color=INK, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(
        49.2,
        threshold * 1.08,
        f"frozen threshold {threshold:g}",
        ha="right",
        va="bottom",
        fontsize=5.5,
        color=INK,
    )
    ax.set_yscale("log")
    ax.set_xlim(-1, 51)
    ax.set_xticks(FRACTIONS)
    ax.set_ylim(8, 12000)
    ax.set_xlabel("Minor-source reads (%)")
    ax.set_ylabel("Intermediate sites per callable Mb")
    ax.set_title(title, loc="left", fontweight="bold", pad=3)
    ax.grid(axis="y", which="major", color=GRID, linewidth=0.45)


def detection_rows(data: pd.DataFrame) -> pd.DataFrame:
    ranked = data.copy()
    ranked["direction_rank"] = direction_rank(ranked)
    ranked["direction_label"] = ranked.apply(
        lambda row: f"{row['pair_id'][-1]}  {row['major']} + {row['minor']}", axis=1
    )
    return ranked.sort_values(["pair_id", "direction_rank", "minor_percent"])


def plot_detection_heatmap(ax: plt.Axes, data: pd.DataFrame) -> None:
    ordered = detection_rows(data)
    labels = ordered[["pair_id", "direction_rank", "direction_label"]].drop_duplicates()["direction_label"].tolist()
    matrix = np.full((len(labels), len(FRACTIONS)), np.nan)
    for row_index, label in enumerate(labels):
        subset = ordered.loc[ordered["direction_label"].eq(label)].set_index("minor_percent")
        for column_index, fraction in enumerate(FRACTIONS):
            matrix[row_index, column_index] = float(subset.loc[fraction, "combined_detection_rate"])
    cmap = mpl.colors.ListedColormap(["#EEEEEE", "#3E7E6B"])
    ax.imshow(matrix, aspect="auto", interpolation="none", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(FRACTIONS)))
    ax.set_xticklabels(FRACTIONS)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Minor-source reads (%)")
    ax.set_title("Combined detection across three read windows", loc="left", fontweight="bold", pad=3)
    ax.tick_params(length=0)
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix[row_index, column_index]
            ax.text(
                column_index,
                row_index,
                f"{int(round(value * 3))}/3",
                ha="center",
                va="center",
                fontsize=5.4,
                color="white" if value >= 0.67 else MID,
                fontweight="bold" if value >= 0.67 else "normal",
            )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.axhline(3.5, color="white", linewidth=2.0)
    ax.text(
        0.0,
        -0.18,
        "Cells show positive read windows/3; thresholds frozen before C/D validation.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.7,
        color=MID,
    )


def stage2_label(row: pd.Series) -> str:
    suffix = "cal" if row["pair_id"] in {"pairA", "pairB"} else "val"
    return f"{row['pair_id'][-1]} {row['major']}\n{int(row['minor_percent'])}% {suffix}"


def plot_stage2(ax: plt.Axes, data: pd.DataFrame) -> None:
    ordered = data.copy()
    ordered["pair_rank"] = ordered["pair_id"].map({item: index for index, item in enumerate(PAIR_ORDER)})
    ordered["direction_rank"] = direction_rank(ordered)
    ordered = ordered.sort_values(["pair_rank", "minor_percent", "direction_rank"]).reset_index(drop=True)
    x = np.arange(len(ordered))
    for route, offset, marker, label in (
        ("strict", -0.13, "o", "Strict route"),
        ("meta", 0.13, "s", "Meta route"),
    ):
        values = ordered[f"{route}_burden_20_80_per_mbp"].astype(float)
        colors = [PAIR_COLORS[item] for item in ordered["pair_id"]]
        ax.scatter(
            x + offset,
            values,
            s=20,
            marker=marker,
            facecolors=colors if route == "strict" else "white",
            edgecolors=colors,
            linewidths=0.8,
            label=label,
            zorder=3,
        )
    ax.axhline(PRIMARY_THRESHOLD, color=INK, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(
        len(ordered) - 0.6,
        PRIMARY_THRESHOLD * 1.08,
        "161.8",
        ha="right",
        va="bottom",
        fontsize=5.5,
        color=INK,
    )
    ax.set_yscale("log")
    ax.set_ylim(8, 12000)
    ax.set_xlim(-0.7, len(ordered) - 0.3)
    ax.set_xticks(x)
    ax.set_xticklabels([stage2_label(row) for _, row in ordered.iterrows()], rotation=45, ha="right")
    ax.set_ylabel("20-80% sites per callable Mb")
    ax.set_title("Complete two-route reconstructions", loc="left", fontweight="bold", pad=3)
    ax.grid(axis="y", which="major", color=GRID, linewidth=0.45)
    ax.legend(loc="upper left", fontsize=5.8, ncol=2, handletextpad=0.4, columnspacing=0.8)
    ax.text(
        1.0,
        -0.25,
        "Seed 202; A/B route-shape checks, C/D untouched operating-point checks",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.7,
        color=MID,
    )


def write_source_data(output: Path, aggregate: pd.DataFrame, stage2: pd.DataFrame, pair_ani: pd.DataFrame) -> None:
    source = output / "source_data"
    source.mkdir(parents=True, exist_ok=True)
    aggregate.to_csv(source / "Figure4_stage1_dilution_aggregate.tsv", sep="\t", index=False)
    stage2.to_csv(source / "Figure4_stage2_selected_two_route.tsv", sep="\t", index=False)
    pair_ani.to_csv(source / "Figure4_source_pair_ani.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    aggregate, stage2, pair_ani = load_inputs(project)
    output = project / "results" / "figures"
    output.mkdir(parents=True, exist_ok=True)
    write_source_data(output, aggregate, stage2, pair_ani)

    set_style()
    width_inches = 183 / 25.4
    height_inches = 170 / 25.4
    figure = plt.figure(figsize=(width_inches, height_inches), facecolor="white")
    grid = figure.add_gridspec(
        2,
        6,
        height_ratios=[1.0, 1.12],
        width_ratios=[0.95, 0.95, 1.25, 1.25, 1.25, 1.25],
        left=0.125,
        right=0.985,
        top=0.955,
        bottom=0.15,
        wspace=1.05,
        hspace=0.68,
    )
    ax_a = figure.add_subplot(grid[0, 0:2])
    ax_b = figure.add_subplot(grid[0, 2:4])
    ax_c = figure.add_subplot(grid[0, 4:6])
    ax_d = figure.add_subplot(grid[1, 0:3])
    ax_e = figure.add_subplot(grid[1, 3:6])

    plot_pair_ani(ax_a, pair_ani)
    plot_dilution_curves(
        ax_b,
        aggregate,
        "burden_20_80_mean",
        "burden_20_80_min",
        "burden_20_80_max",
        PRIMARY_THRESHOLD,
        "Primary 20-80% burden",
    )
    plot_dilution_curves(
        ax_c,
        aggregate,
        "burden_10_90_mean",
        "burden_10_90_min",
        "burden_10_90_max",
        AUXILIARY_THRESHOLD,
        "Auxiliary 10-90% burden",
    )
    plot_detection_heatmap(ax_d, aggregate)
    plot_stage2(ax_e, stage2)

    for axis, label in zip([ax_a, ax_b, ax_c, ax_d, ax_e], "abcde"):
        panel_label(axis, label)

    handles = [
        mpl.lines.Line2D([], [], color=PAIR_COLORS[pair], marker="o", linewidth=1.2, markersize=3, label=PAIR_LABELS[pair])
        for pair in PAIR_ORDER
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.018),
        ncol=4,
        fontsize=6,
        handlelength=1.6,
        columnspacing=1.2,
    )

    stem = output / "Figure4_nearMAC_dilution_benchmark"
    figure.savefig(f"{stem}.svg", facecolor="white")
    figure.savefig(f"{stem}.pdf", facecolor="white")
    figure.savefig(f"{stem}.png", dpi=300, facecolor="white")
    figure.savefig(f"{stem}.tiff", dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(figure)
    print(stem)


if __name__ == "__main__":
    main()
