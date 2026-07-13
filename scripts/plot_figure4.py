#!/usr/bin/env python3
"""Plot the prespecified near-MAC benchmark after all 40 reconstructions complete."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PAIR_ORDER = ["pairA", "pairB", "pairC", "pairD"]
PAIR_SAMPLES = {
    "pairA": ("Mi1", "Ma20"),
    "pairB": ("Mi2", "Mi22"),
    "pairC": ("Mi23", "Mi8"),
    "pairD": ("Mi25", "Mi4"),
}
PAIR_COLORS = {
    "pairA": "#3F6E9A",
    "pairB": "#5B8C72",
    "pairC": "#C28A32",
    "pairD": "#9A6687",
}
PRIMARY_THRESHOLD = 161.8
AUXILIARY_THRESHOLD = 952.42
INK = "#252525"
MID = "#686868"
GRID = "#D8D8D8"
NEGATIVE = "#E8E8E8"
AUXILIARY = "#D9A441"
PRIMARY = "#3D8875"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    return parser.parse_args()


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.titlesize": 7.4,
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


def panel_label(ax: plt.Axes, label: str, x: float = -0.11, y: float = 1.04) -> None:
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


def read_inputs(
    project: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.Series]:
    analysis = project / "analysis_global_mac_upgrade"
    stage1 = pd.read_csv(
        analysis
        / "results/20_nearmac_dilution/stage1_clean_reference/summary_analysis/"
        "nearmac_dilution_aggregate.tsv",
        sep="\t",
    )
    complete = pd.read_csv(
        analysis
        / "results/27_postreview_validation/nearmac_expanded_two_route/"
        "expanded_two_route_summary_complete.tsv",
        sep="\t",
    )
    ani = pd.read_csv(
        analysis / "results/20_nearmac_dilution/source_pair_fastani_corrected.tsv",
        sep="\t",
    )
    dependence = pd.read_csv(
        analysis
        / "results/32_submission_freeze/technical_window_dependence/"
        "technical_window_overlap_summary.tsv",
        sep="\t",
    )
    if complete.shape[0] != 40:
        raise ValueError(f"Figure 4 requires 40 complete conditions; found {complete.shape[0]}")
    if dependence.shape[0] != 1:
        raise ValueError("Figure 4 requires one technical-window dependence summary row")
    return stage1, complete, ani, dependence.iloc[0]


def pair_ani_summary(ani: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for pair_id in PAIR_ORDER:
        left, right = PAIR_SAMPLES[pair_id]
        selected = ani.loc[
            ani["query_sample"].isin([left, right])
            & ani["reference_sample"].isin([left, right])
            & ani["query_sample"].ne(ani["reference_sample"])
        ]
        if selected.shape[0] != 2:
            raise ValueError(f"Expected two reciprocal FastANI rows for {pair_id}")
        rows.append(
            {
                "pair_id": pair_id,
                "samples": f"{left} / {right}",
                "ani_mean": selected["ani"].mean(),
                "ani_min": selected["ani"].min(),
                "ani_max": selected["ani"].max(),
            }
        )
    return pd.DataFrame(rows)


def direction_label(row: pd.Series) -> str:
    letter = row["pair_id"][-1]
    return f"{letter}  {row['major']} + {row['minor']}"


def plot_pair_ani(ax: plt.Axes, data: pd.DataFrame) -> None:
    positions = np.arange(data.shape[0])
    for position, row in data.iterrows():
        color = PAIR_COLORS[row["pair_id"]]
        ax.plot(
            [row["ani_min"], row["ani_max"]],
            [position, position],
            color=color,
            linewidth=1.5,
            solid_capstyle="round",
        )
        ax.scatter(
            row["ani_mean"],
            position,
            s=25,
            color=color,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
        ax.text(
            row["ani_max"] + 0.012,
            position,
            f"{row['ani_mean']:.2f}",
            va="center",
            fontsize=5.8,
            color=INK,
        )
    ax.set_yticks(positions)
    ax.set_yticklabels([f"{row.pair_id[-1]}  {row.samples}" for row in data.itertuples()])
    ax.invert_yaxis()
    ax.set_xlim(98.84, 99.58)
    ax.set_xlabel("Reciprocal FastANI (%)")
    ax.set_title("Challenge-pair similarity", loc="left", fontweight="bold", pad=3)
    ax.grid(axis="x", color=GRID, linewidth=0.45)
    ax.tick_params(axis="y", length=0)
    ax.text(
        0,
        -0.22,
        "A/B cross-lineage; C/D within MP-MIP",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.6,
        color=MID,
    )


def plot_stage1_detection(
    ax: plt.Axes, data: pd.DataFrame, dependence: pd.Series
) -> pd.DataFrame:
    summary = (
        data.groupby("minor_percent", as_index=False)
        .apply(
            lambda group: pd.Series(
                {
                    "technical_windows": int(group["technical_windows"].sum()),
                    "primary_positive": int(
                        round((group["primary_detection_rate"] * group["technical_windows"]).sum())
                    ),
                    "combined_positive": int(
                        round((group["combined_detection_rate"] * group["technical_windows"]).sum())
                    ),
                }
            ),
            include_groups=False,
        )
        .reset_index(drop=True)
    )
    for column, color, marker, label in (
        ("primary_positive", "#3F6E9A", "o", "20-80% rule"),
        ("combined_positive", PRIMARY, "s", "Combined rule"),
    ):
        ax.plot(
            summary["minor_percent"],
            summary[column],
            color=color,
            marker=marker,
            markersize=3.6,
            linewidth=1.35,
            label=label,
        )
    ax.set_xlim(-1, 51)
    ax.set_ylim(-1, 25.5)
    ax.set_xticks([0, 5, 10, 20, 30, 50])
    ax.set_yticks([0, 6, 12, 18, 24])
    ax.set_xlabel("Minor-source reads (%)")
    ax.set_ylabel("Positive dependent windows / 24")
    ax.set_title("Clean-reference dilution stage", loc="left", fontweight="bold", pad=3)
    ax.grid(axis="y", color=GRID, linewidth=0.45)
    ax.legend(loc="lower right", fontsize=5.7, handlelength=1.7)
    for row in summary.itertuples(index=False):
        if row.minor_percent in (5, 10, 20):
            ax.text(
                row.minor_percent,
                row.combined_positive + 0.8,
                f"{row.combined_positive}/24",
                ha="center",
                va="bottom",
                fontsize=5.3,
                color=PRIMARY,
            )
    mean_overlap = 100 * float(dependence["mean_pairwise_overlap_fraction"])
    maximum_overlap = 100 * float(dependence["maximum_pairwise_overlap_fraction"])
    ax.text(
        0,
        -0.23,
        f"Three deterministic windows per source; mean overlap {mean_overlap:.1f}%, max {maximum_overlap:.1f}%",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.3,
        color=MID,
    )
    return summary


def plot_complete_ten_percent(ax: plt.Axes, complete: pd.DataFrame) -> pd.DataFrame:
    selected = complete.loc[
        complete["minor_percent"].eq(10) & complete["total_pairs"].eq(2_000_000)
    ].copy()
    pair_rank = {pair: index for index, pair in enumerate(PAIR_ORDER)}
    selected["pair_rank"] = selected["pair_id"].map(pair_rank)
    selected["direction_rank"] = selected.apply(
        lambda row: 0 if row["major"] == PAIR_SAMPLES[row["pair_id"]][0] else 1,
        axis=1,
    )
    selected = selected.sort_values(["pair_rank", "direction_rank"]).reset_index(drop=True)
    if selected.shape[0] != 8:
        raise ValueError(f"Expected eight complete 10% conditions; found {selected.shape[0]}")
    ratios = np.column_stack(
        [
            selected["strict_burden_20_80_per_mbp"] / PRIMARY_THRESHOLD,
            selected["meta_burden_20_80_per_mbp"] / PRIMARY_THRESHOLD,
            selected["strict_burden_10_90_per_mbp"] / AUXILIARY_THRESHOLD,
            selected["meta_burden_10_90_per_mbp"] / AUXILIARY_THRESHOLD,
        ]
    )
    values = np.log2(np.clip(ratios, 0.125, 8.0))
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "threshold_ratio", ["#496E88", "#F7F7F7", "#B96346"]
    )
    norm = mpl.colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    ax.imshow(values, aspect="auto", cmap=cmap, norm=norm, interpolation="none")
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(
        ["Strict\nprimary", "Meta\nprimary", "Strict\nauxiliary", "Meta\nauxiliary"]
    )
    ax.set_yticks(np.arange(selected.shape[0]))
    ax.set_yticklabels([direction_label(row) for _, row in selected.iterrows()])
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0, pad=2, labelsize=5.7)
    ax.set_title("Complete 10% reconstruction", loc="left", fontweight="bold", pad=3)
    for row_index in range(ratios.shape[0]):
        for column_index in range(ratios.shape[1]):
            value = ratios[row_index, column_index]
            color = "white" if abs(values[row_index, column_index]) > 1.5 else INK
            ax.text(
                column_index,
                row_index,
                f"{value:.1f}x",
                ha="center",
                va="center",
                fontsize=5.3,
                color=color,
            )
    for boundary in (1.5,):
        ax.axvline(boundary, color="white", linewidth=1.7)
    for boundary in (1.5, 3.5, 5.5):
        ax.axhline(boundary, color="white", linewidth=1.5)
    failed_rows = [
        index
        for index, row in selected.iterrows()
        if not truth(row["prespecified_rule_satisfied"])
    ]
    for row_index in failed_rows:
        ax.add_patch(
            mpl.patches.Rectangle(
                (-0.5, row_index - 0.5),
                4,
                1,
                fill=False,
                edgecolor="#B33A3A",
                linewidth=1.15,
                clip_on=False,
            )
        )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        1,
        -0.16,
        "Cell values are burden / prespecified threshold; values >1 are positive",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.6,
        color=MID,
    )
    if failed_rows:
        ax.text(
            0,
            -0.16,
            "Red outline: prespecified expectation not met",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=5.6,
            color="#B33A3A",
        )
    return selected


def truth(value: object) -> bool:
    return value is True or str(value).strip().lower() == "true"


def classify_condition(row: pd.Series) -> int:
    if truth(row["either_route_primary_positive"]):
        return 2
    if truth(row["either_route_auxiliary_positive"]):
        return 1
    return 0


def plot_depth_matrix(ax: plt.Axes, complete: pd.DataFrame) -> pd.DataFrame:
    selected = complete.loc[complete["pair_id"].isin(["pairA", "pairD"])].copy()
    pair_rank = {"pairA": 0, "pairD": 1}
    selected["pair_rank"] = selected["pair_id"].map(pair_rank)
    selected["direction_rank"] = selected.apply(
        lambda row: 0 if row["major"] == PAIR_SAMPLES[row["pair_id"]][0] else 1,
        axis=1,
    )
    row_keys = (
        selected[["pair_id", "major", "minor", "pair_rank", "direction_rank", "minor_percent"]]
        .drop_duplicates()
        .sort_values(["pair_rank", "direction_rank", "minor_percent"])
        .reset_index(drop=True)
    )
    if row_keys.shape[0] != 12:
        raise ValueError(f"Expected 12 depth-series rows; found {row_keys.shape[0]}")
    depths = [1_000_000, 2_000_000, 4_000_000]
    matrix = np.full((row_keys.shape[0], len(depths)), np.nan)
    expectation_met = np.full((row_keys.shape[0], len(depths)), True, dtype=bool)
    for row_index, key in row_keys.iterrows():
        subset = selected.loc[
            selected["pair_id"].eq(key["pair_id"])
            & selected["major"].eq(key["major"])
            & selected["minor"].eq(key["minor"])
            & selected["minor_percent"].eq(key["minor_percent"])
        ].set_index("total_pairs")
        for column_index, depth in enumerate(depths):
            matrix[row_index, column_index] = classify_condition(subset.loc[depth])
            expectation_met[row_index, column_index] = truth(
                subset.loc[depth, "prespecified_rule_satisfied"]
            )
    cmap = mpl.colors.ListedColormap([NEGATIVE, AUXILIARY, PRIMARY])
    ax.imshow(matrix, aspect="auto", interpolation="none", cmap=cmap, vmin=-0.5, vmax=2.5)
    ax.set_xticks(np.arange(len(depths)))
    ax.set_xticklabels(["1 million", "2 million", "4 million"])
    labels = [
        f"{row.pair_id[-1]}  {row.major} + {int(row.minor_percent)}% {row.minor}"
        for row in row_keys.itertuples(index=False)
    ]
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.tick_params(length=0)
    ax.set_title(
        "Nested-depth sensitivity of the complete two-route rule",
        loc="left",
        fontweight="bold",
        pad=3,
    )
    symbols = {0: "-", 1: "A", 2: "P"}
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            state = int(matrix[row_index, column_index])
            ax.text(
                column_index,
                row_index,
                symbols[state],
                ha="center",
                va="center",
                fontsize=6,
                fontweight="bold",
                color="white" if state == 2 else INK,
            )
            if not expectation_met[row_index, column_index]:
                ax.add_patch(
                    mpl.patches.Rectangle(
                        (column_index - 0.5, row_index - 0.5),
                        1,
                        1,
                        fill=False,
                        edgecolor="#B33A3A",
                        linewidth=1.15,
                    )
                )
    for boundary in (2.5, 5.5, 8.5):
        ax.axhline(boundary, color="white", linewidth=1.6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    legend = [
        mpl.patches.Patch(facecolor=NEGATIVE, label="-  neither rule positive"),
        mpl.patches.Patch(facecolor=AUXILIARY, label="A  auxiliary-only positive"),
        mpl.patches.Patch(facecolor=PRIMARY, label="P  primary positive"),
    ]
    if not expectation_met.all():
        legend.append(
            mpl.patches.Patch(
                facecolor="white",
                edgecolor="#B33A3A",
                label="Prespecified expectation not met",
            )
        )
    ax.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.20),
        ncol=len(legend),
        fontsize=5.7,
        handlelength=1.2,
        columnspacing=1.1,
    )
    return row_keys


def write_source_data(
    output: Path,
    stage1: pd.DataFrame,
    complete: pd.DataFrame,
    pair_ani: pd.DataFrame,
    stage1_summary: pd.DataFrame,
    ten_percent: pd.DataFrame,
) -> None:
    source = output / "source_data"
    source.mkdir(parents=True, exist_ok=True)
    for stale_name in (
        "Figure4_stage1_dilution_aggregate.tsv",
        "Figure4_stage2_selected_two_route.tsv",
    ):
        (source / stale_name).unlink(missing_ok=True)
    stage1.to_csv(source / "Figure4_stage1_technical_windows.tsv", sep="\t", index=False)
    stage1_summary.to_csv(
        source / "Figure4_stage1_detection_summary.tsv", sep="\t", index=False
    )
    complete.to_csv(source / "Figure4_complete_40_conditions.tsv", sep="\t", index=False)
    ten_percent.to_csv(
        source / "Figure4_complete_10_percent_conditions.tsv", sep="\t", index=False
    )
    pair_ani.to_csv(source / "Figure4_source_pair_ani.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    stage1, complete, ani, dependence = read_inputs(project)
    pair_ani = pair_ani_summary(ani)
    set_style()
    width_inches = 183 / 25.4
    height_inches = 151 / 25.4
    figure = plt.figure(figsize=(width_inches, height_inches), facecolor="white")
    grid = figure.add_gridspec(
        2,
        8,
        height_ratios=[0.95, 1.28],
        left=0.115,
        right=0.985,
        top=0.965,
        bottom=0.135,
        wspace=1.35,
        hspace=0.56,
    )
    ax_a = figure.add_subplot(grid[0, 0:2])
    ax_b = figure.add_subplot(grid[0, 2:5])
    ax_c = figure.add_subplot(grid[0, 5:8])
    ax_d = figure.add_subplot(grid[1, 1:7])

    plot_pair_ani(ax_a, pair_ani)
    stage1_summary = plot_stage1_detection(ax_b, stage1, dependence)
    ten_percent = plot_complete_ten_percent(ax_c, complete)
    plot_depth_matrix(ax_d, complete)
    for axis, label in zip((ax_a, ax_b, ax_c, ax_d), "abcd"):
        panel_label(axis, label)

    output = (
        project
        / "analysis_global_mac_upgrade/results/26_review_resolution_figures"
    )
    output.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".png", ".svg", ".tiff"):
        (output / f"Figure4_nearMAC_dilution_benchmark{suffix}").unlink(missing_ok=True)
    write_source_data(output, stage1, complete, pair_ani, stage1_summary, ten_percent)
    stem = output / "Figure4_nearMAC_expanded_benchmark"
    figure.savefig(f"{stem}.svg", facecolor="white")
    figure.savefig(f"{stem}.pdf", facecolor="white")
    figure.savefig(f"{stem}.png", dpi=300, facecolor="white")
    figure.savefig(
        f"{stem}.tiff",
        dpi=600,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(figure)
    print(stem)


if __name__ == "__main__":
    main()
