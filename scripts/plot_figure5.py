#!/usr/bin/env python3
"""Plot external and public evaluation of four lineage-associated intervals."""

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
BLOCK_SHORT = {
    "TMI_aromatic_catabolism_associated_block": "TMI aromatic\n(26 families)",
    "TMI_methyltransferase_hydrolase_cupin_block": "TMI methyl/hydr.\n(5 families)",
    "MP_MIP_nitrogen_redox_associated_block": "MP-MIP N/redox\n(3 families)",
    "MP_MIP_oxidoreductase_associated_block": "MP-MIP oxidored.\n(5 families)",
}
TMI = "#416F96"
MP = "#B76555"
INK = "#292929"
MID = "#6E6E6E"
GRID = "#D7D7D7"
ABSENT = "#ECECEC"
WARNING = "#B68A2D"
MISMATCH = "#B33A3A"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    return parser.parse_args()


def truth(value: object) -> bool:
    return value is True or str(value).strip().lower() == "true"


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


def panel_label(ax: plt.Axes, label: str, x: float = -0.08) -> None:
    ax.text(
        x,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        fontweight="bold",
        color=INK,
        clip_on=False,
    )


def load_external(project: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    root = (
        project
        / "analysis_global_mac_upgrade/results/27_postreview_validation/"
        "external_PRJNA983112"
    )
    manifest = pd.read_csv(
        project
        / "analysis_global_mac_upgrade/config/"
        "external_block_validation_PRJNA983112_corrected_20260712.tsv",
        sep="\t",
    )
    lineage = pd.read_csv(root / "external_lineage_ani.tsv", sep="\t")
    qc = pd.read_csv(root / "external_assembly_qc.tsv", sep="\t")
    calls = pd.read_csv(root / "external_block_calls.tsv", sep="\t")
    reads = pd.read_csv(root / "external_block_read_support.tsv", sep="\t")
    expected_samples = set(manifest["sample_id"])
    if len(manifest) != 9 or set(lineage["sample_id"]) != expected_samples:
        raise ValueError("External lineage table does not match the nine-isolate manifest")
    if set(qc["sample_id"]) != expected_samples:
        raise ValueError("External QC table does not match the nine-isolate manifest")
    expected_pairs = {(sample, block) for sample in expected_samples for block in BLOCK_ORDER}
    for name, frame in (("block calls", calls), ("read support", reads)):
        observed = set(zip(frame["sample_id"], frame["block_id"]))
        if observed != expected_pairs:
            raise ValueError(f"External {name} does not contain all 36 sample-block pairs")
    lineage_fields = [
        "sample_id",
        "principal_lineage_by_anchor_ani",
        "principal_lineage_anchor_max_ani",
        "principal_lineage_qualified",
    ]
    qc_fields = ["sample_id", "external_qc_pass", "qc_failure_reasons"]
    annotation = manifest.merge(
        lineage[lineage_fields], on="sample_id", how="left", validate="one_to_one"
    ).merge(qc[qc_fields], on="sample_id", how="left", validate="one_to_one")
    annotation["qc_failure_reasons"] = annotation["qc_failure_reasons"].fillna("")
    required_annotation = lineage_fields[1:] + ["external_qc_pass"]
    if annotation[required_annotation].isna().any().any():
        raise ValueError("External annotation merge produced missing values")
    return annotation, calls, reads, qc


def load_public_prevalence(project: Path) -> pd.DataFrame:
    base = project / "analysis_global_mac_upgrade/results"
    specifications = [
        (
            "Separate public set",
            base
            / "24_public_non_discovery_evaluation/"
            "non_discovery_public_block_statistics.tsv",
            ("tmi_present", "tmi_total", "mp_mip_present", "mp_mip_total"),
        ),
        (
            "Complete genomes",
            base
            / "22_complete_genome_block_validation/"
            "complete_chromosome_block_prevalence.tsv",
            (
                "complete_chromosome_tmi_present",
                "complete_chromosome_tmi_total",
                "complete_chromosome_mp_mip_present",
                "complete_chromosome_mp_mip_total",
            ),
        ),
        (
            "One per project",
            base
            / "23_clone_project_sensitivity/"
            "one_per_bioproject_lineage_block_statistics.tsv",
            ("tmi_present", "tmi_total", "mp_mip_present", "mp_mip_total"),
        ),
        (
            "Near-clone reduced",
            base
            / "23_clone_project_sensitivity/"
            "strict_d20_kmer002_block_statistics.tsv",
            ("tmi_present", "tmi_total", "mp_mip_present", "mp_mip_total"),
        ),
    ]
    rows: list[dict[str, object]] = []
    for dataset, path, fields in specifications:
        frame = pd.read_csv(path, sep="\t").set_index("block_id")
        for block in BLOCK_ORDER:
            values = frame.loc[block]
            for lineage, present_key, total_key in (
                ("TMI", fields[0], fields[1]),
                ("MP-MIP", fields[2], fields[3]),
            ):
                present = int(values[present_key])
                total = int(values[total_key])
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
    return pd.DataFrame(rows)


def load_flanking_context(project: Path) -> pd.DataFrame:
    path = (
        project
        / "analysis_global_mac_upgrade/results/25_block_flanking_context/"
        "block_flanking_context_summary.tsv"
    )
    return pd.read_csv(path, sep="\t")


def plot_external_matrix(
    ax: plt.Axes,
    annotation: pd.DataFrame,
    calls: pd.DataFrame,
    reads: pd.DataFrame,
) -> pd.DataFrame:
    order = (
        annotation.assign(
            lineage_rank=annotation["frozen_expected_lineage"].map({"TMI": 0, "MP_MIP": 1})
        )
        .sort_values(["lineage_rank", "sample_id"])
        .reset_index(drop=True)
    )
    ax.set_xlim(-3.65, 4.55)
    ax.set_ylim(len(order) - 0.5, -1.55)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        -3.55,
        -1.43,
        "External evaluation: 9 records, 5 anchor-qualified (PRJNA983112)",
        fontsize=7.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax.text(-3.55, -0.72, "Isolate", fontsize=6, color=MID, ha="left")
    ax.text(-2.65, -0.72, "Published", fontsize=6, color=MID, ha="left")
    ax.text(-1.23, -0.72, "Anchor-qualified group", fontsize=6, color=MID, ha="center")
    ax.text(-0.52, -0.72, "QC", fontsize=6, color=MID, ha="center")
    for column_index, block in enumerate(BLOCK_ORDER):
        ax.text(
            column_index,
            -0.48,
            BLOCK_SHORT[block],
            ha="center",
            va="center",
            fontsize=5.7,
            fontweight="bold",
            color=INK,
            linespacing=1.05,
        )
    call_index = calls.set_index(["sample_id", "block_id"])
    read_index = reads.set_index(["sample_id", "block_id"])
    for row_index, row in order.iterrows():
        sample = row["sample_id"]
        expected = row["frozen_expected_lineage"]
        anchor = row["principal_lineage_by_anchor_ani"]
        expected_color = TMI if expected == "TMI" else MP
        ax.add_patch(
            Rectangle(
                (-3.64, row_index - 0.42),
                0.08,
                0.84,
                facecolor=expected_color,
                edgecolor="none",
            )
        )
        ax.text(-3.46, row_index, sample, ha="left", va="center", fontsize=6.3, color=INK)
        published = (
            "M. intracellulare"
            if row["published_species"] == "Mycobacterium intracellulare"
            else "M. paraintracellulare"
        )
        ax.text(
            -2.65,
            row_index,
            published,
            ha="left",
            va="center",
            fontsize=5.9,
            fontstyle="italic",
            color=INK,
        )
        anchor_label = {"TMI": "TMI", "MP_MIP": "MP-MIP"}.get(anchor, "outside")
        ax.text(
            -1.23,
            row_index,
            anchor_label,
            ha="center",
            va="center",
            fontsize=5.9,
            color=expected_color if anchor == expected else MISMATCH,
            fontweight="bold" if anchor == expected else "normal",
        )
        qc_pass = truth(row["external_qc_pass"])
        ax.scatter(
            [-0.52],
            [row_index],
            marker="o" if qc_pass else "^",
            s=18 if qc_pass else 24,
            facecolor="#808080" if qc_pass else WARNING,
            edgecolor="white",
            linewidth=0.5,
            zorder=4,
        )
        for column_index, block in enumerate(BLOCK_ORDER):
            call = call_index.loc[(sample, block)]
            read = read_index.loc[(sample, block)]
            present = truth(call["syntenic_block_present"])
            block_color = TMI if block.startswith("TMI_") else MP
            face = block_color if present else ABSENT
            ax.add_patch(
                Rectangle(
                    (column_index - 0.43, row_index - 0.39),
                    0.86,
                    0.78,
                    facecolor=face,
                    edgecolor="white",
                    linewidth=0.8,
                )
            )
            ax.text(
                column_index,
                row_index,
                "●" if present else "–",
                ha="center",
                va="center",
                fontsize=6.4,
                color="white" if present else MID,
            )
            if not truth(read["read_assembly_content_agreement"]):
                ax.add_patch(
                    Rectangle(
                        (column_index - 0.43, row_index - 0.39),
                        0.86,
                        0.78,
                        fill=False,
                        edgecolor=MISMATCH,
                        linewidth=1.1,
                    )
                )
    split = int((order["frozen_expected_lineage"] == "TMI").sum())
    ax.axhline(split - 0.5, color=INK, linewidth=0.7)
    ax.text(
        4.48,
        (split - 1) / 2,
        "TMI",
        ha="right",
        va="center",
        fontsize=6.1,
        color=TMI,
        fontweight="bold",
    )
    ax.text(
        4.48,
        split + (len(order) - split - 1) / 2,
        "MP-MIP",
        ha="right",
        va="center",
        fontsize=6.1,
        color=MP,
        fontweight="bold",
    )
    return order


def plot_public_matrix(ax: plt.Axes, data: pd.DataFrame) -> None:
    datasets = ["Separate public set", "Complete genomes", "One per project", "Near-clone reduced"]
    ax.set_xlim(-0.5, 7.5)
    ax.set_ylim(3.5, -1.88)
    ax.set_xticks([])
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels([BLOCK_LABELS[block] for block in BLOCK_ORDER], ha="right")
    ax.tick_params(axis="y", length=0, pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for dataset_index, dataset in enumerate(datasets):
        x0 = dataset_index * 2
        ax.text(
            x0 + 0.5,
            -1.00,
            dataset.replace(" ", "\n", 1),
            ha="center",
            va="bottom",
            fontsize=5.8,
            fontweight="bold",
            color=INK,
            clip_on=False,
        )
        for lineage_index, lineage in enumerate(("TMI", "MP-MIP")):
            x = x0 + lineage_index
            color = TMI if lineage == "TMI" else MP
            ax.text(
                x,
                -0.49,
                "T" if lineage == "TMI" else "M",
                ha="center",
                va="bottom",
                fontsize=5.7,
                color=color,
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
                ax.add_patch(
                    Rectangle(
                        (x - 0.43, block_index - 0.40),
                        0.86,
                        0.80,
                        facecolor=mpl.colors.to_rgba(color, 0.10 + 0.90 * prevalence),
                        edgecolor="white",
                        linewidth=0.8,
                    )
                )
                ax.text(
                    x,
                    block_index,
                    f"{int(row['present'])}/{int(row['total'])}",
                    ha="center",
                    va="center",
                    fontsize=5.5,
                    color="white" if prevalence >= 0.70 else INK,
                    fontweight="bold" if prevalence >= 0.70 else "normal",
                )
        if dataset_index < len(datasets) - 1:
            ax.axvline(x0 + 1.5, color=GRID, linewidth=0.6)
    ax.text(
        -0.48,
        -1.76,
        "Public and assembly sensitivities",
        ha="left",
        va="top",
        fontsize=7.5,
        fontweight="bold",
        color=INK,
    )
def draw_locus_row(
    ax: plt.Axes,
    y: float,
    label: str,
    block_color: str,
    block_present: bool,
    count_label: str,
) -> None:
    ax.text(0.01, y, label, ha="left", va="center", fontsize=6.0, color=INK)
    for x, width in ((0.18, 0.10), (0.52, 0.10)):
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
                (0.32, y - 0.075),
                0.16,
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
                (0.32, y - 0.06),
                0.16,
                0.12,
                linewidth=0.8,
                edgecolor=MID,
                facecolor="#F4F4F4",
                linestyle=(0, (2, 1.5)),
            )
        )
        ax.text(0.40, y, "absent", ha="center", va="center", fontsize=5.2, color=MID)
    ax.text(0.72, y, count_label, ha="left", va="center", fontsize=5.4, color=INK)


def context_count(
    context: pd.DataFrame,
    block_id: str,
    lineage: str,
    context_class: str,
) -> int:
    match = context.loc[
        context["block_id"].eq(block_id)
        & context["accessory_lineage"].eq(lineage)
        & context["structural_context_class"].eq(context_class),
        "genomes",
    ]
    return int(match.sum())


def lineage_total(context: pd.DataFrame, block_id: str, lineage: str) -> int:
    return int(
        context.loc[
            context["block_id"].eq(block_id)
            & context["accessory_lineage"].eq(lineage),
            "genomes",
        ].sum()
    )


def plot_structural_context(ax: plt.Axes, context: pd.DataFrame) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.01,
        0.99,
        "Flanking context in complete assemblies",
        ha="left",
        va="top",
        fontsize=7.5,
        fontweight="bold",
        color=INK,
    )
    ax.text(0.72, 0.89, "genomes", ha="left", va="center", fontsize=5.2, color=MID)
    mp_block = "MP_MIP_oxidoreductase_associated_block"
    tmi_block = "TMI_methyltransferase_hydrolase_cupin_block"
    mp_present = context_count(
        context, mp_block, "MI_MP_MIP_lineage", "block_and_both_flanks_supported"
    )
    mp_total = lineage_total(context, mp_block, "MI_MP_MIP_lineage")
    tmi_absent = context_count(
        context, mp_block, "MI_TMI_lineage", "both_flanks_supported_block_absent"
    )
    tmi_total = lineage_total(context, mp_block, "MI_TMI_lineage")
    tmi_present = context_count(
        context, tmi_block, "MI_TMI_lineage", "block_and_both_flanks_supported"
    )
    tmi_block_total = lineage_total(context, tmi_block, "MI_TMI_lineage")
    mp_absent = context_count(
        context, tmi_block, "MI_MP_MIP_lineage", "both_flanks_supported_block_absent"
    )
    mp_block_total = lineage_total(context, tmi_block, "MI_MP_MIP_lineage")
    ax.text(0.01, 0.82, "MP-MIP oxidoreductase interval", fontsize=6.1, fontweight="bold", color=MP)
    draw_locus_row(ax, 0.68, "MP-MIP", MP, True, f"{mp_present}/{mp_total}")
    draw_locus_row(ax, 0.51, "TMI", MP, False, f"{tmi_absent}/{tmi_total}")
    ax.text(0.01, 0.38, "TMI methyl/hydrolase interval", fontsize=6.1, fontweight="bold", color=TMI)
    draw_locus_row(
        ax, 0.24, "TMI", TMI, True, f"{tmi_present}/{tmi_block_total}"
    )
    draw_locus_row(
        ax, 0.09, "MP-MIP", TMI, False, f"{mp_absent}/{mp_block_total}"
    )
    ax.text(0.40, 0.89, "candidate interval", ha="center", va="center", fontsize=5.2, color=MID)
    ax.annotate(
        "",
        xy=(0.40, 0.84),
        xytext=(0.40, 0.875),
        arrowprops={"arrowstyle": "-", "color": MID, "lw": 0.6},
    )


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    annotation, calls, reads, qc = load_external(project)
    prevalence = load_public_prevalence(project)
    flanking_context = load_flanking_context(project)
    output = project / "analysis_global_mac_upgrade/results/26_review_resolution_figures"
    source = output / "source_data"
    output.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".png", ".svg", ".tiff"):
        (output / f"Figure5_structural_block_validation{suffix}").unlink(missing_ok=True)
    for stale_name in (
        "Figure5a_prevalence_sensitivity.tsv",
        "Figure5b_flanking_context.tsv",
        "Figure5c_Mi31_route_read_evidence.tsv",
        "Figure5d_minimum_change_summary.tsv",
    ):
        (source / stale_name).unlink(missing_ok=True)
    annotation.to_csv(source / "Figure5a_external_annotations.tsv", sep="\t", index=False)
    qc.to_csv(source / "Figure5a_external_qc.tsv", sep="\t", index=False)
    calls.to_csv(source / "Figure5a_external_interval_calls.tsv", sep="\t", index=False)
    reads.to_csv(source / "Figure5a_external_family_read_support.tsv", sep="\t", index=False)
    prevalence.to_csv(source / "Figure5b_public_sensitivities.tsv", sep="\t", index=False)
    flanking_context.to_csv(source / "Figure5c_flanking_context.tsv", sep="\t", index=False)

    set_style()
    figure = plt.figure(figsize=(183 / 25.4, 142 / 25.4), facecolor="white")
    grid = figure.add_gridspec(
        2,
        10,
        left=0.145,
        right=0.985,
        top=0.95,
        bottom=0.07,
        height_ratios=[1.18, 1.0],
        hspace=0.16,
        wspace=0.90,
    )
    ax_a = figure.add_subplot(grid[0, :])
    ax_b = figure.add_subplot(grid[1, :6])
    ax_c = figure.add_subplot(grid[1, 6:])
    plot_external_matrix(ax_a, annotation, calls, reads)
    plot_public_matrix(ax_b, prevalence)
    plot_structural_context(ax_c, flanking_context)
    for axis, label in zip((ax_a, ax_b, ax_c), "abc"):
        panel_label(axis, label)
    stem = output / "Figure5_external_interval_evaluation"
    figure.savefig(stem.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white", bbox_inches="tight")
    figure.savefig(
        stem.with_suffix(".tiff"),
        dpi=600,
        facecolor="white",
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(figure)
    print(stem)


if __name__ == "__main__":
    main()
