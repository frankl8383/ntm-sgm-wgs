#!/usr/bin/env python3
"""Plot Figures 1-3 for the NTM SGM WGS project."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd


SPECIES_COLORS = {
    "Mycobacterium avium": "#4C78A8",
    "Mycobacterium intracellulare": "#F58518",
    "Mycobacterium paraintracellulare": "#54A24B",
    "Mycobacterium colombiense": "#B279A2",
    "Mycobacterium timonense": "#9C755F",
    "Mycobacterium yongonense": "#72B7B2",
    "mixed_or_conflict": "#ECA82C",
    "pass": "#4C9F70",
    "warn": "#F2C14E",
    "fail": "#D95F5F",
    "excluded": "#8C8C8C",
    "included": "#2F6F9F",
    "other": "#BDBDBD",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flow", required=True, type=Path)
    parser.add_argument("--fig2-input", required=True, type=Path)
    parser.add_argument("--fig3-input", required=True, type=Path)
    parser.add_argument("--figure1-png", required=True, type=Path)
    parser.add_argument("--figure1-pdf", required=True, type=Path)
    parser.add_argument("--figure2-png", required=True, type=Path)
    parser.add_argument("--figure2-pdf", required=True, type=Path)
    parser.add_argument("--figure3-png", required=True, type=Path)
    parser.add_argument("--figure3-pdf", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def savefig(fig: plt.Figure, png: Path, pdf: Path) -> None:
    png.parent.mkdir(parents=True, exist_ok=True)
    pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def species_short(text: str) -> str:
    text = str(text)
    mapping = [
        ("paraintracellulare", "M. para"),
        ("intracellulare", "M. intra"),
        ("colombiense", "M. colomb"),
        ("timonense", "M. timon"),
        ("yongonense", "M. yongo"),
        ("chimaera", "M. chim"),
        ("marseillense", "M. mars"),
        ("avium", "M. avium"),
    ]
    for key, val in mapping:
        if key in text:
            return val
    if not text or text == "NA":
        return "NA"
    return "Other"


def species_color(text: str) -> str:
    text = str(text)
    for species, color in SPECIES_COLORS.items():
        if species in text:
            return color
    if "conflict" in text.lower() or "ambiguous" in text.lower():
        return SPECIES_COLORS["mixed_or_conflict"]
    return SPECIES_COLORS["other"]


def wrap_label(text: str, width: int = 32) -> str:
    return "\n".join(textwrap.wrap(str(text), width=width, break_long_words=False))


def plot_figure1(flow: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.axis("off")
    steps = list(flow["step"])
    counts = dict(zip(flow["step"], flow["n"]))

    main_steps = [
        "Initial presumed SGM-NTM WGS samples",
        "Analysed clean FASTQ QC files parsed",
        "Read-level Mycobacterium-dominant or possible low-purity candidates",
        "Genome-level Mycobacterium candidate basic QC pass",
        "Priority genome-level reassessment set",
        "High-confidence MAC/SGM retained downstream",
    ]
    x_positions = np.linspace(0.06, 0.84, len(main_steps))
    y = 0.62
    box_w = 0.13
    box_h = 0.25
    colors = ["#E8F1F5", "#E8F1F5", "#DDECCB", "#CDE8D5", "#D8E3F3", "#B9D5EA"]

    for i, step in enumerate(main_steps):
        x = x_positions[i]
        rect = patches.FancyBboxPatch(
            (x - box_w / 2, y - box_h / 2),
            box_w,
            box_h,
            boxstyle="round,pad=0.015,rounding_size=0.018",
            facecolor=colors[i],
            edgecolor="#3A3A3A",
            linewidth=1.0,
        )
        ax.add_patch(rect)
        label = wrap_label(step.replace(" WGS samples", ""), width=19)
        ax.text(x, y + 0.025, label, ha="center", va="center", fontsize=8.2)
        ax.text(x, y - 0.085, f"n = {counts.get(step, 'NA')}", ha="center", va="center", fontsize=11, fontweight="bold")
        if i < len(main_steps) - 1:
            ax.annotate(
                "",
                xy=(x_positions[i + 1] - box_w / 2 - 0.01, y),
                xytext=(x + box_w / 2 + 0.01, y),
                arrowprops=dict(arrowstyle="->", lw=1.2, color="#3A3A3A"),
            )

    warning_step = "Species-supported but excluded due genome QC warning"
    ax.add_patch(
        patches.FancyBboxPatch(
            (0.49, 0.15),
            0.28,
            0.19,
            boxstyle="round,pad=0.015,rounding_size=0.018",
            facecolor="#F2E6C9",
            edgecolor="#8C6D31",
            linewidth=1.0,
        )
    )
    ax.text(
        0.63,
        0.25,
        f"Genome-QC warning branch\nn = {counts.get(warning_step, 'NA')}\nMi31: species signal retained, downstream excluded",
        ha="center",
        va="center",
        fontsize=8.6,
    )
    ax.annotate(
        "",
        xy=(0.60, 0.35),
        xytext=(x_positions[4], y - box_h / 2 - 0.015),
        arrowprops=dict(arrowstyle="->", lw=1.0, color="#8C6D31"),
    )

    ax.add_patch(
        patches.FancyBboxPatch(
            (0.79, 0.15),
            0.18,
            0.19,
            boxstyle="round,pad=0.015,rounding_size=0.018",
            facecolor="#EFEFEF",
            edgecolor="#5F5F5F",
            linewidth=1.0,
        )
    )
    ax.text(
        0.88,
        0.25,
        "Downstream modules\npublic context\nphylogeny\nAMR/stress review",
        ha="center",
        va="center",
        fontsize=8.6,
    )
    ax.annotate("", xy=(0.79, 0.25), xytext=(x_positions[-1] + box_w / 2 + 0.01, y - 0.02), arrowprops=dict(arrowstyle="->", lw=1.0, color="#5F5F5F"))

    ax.text(0.02, 0.94, "Figure 1. Contamination-aware WGS reassessment workflow", ha="left", va="center", fontsize=13, fontweight="bold")
    ax.text(
        0.02,
        0.04,
        "Exclusion is treated as an evidence-supported WGS result, not failed data.",
        ha="left",
        va="center",
        fontsize=8.5,
        color="#555555",
    )
    savefig(fig, out_png, out_pdf)


def inclusion_order(value: str) -> int:
    order = {
        "included_high_confidence_MAC_SGM": 0,
        "excluded_species_supported_but_genome_QC_warning": 1,
        "excluded_mixed_contaminated_or_low_quality": 2,
        "excluded_read_level_non_mycobacterium_or_contaminated": 3,
        "excluded_not_prioritized_for_downstream": 4,
    }
    return order.get(value, 9)


def plot_figure2(df: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    data = df.copy()
    data["mycobacterium_fraction"] = pd.to_numeric(data["mycobacterium_fraction"], errors="coerce").fillna(0)
    data["checkm2_contamination"] = pd.to_numeric(data["checkm2_contamination"], errors="coerce")
    data["gunc_pass_or_fail_bool"] = data["gunc_pass_or_fail"].astype(str).str.lower().eq("true")
    data["sort_key"] = data["final_inclusion_class"].map(inclusion_order)
    data = data.sort_values(["sort_key", "mycobacterium_fraction", "sample_id"], ascending=[True, False, True])
    y = np.arange(len(data))

    fig = plt.figure(figsize=(15.2, 12.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.35, 1.35, 0.95, 1.05], wspace=0.30)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharey=ax1)

    myco = data["mycobacterium_fraction"].clip(0, 1)
    non_myco = (1 - myco).clip(0, 1)
    ax1.barh(y, myco, color="#4C78A8", label="Mycobacterium fraction")
    ax1.barh(y, non_myco, left=myco, color="#CFCFCF", label="Non-Mycobacterium/other")
    ax1.set_xlim(0, 1)
    ax1.set_xlabel("Read-level fraction")
    ax1.set_title("A  Read-level genus signal", loc="left", fontweight="bold")
    ax1.set_yticks(y)
    ax1.set_yticklabels(data["sample_id"], fontsize=8)
    ax1.invert_yaxis()
    ax1.legend(loc="lower right", fontsize=7, frameon=False)

    cont = data["checkm2_contamination"]
    colors = ["#4C9F70" if c <= 5 else "#F2C14E" if c <= 10 else "#D95F5F" for c in cont.fillna(999)]
    clipped_cont = cont.fillna(0).clip(upper=20)
    ax2.barh(y, clipped_cont, color=colors)
    ax2.axvline(5, color="#555555", ls="--", lw=1)
    ax2.axvline(10, color="#555555", ls=":", lw=1)
    ax2.set_xlabel("CheckM2 contamination (%)")
    ax2.set_title("B  CheckM2 contamination", loc="left", fontweight="bold")
    ax2.tick_params(axis="y", labelleft=False)
    ax2.set_xlim(0, 22)
    ax2.set_xticks([0, 5, 10, 20])
    for i, c in enumerate(cont):
        if pd.notna(c) and c > 20:
            ax2.scatter([20.7], [i], marker=">", color="#333333", s=18, clip_on=False)
            ax2.text(21.05, i, f"{c:.0f}", va="center", ha="left", fontsize=6.2, color="#333333", clip_on=False)

    qc_matrix = []
    for _, row in data.iterrows():
        checkm_pass = row["checkm2_contamination"] <= 5 if pd.notna(row["checkm2_contamination"]) else False
        gunc_pass = bool(row["gunc_pass_or_fail_bool"])
        included = row["final_inclusion_class"] == "included_high_confidence_MAC_SGM"
        qc_matrix.append([1 if checkm_pass else 0, 1 if gunc_pass else 0, 1 if included else 0])
    qc = np.array(qc_matrix)
    color_grid = np.empty(qc.shape, dtype=object)
    for i in range(qc.shape[0]):
        color_grid[i, 0] = SPECIES_COLORS["pass"] if qc[i, 0] else SPECIES_COLORS["fail"]
        color_grid[i, 1] = SPECIES_COLORS["pass"] if qc[i, 1] else SPECIES_COLORS["fail"]
        color_grid[i, 2] = SPECIES_COLORS["included"] if qc[i, 2] else SPECIES_COLORS["excluded"]
    for i in range(qc.shape[0]):
        for j in range(qc.shape[1]):
            ax3.add_patch(patches.Rectangle((j, i - 0.45), 1, 0.9, facecolor=color_grid[i, j], edgecolor="white", linewidth=1))
            ax3.text(j + 0.5, i, "Y" if qc[i, j] else "N", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
    ax3.set_xlim(0, 3)
    ax3.set_ylim(-0.5, len(data) - 0.5)
    ax3.invert_yaxis()
    ax3.set_xticks([0.5, 1.5, 2.5])
    ax3.set_xticklabels(["CheckM2\n<=5%", "GUNC\npass", "Downstream\nincluded"], fontsize=8)
    ax3.tick_params(axis="y", labelleft=False)
    ax3.set_title("C  QC gates", loc="left", fontweight="bold")
    ax3.set_frame_on(False)

    class_colors = {
        "included_high_confidence_MAC_SGM": "#2F6F9F",
        "excluded_species_supported_but_genome_QC_warning": "#B8860B",
        "excluded_mixed_contaminated_or_low_quality": "#D95F5F",
        "excluded_read_level_non_mycobacterium_or_contaminated": "#7F7F7F",
        "excluded_not_prioritized_for_downstream": "#BDBDBD",
    }
    class_labels = {
        "included_high_confidence_MAC_SGM": "Included high-confidence MAC/SGM",
        "excluded_species_supported_but_genome_QC_warning": "Species signal, QC warning",
        "excluded_mixed_contaminated_or_low_quality": "Mixed/contaminated/low quality",
        "excluded_read_level_non_mycobacterium_or_contaminated": "Read-level non-Mycobacterium",
        "excluded_not_prioritized_for_downstream": "Not prioritized downstream",
    }
    for i, (_, row) in enumerate(data.iterrows()):
        cls = row["final_inclusion_class"]
        ax4.add_patch(patches.Rectangle((0, i - 0.45), 1, 0.9, facecolor=class_colors.get(cls, "#BDBDBD"), edgecolor="white", linewidth=1))
    ax4.set_xlim(0, 1)
    ax4.set_ylim(-0.5, len(data) - 0.5)
    ax4.invert_yaxis()
    ax4.axis("off")
    ax4.set_title("D  Final inclusion class", loc="left", fontweight="bold")
    handles = [patches.Patch(color=color, label=class_labels.get(cls, cls)) for cls, color in class_colors.items()]
    ax4.legend(handles=handles, loc="center left", bbox_to_anchor=(1.03, 0.5), frameon=False, fontsize=7)

    fig.suptitle("Figure 2. Read-level and assembly-level evidence for contamination-aware sample inclusion", x=0.02, ha="left", fontsize=13, fontweight="bold")
    savefig(fig, out_png, out_pdf)


def evidence_state(sample: pd.Series, evidence_col: str) -> tuple[str, str, str]:
    current = sample["current_final_norm"]
    if evidence_col == "public":
        species = sample["public_final_best_species"]
        if species == current:
            return "support", species_short(species), species_color(species)
        return "conflict", species_short(species), SPECIES_COLORS["mixed_or_conflict"]
    if evidence_col == "type":
        species = sample["type_panel_v2_best_species"]
        if species == current:
            return "support", species_short(species), species_color(species)
        return "conflict", species_short(species), SPECIES_COLORS["mixed_or_conflict"]
    if evidence_col == "ntm":
        species = sample["ntm_profiler_predicted_species"]
        if species == current:
            return "support", species_short(species), species_color(species)
        return "conflict", species_short(species), SPECIES_COLORS["mixed_or_conflict"]
    if evidence_col == "markers":
        marker_set = sample["independent_marker_top1_species_set"]
        if current in marker_set:
            return "support", "marker\nsupports", SPECIES_COLORS["pass"]
        return "conflict", species_short(marker_set), SPECIES_COLORS["mixed_or_conflict"]
    if evidence_col == "qc":
        gunc = str(sample["gunc_status"]).lower() == "true"
        contam = pd.to_numeric(sample["checkm2_contamination"], errors="coerce")
        if gunc and pd.notna(contam) and contam <= 5:
            return "support", "pass", SPECIES_COLORS["pass"]
        return "fail", "QC\nwarn", SPECIES_COLORS["fail"]
    if evidence_col == "integrated":
        category = str(sample.get("integrated_support_category", "not_run"))
        if category == "integrated_type_public_ANI_supports_analysis_clade":
            return "support", "type+\npublic", SPECIES_COLORS["pass"]
        if category == "integrated_ANI_supports_neighbor_clade_not_current_analysis_label":
            return "conflict", "neighbor\nclade", SPECIES_COLORS["warn"]
        if category == "integrated_type_public_ANI_conflict_MAC_boundary":
            return "conflict", "MAC\nboundary", SPECIES_COLORS["mixed_or_conflict"]
        if category == "QC_warning_excluded_despite_integrated_ANI_support":
            return "fail", "QC\nexcluded", SPECIES_COLORS["excluded"]
        return "other", "not\nrun", SPECIES_COLORS["other"]
    if evidence_col == "tier":
        tier = sample["species_level_confidence_tier"]
        if tier in {"MAC_complex_high_species_clade_supported", "MAC_complex_high_integrated_type_public_support"}:
            return "support", "type/\npublic", SPECIES_COLORS["pass"]
        if "QC_warning" in tier:
            return "fail", "QC\nexcluded", SPECIES_COLORS["excluded"]
        if "neighbor_clade" in tier:
            return "conflict", "neighbor\nclade", SPECIES_COLORS["warn"]
        return "conflict", "clade\ncaution", SPECIES_COLORS["mixed_or_conflict"]
    if evidence_col == "include":
        incl = "included_high_confidence" in str(sample["downstream_inclusion_status_revised"])
        return ("support", "include", SPECIES_COLORS["included"]) if incl else ("fail", "exclude", SPECIES_COLORS["excluded"])
    return "other", "NA", SPECIES_COLORS["other"]


def plot_figure3(df: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    data = df.copy()
    data["included"] = data["downstream_inclusion_status_revised"].astype(str).str.contains("included_high_confidence")
    data["order"] = data["included"].map({True: 0, False: 1})
    data = data.sort_values(["order", "current_final_norm", "sample_id"])

    cols = [
        ("public", "Public\nnearest"),
        ("type", "Type-panel\nANI"),
        ("integrated", "Integrated\nANI"),
        ("ntm", "NTM-\nProfiler"),
        ("markers", "Independent\nmarkers"),
        ("qc", "CheckM2+\nGUNC"),
        ("tier", "Species\nconfidence"),
        ("include", "Downstream"),
    ]
    nrow, ncol = len(data), len(cols)
    fig, ax = plt.subplots(figsize=(14.3, 7.7))
    for i, (_, sample) in enumerate(data.iterrows()):
        is_neighbor_warning = str(sample.get("integrated_support_category", "")) == "integrated_ANI_supports_neighbor_clade_not_current_analysis_label"
        if is_neighbor_warning:
            ax.add_patch(
                patches.Rectangle(
                    (-1.56, i + 0.03),
                    ncol + 1.52,
                    0.94,
                    facecolor="none",
                    edgecolor="#8A5A00",
                    linewidth=2.0,
                    zorder=8,
                )
            )
        ax.text(-0.78, i + 0.5, sample["sample_id"], ha="right", va="center", fontsize=8, fontweight="bold")
        ax.text(-0.03, i + 0.5, species_short(sample["current_final_norm"]), ha="right", va="center", fontsize=7.5)
        for j, (key, _) in enumerate(cols):
            _, label, color = evidence_state(sample, key)
            ax.add_patch(patches.Rectangle((j, i), 1, 1, facecolor=color, edgecolor="white", linewidth=1.4))
            ax.text(j + 0.5, i + 0.5, label, ha="center", va="center", fontsize=7, color="white" if color not in {"#F2C14E", "#ECA82C"} else "#222222", fontweight="bold")

    ax.set_xlim(-1.6, ncol)
    ax.set_ylim(nrow, 0)
    ax.set_xticks([j + 0.5 for j in range(ncol)])
    ax.set_xticklabels([label for _, label in cols], fontsize=8)
    ax.set_yticks([])
    ax.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.set_frame_on(False)
    ax.text(-1.52, -0.35, "Sample", ha="left", va="center", fontsize=8, fontweight="bold")
    ax.text(-0.62, -0.35, "Analysis\ncall", ha="left", va="center", fontsize=8, fontweight="bold")
    ax.set_title("Figure 3. Multi-evidence WGS reassessment matrix with explicit MAC species-level confidence tiers", loc="left", fontsize=13, fontweight="bold", pad=34)

    legend_items = [
        ("type/public support; unqualified wording only where tier permits", SPECIES_COLORS["pass"]),
        ("MAC boundary conflict / clade wording", SPECIES_COLORS["mixed_or_conflict"]),
        ("neighbor-clade warning", SPECIES_COLORS["warn"]),
        ("QC fail/excluded", SPECIES_COLORS["fail"]),
    ]
    x0, y0 = 0.02, -0.09
    for k, (text, color) in enumerate(legend_items):
        ax.add_patch(patches.Rectangle((x0 + k * 0.235, y0), 0.025, 0.035, transform=ax.transAxes, facecolor=color, clip_on=False))
        ax.text(x0 + k * 0.235 + 0.03, y0 + 0.018, text, transform=ax.transAxes, va="center", fontsize=7.3)
    ax.text(
        0.02,
        -0.165,
        "Abbreviations: M. intra, M. intracellulare; M. para, M. paraintracellulare; M. timon, M. timonense; M. colomb, M. colombiense. Wording tiers are interpretation safeguards, not formal taxonomic acts. Brown row border marks neighbor-clade warning.",
        transform=ax.transAxes,
        va="center",
        fontsize=7.1,
        color="#444444",
    )

    savefig(fig, out_png, out_pdf)


def main() -> int:
    args = parse_args()
    flow = read_tsv(args.flow)
    fig2 = read_tsv(args.fig2_input)
    fig3 = read_tsv(args.fig3_input)

    plot_figure1(flow, args.figure1_png, args.figure1_pdf)
    plot_figure2(fig2, args.figure2_png, args.figure2_pdf)
    plot_figure3(fig3, args.figure3_png, args.figure3_pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
