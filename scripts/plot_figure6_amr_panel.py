#!/usr/bin/env python3
"""Build draft Figure 6 AMR panel from AMRFinderPlus and NTM locus review tables."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle


SPECIES_ORDER = [
    "Mycobacterium avium",
    "Mycobacterium colombiense",
    "Mycobacterium intracellulare",
    "Mycobacterium paraintracellulare",
]

SPECIES_COLORS = {
    "Mycobacterium avium": "#0072B2",
    "Mycobacterium colombiense": "#CC79A7",
    "Mycobacterium intracellulare": "#E69F00",
    "Mycobacterium paraintracellulare": "#009E73",
}

STATUS_TO_CODE = {
    "not_detected_or_wildtype": 0,
    "amrfinderplus_stress_hit": 1,
    "low_confidence_amr_gene_hit": 2,
    "coding_variant_review_flag": 3,
    "technical_or_reference_warning": 4,
    "known_or_confident_amr_signal": 5,
}

STATUS_COLORS = {
    0: "#f1f1f1",
    1: "#59A14F",
    2: "#F28E2B",
    3: "#B07AA1",
    4: "#9C755F",
    5: "#D62728",
}

STATUS_LABELS = {
    0: "Not detected / no reviewed hotspot alternate",
    1: "AMRFinderPlus stress hit",
    2: "Low-confidence AMR gene hit",
    3: "Coding variant review flag",
    4: "Technical/reference warning",
    5: "Known hotspot or confident AMR signal",
}

AMR_COLUMNS = [
    ("AMRFinderPlus", "AMR:aph(3')-IIa", "aph(3')-IIa"),
    ("AMRFinderPlus", "STRESS:arsN1", "arsN1"),
]

LOCUS_COLUMNS = [
    ("NTM locus review", "rrl", "rrl"),
    ("NTM locus review", "rrs", "rrs"),
    ("NTM locus review", "erm", "erm"),
    ("NTM locus review", "gyrA", "gyrA"),
    ("NTM locus review", "gyrB", "gyrB"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amrfinder-matrix", required=True, type=Path)
    parser.add_argument("--amrfinder-summary", required=True, type=Path)
    parser.add_argument("--ntm-review", required=True, type=Path)
    parser.add_argument("--sample-table", required=True, type=Path)
    parser.add_argument("--out-png", required=True, type=Path)
    parser.add_argument("--out-pdf", required=True, type=Path)
    parser.add_argument("--matrix-out", required=True, type=Path)
    parser.add_argument("--sample-meta-out", required=True, type=Path)
    parser.add_argument("--column-meta-out", required=True, type=Path)
    return parser.parse_args()


def sample_order(sample_table: pd.DataFrame) -> pd.DataFrame:
    samples = sample_table[["sample_id", "amr_module_species", "checkm2_contamination", "evidence_conflict_flag"]].copy()
    species_rank = {species: i for i, species in enumerate(SPECIES_ORDER)}
    samples["species_rank"] = samples["amr_module_species"].map(species_rank).fillna(99).astype(int)
    samples = samples.sort_values(["species_rank", "sample_id"]).reset_index(drop=True)
    return samples


def amrfinder_cell(
    sample_id: str,
    element_label: str,
    amr_matrix: pd.DataFrame,
    amr_summary: pd.DataFrame,
) -> tuple[int, str, str]:
    if element_label not in amr_matrix.columns:
        return STATUS_TO_CODE["not_detected_or_wildtype"], "", "Element absent from AMRFinderPlus matrix."
    present = int(amr_matrix.loc[sample_id, element_label]) if sample_id in amr_matrix.index else 0
    if present == 0:
        return STATUS_TO_CODE["not_detected_or_wildtype"], "", "No AMRFinderPlus hit."
    hits = amr_summary[(amr_summary["sample_id"] == sample_id) & (amr_summary["element_label"] == element_label)]
    warning = "none"
    if not hits.empty and "hit_warning" in hits.columns:
        warning = ";".join(sorted(set(hits["hit_warning"].fillna("none").astype(str))))
    if element_label.startswith("STRESS:"):
        return STATUS_TO_CODE["amrfinderplus_stress_hit"], "+", "AMRFinderPlus --plus stress feature present."
    if warning and warning != "none":
        return STATUS_TO_CODE["low_confidence_amr_gene_hit"], "LC", f"AMRFinderPlus hit warning: {warning}."
    return STATUS_TO_CODE["known_or_confident_amr_signal"], "+", "AMRFinderPlus AMR feature present."


def summarize_locus_cell(ntm_review: pd.DataFrame, sample_id: str, locus: str) -> tuple[int, str, str]:
    sub = ntm_review[(ntm_review["sample_id"] == sample_id) & (ntm_review["locus"] == locus)].copy()
    if sub.empty:
        return STATUS_TO_CODE["technical_or_reference_warning"], "NA", "No review row was generated for this locus."

    statuses = set(sub["review_status"].fillna("").astype(str))
    if any("known_resistance_associated_alt_detected" in status for status in statuses):
        changes = ";".join(sub["mutation_or_change"].dropna().astype(str).unique()[:3])
        return STATUS_TO_CODE["known_or_confident_amr_signal"], changes or "+", "Known curated rRNA hotspot alternate base detected."
    if "confident_erm_like_hit" in statuses:
        changes = ";".join(sub["mutation_or_change"].dropna().astype(str).unique()[:3])
        return STATUS_TO_CODE["known_or_confident_amr_signal"], changes or "+", "Confident Erm-like protein hit detected."
    warning_statuses = {
        "low_identity_reference_for_manual_review_only",
        "uncalled_or_indel_codon",
        "locus_missing_or_below_threshold",
        "locus_not_evaluated",
        "not_evaluated",
    }
    if statuses & warning_statuses:
        changes = ";".join(
            str(x)
            for x in sub.loc[sub["review_status"].isin(warning_statuses), "mutation_or_change"].dropna().unique()[:2]
            if str(x) and str(x) != "nan"
        )
        label = "warn" if not changes else "warn"
        return STATUS_TO_CODE["technical_or_reference_warning"], label, "Technical or low-reference-identity warning; manual review required."
    if "nonsynonymous" in statuses:
        changes = [str(x) for x in sub.loc[sub["review_status"].eq("nonsynonymous"), "mutation_or_change"].dropna().unique()]
        label = changes[0] if len(changes) == 1 and len(changes[0]) <= 7 else f"{len(changes)} var"
        return STATUS_TO_CODE["coding_variant_review_flag"], label, "Amino-acid difference relative to selected public/reference coding locus."
    if locus in {"rrl", "rrs"} and "all_defined_hotspots_wildtype" in statuses:
        return STATUS_TO_CODE["not_detected_or_wildtype"], "no alt", "No reviewed NTM-Profiler rRNA hotspot alternate allele was detected."
    if locus == "erm" and "no_confident_erm_hit" in statuses:
        return STATUS_TO_CODE["not_detected_or_wildtype"], "-", "No confident Erm-like hit."
    if "no_nonsynonymous_changes_detected" in statuses:
        return STATUS_TO_CODE["not_detected_or_wildtype"], "-", "No nonsynonymous coding difference relative to selected reference."
    return STATUS_TO_CODE["technical_or_reference_warning"], "check", "Unexpected review status; inspect source table."


def build_plot_matrix(
    samples: pd.DataFrame,
    amr_matrix: pd.DataFrame,
    amr_summary: pd.DataFrame,
    ntm_review: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = AMR_COLUMNS + LOCUS_COLUMNS
    matrix_rows: list[dict[str, object]] = []
    cell_rows: list[dict[str, object]] = []
    for sample in samples.itertuples(index=False):
        sample_id = str(sample.sample_id)
        row: dict[str, object] = {"sample_id": sample_id}
        for group, key, label in columns:
            if group == "AMRFinderPlus":
                code, cell_label, note = amrfinder_cell(sample_id, key, amr_matrix, amr_summary)
            else:
                code, cell_label, note = summarize_locus_cell(ntm_review, sample_id, key)
            row[label] = code
            cell_rows.append(
                {
                    "sample_id": sample_id,
                    "column_group": group,
                    "feature": label,
                    "source_key": key,
                    "status_code": code,
                    "status_label": STATUS_LABELS[code],
                    "cell_label": cell_label,
                    "manual_review_note": note,
                }
            )
        matrix_rows.append(row)
    matrix = pd.DataFrame(matrix_rows).set_index("sample_id")
    cell_meta = pd.DataFrame(cell_rows)
    return matrix, cell_meta


def abbreviate_species(species: str) -> str:
    return (
        species.replace("Mycobacterium ", "M. ")
        .replace("intracellulare", "intra.")
        .replace("paraintracellulare", "paraintra.")
        .replace("colombiense", "colomb.")
    )


def wrapped_tick(label: str) -> str:
    if label == "aph(3')-IIa":
        return "aph(3')\n-IIa"
    return "\n".join(textwrap.wrap(label, width=8, break_long_words=False)) or label


def draw_figure(
    matrix: pd.DataFrame,
    cell_meta: pd.DataFrame,
    samples: pd.DataFrame,
    out_png: Path,
    out_pdf: Path,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    values = matrix.to_numpy(dtype=int)
    n_samples, n_features = values.shape
    fig_width = max(8.0, 1.02 * n_features + 3.8)
    fig_height = max(5.4, 0.39 * n_samples + 2.2)
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(
        nrows=1,
        ncols=3,
        width_ratios=[0.34, 1.0, 0.48],
        left=0.10,
        right=0.78,
        bottom=0.18,
        top=0.82,
        wspace=0.04,
    )
    ax_species = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_counts = fig.add_subplot(gs[0, 2], sharey=ax)

    cmap = ListedColormap([STATUS_COLORS[i] for i in range(6)])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5], cmap.N)
    ax.imshow(values, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")

    ax.set_yticks(range(n_samples))
    ax.set_yticklabels(samples["sample_id"], fontsize=8, fontweight="bold")
    ax.set_xticks(range(n_features))
    ax.set_xticklabels([wrapped_tick(c) for c in matrix.columns], fontsize=8, rotation=0)
    ax.tick_params(axis="both", length=0)
    ax.set_xlim(-0.5, n_features - 0.5)
    ax.set_ylim(n_samples - 0.5, -0.5)
    ax.set_xticks([x - 0.5 for x in range(1, n_features)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, n_samples)], minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    labels = cell_meta.pivot(index="sample_id", columns="feature", values="cell_label").reindex(matrix.index)[matrix.columns]
    for y, sample_id in enumerate(matrix.index):
        for x, feature in enumerate(matrix.columns):
            label = str(labels.loc[sample_id, feature])
            if label and label != "nan":
                color = "white" if values[y, x] in {1, 2, 3, 4, 5} else "#444444"
                ax.text(x, y, label, ha="center", va="center", fontsize=6.4, color=color)

    # Column group brackets.
    group_bounds = [(0, len(AMR_COLUMNS) - 1, "AMRFinderPlus gene screen"), (len(AMR_COLUMNS), n_features - 1, "Curated NTM locus review")]
    for start, end, label in group_bounds:
        ax.plot([start - 0.45, end + 0.45], [-1.15, -1.15], color="#333333", linewidth=1.0, clip_on=False)
        ax.plot([start - 0.45, start - 0.45], [-1.15, -0.83], color="#333333", linewidth=1.0, clip_on=False)
        ax.plot([end + 0.45, end + 0.45], [-1.15, -0.83], color="#333333", linewidth=1.0, clip_on=False)
        ax.text((start + end) / 2, -1.42, label, ha="center", va="bottom", fontsize=8, weight="bold", clip_on=False)
    ax.axvline(len(AMR_COLUMNS) - 0.5, color="#333333", linewidth=1.2)

    # Species color strip.
    ax_species.set_xlim(-1.35, 1)
    ax_species.set_ylim(n_samples - 0.5, -0.5)
    ax_species.axis("off")
    for y, species in enumerate(samples["amr_module_species"]):
        ax_species.text(-0.12, y, samples.iloc[y]["sample_id"], ha="right", va="center", fontsize=8, fontweight="bold")
        ax_species.add_patch(
            Rectangle((0.25, y - 0.38), 0.5, 0.76, color=SPECIES_COLORS.get(species, "#777777"), ec="white", lw=0.4)
        )
    ax_species.text(-0.12, -0.78, "Sample", ha="right", va="bottom", fontsize=8, weight="bold", clip_on=False)
    ax_species.text(0.5, -0.78, "Species", ha="center", va="bottom", fontsize=8, weight="bold", clip_on=False)

    # Review burden bars.
    review_counts = []
    amr_gene_counts = []
    stress_counts = []
    for sample_id in matrix.index:
        row_meta = cell_meta[cell_meta["sample_id"] == sample_id]
        review_counts.append(
            int(row_meta["status_code"].isin([STATUS_TO_CODE["coding_variant_review_flag"], STATUS_TO_CODE["technical_or_reference_warning"], STATUS_TO_CODE["known_or_confident_amr_signal"]]).sum())
        )
        amr_gene_counts.append(int(matrix.loc[sample_id, "aph(3')-IIa"] > 0))
        stress_counts.append(int(matrix.loc[sample_id, "arsN1"] == STATUS_TO_CODE["amrfinderplus_stress_hit"]))
    ax_counts.barh(range(n_samples), stress_counts, color=STATUS_COLORS[1], height=0.55, label="stress")
    ax_counts.barh(range(n_samples), amr_gene_counts, left=stress_counts, color=STATUS_COLORS[2], height=0.55, label="low-conf. AMR gene")
    left = [a + b for a, b in zip(stress_counts, amr_gene_counts)]
    ax_counts.barh(range(n_samples), review_counts, left=left, color=STATUS_COLORS[3], height=0.55, label="locus review flag")
    ax_counts.set_xlim(0, max(4, max([a + b + c for a, b, c in zip(stress_counts, amr_gene_counts, review_counts)]) + 0.6))
    ax_counts.set_xlabel("Feature count", fontsize=8)
    ax_counts.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, length=0)
    ax_counts.set_yticklabels([])
    ax_counts.tick_params(axis="x", labelsize=7)
    ax_counts.grid(axis="x", color="#dddddd", linewidth=0.5)
    for spine in ["top", "right", "left"]:
        ax_counts.spines[spine].set_visible(False)

    # Species separators.
    species_values = samples["amr_module_species"].tolist()
    for i in range(1, n_samples):
        if species_values[i] != species_values[i - 1]:
            for target_ax in [ax, ax_counts]:
                target_ax.axhline(i - 0.5, color="#333333", linewidth=1.1)
            ax_species.axhline(i - 0.5, color="#333333", linewidth=1.1)

    species_handles = [
        Patch(facecolor=color, edgecolor="none", label=abbreviate_species(species))
        for species, color in SPECIES_COLORS.items()
        if species in set(species_values)
    ]
    status_handles = [Patch(facecolor=STATUS_COLORS[i], edgecolor="none", label=STATUS_LABELS[i]) for i in range(6)]
    leg1 = fig.legend(handles=species_handles, loc="upper left", bbox_to_anchor=(0.80, 0.80), frameon=False, fontsize=8, title="Species", title_fontsize=8)
    fig.add_artist(leg1)
    fig.legend(handles=status_handles, loc="upper left", bbox_to_anchor=(0.80, 0.52), frameon=False, fontsize=8, title="Cell status", title_fontsize=8)

    fig.text(0.08, 0.93, "Figure 6. AMR/stress feature screening and curated NTM resistance-locus review", fontsize=12, weight="bold")
    fig.text(
        0.08,
        0.885,
        "Genome-only review; not clinical resistance prediction. High-confidence local MAC/SGM assemblies only (n=13). Low-confidence and coding-variant cells are manual-review prompts.",
        fontsize=8.5,
        color="#333333",
    )
    fig.text(
        0.08,
        0.055,
        "no alt: no reviewed rRNA hotspot alternate allele detected; this is not phenotypic susceptibility. LC: low-confidence AMRFinderPlus AMR hit; var: amino-acid difference relative to selected public/reference coding locus.",
        fontsize=7.5,
        color="#444444",
    )

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    amr_matrix = pd.read_csv(args.amrfinder_matrix, sep="\t")
    amr_summary = pd.read_csv(args.amrfinder_summary, sep="\t")
    ntm_review = pd.read_csv(args.ntm_review, sep="\t")
    sample_table = pd.read_csv(args.sample_table, sep="\t")

    samples = sample_order(sample_table)
    amr_matrix = amr_matrix.set_index("sample_id")
    matrix, cell_meta = build_plot_matrix(samples, amr_matrix, amr_summary, ntm_review)

    args.matrix_out.parent.mkdir(parents=True, exist_ok=True)
    args.sample_meta_out.parent.mkdir(parents=True, exist_ok=True)
    args.column_meta_out.parent.mkdir(parents=True, exist_ok=True)

    matrix_out = matrix.reset_index()
    matrix_out.to_csv(args.matrix_out, sep="\t", index=False)
    samples.drop(columns=["species_rank"]).to_csv(args.sample_meta_out, sep="\t", index=False)
    cell_meta.to_csv(args.column_meta_out, sep="\t", index=False)

    draw_figure(matrix, cell_meta, samples, args.out_png, args.out_pdf)


if __name__ == "__main__":
    main()
