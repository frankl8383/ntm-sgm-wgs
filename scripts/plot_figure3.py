#!/usr/bin/env python3
"""Plot Figure 3: MI phylogeny and candidate syntenic intervals."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.patches import Polygon, Rectangle


MM = 1 / 25.4
TMI = "MI_TMI_lineage"
MP = "MI_MP_MIP_lineage"
CORE_BLOCKS = [
    "TMI_aromatic_catabolism_associated_block",
    "TMI_methyltransferase_hydrolase_cupin_block",
    "MP_MIP_nitrogen_redox_associated_block",
    "MP_MIP_oxidoreductase_associated_block",
]
BLOCK_LABELS = {
    CORE_BLOCKS[0]: "Aromatic\ncatabolism",
    CORE_BLOCKS[1]: "Methyl./hydrolase/\ncupin",
    CORE_BLOCKS[2]: "Nitrogen/\nredox",
    CORE_BLOCKS[3]: "Oxidoreductase",
}
ROW_LABELS = {
    CORE_BLOCKS[0]: "Aromatic",
    CORE_BLOCKS[1]: "Methyl./cupin",
    CORE_BLOCKS[2]: "N-redox",
    CORE_BLOCKS[3]: "Oxidoreductase",
}
HEAT_LABELS = {
    CORE_BLOCKS[0]: "Aromatic",
    CORE_BLOCKS[1]: "Methyl.",
    CORE_BLOCKS[2]: "N/redox",
    CORE_BLOCKS[3]: "Oxidored.",
}
COLORS = {
    "ink": "#242424",
    "muted": "#777777",
    "branch": "#B7B7B7",
    "absent": "#EEEEEE",
    "tmi": "#3E6FA3",
    "tmi_light": "#B9CCE2",
    "mp": "#B85C55",
    "mp_light": "#E2BBB7",
    "other": "#B6B6B6",
    "local": "#202020",
    "intervening": "#D9D9D9",
}
GENE_SUFFIX = re.compile(r"_(\d+)$")


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.6,
            "axes.titlesize": 7.4,
            "axes.labelsize": 6.8,
            "xtick.labelsize": 6.1,
            "ytick.labelsize": 6.1,
            "axes.linewidth": 0.7,
            "legend.fontsize": 6.1,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.06) -> None:
    ax.text(
        x,
        1.025,
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


def plot_tree_and_blocks(
    ax_tree: plt.Axes,
    ax_labels: plt.Axes,
    ax_heat: plt.Axes,
    tree_path: str,
    panel: pd.DataFrame,
    per_genome: pd.DataFrame,
) -> pd.DataFrame:
    tree = Phylo.read(tree_path, "newick")
    anchor_name = panel.loc[panel.source.eq("anchor"), "tree_id"].iat[0]
    anchor = next(tip for tip in tree.get_terminals() if tip.name == anchor_name)
    tree.root_with_outgroup(anchor)
    comparison_set = set(
        panel[
            panel.primary_tmi_vs_mp_mip_include.astype(str).str.lower().eq("true")
        ].tree_id
    )
    for tip in list(tree.get_terminals()):
        if tip.name not in comparison_set:
            tree.prune(tip)
    tree.ladderize(reverse=True)
    x, y = clade_coordinates(tree)
    terminals = sorted(tree.get_terminals(), key=lambda tip: y[tip])
    panel_index = panel.set_index("tree_id")
    max_x = max(x.values())

    for clade in tree.find_clades(order="preorder"):
        if clade.clades:
            child_y = [y[child] for child in clade.clades]
            ax_tree.plot(
                [x[clade], x[clade]],
                [min(child_y), max(child_y)],
                color=COLORS["branch"],
                linewidth=0.45,
                zorder=1,
            )
        if clade is not tree.root:
            parent = next(
                parent
                for parent in tree.find_clades(order="preorder")
                if clade in parent.clades
            )
            ax_tree.plot(
                [x[parent], x[clade]],
                [y[clade], y[clade]],
                color=COLORS["branch"],
                linewidth=0.45,
                zorder=1,
            )

    per_genome = per_genome[per_genome.block_id.isin(CORE_BLOCKS)].copy()
    per_genome["syntenic_block_present"] = (
        per_genome.syntenic_block_present.astype(str).str.lower().eq("true")
    )
    matrix = per_genome.pivot(
        index="tree_id", columns="block_id", values="syntenic_block_present"
    )
    source_rows: list[dict[str, object]] = []
    for tip in terminals:
        row = panel_index.loc[tip.name]
        lineage = row.accessory_lineage
        color = COLORS["tmi"] if lineage == TMI else COLORS["mp"] if lineage == MP else COLORS["other"]
        marker = "*" if row.source == "local" else "D" if row.source == "anchor" else "o"
        size = 19 if row.source == "local" else 10 if row.source == "anchor" else 5
        ax_tree.scatter(
            x[tip],
            y[tip],
            marker=marker,
            s=size,
            facecolor=color,
            edgecolor="white" if row.source == "local" else "none",
            linewidth=0.35,
            zorder=3,
        )
        if row.source == "local":
            ax_labels.text(
                0.02,
                y[tip],
                row.sample_id,
                fontsize=5.2,
                fontweight="bold",
                ha="left",
                va="center",
                color=color,
            )
        values = []
        for column_index, block in enumerate(CORE_BLOCKS):
            present = bool(matrix.loc[tip.name, block]) if tip.name in matrix.index else False
            values.append(int(present))
            if present:
                face = COLORS["tmi"] if block.startswith("TMI_") else COLORS["mp"]
            else:
                face = COLORS["absent"]
            ax_heat.add_patch(
                Rectangle(
                    (column_index - 0.43, y[tip] - 0.43),
                    0.86,
                    0.86,
                    facecolor=face,
                    edgecolor="white",
                    linewidth=0.25,
                )
            )
        source_rows.append(
            {
                "tree_id": tip.name,
                "source": row.source,
                "sample_id": row.sample_id,
                "accessory_lineage": lineage,
                **{block: value for block, value in zip(CORE_BLOCKS, values, strict=True)},
            }
        )

    y_limits = (-1, len(terminals))
    ax_tree.set_ylim(*y_limits)
    ax_labels.set_ylim(*y_limits)
    ax_heat.set_ylim(*y_limits)
    ax_tree.set_xlim(-0.002, max_x * 1.03)
    ax_tree.set_yticks([])
    ax_tree.set_xlabel("Substitutions per core SNP site")
    ax_tree.spines[["left", "right", "top"]].set_visible(False)
    ax_tree.tick_params(axis="x", length=2.5)
    ax_tree.set_title(
        f"TMI/MP-MIP complete-core SNP phylogeny (n = {len(comparison_set)})",
        loc="left",
        pad=8,
    )
    panel_label(ax_tree, "a")

    ax_labels.set_xlim(0, 1)
    ax_labels.axis("off")
    ax_heat.set_xlim(-0.5, len(CORE_BLOCKS) - 0.5)
    ax_heat.set_xticks(range(len(CORE_BLOCKS)))
    ax_heat.set_xticklabels(
        [HEAT_LABELS[block] for block in CORE_BLOCKS],
        rotation=90,
        ha="left",
        va="center",
        rotation_mode="anchor",
        fontsize=4.9,
    )
    ax_heat.xaxis.tick_top()
    ax_heat.tick_params(axis="x", length=0, pad=3)
    ax_heat.set_yticks([])
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["tmi"], markeredgecolor="none", label="TMI"),
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor=COLORS["mp"], markeredgecolor="none", label="MP-MIP"),
        mpl.lines.Line2D([], [], marker="*", linestyle="", markerfacecolor="#555555", markeredgecolor="none", label="Local genome"),
    ]
    ax_tree.legend(
        handles=handles,
        loc="upper left",
        ncol=3,
        handletextpad=0.25,
        columnspacing=0.7,
        borderaxespad=0.35,
    )
    return pd.DataFrame(source_rows)


def plot_prevalence(ax: plt.Axes, validation: pd.DataFrame) -> pd.DataFrame:
    panel_label(ax, "b", x=-0.12)
    ax.set_title("Full-atlas prevalence and local directional concordance", loc="left", pad=8)
    data = validation[validation.block_id.isin(CORE_BLOCKS)].copy()
    data["order"] = data.block_id.map({block: index for index, block in enumerate(CORE_BLOCKS)})
    data = data.sort_values("order")
    y = np.arange(len(data))[::-1]
    tmi = data.full_public_tmi_prevalence.astype(float).to_numpy() * 100
    mp = data.full_public_mp_mip_prevalence.astype(float).to_numpy() * 100
    for index, row in enumerate(data.itertuples(index=False)):
        yy = y[index]
        ax.plot([tmi[index], mp[index]], [yy, yy], color="#C6C6C6", linewidth=1)
        ax.scatter(tmi[index], yy, s=25, color=COLORS["tmi"], zorder=3)
        ax.scatter(mp[index], yy, s=25, color=COLORS["mp"], zorder=3)
        ax.text(tmi[index], yy + 0.18, f"{tmi[index]:.0f}", color=COLORS["tmi"], ha="center", fontsize=5.4)
        ax.text(mp[index], yy - 0.21, f"{mp[index]:.0f}", color=COLORS["mp"], ha="center", va="top", fontsize=5.4)
        local_tmi = f"{int(row.full_local_tmi_present)}/{int(row.full_local_tmi_total)}"
        local_mp = f"{int(row.full_local_mp_mip_present)}/{int(row.full_local_mp_mip_total)}"
        ax.text(
            103,
            yy,
            f"local {local_tmi} | {local_mp}",
            fontsize=5.2,
            color=COLORS["muted"],
            ha="left",
            va="center",
        )
    ax.set_xlim(0, 142)
    ax.set_ylim(-0.6, len(data) - 0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([ROW_LABELS[block] for block in data.block_id])
    ax.set_xlabel("Public genomes with candidate interval (%)")
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(22, len(data) - 0.48, "TMI (n=64)", color=COLORS["tmi"], fontsize=5.7, ha="center")
    ax.text(78, len(data) - 0.48, "MP-MIP (n=43)", color=COLORS["mp"], fontsize=5.7, ha="center")
    return data


def plot_project_forest(ax: plt.Axes, stats: pd.DataFrame) -> pd.DataFrame:
    panel_label(ax, "c", x=-0.12)
    ax.set_title("Association within mixed BioProjects", loc="left", pad=8)
    data = stats[stats.block_id.isin(CORE_BLOCKS)].copy()
    data["order"] = data.block_id.map({block: index for index, block in enumerate(CORE_BLOCKS)})
    data = data.sort_values("order")
    y = np.arange(len(data))[::-1]
    odds = data.mantel_haenszel_common_odds_ratio_expected_direction.astype(float).to_numpy()
    lower = data.mantel_haenszel_odds_ratio_95ci_lower.astype(float).to_numpy()
    upper = data.mantel_haenszel_odds_ratio_95ci_upper.astype(float).to_numpy()
    colors = [COLORS["tmi"] if direction == "TMI_enriched" else COLORS["mp"] for direction in data.association_direction]
    for index in range(len(data)):
        ax.plot([lower[index], upper[index]], [y[index], y[index]], color=colors[index], linewidth=1.2)
        ax.scatter(odds[index], y[index], s=26, color=colors[index], zorder=3)
        informative_column = (
            "nonzero_prevalence_gap_bioprojects"
            if "nonzero_prevalence_gap_bioprojects" in data.columns
            else "informative_bioprojects"
        )
        ax.text(
            305,
            y[index],
            f"{int(data.iloc[index][informative_column])}/5",
            ha="right",
            va="center",
            fontsize=5.0,
            color=COLORS["muted"],
        )
    ax.axvline(1, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlim(1, 320)
    ax.set_ylim(-0.6, len(data) - 0.4)
    ax.set_yticks(y)
    ax.set_yticklabels([ROW_LABELS[block] for block in data.block_id])
    ax.set_xlabel("Mantel-Haenszel odds ratio for expected lineage (95% CI)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(
        0.99,
        0.03,
        "Informative/mixed-lineage BioProjects shown by row",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color=COLORS["muted"],
        fontsize=5.4,
    )
    return data


def parse_gff(path: Path) -> pd.DataFrame:
    rows = []
    contig_counts: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            contig = fields[0]
            contig_counts[contig] = contig_counts.get(contig, 0) + 1
            gene_index = contig_counts[contig]
            rows.append(
                {
                    "contig": contig,
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "original_prodigal_id": f"{contig}_{gene_index}",
                    "gene_index_on_contig": gene_index,
                }
            )
    return pd.DataFrame(rows)


def arrow_polygon(x0: float, x1: float, y: float, strand: str) -> np.ndarray:
    height = 0.20
    head = min(0.035, (x1 - x0) * 0.42)
    if strand == "+":
        return np.asarray(
            [[x0, y - height], [x1 - head, y - height], [x1, y], [x1 - head, y + height], [x0, y + height]]
        )
    return np.asarray(
        [[x1, y - height], [x0 + head, y - height], [x0, y], [x0 + head, y + height], [x1, y + height]]
    )


def neighborhood_table(
    block_members: pd.DataFrame, sequence_meta: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    meta = sequence_meta.set_index("sequence_id")
    for block_id in CORE_BLOCKS:
        members = block_members[block_members.block_id.eq(block_id)].copy()
        representative_ids = list(members.cluster_id)
        representative = representative_ids[0].split("__gene")[0]
        member_meta = meta.loc[representative_ids].reset_index()
        contigs = member_meta.original_prodigal_id.map(lambda value: GENE_SUFFIX.sub("", value))
        contig = contigs.mode().iat[0]
        gff_path = Path(member_meta.per_genome_faa.iat[0]).with_suffix(".gff")
        gff = parse_gff(gff_path)
        selected = set(member_meta.original_prodigal_id)
        selected_rows = gff[gff.original_prodigal_id.isin(selected)]
        min_index = int(selected_rows.gene_index_on_contig.min())
        max_index = int(selected_rows.gene_index_on_contig.max())
        context = gff[
            gff.contig.eq(contig)
            & gff.gene_index_on_contig.between(min_index, max_index)
        ].copy()
        products = members.set_index("cluster_id").conservative_product.to_dict()
        original_to_cluster = member_meta.set_index("original_prodigal_id").sequence_id.to_dict()
        for row in context.itertuples(index=False):
            cluster_id = original_to_cluster.get(row.original_prodigal_id, "")
            rows.append(
                {
                    "block_id": block_id,
                    "representative_genome": representative,
                    "contig": row.contig,
                    "start": row.start,
                    "end": row.end,
                    "strand": row.strand,
                    "original_prodigal_id": row.original_prodigal_id,
                    "stable_family": bool(cluster_id),
                    "cluster_id": cluster_id,
                    "conservative_product": products.get(cluster_id, "intervening ORF"),
                }
            )
    return pd.DataFrame(rows)


def plot_neighborhoods(ax: plt.Axes, genes: pd.DataFrame) -> None:
    panel_label(ax, "d", x=-0.12)
    ax.set_title("Representative local gene order", loc="left", pad=8)
    y_positions = {block: 3 - index for index, block in enumerate(CORE_BLOCKS)}
    for block in CORE_BLOCKS:
        frame = genes[genes.block_id.eq(block)].copy().sort_values("start")
        start = frame.start.min()
        end = frame.end.max()
        span = max(1, end - start)
        y = y_positions[block]
        direction_color = COLORS["tmi"] if block.startswith("TMI_") else COLORS["mp"]
        for row in frame.itertuples(index=False):
            x0 = 0.36 + 0.61 * (row.start - start) / span
            x1 = 0.36 + 0.61 * (row.end - start) / span
            face = direction_color if row.stable_family else COLORS["intervening"]
            ax.add_patch(
                Polygon(
                    arrow_polygon(x0, x1, y, row.strand),
                    closed=True,
                    facecolor=face,
                    edgecolor="white",
                    linewidth=0.3,
                )
            )
        ax.text(
            0.0,
            y,
            f"{ROW_LABELS[block]} ({int(frame.stable_family.sum())})",
            ha="left",
            va="center",
            fontsize=5.6,
        )
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.55, 3.55)
    ax.axis("off")
    ax.text(
        0.36,
        -0.42,
        "coloured, stable family; grey, intervening ORF; row scales are independent",
        fontsize=5.0,
        color=COLORS["muted"],
        ha="left",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--panel-manifest", required=True)
    parser.add_argument("--block-per-genome", required=True)
    parser.add_argument("--block-validation", required=True)
    parser.add_argument("--project-statistics", required=True)
    parser.add_argument("--block-members", required=True)
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    configure()

    panel = pd.read_csv(args.panel_manifest, sep="\t", dtype=str).fillna("")
    per_genome = pd.read_csv(args.block_per_genome, sep="\t", dtype=str).fillna("")
    validation = pd.read_csv(args.block_validation, sep="\t")
    project_stats = pd.read_csv(args.project_statistics, sep="\t")
    block_members = pd.read_csv(args.block_members, sep="\t", dtype=str).fillna("")
    sequence_meta = pd.read_csv(args.sequence_metadata, sep="\t", dtype=str).fillna("")
    genes = neighborhood_table(block_members, sequence_meta)

    output_dir = Path(args.output_dir)
    source_dir = output_dir / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(183 * MM, 220 * MM))
    outer = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.72, 1.0],
        left=0.055,
        right=0.985,
        top=0.94,
        bottom=0.06,
        wspace=0.30,
    )
    left = outer[0].subgridspec(1, 3, width_ratios=[3.4, 0.72, 1.05], wspace=0.01)
    right = outer[1].subgridspec(3, 1, height_ratios=[1.0, 0.9, 1.08], hspace=0.48)
    ax_tree = fig.add_subplot(left[0, 0])
    ax_labels = fig.add_subplot(left[0, 1])
    ax_heat = fig.add_subplot(left[0, 2])
    ax_prev = fig.add_subplot(right[0, 0])
    ax_forest = fig.add_subplot(right[1, 0])
    ax_genes = fig.add_subplot(right[2, 0])

    tree_source = plot_tree_and_blocks(
        ax_tree, ax_labels, ax_heat, args.tree, panel, per_genome
    )
    prevalence_source = plot_prevalence(ax_prev, validation)
    project_source = plot_project_forest(ax_forest, project_stats)
    plot_neighborhoods(ax_genes, genes)

    tree_source.to_csv(
        source_dir / "Figure3a_tree_aligned_interval_presence.tsv", sep="\t", index=False
    )
    prevalence_source.to_csv(
        source_dir / "Figure3b_full_atlas_interval_prevalence.tsv", sep="\t", index=False
    )
    project_source.to_csv(
        source_dir / "Figure3c_bioproject_stratified_statistics.tsv", sep="\t", index=False
    )
    genes.to_csv(
        source_dir / "Figure3d_representative_gene_order.tsv", sep="\t", index=False
    )

    stem = output_dir / "Figure3_MI_lineage_syntenic_intervals"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(stem)


if __name__ == "__main__":
    main()
