#!/usr/bin/env python3
"""Parse integrated-panel FastANI and plot ANI network + clustering support figures."""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


try:
    import networkx as nx
except Exception:  # pragma: no cover - optional runtime dependency
    nx = None


SPECIES_COLORS = {
    "M. avium": "#4C78A8",
    "M. intracellulare": "#F58518",
    "M. paraintracellulare": "#54A24B",
    "M. colombiense": "#B279A2",
    "M. timonense": "#9C755F",
    "M. bouchedurhonense": "#A0CBE8",
    "M. yongonense": "#72B7B2",
    "M. chimaera": "#BAB0AC",
    "M. marseillense": "#E45756",
    "M. arosiense": "#FF9DA6",
    "M. xenopi": "#8CD17D",
    "Other": "#BDBDBD",
}
ROLE_MARKERS = {
    "local_downstream": "*",
    "local_qc_warning": "X",
    "type_or_representative": "^",
    "public_near_neighbor": "o",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-metadata", required=True, type=Path)
    parser.add_argument("--fastani", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--network-ani-threshold", type=float, default=96.0)
    parser.add_argument("--network-af-threshold", type=float, default=0.50)
    parser.add_argument("--top-neighbors", type=int, default=12)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def species_color(species_short: str) -> str:
    for key, color in SPECIES_COLORS.items():
        if str(species_short).startswith(key):
            return color
    return SPECIES_COLORS["Other"]


def canonical_species(value: object) -> str:
    text = str(value)
    species_names = [
        "Mycobacterium paraintracellulare",
        "Mycobacterium intracellulare",
        "Mycobacterium colombiense",
        "Mycobacterium timonense",
        "Mycobacterium bouchedurhonense",
        "Mycobacterium yongonense",
        "Mycobacterium chimaera",
        "Mycobacterium marseillense",
        "Mycobacterium arosiense",
        "Mycobacterium avium",
        "Mycobacterium xenopi",
        "Mycobacterium indicus pranii",
        "Mycobacterium gordonae",
        "Mycobacterium kansasii",
        "Mycobacterium simiae",
        "Mycobacterium malmoense",
        "Mycobacterium marinum",
        "Mycobacterium ulcerans",
    ]
    for species in species_names:
        if species in text:
            return species
    return text


def parse_fastani(path: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    fasta_to_node = dict(zip(metadata["fasta"], metadata["node_id"]))
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            q, r, ani, matched, total = parts[:5]
            qid = fasta_to_node.get(q)
            rid = fasta_to_node.get(r)
            if qid is None or rid is None:
                continue
            total_i = int(total)
            rows.append(
                {
                    "query_node": qid,
                    "reference_node": rid,
                    "query_fasta": q,
                    "reference_fasta": r,
                    "ani": float(ani),
                    "matched_fragments": int(matched),
                    "total_fragments": total_i,
                    "alignment_fraction": int(matched) / total_i if total_i else np.nan,
                }
            )
    return pd.DataFrame(rows)


def symmetric_pairs(all_hits: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    by_pair: dict[tuple[str, str], list[pd.Series]] = defaultdict(list)
    for _, row in all_hits.iterrows():
        a, b = sorted([row["query_node"], row["reference_node"]])
        if a == b:
            continue
        by_pair[(a, b)].append(row)

    meta = metadata.set_index("node_id")
    rows: list[dict[str, object]] = []
    for (a, b), vals in by_pair.items():
        ani = float(np.mean([v["ani"] for v in vals]))
        af = float(np.mean([v["alignment_fraction"] for v in vals]))
        rows.append(
            {
                "node_a": a,
                "node_b": b,
                "ani_mean": ani,
                "alignment_fraction_mean": af,
                "node_a_role": meta.loc[a, "node_role"],
                "node_b_role": meta.loc[b, "node_role"],
                "node_a_species": meta.loc[a, "species"],
                "node_b_species": meta.loc[b, "species"],
                "node_a_label": meta.loc[a, "display_label"],
                "node_b_label": meta.loc[b, "display_label"],
                "node_a_source": meta.loc[a, "source"],
                "node_b_source": meta.loc[b, "source"],
            }
        )
    return pd.DataFrame(rows).sort_values(["ani_mean", "alignment_fraction_mean"], ascending=[False, False])


def build_ani_matrix(pairs: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    ids = list(metadata["node_id"])
    mat = pd.DataFrame(np.nan, index=ids, columns=ids)
    np.fill_diagonal(mat.values, 100.0)
    for _, row in pairs.iterrows():
        mat.loc[row["node_a"], row["node_b"]] = row["ani_mean"]
        mat.loc[row["node_b"], row["node_a"]] = row["ani_mean"]
    return mat


def connected_components(nodes: Iterable[str], pairs: pd.DataFrame, threshold: float) -> dict[str, int]:
    graph: dict[str, set[str]] = {n: set() for n in nodes}
    for _, row in pairs.iterrows():
        if row["ani_mean"] >= threshold:
            graph[row["node_a"]].add(row["node_b"])
            graph[row["node_b"]].add(row["node_a"])
    comp_id: dict[str, int] = {}
    comp = 0
    for node in nodes:
        if node in comp_id:
            continue
        comp += 1
        queue = deque([node])
        comp_id[node] = comp
        while queue:
            cur = queue.popleft()
            for nxt in graph[cur]:
                if nxt not in comp_id:
                    comp_id[nxt] = comp
                    queue.append(nxt)
    return comp_id


def local_neighbor_tables(
    metadata: pd.DataFrame,
    pairs: pd.DataFrame,
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = metadata.set_index("node_id")
    local_nodes = metadata[metadata["node_role"].str.startswith("local")]["node_id"].tolist()
    rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for node in local_nodes:
        sub = pairs[(pairs["node_a"].eq(node)) | (pairs["node_b"].eq(node))].copy()
        def other(row: pd.Series) -> str:
            return row["node_b"] if row["node_a"] == node else row["node_a"]
        sub["neighbor_node"] = sub.apply(other, axis=1)
        sub["neighbor_role"] = sub["neighbor_node"].map(meta["node_role"])
        sub["neighbor_species"] = sub["neighbor_node"].map(meta["species"])
        sub["neighbor_species_short"] = sub["neighbor_node"].map(meta["species_short"])
        sub["neighbor_label"] = sub["neighbor_node"].map(meta["display_label"])
        sub = sub.sort_values(["ani_mean", "alignment_fraction_mean"], ascending=[False, False])
        for rank, (_, row) in enumerate(sub.head(top_n).iterrows(), start=1):
            rows.append(
                {
                    "sample_id": meta.loc[node, "sample_id"],
                    "local_node": node,
                    "local_analysis_group": meta.loc[node, "analysis_group"],
                    "local_species_confidence_tier": meta.loc[node, "species_confidence_tier"],
                    "rank": rank,
                    "neighbor_node": row["neighbor_node"],
                    "neighbor_role": row["neighbor_role"],
                    "neighbor_species": row["neighbor_species"],
                    "neighbor_label": row["neighbor_label"],
                    "ani_mean": row["ani_mean"],
                    "alignment_fraction_mean": row["alignment_fraction_mean"],
                    "above_95": row["ani_mean"] >= 95,
                    "above_96": row["ani_mean"] >= 96,
                    "above_98": row["ani_mean"] >= 98,
                }
            )
        best_type = sub[sub["neighbor_role"].eq("type_or_representative")].head(1)
        best_public = sub[sub["neighbor_role"].eq("public_near_neighbor")].head(1)
        best_nonlocal = sub[~sub["neighbor_role"].str.startswith("local")].head(1)
        best_any = sub.head(1)
        def get_best(df: pd.DataFrame, field: str) -> object:
            return df.iloc[0][field] if not df.empty else "NA"
        summary_rows.append(
            {
                "sample_id": meta.loc[node, "sample_id"],
                "node_role": meta.loc[node, "node_role"],
                "analysis_group": meta.loc[node, "analysis_group"],
                "analysis_group_canonical": canonical_species(meta.loc[node, "analysis_group"]),
                "species_confidence_tier": meta.loc[node, "species_confidence_tier"],
                "best_any_neighbor": get_best(best_any, "neighbor_label"),
                "best_any_species": get_best(best_any, "neighbor_species"),
                "best_any_ani": get_best(best_any, "ani_mean"),
                "best_nonlocal_neighbor": get_best(best_nonlocal, "neighbor_label"),
                "best_nonlocal_species": get_best(best_nonlocal, "neighbor_species"),
                "best_nonlocal_ani": get_best(best_nonlocal, "ani_mean"),
                "best_type_neighbor": get_best(best_type, "neighbor_label"),
                "best_type_species": get_best(best_type, "neighbor_species"),
                "best_type_canonical": canonical_species(get_best(best_type, "neighbor_species")),
                "best_type_ani": get_best(best_type, "ani_mean"),
                "best_public_neighbor": get_best(best_public, "neighbor_label"),
                "best_public_species": get_best(best_public, "neighbor_species"),
                "best_public_canonical": canonical_species(get_best(best_public, "neighbor_species")),
                "best_public_ani": get_best(best_public, "ani_mean"),
                "type_and_public_species_match": canonical_species(get_best(best_type, "neighbor_species")) == canonical_species(get_best(best_public, "neighbor_species")),
                "strict_species_wording_safe": meta.loc[node, "strict_species_claim_allowed"] in {True, "True", "true"},
            }
        )
    summary = pd.DataFrame(summary_rows)
    def category(row: pd.Series) -> tuple[str, bool, str]:
        if row["node_role"] == "local_qc_warning":
            return (
                "QC_warning_excluded_despite_integrated_ANI_support",
                False,
                "Keep as species-supported QC warning; do not include in downstream trees or AMR interpretation.",
            )
        analysis = row["analysis_group_canonical"]
        type_sp = row["best_type_canonical"]
        public_sp = row["best_public_canonical"]
        if analysis == type_sp == public_sp:
            return (
                "integrated_type_public_ANI_supports_analysis_clade",
                True,
                "Strict species-clade wording is supported by integrated ANI, while retaining usual MAC database caveats.",
            )
        if row["type_and_public_species_match"] and analysis != type_sp:
            return (
                "integrated_ANI_supports_neighbor_clade_not_current_analysis_label",
                False,
                "Describe as MAC boundary/near-neighbor clade and review final analysis label before any unqualified species wording.",
            )
        return (
            "integrated_type_public_ANI_conflict_MAC_boundary",
            False,
            "Use public-context clade wording; explicitly discuss MAC type/public/database granularity conflict.",
        )
    cats = summary.apply(category, axis=1, result_type="expand")
    cats.columns = ["integrated_support_category", "integrated_strict_species_wording_safe", "integrated_wording_recommendation"]
    summary = pd.concat([summary, cats], axis=1)
    return pd.DataFrame(rows), summary


def plot_network(metadata: pd.DataFrame, pairs: pd.DataFrame, threshold: float, af_threshold: float, out_png: Path, out_pdf: Path) -> None:
    edges = pairs[(pairs["ani_mean"] >= threshold) & (pairs["alignment_fraction_mean"] >= af_threshold)].copy()
    if nx is None:
        raise RuntimeError("networkx is required for network plotting")
    graph = nx.Graph()
    meta = metadata.set_index("node_id")
    for node in metadata["node_id"]:
        graph.add_node(node)
    for _, row in edges.iterrows():
        graph.add_edge(row["node_a"], row["node_b"], weight=max(0.1, row["ani_mean"] - threshold + 0.2))
    pos = nx.spring_layout(graph, seed=42, k=0.72, iterations=250, weight="weight")

    fig, ax = plt.subplots(figsize=(13, 10))
    for lo, hi, color, width, alpha in [
        (98, 101, "#2B6CB0", 1.2, 0.25),
        (96, 98, "#7F7F7F", 0.7, 0.16),
    ]:
        sub_edges = [(u, v) for u, v, d in graph.edges(data=True) if lo <= pairs_key_ani(pairs, u, v) < hi]
        if sub_edges:
            nx.draw_networkx_edges(graph, pos, edgelist=sub_edges, edge_color=color, width=width, alpha=alpha, ax=ax)

    for role, marker in ROLE_MARKERS.items():
        nodes = metadata[metadata["node_role"].eq(role)]["node_id"].tolist()
        if not nodes:
            continue
        colors = [species_color(meta.loc[n, "species_short"]) for n in nodes]
        sizes = [480 if role == "local_downstream" else 520 if role == "local_qc_warning" else 165 if role == "type_or_representative" else 95 for _ in nodes]
        nx.draw_networkx_nodes(
            graph,
            pos,
            nodelist=nodes,
            node_color=colors,
            node_shape=marker,
            node_size=sizes,
            edgecolors="#222222",
            linewidths=0.7,
            alpha=0.95,
            ax=ax,
        )

    # Dense MAC near-neighbor graphs become unreadable when every local/type
    # anchor is labeled. Keep this panel as an overview and put identities in
    # the local type/public summary panel and TSV tables.

    ax.set_title(
        f"Integrated MAC/SGM ANI network overview: local priority isolates, type anchors, and public near-neighbors\nEdges shown at ANI >= {threshold:.0f}% and alignment fraction >= {af_threshold:.2f}; node labels omitted to avoid crowding",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.axis("off")
    legend_handles = [
        plt.Line2D([0], [0], marker="*", color="w", markerfacecolor="#FFFFFF", markeredgecolor="#222222", markersize=12, label="Local downstream"),
        plt.Line2D([0], [0], marker="X", color="w", markerfacecolor="#FFFFFF", markeredgecolor="#222222", markersize=10, label="Local QC warning"),
        plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="#FFFFFF", markeredgecolor="#222222", markersize=8, label="Type/representative anchor"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#FFFFFF", markeredgecolor="#222222", markersize=7, label="Public near-neighbor"),
        plt.Line2D([0], [0], color="#2B6CB0", lw=2, label="ANI >=98%"),
        plt.Line2D([0], [0], color="#7F7F7F", lw=2, label="ANI 96-98%"),
    ]
    ax.legend(handles=legend_handles, loc="lower left", fontsize=8, frameon=False)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def pairs_key_ani(pairs: pd.DataFrame, u: str, v: str) -> float:
    a, b = sorted([u, v])
    row = pairs[(pairs["node_a"].eq(a)) & (pairs["node_b"].eq(b))]
    if row.empty:
        return np.nan
    return float(row.iloc[0]["ani_mean"])


def plot_local_heatmap_and_tree(metadata: pd.DataFrame, ani_mat: pd.DataFrame, local_summary: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    meta = metadata.set_index("node_id")
    # Focus heatmap on local isolates plus their best type and public neighbors.
    focus_nodes: set[str] = set(metadata[metadata["node_role"].str.startswith("local")]["node_id"])
    for _, row in local_summary.iterrows():
        for col in ["best_type_neighbor", "best_public_neighbor"]:
            label = row[col]
            matches = metadata[metadata["display_label"].eq(label)]["node_id"].tolist()
            focus_nodes.update(matches[:1])
    focus = [n for n in metadata["node_id"] if n in focus_nodes]
    sub = ani_mat.loc[focus, focus].copy()
    # Fill distant/missing comparisons conservatively for clustering only.
    dist = 1 - (sub.fillna(80.0) / 100.0)
    np.fill_diagonal(dist.values, 0)
    condensed = squareform(dist.values, checks=False)
    z = linkage(condensed, method="average")

    fig = plt.figure(figsize=(18, max(10, 0.30 * len(focus) + 4)))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.75], wspace=0.42)
    ax_tree = fig.add_subplot(gs[0, 0])
    dendro = dendrogram(
        z,
        orientation="left",
        labels=[meta.loc[n, "display_label"] for n in focus],
        ax=ax_tree,
        color_threshold=0.04,
        leaf_font_size=7,
    )
    ax_tree.set_xlabel("1 - ANI/100 (average-linkage; missing filled as 80% ANI)")
    ax_tree.set_title("A  ANI clustering of local isolates and nearest anchors", loc="left", fontweight="bold")

    ordered_labels = dendro["ivl"]
    label_to_node = dict(zip([meta.loc[n, "display_label"] for n in focus], focus))
    ordered_nodes = [label_to_node[l] for l in ordered_labels]
    heat = sub.loc[ordered_nodes, ordered_nodes]
    ax_heat = fig.add_subplot(gs[0, 1])
    im = ax_heat.imshow(heat.values, vmin=94, vmax=100, cmap="viridis")
    ax_heat.set_xticks(np.arange(len(ordered_nodes)))
    ax_heat.set_yticks([])
    ax_heat.set_xticklabels([meta.loc[n, "display_label"] for n in ordered_nodes], rotation=90, fontsize=6)
    ax_heat.set_title("B  Pairwise FastANI among focused anchors", loc="left", fontweight="bold")
    for i, n in enumerate(ordered_nodes):
        role = meta.loc[n, "node_role"]
        color = "#000000" if role.startswith("local") else species_color(meta.loc[n, "species_short"])
        ax_heat.add_patch(patches.Rectangle((-0.5, i - 0.5), 0.25, 1, facecolor=color, clip_on=False))
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.02)
    cbar.set_label("ANI (%)")
    fig.suptitle("Integrated ANI tree/heatmap for MAC species-granularity review", x=0.02, ha="left", fontsize=13, fontweight="bold")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def plot_local_type_public_summary(local_summary: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    data = local_summary.copy()
    order_map = {
        "integrated_type_public_ANI_supports_analysis_clade": 0,
        "integrated_type_public_ANI_conflict_MAC_boundary": 1,
        "integrated_ANI_supports_neighbor_clade_not_current_analysis_label": 2,
        "QC_warning_excluded_despite_integrated_ANI_support": 3,
    }
    data["order"] = data["integrated_support_category"].map(order_map).fillna(9)
    data = data.sort_values(["order", "analysis_group_canonical", "sample_id"]).reset_index(drop=True)
    y = np.arange(len(data))

    category_colors = {
        "integrated_type_public_ANI_supports_analysis_clade": "#DDECCB",
        "integrated_type_public_ANI_conflict_MAC_boundary": "#F8E3A3",
        "integrated_ANI_supports_neighbor_clade_not_current_analysis_label": "#F2C6B6",
        "QC_warning_excluded_despite_integrated_ANI_support": "#E0E0E0",
    }
    fig, ax = plt.subplots(figsize=(13.8, 7.4))
    category_labels = {
        "integrated_type_public_ANI_supports_analysis_clade": "type+public support",
        "integrated_type_public_ANI_conflict_MAC_boundary": "MAC boundary conflict",
        "integrated_ANI_supports_neighbor_clade_not_current_analysis_label": "neighbor-clade warning",
        "QC_warning_excluded_despite_integrated_ANI_support": "QC excluded",
    }
    for i, row in data.iterrows():
        ax.axhspan(i - 0.46, i + 0.46, color=category_colors.get(row["integrated_support_category"], "#FFFFFF"), alpha=0.95)
        type_ani = pd.to_numeric(row["best_type_ani"], errors="coerce")
        pub_ani = pd.to_numeric(row["best_public_ani"], errors="coerce")
        if pd.notna(type_ani) and pd.notna(pub_ani):
            ax.plot([type_ani, pub_ani], [i, i], color="#555555", lw=1.0, alpha=0.65)
        if pd.notna(type_ani):
            ax.scatter(type_ani, i, marker="^", s=95, color=species_color(short_species(row["best_type_canonical"])), edgecolor="#222222", zorder=3, label="Best type anchor" if i == 0 else None)
        if pd.notna(pub_ani):
            ax.scatter(pub_ani, i, marker="o", s=80, color=species_color(short_species(row["best_public_canonical"])), edgecolor="#222222", zorder=3, label="Best public neighbor" if i == 0 else None)
        category = category_labels.get(row["integrated_support_category"], row["integrated_support_category"])
        ax.text(96.16, i, category, ha="left", va="center", fontsize=7.2, color="#333333")
        type_label = f"type: {short_species(row['best_type_canonical'])}"
        public_label = f"public: {short_species(row['best_public_canonical'])}"
        if pd.notna(type_ani):
            type_label += f" ({type_ani:.2f})"
        if pd.notna(pub_ani):
            public_label += f" ({pub_ani:.2f})"
        ax.text(100.18, i - 0.14, type_label, ha="left", va="center", fontsize=7.1, color="#333333")
        ax.text(100.18, i + 0.14, public_label, ha="left", va="center", fontsize=7.1, color="#333333")

    ax.axvline(95, color="#555555", ls="--", lw=1, alpha=0.75)
    ax.axvline(96, color="#555555", ls=":", lw=1, alpha=0.75)
    ax.set_xlim(95.0, 101.58)
    ax.set_ylim(len(data) - 0.5, -0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.sample_id}  ({short_species(r.analysis_group_canonical)})" for r in data.itertuples()], fontsize=8)
    ax.set_xticks([95, 96, 97, 98, 99, 100])
    ax.set_xlabel("FastANI to best integrated-panel anchor (%)")
    ax.set_title(
        "Integrated local ANI support: best type anchor versus best public near-neighbor",
        loc="left",
        fontsize=12,
        fontweight="bold",
        pad=10,
    )
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E5E5E5", lw=0.6)
    fig.subplots_adjust(bottom=0.16)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def short_species(value: object) -> str:
    text = canonical_species(value)
    return text.replace("Mycobacterium ", "M. ")


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    metadata = read_tsv(args.panel_metadata)
    all_hits = parse_fastani(args.fastani, metadata)
    pairs = symmetric_pairs(all_hits, metadata)
    ani_mat = build_ani_matrix(pairs, metadata)

    all_hits.to_csv(args.outdir / "integrated_panel_fastani_all_hits_annotated.tsv", sep="\t", index=False)
    pairs.to_csv(args.outdir / "integrated_panel_fastani_symmetric_pairs.tsv", sep="\t", index=False)
    ani_mat.to_csv(args.outdir / "integrated_panel_ani_matrix.tsv", sep="\t")

    edges_95 = pairs[pairs["ani_mean"] >= 95.0].copy()
    edges_96 = pairs[pairs["ani_mean"] >= 96.0].copy()
    edges_98 = pairs[pairs["ani_mean"] >= 98.0].copy()
    edges_95.to_csv(args.outdir / "integrated_panel_edges_ani95.tsv", sep="\t", index=False)
    edges_96.to_csv(args.outdir / "integrated_panel_edges_ani96.tsv", sep="\t", index=False)
    edges_98.to_csv(args.outdir / "integrated_panel_edges_ani98.tsv", sep="\t", index=False)

    local_neighbors, local_summary = local_neighbor_tables(metadata, pairs, args.top_neighbors)
    comp96 = connected_components(metadata["node_id"], pairs[pairs["alignment_fraction_mean"] >= args.network_af_threshold], 96.0)
    comp98 = connected_components(metadata["node_id"], pairs[pairs["alignment_fraction_mean"] >= args.network_af_threshold], 98.0)
    metadata2 = metadata.copy()
    metadata2["ani96_component"] = metadata2["node_id"].map(comp96)
    metadata2["ani98_component"] = metadata2["node_id"].map(comp98)
    metadata2.to_csv(args.outdir / "integrated_ani_panel_metadata_with_components.tsv", sep="\t", index=False)
    local_neighbors.to_csv(args.outdir / "local_priority14_integrated_top_neighbors.tsv", sep="\t", index=False)
    local_summary.to_csv(args.outdir / "local_priority14_integrated_neighbor_summary.tsv", sep="\t", index=False)

    plot_network(
        metadata2,
        pairs,
        args.network_ani_threshold,
        args.network_af_threshold,
        args.outdir / "integrated_ani_network.png",
        args.outdir / "integrated_ani_network.pdf",
    )
    plot_local_heatmap_and_tree(
        metadata2,
        ani_mat,
        local_summary,
        args.outdir / "integrated_ani_tree_heatmap.png",
        args.outdir / "integrated_ani_tree_heatmap.pdf",
    )
    plot_local_type_public_summary(
        local_summary,
        args.outdir / "integrated_local_type_public_ani_summary.png",
        args.outdir / "integrated_local_type_public_ani_summary.pdf",
    )

    with (args.outdir / "integrated_ani_taxonomy_notes.md").open("w") as handle:
        handle.write(
            "# Integrated ANI Taxonomy Notes\n\n"
            "This analysis combines local priority isolates, MAC/SGM type-strain or representative anchors, "
            "and per-isolate public near-neighbors. FastANI pairwise values are used for taxonomy-context "
            "support, not transmission inference. Network edges are shown at ANI >=96% with alignment "
            f"fraction >= {args.network_af_threshold:.2f}; edge tables are also exported at 95%, 96%, and 98%.\n\n"
            "Interpret exact MAC species names conservatively when type-panel, public-neighbor, marker-gene, "
            "or NTM-Profiler evidence disagree. This panel is intended to support public-context clade wording "
            "and highlight where unqualified species wording requires caution.\n\n"
            "Visualization note: the integrated network deliberately omits node labels because dense MAC "
            "near-neighbor clusters become unreadable when all local, public, and type-anchor labels are shown. "
            "Use the local type/public ANI support panel and exported TSV tables for exact sample/reference "
            "identities. The focused tree/heatmap is enlarged and suppresses duplicate heatmap row labels to "
            "avoid label collisions.\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
