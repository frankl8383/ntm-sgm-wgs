#!/usr/bin/env python3
"""Group threshold-stable accessory families by co-occurrence and gene adjacency."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


GENE_PATTERN = re.compile(r"__gene(\d+)$")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotated-stable", required=True)
    parser.add_argument("--primary-clusters", required=True)
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--panel-manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--max-intervening-genes", type=int, default=5)
    parser.add_argument("--minimum-pair-genomes", type=int, default=10)
    parser.add_argument("--minimum-presence-jaccard", type=float, default=0.80)
    parser.add_argument("--minimum-adjacency-fraction", type=float, default=0.80)
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.annotated_stable, sep="\t", dtype=str).fillna("")
    candidate_ids = set(candidates.primary_cluster_id)
    clusters = pd.read_csv(
        args.primary_clusters,
        sep="\t",
        names=["cluster_id", "sequence_id"],
        dtype=str,
    )
    clusters = clusters[clusters.cluster_id.isin(candidate_ids)].copy()
    metadata = pd.read_csv(args.sequence_metadata, sep="\t", dtype=str).fillna("")
    members = clusters.merge(
        metadata[["sequence_id", "tree_id"]],
        on="sequence_id",
        how="left",
        validate="many_to_one",
    )
    members["gene_index"] = members.sequence_id.map(
        lambda value: int(GENE_PATTERN.search(value).group(1))
        if GENE_PATTERN.search(value)
        else -1
    )
    panel = pd.read_csv(args.panel_manifest, sep="\t", dtype=str).fillna("")
    panel_meta = panel.set_index("tree_id")
    candidate_meta = candidates.set_index("primary_cluster_id")

    genomes_by_cluster = (
        members.groupby("cluster_id").tree_id.apply(lambda values: frozenset(values)).to_dict()
    )
    positions_by_cluster: dict[str, dict[str, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in members.itertuples(index=False):
        if int(row.gene_index) >= 0:
            positions_by_cluster[row.cluster_id][row.tree_id].append(int(row.gene_index))
    signature_groups: dict[tuple[str, frozenset[str]], list[str]] = defaultdict(list)
    for cluster_id, genomes in genomes_by_cluster.items():
        direction = candidate_meta.loc[cluster_id, "association_direction"]
        signature_groups[(direction, genomes)].append(cluster_id)

    module_rows: list[dict[str, object]] = []
    module_member_rows: list[dict[str, object]] = []
    direction_index = defaultdict(int)
    for (direction, genomes), cluster_ids in sorted(
        signature_groups.items(), key=lambda item: (-len(item[1]), item[0][0])
    ):
        if len(cluster_ids) < 2:
            continue
        direction_index[direction] += 1
        prefix = "MP" if direction == "MP_MIP_enriched" else "TMI"
        module_id = f"{prefix}_module_{direction_index[direction]:02d}"
        positions: dict[str, list[int]] = defaultdict(list)
        for cluster_id in cluster_ids:
            subset = members[members.cluster_id.eq(cluster_id)]
            for tree_id, frame in subset.groupby("tree_id"):
                valid = [value for value in frame.gene_index.astype(int) if value >= 0]
                if valid:
                    positions[tree_id].append(min(valid))
            module_member_rows.append(
                {
                    "module_id": module_id,
                    "association_direction": direction,
                    "cluster_id": cluster_id,
                    "conservative_product": candidate_meta.loc[
                        cluster_id, "conservative_product"
                    ],
                    "representative_aa_length": candidate_meta.loc[
                        cluster_id, "representative_aa_length"
                    ],
                }
            )

        adjacent = 0
        complete_with_positions = 0
        max_span = 0
        for tree_id in genomes:
            gene_positions = positions.get(tree_id, [])
            if len(gene_positions) != len(cluster_ids):
                continue
            complete_with_positions += 1
            span = max(gene_positions) - min(gene_positions) + 1
            max_span = max(max_span, span)
            if span <= len(cluster_ids) + args.max_intervening_genes:
                adjacent += 1

        public_tmi = sum(
            panel_meta.loc[tree_id, "source"] == "public"
            and panel_meta.loc[tree_id, "accessory_lineage"] == "MI_TMI_lineage"
            for tree_id in genomes
        )
        public_mp = sum(
            panel_meta.loc[tree_id, "source"] == "public"
            and panel_meta.loc[tree_id, "accessory_lineage"] == "MI_MP_MIP_lineage"
            for tree_id in genomes
        )
        local_tmi = sum(
            panel_meta.loc[tree_id, "source"] == "local"
            and panel_meta.loc[tree_id, "accessory_lineage"] == "MI_TMI_lineage"
            for tree_id in genomes
        )
        local_mp = sum(
            panel_meta.loc[tree_id, "source"] == "local"
            and panel_meta.loc[tree_id, "accessory_lineage"] == "MI_MP_MIP_lineage"
            for tree_id in genomes
        )
        products = [
            candidate_meta.loc[cluster_id, "conservative_product"]
            for cluster_id in cluster_ids
        ]
        module_rows.append(
            {
                "module_id": module_id,
                "association_direction": direction,
                "protein_families": len(cluster_ids),
                "genomes_with_complete_cooccurrence": len(genomes),
                "public_tmi_present": public_tmi,
                "public_mp_present": public_mp,
                "local_tmi_present": local_tmi,
                "local_mp_present": local_mp,
                "genomes_with_all_gene_positions": complete_with_positions,
                "genomes_with_adjacent_module": adjacent,
                "adjacent_fraction": (
                    adjacent / complete_with_positions if complete_with_positions else 0.0
                ),
                "maximum_observed_gene_span": max_span,
                "products": "; ".join(products),
                "interpretation": "co-occurring lineage-associated protein-family module; adjacency is based on uniform Prodigal gene order",
            }
        )

    module_frame = pd.DataFrame(module_rows)
    if not module_frame.empty:
        module_frame = module_frame.sort_values(
            ["protein_families", "adjacent_fraction"], ascending=[False, False]
        )
        module_frame.to_csv(output_dir / "stable_accessory_modules.tsv", sep="\t", index=False)
        pd.DataFrame(module_member_rows).to_csv(
            output_dir / "stable_accessory_module_members.tsv", sep="\t", index=False
        )
    else:
        pd.DataFrame(columns=["module_id"]).to_csv(
            output_dir / "stable_accessory_modules.tsv", sep="\t", index=False
        )
        pd.DataFrame(columns=["module_id"]).to_csv(
            output_dir / "stable_accessory_module_members.tsv", sep="\t", index=False
        )

    syntenic = (
        module_frame[
            (module_frame.protein_families >= 3) & (module_frame.adjacent_fraction >= 0.80)
        ]
        if not module_frame.empty
        else module_frame
    )
    if not syntenic.empty:
        syntenic.to_csv(output_dir / "stable_syntenic_accessory_modules.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=module_frame.columns).to_csv(
            output_dir / "stable_syntenic_accessory_modules.tsv", sep="\t", index=False
        )

    candidate_list = sorted(candidate_ids)
    graph: dict[str, set[str]] = {cluster_id: set() for cluster_id in candidate_list}
    edge_rows: list[dict[str, object]] = []
    for left_index, left in enumerate(candidate_list):
        left_direction = candidate_meta.loc[left, "association_direction"]
        left_genomes = genomes_by_cluster[left]
        for right in candidate_list[left_index + 1 :]:
            if candidate_meta.loc[right, "association_direction"] != left_direction:
                continue
            right_genomes = genomes_by_cluster[right]
            intersection = left_genomes & right_genomes
            union = left_genomes | right_genomes
            if len(intersection) < args.minimum_pair_genomes:
                continue
            jaccard = len(intersection) / len(union)
            if jaccard < args.minimum_presence_jaccard:
                continue
            valid = 0
            adjacent_pair = 0
            for tree_id in intersection:
                left_positions = positions_by_cluster[left].get(tree_id, [])
                right_positions = positions_by_cluster[right].get(tree_id, [])
                if not left_positions or not right_positions:
                    continue
                valid += 1
                minimum_distance = min(
                    abs(left_position - right_position)
                    for left_position in left_positions
                    for right_position in right_positions
                )
                if minimum_distance <= args.max_intervening_genes + 1:
                    adjacent_pair += 1
            adjacency_fraction = adjacent_pair / valid if valid else 0.0
            if adjacency_fraction < args.minimum_adjacency_fraction:
                continue
            graph[left].add(right)
            graph[right].add(left)
            edge_rows.append(
                {
                    "cluster_a": left,
                    "cluster_b": right,
                    "association_direction": left_direction,
                    "co_present_genomes": len(intersection),
                    "presence_jaccard": jaccard,
                    "genomes_with_gene_positions": valid,
                    "adjacent_genomes": adjacent_pair,
                    "adjacency_fraction": adjacency_fraction,
                }
            )

    visited: set[str] = set()
    components: list[list[str]] = []
    for start in candidate_list:
        if start in visited or not graph[start]:
            continue
        stack = [start]
        component: list[str] = []
        visited.add(start)
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbour in graph[node]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    stack.append(neighbour)
        if len(component) >= 2:
            components.append(sorted(component))

    component_rows: list[dict[str, object]] = []
    component_member_rows: list[dict[str, object]] = []
    component_direction_index = defaultdict(int)
    for component in sorted(components, key=lambda values: (-len(values), values[0])):
        direction = candidate_meta.loc[component[0], "association_direction"]
        component_direction_index[direction] += 1
        prefix = "MP" if direction == "MP_MIP_enriched" else "TMI"
        component_id = f"{prefix}_syntenic_component_{component_direction_index[direction]:02d}"
        complete_genomes = set.intersection(
            *(set(genomes_by_cluster[cluster_id]) for cluster_id in component)
        )
        component_edges = [
            row
            for row in edge_rows
            if row["cluster_a"] in component and row["cluster_b"] in component
        ]
        products = [
            candidate_meta.loc[cluster_id, "conservative_product"]
            for cluster_id in component
        ]
        component_rows.append(
            {
                "component_id": component_id,
                "association_direction": direction,
                "protein_families": len(component),
                "complete_module_genomes": len(complete_genomes),
                "supporting_adjacency_edges": len(component_edges),
                "minimum_edge_presence_jaccard": min(
                    float(row["presence_jaccard"]) for row in component_edges
                ),
                "minimum_edge_adjacency_fraction": min(
                    float(row["adjacency_fraction"]) for row in component_edges
                ),
                "products": "; ".join(products),
                "interpretation": "syntenic lineage-associated protein-family component defined by pairwise co-presence and conserved gene-order adjacency",
            }
        )
        for cluster_id in component:
            component_member_rows.append(
                {
                    "component_id": component_id,
                    "association_direction": direction,
                    "cluster_id": cluster_id,
                    "conservative_product": candidate_meta.loc[
                        cluster_id, "conservative_product"
                    ],
                    "representative_aa_length": candidate_meta.loc[
                        cluster_id, "representative_aa_length"
                    ],
                }
            )
    pd.DataFrame(edge_rows).to_csv(
        output_dir / "stable_syntenic_component_edges.tsv", sep="\t", index=False
    )
    pd.DataFrame(component_rows).to_csv(
        output_dir / "stable_syntenic_components.tsv", sep="\t", index=False
    )
    pd.DataFrame(component_member_rows).to_csv(
        output_dir / "stable_syntenic_component_members.tsv", sep="\t", index=False
    )
    report = [
        "# Stable accessory co-occurrence modules",
        "",
        f"Threshold-stable protein families: {len(candidate_ids)}",
        f"Exact co-occurrence modules with at least two families: {len(module_rows)}",
        f"Syntenic modules with at least three families and at least 80% adjacency: {syntenic.shape[0]}",
        f"Pairwise adjacency-supported syntenic components: {len(component_rows)}",
        f"Largest adjacency-supported component: {max((len(component) for component in components), default=0)} families",
        "",
        "Modules are descriptive lineage-associated gene-content blocks. They do not establish phenotype, selective advantage, horizontal transfer or diagnostic performance.",
    ]
    (output_dir / "stable_accessory_module_report.md").write_text(
        "\n".join(report) + "\n"
    )
    print(
        f"stable_families={len(candidate_ids)} modules={len(module_rows)} "
        f"syntenic_modules={syntenic.shape[0]}"
    )


if __name__ == "__main__":
    main()
