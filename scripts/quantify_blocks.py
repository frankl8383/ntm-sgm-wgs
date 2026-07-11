#!/usr/bin/env python3
"""Quantify selected stable MI lineage-associated blocks across panel strata."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


GENE_PATTERN = re.compile(r"__gene(\d+)$")
TMI = "MI_TMI_lineage"
MP = "MI_MP_MIP_lineage"


def benjamini_hochberg(values: np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted P values."""
    n = len(values)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for index in range(n - 1, -1, -1):
        running = min(running, ranked[index] * n / (index + 1))
        adjusted[index] = running
    output = np.empty(n, dtype=float)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotated-stable", required=True)
    parser.add_argument("--primary-clusters", required=True)
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--panel-manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stable = pd.read_csv(args.annotated_stable, sep="\t", dtype=str).fillna("")
    stable_ids = set(stable.primary_cluster_id)
    stable_meta = stable.set_index("primary_cluster_id")

    def in_gene_range(prefix: str, start: int, end: int) -> list[str]:
        selected: list[str] = []
        for cluster_id in stable_ids:
            if not cluster_id.startswith(prefix):
                continue
            match = GENE_PATTERN.search(cluster_id)
            if match and start <= int(match.group(1)) <= end:
                selected.append(cluster_id)
        return sorted(selected)

    definitions = {
        "TMI_aromatic_catabolism_associated_block": {
            "clusters": in_gene_range("LOCAL_Mi1__", 14, 39),
            "expected_direction": "TMI_enriched",
            "maximum_span": 35,
        },
        "TMI_methyltransferase_hydrolase_cupin_block": {
            "clusters": in_gene_range("LOCAL_Mi1__", 389, 393),
            "expected_direction": "TMI_enriched",
            "maximum_span": 10,
        },
        "MP_MIP_nitrogen_redox_associated_block": {
            "clusters": [
                "LOCAL_Ma20__gene01170",
                "LOCAL_Ma20__gene01177",
                "LOCAL_Ma20__gene01178",
            ],
            "expected_direction": "MP_MIP_enriched",
            "maximum_span": 15,
        },
        "MP_MIP_oxidoreductase_associated_block": {
            "clusters": in_gene_range("LOCAL_Ma20__", 2160, 2164),
            "expected_direction": "MP_MIP_enriched",
            "maximum_span": 10,
        },
        "MP_MIP_ArsN1_family_B_signal": {
            "clusters": ["LOCAL_Mi32__gene00261"],
            "expected_direction": "MP_MIP_enriched",
            "maximum_span": 1,
        },
    }
    for name, definition in definitions.items():
        missing = set(definition["clusters"]) - stable_ids
        if missing or not definition["clusters"]:
            raise SystemExit(f"Curated block {name} missing stable clusters: {sorted(missing)}")
        observed = {
            stable_meta.loc[cluster_id, "association_direction"]
            for cluster_id in definition["clusters"]
        }
        if observed != {definition["expected_direction"]}:
            raise SystemExit(f"Curated block {name} has mixed directions: {observed}")

    clusters = pd.read_csv(
        args.primary_clusters,
        sep="\t",
        names=["cluster_id", "sequence_id"],
        dtype=str,
    )
    selected_ids = {
        cluster_id
        for definition in definitions.values()
        for cluster_id in definition["clusters"]
    }
    clusters = clusters[clusters.cluster_id.isin(selected_ids)].copy()
    sequence_meta = pd.read_csv(args.sequence_metadata, sep="\t", dtype=str).fillna("")
    members = clusters.merge(
        sequence_meta[["sequence_id", "tree_id"]],
        on="sequence_id",
        how="left",
        validate="many_to_one",
    )
    members["gene_index"] = members.sequence_id.map(
        lambda value: int(GENE_PATTERN.search(value).group(1))
        if GENE_PATTERN.search(value)
        else -1
    )
    presence = (
        members.groupby("cluster_id").tree_id.apply(lambda values: set(values)).to_dict()
    )
    positions = (
        members[members.gene_index.ge(0)]
        .groupby(["cluster_id", "tree_id"])
        .gene_index.apply(lambda values: list(map(int, values)))
        .to_dict()
    )
    panel = pd.read_csv(args.panel_manifest, sep="\t", dtype=str).fillna("")
    groups = {
        "public_tmi": set(
            panel[(panel.source == "public") & (panel.accessory_lineage == TMI)].tree_id
        ),
        "public_mp_mip": set(
            panel[(panel.source == "public") & (panel.accessory_lineage == MP)].tree_id
        ),
        "local_tmi": set(
            panel[(panel.source == "local") & (panel.accessory_lineage == TMI)].tree_id
        ),
        "local_mp_mip": set(
            panel[(panel.source == "local") & (panel.accessory_lineage == MP)].tree_id
        ),
    }

    rows: list[dict[str, object]] = []
    member_rows: list[dict[str, object]] = []
    for name, definition in definitions.items():
        cluster_ids = list(definition["clusters"])
        complete = set.intersection(*(presence[cluster_id] for cluster_id in cluster_ids))
        adjacent = 0
        with_positions = 0
        observed_spans: list[int] = []
        for tree_id in complete:
            gene_positions: list[int] = []
            valid = True
            for cluster_id in cluster_ids:
                values = positions.get((cluster_id, tree_id), [])
                if not values:
                    valid = False
                    break
                gene_positions.append(min(values))
            if not valid:
                continue
            with_positions += 1
            span = max(gene_positions) - min(gene_positions) + 1
            observed_spans.append(span)
            if span <= int(definition["maximum_span"]):
                adjacent += 1
        row: dict[str, object] = {
            "block_id": name,
            "association_direction": definition["expected_direction"],
            "stable_protein_families": len(cluster_ids),
            "complete_block_genomes": len(complete),
            "genomes_with_gene_positions": with_positions,
            "genomes_with_expected_span": adjacent,
            "expected_span_fraction": adjacent / with_positions if with_positions else 0.0,
            "median_observed_gene_span": (
                float(pd.Series(observed_spans).median()) if observed_spans else 0.0
            ),
            "maximum_allowed_gene_span": definition["maximum_span"],
            "products": "; ".join(
                stable_meta.loc[cluster_id, "conservative_product"]
                for cluster_id in cluster_ids
            ),
            "interpretation": "stable lineage-associated block; functional names are conservative homology annotations, not phenotype or mechanism",
        }
        for group_name, genomes in groups.items():
            present = len(complete & genomes)
            row[f"{group_name}_complete"] = present
            row[f"{group_name}_total"] = len(genomes)
            row[f"{group_name}_prevalence"] = present / len(genomes) if genomes else 0.0
        rows.append(row)
        for cluster_id in cluster_ids:
            member_rows.append(
                {
                    "block_id": name,
                    "cluster_id": cluster_id,
                    "conservative_product": stable_meta.loc[
                        cluster_id, "conservative_product"
                    ],
                    "association_direction": stable_meta.loc[
                        cluster_id, "association_direction"
                    ],
                }
            )

    frame = pd.DataFrame(rows)
    p_values: list[float] = []
    odds_ratios: list[float] = []
    expected_differences: list[float] = []
    for row in rows:
        odds_ratio, p_value = fisher_exact(
            [
                [
                    int(row["public_mp_mip_complete"]),
                    int(row["public_mp_mip_total"])
                    - int(row["public_mp_mip_complete"]),
                ],
                [
                    int(row["public_tmi_complete"]),
                    int(row["public_tmi_total"])
                    - int(row["public_tmi_complete"]),
                ],
            ]
        )
        odds_ratios.append(float(odds_ratio))
        p_values.append(float(p_value))
        if row["association_direction"] == "TMI_enriched":
            expected_differences.append(
                float(row["public_tmi_prevalence"])
                - float(row["public_mp_mip_prevalence"])
            )
        else:
            expected_differences.append(
                float(row["public_mp_mip_prevalence"])
                - float(row["public_tmi_prevalence"])
            )
    frame["public_fisher_odds_ratio_mp_vs_tmi"] = odds_ratios
    frame["public_fisher_p_value"] = p_values
    frame["public_fisher_fdr"] = benjamini_hochberg(np.asarray(p_values))
    frame["public_prevalence_gap_expected_direction"] = expected_differences
    frame.to_csv(output_dir / "curated_accessory_blocks.tsv", sep="\t", index=False)
    pd.DataFrame(member_rows).to_csv(
        output_dir / "curated_accessory_block_members.tsv", sep="\t", index=False
    )
    report = [
        "# Curated MI lineage-associated accessory blocks",
        "",
        "| Block | Families | Public TMI | Public MP-MIP | Local TMI | Local MP-MIP | Expected-span fraction | Fisher FDR |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in frame.to_dict(orient="records"):
        report.append(
            "| {block_id} | {stable_protein_families} | "
            "{public_tmi_complete}/{public_tmi_total} | "
            "{public_mp_mip_complete}/{public_mp_mip_total} | "
            "{local_tmi_complete}/{local_tmi_total} | "
            "{local_mp_mip_complete}/{local_mp_mip_total} | "
            "{expected_span_fraction:.3f} | {public_fisher_fdr:.3g} |".format(**row)
        )
    report.extend(
        [
            "",
            "Blocks were selected from threshold-stable, public-discovered families and quantified without phenotype claims. ArsN1 family B is reported as a stress-associated family signal rather than a clinical resistance determinant.",
        ]
    )
    (output_dir / "curated_accessory_blocks_report.md").write_text(
        "\n".join(report) + "\n"
    )
    print(f"blocks={len(rows)} families={len(selected_ids)}")


if __name__ == "__main__":
    main()
