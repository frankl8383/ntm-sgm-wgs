#!/usr/bin/env python3
"""Discover public MI lineage-associated clusters and validate them locally."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


TMI = "MI_TMI_lineage"
MP = "MI_MP_MIP_lineage"


def benjamini_hochberg(values: np.ndarray) -> np.ndarray:
    n = len(values)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for index in range(n - 1, -1, -1):
        rank = index + 1
        running = min(running, ranked[index] * n / rank)
        adjusted[index] = running
    output = np.empty(n, dtype=float)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def parse_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    header = ""
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if header:
            sequences[header] = "".join(parts)
    return sequences


def write_fasta(path: Path, ids: list[str], sequences: dict[str, str]) -> None:
    with path.open("w") as handle:
        for sequence_id in ids:
            sequence = sequences.get(sequence_id)
            if not sequence:
                continue
            handle.write(f">{sequence_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-tsv", required=True)
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--panel-manifest", required=True)
    parser.add_argument("--representative-fasta", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--run-label", default="primary70_cov80")
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.sequence_metadata, sep="\t", dtype=str).fillna("")
    panel = pd.read_csv(args.panel_manifest, sep="\t", dtype=str).fillna("")
    clusters = pd.read_csv(
        args.cluster_tsv,
        sep="\t",
        names=["cluster_id", "sequence_id"],
        dtype=str,
    )
    merged = clusters.merge(
        metadata[
            [
                "sequence_id",
                "tree_id",
                "source",
                "accessory_lineage",
                "cohort_source",
                "aa_length",
            ]
        ],
        on="sequence_id",
        how="left",
        validate="many_to_one",
    )
    if merged["tree_id"].isna().any():
        raise SystemExit("MMseqs members missing from sequence metadata")
    presence = merged[["cluster_id", "tree_id"]].drop_duplicates()

    panel_by_tree = panel.set_index("tree_id", drop=False)
    groups: dict[str, set[str]] = {
        "public_tmi": set(
            panel[(panel.source == "public") & (panel.accessory_lineage == TMI)].tree_id
        ),
        "public_mp": set(
            panel[(panel.source == "public") & (panel.accessory_lineage == MP)].tree_id
        ),
        "local_tmi": set(
            panel[(panel.source == "local") & (panel.accessory_lineage == TMI)].tree_id
        ),
        "local_mp": set(
            panel[(panel.source == "local") & (panel.accessory_lineage == MP)].tree_id
        ),
    }
    groups["local_tmi_direct"] = {
        tree
        for tree in groups["local_tmi"]
        if panel_by_tree.loc[tree, "cohort_source"] == "directly_retained_original"
    }
    groups["local_tmi_rescued"] = groups["local_tmi"] - groups["local_tmi_direct"]
    groups["local_mp_direct"] = {
        tree
        for tree in groups["local_mp"]
        if panel_by_tree.loc[tree, "cohort_source"] == "directly_retained_original"
    }
    groups["local_mp_rescued"] = groups["local_mp"] - groups["local_mp_direct"]
    for required in ("public_tmi", "public_mp", "local_tmi", "local_mp"):
        if not groups[required]:
            raise SystemExit(f"Empty required group: {required}")

    members_by_cluster = presence.groupby("cluster_id").tree_id.apply(set)
    rep_lengths = (
        merged.drop_duplicates("cluster_id")
        .set_index("cluster_id")["aa_length"]
        .astype(int)
        .to_dict()
    )
    rows: list[dict[str, object]] = []
    p_values: list[float] = []
    for cluster_id, genomes in members_by_cluster.items():
        counts = {name: len(genomes & members) for name, members in groups.items()}
        pub_tmi_total = len(groups["public_tmi"])
        pub_mp_total = len(groups["public_mp"])
        odds, p_value = fisher_exact(
            [
                [counts["public_mp"], pub_mp_total - counts["public_mp"]],
                [counts["public_tmi"], pub_tmi_total - counts["public_tmi"]],
            ],
            alternative="two-sided",
        )
        p_values.append(float(p_value))
        row: dict[str, object] = {
            "cluster_id": cluster_id,
            "representative_aa_length": rep_lengths.get(cluster_id, 0),
            "total_genomes_present": len(genomes),
            "public_tmi_present": counts["public_tmi"],
            "public_tmi_total": pub_tmi_total,
            "public_tmi_prevalence": counts["public_tmi"] / pub_tmi_total,
            "public_mp_present": counts["public_mp"],
            "public_mp_total": pub_mp_total,
            "public_mp_prevalence": counts["public_mp"] / pub_mp_total,
            "public_mp_minus_tmi_prevalence": (
                counts["public_mp"] / pub_mp_total - counts["public_tmi"] / pub_tmi_total
            ),
            "fisher_odds_ratio_mp_vs_tmi": odds,
            "fisher_p_value": p_value,
            "local_tmi_present": counts["local_tmi"],
            "local_tmi_total": len(groups["local_tmi"]),
            "local_tmi_prevalence": counts["local_tmi"] / len(groups["local_tmi"]),
            "local_mp_present": counts["local_mp"],
            "local_mp_total": len(groups["local_mp"]),
            "local_mp_prevalence": counts["local_mp"] / len(groups["local_mp"]),
        }
        for name in (
            "local_tmi_direct",
            "local_tmi_rescued",
            "local_mp_direct",
            "local_mp_rescued",
        ):
            total = len(groups[name])
            row[f"{name}_present"] = counts[name]
            row[f"{name}_total"] = total
            row[f"{name}_prevalence"] = counts[name] / total if total else np.nan
        rows.append(row)

    q_values = benjamini_hochberg(np.asarray(p_values, dtype=float))
    for row, q_value in zip(rows, q_values, strict=True):
        row["fisher_fdr"] = q_value
        difference = float(row["public_mp_minus_tmi_prevalence"])
        if difference >= 0.60:
            direction = "MP_MIP_enriched"
        elif difference <= -0.60:
            direction = "TMI_enriched"
        else:
            direction = "not_strongly_lineage_associated"
        row["public_association_direction"] = direction
        public_candidate = (
            float(row["fisher_fdr"]) <= 0.01
            and (
                (float(row["public_mp_prevalence"]) >= 0.80 and float(row["public_tmi_prevalence"]) <= 0.20)
                or (float(row["public_tmi_prevalence"]) >= 0.80 and float(row["public_mp_prevalence"]) <= 0.20)
            )
        )
        row["public_discovery_candidate"] = public_candidate
        if direction == "MP_MIP_enriched":
            directionally_concordant = (
                float(row["local_mp_prevalence"]) >= 0.75
                and float(row["local_tmi_prevalence"]) <= 0.25
            )
            route_consistent = (
                float(row["local_mp_direct_prevalence"]) >= 0.60
                and float(row["local_mp_rescued_prevalence"]) >= 0.60
                and float(row["local_tmi_direct_prevalence"]) <= 0.50
                and float(row["local_tmi_rescued_prevalence"]) <= 0.50
            )
        elif direction == "TMI_enriched":
            directionally_concordant = (
                float(row["local_tmi_prevalence"]) >= 0.75
                and float(row["local_mp_prevalence"]) <= 0.25
            )
            route_consistent = (
                float(row["local_tmi_direct_prevalence"]) >= 0.50
                and float(row["local_tmi_rescued_prevalence"]) >= 0.50
                and float(row["local_mp_direct_prevalence"]) <= 0.40
                and float(row["local_mp_rescued_prevalence"]) <= 0.40
            )
        else:
            directionally_concordant = False
            route_consistent = False
        row["local_directional_concordance"] = bool(
            public_candidate and directionally_concordant
        )
        row["local_route_consistent"] = bool(
            public_candidate and directionally_concordant and route_consistent
        )

    result = pd.DataFrame(rows)
    result = result.sort_values(
        ["public_discovery_candidate", "local_directional_concordance", "fisher_fdr"],
        ascending=[False, False, True],
    )
    result.to_csv(output_dir / "all_cluster_lineage_statistics.tsv", sep="\t", index=False)
    candidates = result[result.public_discovery_candidate].copy()
    candidates.to_csv(output_dir / "public_discovery_candidates.tsv", sep="\t", index=False)
    concordant = candidates[candidates.local_directional_concordance].copy()
    concordant.to_csv(
        output_dir / "public_candidates_with_local_directional_concordance.tsv",
        sep="\t",
        index=False,
    )
    robust = concordant[concordant.local_route_consistent].copy()
    robust.to_csv(output_dir / "public_candidates_route_robust.tsv", sep="\t", index=False)

    sequences = parse_fasta(Path(args.representative_fasta))
    write_fasta(
        output_dir / "public_candidates_with_local_directional_concordance.faa",
        list(concordant.cluster_id),
        sequences,
    )
    write_fasta(
        output_dir / "public_candidates_route_robust.faa",
        list(robust.cluster_id),
        sequences,
    )

    summary = [
        {
            "run_label": args.run_label,
            "n_genomes_total": panel.shape[0],
            "public_tmi_n": len(groups["public_tmi"]),
            "public_mp_mip_n": len(groups["public_mp"]),
            "local_tmi_n": len(groups["local_tmi"]),
            "local_mp_mip_n": len(groups["local_mp"]),
            "protein_clusters": result.shape[0],
            "public_discovery_candidates": candidates.shape[0],
            "local_directionally_concordant_candidates": concordant.shape[0],
            "route_robust_candidates": robust.shape[0],
            "interpretation": "lineage-associated protein families; not phenotype, virulence, transmission, HGT or diagnostic markers",
        }
    ]
    pd.DataFrame(summary).to_csv(output_dir / "accessory_analysis_summary.tsv", sep="\t", index=False)
    report = [
        "# MI lineage-associated accessory protein analysis",
        "",
        f"Run: {args.run_label}",
        f"Protein clusters: {result.shape[0]}",
        f"Public discovery candidates: {candidates.shape[0]}",
        f"Candidates with the same direction in local TMI versus MP-MIP genomes: {concordant.shape[0]}",
        f"Candidates also consistent across direct and rescued local assemblies: {robust.shape[0]}",
        "",
        "Candidates were discovered in public genomes and evaluated in local genomes. They are lineage-associated protein families, not phenotype, virulence, transmission, HGT or diagnostic-marker claims.",
    ]
    (output_dir / "accessory_analysis_report.md").write_text("\n".join(report) + "\n")
    print(
        f"clusters={result.shape[0]} public_candidates={candidates.shape[0]} "
        f"local_concordant={concordant.shape[0]} route_robust={robust.shape[0]}"
    )


if __name__ == "__main__":
    main()
