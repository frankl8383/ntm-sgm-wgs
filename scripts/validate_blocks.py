#!/usr/bin/env python3
"""Validate curated MI accessory blocks in the full current public panel."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable


TMI = "MI_TMI_lineage"
MP = "MI_MP_MIP_lineage"
GENE_SUFFIX = re.compile(r"_(\d+)$")


def benjamini_hochberg(values: np.ndarray) -> np.ndarray:
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


def prodigal_location(original_id: str) -> tuple[str, int]:
    """Split a Prodigal protein ID into contig and within-contig gene index."""
    match = GENE_SUFFIX.search(original_id)
    if not match:
        return original_id, -1
    return original_id[: match.start()], int(match.group(1))


def minimum_covering_span(
    hits: pd.DataFrame, required_families: set[str]
) -> tuple[str, int] | None:
    """Find the shortest same-contig gene window containing every family."""
    best: tuple[str, int] | None = None
    for contig, frame in hits.groupby("contig", sort=False):
        if set(frame["query"]) != required_families:
            continue
        points = sorted(
            {(int(row.gene_index), str(row.query)) for row in frame.itertuples(index=False)}
        )
        counts: defaultdict[str, int] = defaultdict(int)
        covered = 0
        left = 0
        for right, (right_pos, right_family) in enumerate(points):
            if counts[right_family] == 0:
                covered += 1
            counts[right_family] += 1
            while covered == len(required_families) and left <= right:
                left_pos, left_family = points[left]
                span = right_pos - left_pos + 1
                if best is None or span < best[1]:
                    best = (str(contig), int(span))
                counts[left_family] -= 1
                if counts[left_family] == 0:
                    covered -= 1
                left += 1
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--block-summary", required=True)
    parser.add_argument("--block-members", required=True)
    parser.add_argument("--search-tsv", required=True)
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--panel-manifest", required=True)
    parser.add_argument("--discovery-panel-manifest")
    parser.add_argument("--atlas-metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    blocks = pd.read_csv(args.block_summary, sep="\t", dtype=str).fillna("")
    members = pd.read_csv(args.block_members, sep="\t", dtype=str).fillna("")
    panel = pd.read_csv(args.panel_manifest, sep="\t", dtype=str).fillna("")
    sequence_meta = pd.read_csv(args.sequence_metadata, sep="\t", dtype=str).fillna("")
    atlas = pd.read_csv(args.atlas_metadata, sep="\t", dtype=str).fillna("")
    hits = pd.read_csv(
        args.search_tsv,
        sep="\t",
        names=["query", "target", "pident", "alnlen", "qcov", "tcov", "evalue", "bits"],
        dtype=str,
    )

    wanted = set(members.cluster_id)
    hits = hits[hits["query"].isin(wanted)].copy()
    hits = hits.merge(
        sequence_meta[["sequence_id", "tree_id", "original_prodigal_id"]],
        left_on="target",
        right_on="sequence_id",
        how="left",
        validate="many_to_one",
    )
    if hits.tree_id.isna().any():
        raise SystemExit("Search targets are missing sequence metadata")
    locations = hits.original_prodigal_id.map(prodigal_location)
    hits["contig"] = locations.map(lambda item: item[0])
    hits["gene_index"] = locations.map(lambda item: item[1])

    atlas_index = atlas.set_index("assembly_accession")
    public_project = atlas_index.bioproject.to_dict()
    panel["bioproject"] = panel.apply(
        lambda row: public_project.get(row.tree_id.removeprefix("PUB_"), "")
        if row.source == "public"
        else row.source,
        axis=1,
    )
    for column in [
        "checkm_completeness",
        "checkm_contamination",
        "number_of_contigs",
        "contig_n50",
    ]:
        values = pd.to_numeric(atlas_index[column], errors="coerce").to_dict()
        panel[column] = panel.apply(
            lambda row: values.get(row.tree_id.removeprefix("PUB_"), np.nan)
            if row.source == "public"
            else np.nan,
            axis=1,
        )
    panel["high_qc_public_sensitivity_include"] = (
        panel.source.eq("public")
        & panel.checkm_completeness.ge(95.0)
        & panel.checkm_contamination.le(2.0)
        & panel.number_of_contigs.le(200)
        & panel.contig_n50.ge(20000)
    )
    panel.to_csv(output_dir / "full_mi_panel_with_bioproject.tsv", sep="\t", index=False)

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

    per_genome_rows: list[dict[str, object]] = []
    block_presence: dict[str, set[str]] = {}
    block_family_presence: dict[str, set[str]] = {}
    for block in blocks.itertuples(index=False):
        block_id = str(block.block_id)
        required = set(members[members.block_id == block_id].cluster_id)
        maximum_span = int(float(block.maximum_allowed_gene_span))
        present: set[str] = set()
        family_complete_genomes: set[str] = set()
        for genome in panel.tree_id:
            genome_hits = hits[
                (hits.tree_id == genome) & (hits["query"].isin(required))
            ]
            observed = set(genome_hits["query"])
            complete = observed == required
            if complete:
                family_complete_genomes.add(genome)
            window = minimum_covering_span(genome_hits, required) if complete else None
            same_contig = window is not None
            minimum_span = window[1] if window else np.nan
            expected_span = bool(window is not None and window[1] <= maximum_span)
            if expected_span:
                present.add(genome)
            per_genome_rows.append(
                {
                    "block_id": block_id,
                    "association_direction": block.association_direction,
                    "tree_id": genome,
                    "source": panel.loc[panel.tree_id == genome, "source"].iat[0],
                    "accessory_lineage": panel.loc[
                        panel.tree_id == genome, "accessory_lineage"
                    ].iat[0],
                    "bioproject": panel.loc[panel.tree_id == genome, "bioproject"].iat[0],
                    "required_families": len(required),
                    "observed_families": len(observed),
                    "all_families_present": complete,
                    "all_families_on_one_contig": same_contig,
                    "minimum_gene_span": minimum_span,
                    "maximum_allowed_gene_span": maximum_span,
                    "syntenic_block_present": expected_span,
                }
            )
        block_presence[block_id] = present
        block_family_presence[block_id] = family_complete_genomes

    per_genome = pd.DataFrame(per_genome_rows)
    per_genome.to_csv(
        output_dir / "curated_blocks_full_atlas_per_genome.tsv", sep="\t", index=False
    )

    if args.discovery_panel_manifest:
        discovery_panel = pd.read_csv(
            args.discovery_panel_manifest, sep="\t", dtype=str
        ).fillna("")
        discovery_ids = set(
            discovery_panel[
                discovery_panel.source.eq("public")
                & discovery_panel.accessory_lineage.isin([TMI, MP])
            ].tree_id
        )
        evaluation_panel = panel[
            panel.source.eq("public")
            & panel.accessory_lineage.isin([TMI, MP])
            & ~panel.tree_id.isin(discovery_ids)
        ].copy()
        evaluation_panel[["tree_id", "accessory_lineage", "bioproject"]].to_csv(
            output_dir / "non_discovery_public_membership.tsv", sep="\t", index=False
        )
        evaluation_tmi = set(
            evaluation_panel[evaluation_panel.accessory_lineage.eq(TMI)].tree_id
        )
        evaluation_mp = set(
            evaluation_panel[evaluation_panel.accessory_lineage.eq(MP)].tree_id
        )
        evaluation_rows: list[dict[str, object]] = []
        evaluation_p_values: list[float] = []
        for block in blocks.itertuples(index=False):
            block_id = str(block.block_id)
            present = block_presence[block_id]
            tmi_present = len(present & evaluation_tmi)
            mp_present = len(present & evaluation_mp)
            if str(block.association_direction) == "TMI_enriched":
                table = [
                    [tmi_present, len(evaluation_tmi) - tmi_present],
                    [mp_present, len(evaluation_mp) - mp_present],
                ]
                expected_gap = (
                    tmi_present / len(evaluation_tmi)
                    - mp_present / len(evaluation_mp)
                )
            else:
                table = [
                    [mp_present, len(evaluation_mp) - mp_present],
                    [tmi_present, len(evaluation_tmi) - tmi_present],
                ]
                expected_gap = (
                    mp_present / len(evaluation_mp)
                    - tmi_present / len(evaluation_tmi)
                )
            odds_ratio, p_value = fisher_exact(table)
            evaluation_p_values.append(float(p_value))
            evaluation_rows.append(
                {
                    "block_id": block_id,
                    "association_direction": block.association_direction,
                    "tmi_present": tmi_present,
                    "tmi_total": len(evaluation_tmi),
                    "mp_mip_present": mp_present,
                    "mp_mip_total": len(evaluation_mp),
                    "expected_direction_prevalence_gap": expected_gap,
                    "expected_direction_odds_ratio": float(odds_ratio),
                    "fisher_p_value": float(p_value),
                }
            )
        evaluation = pd.DataFrame(evaluation_rows)
        evaluation["fisher_fdr"] = benjamini_hochberg(
            np.asarray(evaluation_p_values)
        )
        evaluation.to_csv(
            output_dir / "non_discovery_public_block_statistics.tsv",
            sep="\t",
            index=False,
        )

    summary_rows: list[dict[str, object]] = []
    p_values: list[float] = []
    high_qc_p_values: list[float] = []
    high_qc_public_tmi = set(
        panel[
            panel.high_qc_public_sensitivity_include
            & panel.accessory_lineage.eq(TMI)
        ].tree_id
    )
    high_qc_public_mp = set(
        panel[
            panel.high_qc_public_sensitivity_include
            & panel.accessory_lineage.eq(MP)
        ].tree_id
    )
    for block in blocks.itertuples(index=False):
        block_id = str(block.block_id)
        present = block_presence[block_id]
        family_complete = block_family_presence[block_id]
        counts = {name: len(present & genomes) for name, genomes in groups.items()}
        tmi_total = len(groups["public_tmi"])
        mp_total = len(groups["public_mp_mip"])
        odds, p_value = fisher_exact(
            [
                [counts["public_mp_mip"], mp_total - counts["public_mp_mip"]],
                [counts["public_tmi"], tmi_total - counts["public_tmi"]],
            ]
        )
        p_values.append(float(p_value))
        high_qc_tmi_present = len(present & high_qc_public_tmi)
        high_qc_mp_present = len(present & high_qc_public_mp)
        high_qc_odds, high_qc_p_value = fisher_exact(
            [
                [
                    high_qc_mp_present,
                    len(high_qc_public_mp) - high_qc_mp_present,
                ],
                [
                    high_qc_tmi_present,
                    len(high_qc_public_tmi) - high_qc_tmi_present,
                ],
            ]
        )
        high_qc_p_values.append(float(high_qc_p_value))
        direction = str(block.association_direction)
        tmi_prev = counts["public_tmi"] / tmi_total
        mp_prev = counts["public_mp_mip"] / mp_total
        expected_gap = tmi_prev - mp_prev if direction == "TMI_enriched" else mp_prev - tmi_prev
        summary_rows.append(
            {
                "block_id": block_id,
                "association_direction": direction,
                "stable_protein_families": int(block.stable_protein_families),
                "full_public_tmi_family_complete": len(
                    family_complete & groups["public_tmi"]
                ),
                "full_public_tmi_present": counts["public_tmi"],
                "full_public_tmi_total": tmi_total,
                "full_public_tmi_prevalence": tmi_prev,
                "full_public_mp_mip_present": counts["public_mp_mip"],
                "full_public_mp_mip_family_complete": len(
                    family_complete & groups["public_mp_mip"]
                ),
                "full_public_mp_mip_total": mp_total,
                "full_public_mp_mip_prevalence": mp_prev,
                "full_local_tmi_present": counts["local_tmi"],
                "full_local_tmi_family_complete": len(
                    family_complete & groups["local_tmi"]
                ),
                "full_local_tmi_total": len(groups["local_tmi"]),
                "full_local_mp_mip_present": counts["local_mp_mip"],
                "full_local_mp_mip_family_complete": len(
                    family_complete & groups["local_mp_mip"]
                ),
                "full_local_mp_mip_total": len(groups["local_mp_mip"]),
                "expected_direction_prevalence_gap": expected_gap,
                "full_public_fisher_odds_ratio_mp_vs_tmi": float(odds),
                "full_public_fisher_p_value": float(p_value),
                "high_qc_public_tmi_present": high_qc_tmi_present,
                "high_qc_public_tmi_total": len(high_qc_public_tmi),
                "high_qc_public_tmi_prevalence": high_qc_tmi_present
                / len(high_qc_public_tmi),
                "high_qc_public_mp_mip_present": high_qc_mp_present,
                "high_qc_public_mp_mip_total": len(high_qc_public_mp),
                "high_qc_public_mp_mip_prevalence": high_qc_mp_present
                / len(high_qc_public_mp),
                "high_qc_fisher_odds_ratio_mp_vs_tmi": float(high_qc_odds),
                "high_qc_fisher_p_value": float(high_qc_p_value),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary["full_public_fisher_fdr"] = benjamini_hochberg(np.asarray(p_values))
    summary["high_qc_fisher_fdr"] = benjamini_hochberg(
        np.asarray(high_qc_p_values)
    )

    public_panel = panel[
        (panel.source == "public") & panel.accessory_lineage.isin([TMI, MP])
    ].copy()
    mixed_projects = []
    for project, frame in public_panel.groupby("bioproject"):
        lineages = set(frame.accessory_lineage)
        if project and {TMI, MP}.issubset(lineages):
            mixed_projects.append(project)

    project_rows: list[dict[str, object]] = []
    cmh_rows: list[dict[str, object]] = []
    leave_one_out_rows: list[dict[str, object]] = []
    all_projects = sorted(set(public_panel.bioproject) - {""})
    for block in blocks.itertuples(index=False):
        block_id = str(block.block_id)
        direction = str(block.association_direction)
        present = block_presence[block_id]
        strata: list[np.ndarray] = []
        concordant = 0
        informative = 0
        for project in sorted(mixed_projects):
            stratum = public_panel[public_panel.bioproject == project]
            tmi = set(stratum[stratum.accessory_lineage == TMI].tree_id)
            mp = set(stratum[stratum.accessory_lineage == MP].tree_id)
            tmi_present = len(tmi & present)
            mp_present = len(mp & present)
            tmi_prev = tmi_present / len(tmi)
            mp_prev = mp_present / len(mp)
            gap = tmi_prev - mp_prev if direction == "TMI_enriched" else mp_prev - tmi_prev
            if gap != 0:
                informative += 1
                concordant += int(gap > 0)
            project_rows.append(
                {
                    "block_id": block_id,
                    "association_direction": direction,
                    "bioproject": project,
                    "public_tmi_present": tmi_present,
                    "public_tmi_total": len(tmi),
                    "public_tmi_prevalence": tmi_prev,
                    "public_mp_mip_present": mp_present,
                    "public_mp_mip_total": len(mp),
                    "public_mp_mip_prevalence": mp_prev,
                    "expected_direction_prevalence_gap": gap,
                    "direction_concordant": gap > 0,
                }
            )
            if direction == "TMI_enriched":
                enriched_present, enriched_total = tmi_present, len(tmi)
                other_present, other_total = mp_present, len(mp)
            else:
                enriched_present, enriched_total = mp_present, len(mp)
                other_present, other_total = tmi_present, len(tmi)
            strata.append(
                np.asarray(
                    [
                        [enriched_present, enriched_total - enriched_present],
                        [other_present, other_total - other_present],
                    ],
                    dtype=float,
                )
            )
        if strata:
            table = StratifiedTable(np.stack(strata, axis=2), shift_zeros=True)
            test = table.test_null_odds()
            odds_lcl, odds_ucl = table.oddsratio_pooled_confint()
            cmh_rows.append(
                {
                    "block_id": block_id,
                    "association_direction": direction,
                    "mixed_bioprojects": len(strata),
                    "informative_bioprojects": informative,
                    "direction_concordant_bioprojects": concordant,
                    "mantel_haenszel_common_odds_ratio_expected_direction": float(
                        table.oddsratio_pooled
                    ),
                    "mantel_haenszel_odds_ratio_95ci_lower": float(odds_lcl),
                    "mantel_haenszel_odds_ratio_95ci_upper": float(odds_ucl),
                    "mantel_haenszel_p_value": float(test.pvalue),
                }
            )

        for omitted in all_projects:
            remaining = public_panel[public_panel.bioproject != omitted]
            tmi = set(remaining[remaining.accessory_lineage == TMI].tree_id)
            mp = set(remaining[remaining.accessory_lineage == MP].tree_id)
            tmi_prev = len(tmi & present) / len(tmi)
            mp_prev = len(mp & present) / len(mp)
            gap = tmi_prev - mp_prev if direction == "TMI_enriched" else mp_prev - tmi_prev
            leave_one_out_rows.append(
                {
                    "block_id": block_id,
                    "association_direction": direction,
                    "omitted_bioproject": omitted,
                    "remaining_public_tmi": len(tmi),
                    "remaining_public_mp_mip": len(mp),
                    "expected_direction_prevalence_gap": gap,
                    "direction_retained": gap > 0,
                }
            )

    project_frame = pd.DataFrame(project_rows)
    project_frame.to_csv(
        output_dir / "curated_blocks_mixed_bioproject_strata.tsv", sep="\t", index=False
    )
    cmh = pd.DataFrame(cmh_rows)
    if not cmh.empty:
        cmh["mantel_haenszel_fdr"] = benjamini_hochberg(
            cmh.mantel_haenszel_p_value.to_numpy(dtype=float)
        )
    cmh.to_csv(
        output_dir / "curated_blocks_bioproject_stratified_statistics.tsv",
        sep="\t",
        index=False,
    )
    leave_one_out = pd.DataFrame(leave_one_out_rows)
    leave_one_out.to_csv(
        output_dir / "curated_blocks_leave_one_bioproject_out.tsv", sep="\t", index=False
    )
    robustness = (
        leave_one_out.groupby(["block_id", "association_direction"], as_index=False)
        .agg(
            omitted_bioprojects=("omitted_bioproject", "nunique"),
            all_omissions_retain_direction=("direction_retained", "all"),
            minimum_expected_direction_prevalence_gap=(
                "expected_direction_prevalence_gap",
                "min",
            ),
            maximum_expected_direction_prevalence_gap=(
                "expected_direction_prevalence_gap",
                "max",
            ),
        )
    )
    robustness.to_csv(
        output_dir / "curated_blocks_project_robustness_summary.tsv", sep="\t", index=False
    )

    summary = summary.merge(cmh, on=["block_id", "association_direction"], how="left")
    summary = summary.merge(
        robustness, on=["block_id", "association_direction"], how="left"
    )
    summary.to_csv(
        output_dir / "curated_blocks_full_atlas_validation.tsv", sep="\t", index=False
    )

    report = [
        "# Curated accessory blocks in the full current MI atlas",
        "",
        (
            f"The sensitivity panel contained {len(groups['public_tmi'])} public TMI, "
            f"{len(groups['public_mp_mip'])} public MP-MIP, {len(groups['local_tmi'])} local TMI "
            f"and {len(groups['local_mp_mip'])} local MP-MIP genomes."
        ),
        (
            f"The stricter assembly-quality sensitivity subset contained "
            f"{len(high_qc_public_tmi)} TMI and {len(high_qc_public_mp)} MP-MIP genomes "
            "(CheckM2 completeness >=95%, contamination <=2%, contigs <=200 and N50 >=20 kb)."
        ),
        "",
        "| Block | Public TMI | Public MP-MIP | Local TMI | Local MP-MIP | Fisher FDR | High-QC FDR | Mixed projects concordant | Leave-one-project-out |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary.to_dict(orient="records"):
        report.append(
            "| {block_id} | {full_public_tmi_present}/{full_public_tmi_total} | "
            "{full_public_mp_mip_present}/{full_public_mp_mip_total} | "
            "{full_local_tmi_present}/{full_local_tmi_total} | "
            "{full_local_mp_mip_present}/{full_local_mp_mip_total} | "
            "{full_public_fisher_fdr:.3g} | {high_qc_fisher_fdr:.3g} | "
            "{direction_concordant_bioprojects:.0f}/{informative_bioprojects:.0f} | "
            "{all_omissions_retain_direction} |".format(**row)
        )
    report.extend(
        [
            "",
            "A block was counted only when all member families occurred on one contig within the predefined gene-span limit. BioProject-stratified and leave-one-project-out analyses test study-source robustness; they do not imply phenotype, adaptation, horizontal transfer or causality.",
        ]
    )
    (output_dir / "curated_blocks_full_atlas_validation_report.md").write_text(
        "\n".join(report) + "\n"
    )
    print(
        f"blocks={summary.shape[0]} mixed_projects={len(mixed_projects)} "
        f"all_leave_one_out_robust={int(robustness.all_omissions_retain_direction.sum())}"
    )


if __name__ == "__main__":
    main()
