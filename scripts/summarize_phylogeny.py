#!/usr/bin/env python3
"""Audit current SKA/IQ-TREE panels and finalize local analysis labels."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from Bio import Phylo, SeqIO


PANELS = (
    "intracellulare_complex_context",
    "avium_timonense_context",
    "colombiense_context",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def support_value(clade) -> str:
    if clade.confidence is None:
        return ""
    return f"{float(clade.confidence):.2f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-root", required=True)
    parser.add_argument("--pre-phylogeny-labels", required=True)
    parser.add_argument("--public-anchor-proximity", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    panel_root = Path(args.panel_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    public_anchor = {
        row["public_reference_id"]: row
        for row in read_tsv(Path(args.public_anchor_proximity))
    }
    label_rows = read_tsv(Path(args.pre_phylogeny_labels))
    label_by_sample = {row["sample_id"]: row for row in label_rows}

    panel_stats: list[dict[str, object]] = []
    nearest_rows: list[dict[str, object]] = []
    monophyly_rows: list[dict[str, object]] = []

    for panel in PANELS:
        panel_dir = panel_root / panel
        manifest = read_tsv(panel_dir / "panel_manifest.tsv")
        meta_by_id = {row["tree_id"]: row for row in manifest}
        alignment = list(SeqIO.parse(panel_dir / "core_k31_complete_snps.fasta", "fasta"))
        lengths = {len(record.seq) for record in alignment}
        if len(lengths) != 1:
            raise SystemExit(f"Non-rectangular alignment for {panel}")
        alignment_length = next(iter(lengths))
        missing = sum(
            str(record.seq).upper().count("N") + str(record.seq).count("-")
            for record in alignment
        )
        tree = Phylo.read(panel_dir / "iqtree_gtr_asc.treefile", "newick")
        terminals = tree.get_terminals()
        terminal_by_name = {tip.name: tip for tip in terminals}
        if set(terminal_by_name) != set(meta_by_id):
            raise SystemExit(
                f"Tree/manifest mismatch for {panel}: tree={len(terminal_by_name)} "
                f"manifest={len(meta_by_id)}"
            )
        panel_stats.append(
            {
                "panel": panel,
                "n_genomes": len(manifest),
                "n_local": sum(row["source"] == "local" for row in manifest),
                "n_public": sum(row["source"] == "public" for row in manifest),
                "n_added_anchors": sum(row["source"] == "anchor" for row in manifest),
                "snp_alignment_sites": alignment_length,
                "missing_or_gap_cells": missing,
                "model": "GTR+F+ASC",
                "ultrafast_bootstrap_replicates": 1000,
            }
        )

        public_tips = [tip for tip in terminals if meta_by_id[tip.name]["source"] == "public"]
        local_tips = [tip for tip in terminals if meta_by_id[tip.name]["source"] == "local"]
        for local_tip in local_tips:
            nearest = min(public_tips, key=lambda tip: tree.distance(local_tip, tip))
            local_meta = meta_by_id[local_tip.name]
            public_meta = meta_by_id[nearest.name]
            accession = nearest.name.removeprefix("PUB_")
            anchor = public_anchor[accession]
            nearest_rows.append(
                {
                    "sample_id": local_meta["sample_id"],
                    "panel": panel,
                    "provisional_genomic_lineage": local_meta["analysis_lineage"],
                    "nearest_public_accession": accession,
                    "nearest_public_reporting_label": public_meta["reporting_label"],
                    "nearest_public_fine_label": public_meta["fine_label"],
                    "nearest_public_anchor_proximity_group": anchor["best_type_broad_group"],
                    "nearest_public_anchor_ani": anchor["best_type_ani"],
                    "patristic_distance": f"{tree.distance(local_tip, nearest):.10f}",
                }
            )

        by_lineage: dict[str, list] = {}
        for tip in local_tips:
            lineage = meta_by_id[tip.name]["analysis_lineage"]
            by_lineage.setdefault(lineage, []).append(tip)
        for lineage, tips in sorted(by_lineage.items()):
            if len(tips) == 1:
                monophyly_rows.append(
                    {
                        "panel": panel,
                        "analysis_lineage": lineage,
                        "n_local": 1,
                        "local_monophyly_testable": "false",
                        "local_monophyletic_relative_to_other_local_labels": "not_applicable_singleton",
                        "mrca_ultrafast_bootstrap": "",
                        "public_labels_inside_mrca": "",
                    }
                )
                continue
            mrca = tree.common_ancestor(tips)
            mrca_local = {
                tip.name
                for tip in mrca.get_terminals()
                if meta_by_id[tip.name]["source"] == "local"
            }
            expected = {tip.name for tip in tips}
            public_labels = Counter(
                meta_by_id[tip.name]["fine_label"]
                for tip in mrca.get_terminals()
                if meta_by_id[tip.name]["source"] == "public"
            )
            monophyly_rows.append(
                {
                    "panel": panel,
                    "analysis_lineage": lineage,
                    "n_local": len(tips),
                    "local_monophyly_testable": "true",
                    "local_monophyletic_relative_to_other_local_labels": (
                        "true" if mrca_local == expected else "false"
                    ),
                    "mrca_ultrafast_bootstrap": support_value(mrca),
                    "public_labels_inside_mrca": ";".join(
                        f"{label}:{count}" for label, count in public_labels.most_common()
                    ),
                }
            )

    write_tsv(output_dir / "current_ska_iqtree_panel_stats.tsv", panel_stats)
    write_tsv(output_dir / "local21_iqtree_nearest_public.tsv", nearest_rows)
    write_tsv(output_dir / "local_lineage_monophyly.tsv", monophyly_rows)

    nearest_by_sample = {row["sample_id"]: row for row in nearest_rows}
    mono_by_lineage = {row["analysis_lineage"]: row for row in monophyly_rows}
    final_rows: list[dict[str, object]] = []
    for row in label_rows:
        sample = row["sample_id"]
        nearest = nearest_by_sample[sample]
        mono = mono_by_lineage[row["provisional_genomic_lineage"]]
        public_label = nearest["nearest_public_reporting_label"]
        consistency = "context_supported"
        if row["broad_analysis_panel"] == "M_colombiense" and public_label != "Mycobacterium colombiense":
            consistency = "context_conflict"
        elif row["broad_analysis_panel"] == "M_avium_timonense_boundary" and public_label != "Mycobacterium avium":
            consistency = "context_conflict"
        elif row["broad_analysis_panel"] == "M_intracellulare_complex" and public_label not in {
            "Mycobacterium intracellulare",
            "Mycobacterium paraintracellulare",
            "Mycobacterium marseillense",
        }:
            consistency = "context_conflict"
        if mono["local_monophyletic_relative_to_other_local_labels"] == "false":
            consistency = "lineage_not_monophyletic_among_local_labels"
        final_rows.append(
            {
                **row,
                "final_manuscript_wording": (
                    "Mycobacterium avium public-context group"
                    if row["broad_analysis_panel"] == "M_avium_timonense_boundary"
                    else row["manuscript_wording"]
                ),
                "nearest_tree_public_accession": nearest["nearest_public_accession"],
                "nearest_tree_public_reporting_label": public_label,
                "nearest_tree_public_fine_label": nearest["nearest_public_fine_label"],
                "nearest_tree_public_anchor_proximity_group": nearest[
                    "nearest_public_anchor_proximity_group"
                ],
                "tree_context_consistency": consistency,
                "lineage_mrca_ultrafast_bootstrap": mono["mrca_ultrafast_bootstrap"],
                "label_status": "final_analysis_label_after_context_tree_review",
            }
        )
    write_tsv(output_dir / "local21_analysis_labels_final.tsv", final_rows)

    report = [
        "# Current SKA/IQ-TREE validation of the updated 21-genome cohort",
        "",
        "## Panels",
        "",
        "| Panel | Genomes | Local | Public | Anchors | Complete-core SNP sites | Missing/gap cells |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in panel_stats:
        report.append(
            "| {panel} | {n_genomes} | {n_local} | {n_public} | {n_added_anchors} | "
            "{snp_alignment_sites} | {missing_or_gap_cells} |".format(**row)
        )
    report.extend(
        [
            "",
            "## Local-lineage topology",
            "",
            "| Panel | Lineage | Local n | Monophyletic relative to other local labels | MRCA UFBoot | Public context inside MRCA |",
            "|---|---|---:|---|---:|---|",
        ]
    )
    for row in monophyly_rows:
        report.append(
            "| {panel} | {analysis_lineage} | {n_local} | "
            "{local_monophyletic_relative_to_other_local_labels} | "
            "{mrca_ultrafast_bootstrap} | {public_labels_inside_mrca} |".format(**row)
        )
    report.extend(
        [
            "",
            "All trees are public-context phylogenies. They are not interpreted as evidence of transmission, outbreak linkage or formal taxonomic acts.",
        ]
    )
    (output_dir / "current_ska_iqtree_validation_report.md").write_text(
        "\n".join(report) + "\n"
    )
    print(
        f"panels={len(panel_stats)} local={len(nearest_rows)} "
        f"lineages={len(monophyly_rows)}"
    )


if __name__ == "__main__":
    main()
