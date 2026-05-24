#!/usr/bin/env python3
"""Curate downloaded M. paraintracellulare public genomes against a type panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-raw", required=True, type=Path)
    parser.add_argument("--datasets-root", required=True, type=Path)
    parser.add_argument("--fastani-best", required=True, type=Path)
    parser.add_argument("--evaluated-output", required=True, type=Path)
    parser.add_argument("--curated-output", required=True, type=Path)
    parser.add_argument("--reference-list-output", required=True, type=Path)
    return parser.parse_args()


def accession_core(accession: str) -> str:
    if "_" in accession:
        return accession.split("_", 1)[1]
    return accession


def attr_value(row: dict, name: str) -> str:
    biosample = (row.get("assemblyInfo") or {}).get("biosample") or {}
    for attr in biosample.get("attributes") or []:
        if attr.get("name") == name:
            return attr.get("value", "NA")
    return "NA"


def norm_taxon(name: str) -> str:
    low = str(name).lower()
    if "paraintracellulare" in low:
        return "M_paraintracellulare"
    if "yongonense" in low:
        return "M_yongonense"
    if "chimaera" in low:
        return "M_chimaera"
    if "intracellulare" in low:
        return "M_intracellulare"
    if "avium" in low:
        return "M_avium"
    if "colombiense" in low:
        return "M_colombiense"
    return "other_or_unresolved"


def assembly_level_rank(level: str) -> int:
    return {
        "Complete Genome": 0,
        "Chromosome": 1,
        "Scaffold": 2,
        "Contig": 3,
    }.get(str(level), 4)


def parse_reports(root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for report in sorted(root.glob("*/ncbi_dataset/data/assembly_data_report.jsonl")):
        for line in report.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            accession = row.get("accession") or row.get("currentAccession")
            if not accession:
                continue
            stats = row.get("assemblyStats") or {}
            biosample = ((row.get("assemblyInfo") or {}).get("biosample") or {})
            rows.append(
                {
                    "accession": accession,
                    "paired_accession": row.get("pairedAccession", "NA"),
                    "assembly_size": stats.get("totalSequenceLength", "NA"),
                    "gc_percent": stats.get("gcPercent", "NA"),
                    "num_contigs": stats.get("numberOfContigs", "NA"),
                    "contig_n50": stats.get("contigN50", "NA"),
                    "host": biosample.get("host", attr_value(row, "host")),
                    "isolation_source": biosample.get("isolationSource", attr_value(row, "isolation_source")),
                    "collection_date": biosample.get("collectionDate", attr_value(row, "collection_date")),
                    "geo_loc_name": biosample.get("geoLocName", attr_value(row, "geo_loc_name")),
                }
            )
    return pd.DataFrame(rows).drop_duplicates("accession")


def main() -> int:
    args = parse_args()
    meta = pd.read_csv(args.metadata_raw, sep="\t", keep_default_na=False)
    stats = parse_reports(args.datasets_root)
    best = pd.read_csv(args.fastani_best, sep="\t", keep_default_na=False)

    df = meta.merge(stats, on="accession", how="left")
    best_subset = best[
        [
            "query_fasta",
            "fastani_best_hit_accession",
            "fastani_best_hit_species",
            "fastani_best_hit_strain",
            "fastani_best_hit_ani",
            "fastani_best_hit_alignment_fraction",
        ]
    ].rename(columns={"query_fasta": "reference_fasta"})
    df = df.merge(best_subset, on="reference_fasta", how="left")
    df["best_type_panel_norm"] = df["fastani_best_hit_species"].map(norm_taxon)

    for col in [
        "checkm_completeness",
        "checkm_contamination",
        "fastani_best_hit_ani",
        "fastani_best_hit_alignment_fraction",
        "assembly_size",
        "gc_percent",
        "num_contigs",
        "contig_n50",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["passes_basic_qc"] = (
        (df["checkm_completeness"] >= 90)
        & (df["checkm_contamination"] <= 5)
        & (df["num_contigs"] <= 500)
        & (df["contig_n50"] >= 20_000)
        & (df["assembly_size"].between(4_000_000, 7_000_000))
        & (df["gc_percent"].between(60, 72))
    )
    df["passes_paraintracellulare_screen"] = (
        (df["best_type_panel_norm"] == "M_paraintracellulare")
        & (df["fastani_best_hit_ani"] >= 98.0)
        & (df["fastani_best_hit_alignment_fraction"] >= 0.80)
    )
    df["curation_status"] = "exclude_or_manual_review"
    strict = df["passes_basic_qc"] & df["passes_paraintracellulare_screen"] & (df["fastani_best_hit_ani"] >= 98.5)
    review = df["passes_basic_qc"] & df["passes_paraintracellulare_screen"] & (df["fastani_best_hit_ani"] < 98.5)
    df.loc[strict, "curation_status"] = "tier2_curated_strict"
    df.loc[review, "curation_status"] = "tier2_candidate_review_98_0_to_98_5"

    df["paired_group"] = df["accession"].map(accession_core)
    df["source_rank"] = df["source_database"].map({"SOURCE_DATABASE_REFSEQ": 0, "SOURCE_DATABASE_GENBANK": 1}).fillna(2)
    df["assembly_level_rank"] = df["assembly_level"].map(assembly_level_rank)
    df = df.sort_values(
        [
            "paired_group",
            "curation_status",
            "source_rank",
            "assembly_level_rank",
            "checkm_contamination",
            "contig_n50",
        ],
        ascending=[True, True, True, True, True, False],
    )

    curated = df[df["curation_status"].str.startswith("tier2_")].copy()
    curated = curated.sort_values(
        ["paired_group", "source_rank", "assembly_level_rank", "checkm_contamination", "contig_n50"],
        ascending=[True, True, True, True, False],
    ).drop_duplicates("paired_group", keep="first")

    args.evaluated_output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.evaluated_output, sep="\t", index=False)
    curated.to_csv(args.curated_output, sep="\t", index=False)
    args.reference_list_output.write_text(
        "\n".join(curated["reference_fasta"].astype(str)) + ("\n" if not curated.empty else "")
    )
    print(f"evaluated={len(df)} curated_deduplicated={len(curated)}")
    print(df["curation_status"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
