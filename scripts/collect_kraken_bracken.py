#!/usr/bin/env python3
"""Collect fastp, Kraken2, and Bracken read-level taxonomy summaries."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


SUMMARY_FIELDS = [
    "sample_id",
    "total_reads",
    "clean_reads",
    "q30_rate",
    "top_kingdom",
    "top_genus",
    "top_species",
    "mycobacterium_fraction",
    "top_species_fraction",
    "second_species",
    "second_species_fraction",
    "suspected_mixed_or_contaminated",
]

LONG_FIELDS = [
    "sample_id",
    "species",
    "taxonomy_id",
    "taxonomy_level",
    "kraken_assigned_reads",
    "added_reads",
    "estimated_reads",
    "fraction_total_reads",
]


@dataclass
class TaxonRow:
    name: str
    taxonomy_id: str
    taxonomy_level: str
    kraken_assigned_reads: float
    added_reads: float
    estimated_reads: float
    fraction_total_reads: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--fastp-dir", required=True)
    parser.add_argument("--kraken-dir", required=True)
    parser.add_argument("--bracken-dir", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--species-long-output", required=True)
    return parser.parse_args()


def read_samples(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row["sample_id"] for row in reader if row.get("sample_id")]


def get_nested(data: dict[str, Any], *keys: str, default: Any = None) -> Any:
    value: Any = data
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def parse_fastp_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    before = get_nested(data, "summary", "before_filtering", default={}) or {}
    after = get_nested(data, "summary", "after_filtering", default={}) or {}
    total_reads = before.get("total_reads", "NA")
    clean_reads = after.get("total_reads", "NA")
    q30_rate = after.get("q30_rate", before.get("q30_rate", "NA"))
    return {
        "total_reads": total_reads,
        "clean_reads": clean_reads,
        "q30_rate": q30_rate,
    }


def parse_kraken_report(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                pct = float(parts[0].strip())
                reads_clade = int(parts[1].strip())
                reads_taxon = int(parts[2].strip())
            except ValueError:
                continue
            rows.append(
                {
                    "pct": pct,
                    "reads_clade": reads_clade,
                    "reads_taxon": reads_taxon,
                    "rank": parts[3].strip(),
                    "taxid": parts[4].strip(),
                    "name": parts[5].strip(),
                }
            )
    return rows


def top_kraken_name(rows: list[dict[str, Any]], ranks: set[str]) -> str:
    candidates = [row for row in rows if row["rank"] in ranks]
    if not candidates:
        return "NA"
    return max(candidates, key=lambda row: row["pct"])["name"]


def parse_float(value: str | None) -> float:
    if value is None or value == "":
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def parse_bracken(path: Path) -> list[TaxonRow]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[TaxonRow] = []
        for row in reader:
            name = (row.get("name") or row.get("Name") or "").strip()
            if not name:
                continue
            rows.append(
                TaxonRow(
                    name=name,
                    taxonomy_id=(row.get("taxonomy_id") or row.get("taxonomy_lvl_id") or row.get("taxid") or "").strip(),
                    taxonomy_level=(row.get("taxonomy_lvl") or row.get("level") or "").strip(),
                    kraken_assigned_reads=parse_float(row.get("kraken_assigned_reads")),
                    added_reads=parse_float(row.get("added_reads")),
                    estimated_reads=parse_float(row.get("new_est_reads") or row.get("estimated_reads")),
                    fraction_total_reads=parse_float(row.get("fraction_total_reads")),
                )
            )
    rows.sort(key=lambda item: item.fraction_total_reads, reverse=True)
    return rows


def fraction_for_name(rows: list[TaxonRow], name: str) -> float:
    target = name.lower()
    for row in rows:
        if row.name.lower() == target:
            return row.fraction_total_reads
    return 0.0


def format_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def main() -> int:
    args = parse_args()
    samples = read_samples(Path(args.samplesheet))
    fastp_dir = Path(args.fastp_dir)
    kraken_dir = Path(args.kraken_dir)
    bracken_dir = Path(args.bracken_dir)

    with Path(args.thresholds).open("r", encoding="utf-8") as handle:
        thresholds = yaml.safe_load(handle) or {}
    tax_thresholds = thresholds.get("taxonomy_reads", {})
    min_mycobacterium = float(tax_thresholds.get("min_mycobacterium_fraction_for_ntm_candidate", 0.70))
    min_top_species = float(tax_thresholds.get("min_top_species_fraction_for_single_species", 0.60))
    mixed_warning = float(tax_thresholds.get("mixed_species_warning_fraction", 0.10))

    summary_rows: list[dict[str, Any]] = []
    species_long_rows: list[dict[str, Any]] = []

    for sample_id in samples:
        fastp_stats = parse_fastp_json(fastp_dir / f"{sample_id}.fastp.json")
        kraken_rows = parse_kraken_report(kraken_dir / f"{sample_id}.kraken2.report.txt")
        species_rows = parse_bracken(bracken_dir / f"{sample_id}.species.bracken.tsv")
        genus_rows = parse_bracken(bracken_dir / f"{sample_id}.genus.bracken.tsv")

        top_kingdom = top_kraken_name(kraken_rows, {"D", "K"})
        if genus_rows:
            top_genus = genus_rows[0].name
        elif species_rows:
            top_genus = species_rows[0].name.split()[0]
        else:
            top_genus = top_kraken_name(kraken_rows, {"G"})
        top_species = species_rows[0].name if species_rows else top_kraken_name(kraken_rows, {"S"})
        top_species_fraction = species_rows[0].fraction_total_reads if species_rows else 0.0
        second_species = species_rows[1].name if len(species_rows) > 1 else "NA"
        second_species_fraction = species_rows[1].fraction_total_reads if len(species_rows) > 1 else 0.0
        mycobacterium_fraction = fraction_for_name(genus_rows, "Mycobacterium")
        if mycobacterium_fraction == 0.0:
            for row in kraken_rows:
                if row["rank"] == "G" and row["name"].lower() == "mycobacterium":
                    mycobacterium_fraction = row["pct"] / 100.0
                    break

        suspected = (
            mycobacterium_fraction < min_mycobacterium
            or second_species_fraction >= mixed_warning
            or (top_species_fraction > 0 and top_species_fraction < min_top_species)
        )

        summary_row: dict[str, Any] = {
            "sample_id": sample_id,
            **fastp_stats,
            "top_kingdom": top_kingdom,
            "top_genus": top_genus,
            "top_species": top_species,
            "mycobacterium_fraction": mycobacterium_fraction,
            "top_species_fraction": top_species_fraction,
            "second_species": second_species,
            "second_species_fraction": second_species_fraction,
            "suspected_mixed_or_contaminated": "TRUE" if suspected else "FALSE",
        }
        summary_rows.append(summary_row)

        for row in species_rows:
            species_long_rows.append(
                {
                    "sample_id": sample_id,
                    "species": row.name,
                    "taxonomy_id": row.taxonomy_id,
                    "taxonomy_level": row.taxonomy_level,
                    "kraken_assigned_reads": row.kraken_assigned_reads,
                    "added_reads": row.added_reads,
                    "estimated_reads": row.estimated_reads,
                    "fraction_total_reads": row.fraction_total_reads,
                }
            )

    summary_output = Path(args.summary_output)
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    with summary_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({field: format_value(row.get(field, "NA")) for field in SUMMARY_FIELDS})

    species_output = Path(args.species_long_output)
    species_output.parent.mkdir(parents=True, exist_ok=True)
    with species_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LONG_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in species_long_rows:
            writer.writerow({field: format_value(row.get(field, "NA")) for field in LONG_FIELDS})

    print(f"Wrote summary for {len(summary_rows)} samples: {summary_output}")
    print(f"Wrote {len(species_long_rows)} species abundance rows: {species_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
