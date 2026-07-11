#!/usr/bin/env python3
"""Flatten, quality-filter and deduplicate a current NCBI MAC genome snapshot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


MIN_LENGTH = 4_000_000
MAX_LENGTH = 7_500_000
MIN_GC = 66.0
MAX_GC = 71.0
MAX_CONTIGS = 1_000
MIN_N50 = 20_000
MIN_COMPLETENESS = 90.0
MAX_CONTAMINATION = 5.0


def text(value: object) -> str:
    return "" if value is None else str(value)


def fine_label(organism: str, reporting_species: str) -> str:
    lower = organism.lower()
    for token, label in (
        ("chimaera", "Mycobacterium chimaera legacy subgroup"),
        ("yongonense", "Mycobacterium yongonense legacy subgroup"),
        ("hominissuis", "Mycobacterium avium subsp. hominissuis"),
        ("paratuberculosis", "Mycobacterium avium subsp. paratuberculosis"),
        ("subsp. avium", "Mycobacterium avium subsp. avium"),
        ("silvaticum", "Mycobacterium avium subsp. silvaticum"),
    ):
        if token in lower:
            return label
    return reporting_species or organism


def sample_attribute(biosample: dict, name: str) -> str:
    for item in biosample.get("attributes", []) or []:
        if text(item.get("name")).lower() == name.lower():
            return text(item.get("value"))
    return ""


def assembly_score(row: dict[str, object]) -> float:
    score = 100.0 if row["source_database"] == "RefSeq" else 0.0
    category = str(row["refseq_category"]).lower()
    if category == "reference genome":
        score += 50
    elif category == "representative genome":
        score += 35
    if row["type_material"] == "true":
        score += 30
    score += {
        "complete genome": 25,
        "chromosome": 20,
        "scaffold": 10,
        "contig": 5,
    }.get(str(row["assembly_level"]).lower(), 0)
    completeness = float(row["checkm_completeness"] or 0)
    contamination = float(row["checkm_contamination"] or 100)
    score += completeness - 2 * contamination
    score += min(20.0, math.log10(max(1, int(row["contig_n50"]))) * 3)
    return score


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonl", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir)

    flattened: list[dict[str, object]] = []
    with Path(args.jsonl).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            assembly = record.get("assembly_info", {}) or {}
            stats = record.get("assembly_stats", {}) or {}
            organism = record.get("organism", {}) or {}
            ani = record.get("average_nucleotide_identity", {}) or {}
            checkm = record.get("checkm_info", {}) or {}
            biosample = assembly.get("biosample", {}) or {}
            accession = text(record.get("current_accession") or record.get("accession"))
            paired = text(record.get("paired_accession") or (assembly.get("paired_assembly", {}) or {}).get("accession"))
            source = "RefSeq" if accession.startswith("GCF_") else "GenBank"
            reporting_species = text(ani.get("submitted_species"))
            organism_name = text(organism.get("organism_name"))
            if not reporting_species:
                reporting_species = organism_name
            pair_key = "|".join(sorted(x for x in {accession, paired} if x))
            bioproject = text(assembly.get("bioproject_accession"))
            if not bioproject:
                projects = biosample.get("bioprojects", []) or []
                bioproject = text(projects[0].get("accession")) if projects else ""
            row: dict[str, object] = {
                "assembly_accession": accession,
                "paired_accession": paired,
                "pair_key": pair_key or accession,
                "source_database": source,
                "organism_name": organism_name,
                "reporting_species": reporting_species,
                "fine_label": fine_label(organism_name, reporting_species),
                "tax_id": text(organism.get("tax_id")),
                "strain": text((organism.get("infraspecific_names", {}) or {}).get("strain") or biosample.get("strain")),
                "biosample": text(biosample.get("accession")),
                "bioproject": bioproject,
                "host": text(biosample.get("host") or sample_attribute(biosample, "host")),
                "geo_loc_name": text(biosample.get("geo_loc_name") or sample_attribute(biosample, "geo_loc_name")),
                "isolation_source": text(biosample.get("isolation_source") or sample_attribute(biosample, "isolation_source")),
                "assembly_level": text(assembly.get("assembly_level")),
                "assembly_name": text(assembly.get("assembly_name")),
                "refseq_category": text(assembly.get("refseq_category")),
                "release_date": text(assembly.get("release_date")),
                "type_material": "true" if record.get("type_material") else "false",
                "total_length": int(stats.get("total_sequence_length") or 0),
                "gc_percent": float(stats.get("gc_percent") or 0),
                "number_of_contigs": int(stats.get("number_of_contigs") or 0),
                "contig_n50": int(stats.get("contig_n50") or 0),
                "checkm_completeness": text(checkm.get("completeness")),
                "checkm_contamination": text(checkm.get("contamination")),
                "taxonomy_check_status": text(ani.get("taxonomy_check_status")),
            }
            row["selection_score"] = f"{assembly_score(row):.6f}"
            flattened.append(row)

    by_pair: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in flattened:
        by_pair[str(row["pair_key"])].append(row)
    pair_selected: list[dict[str, object]] = []
    for group in by_pair.values():
        pair_selected.append(max(group, key=lambda row: float(row["selection_score"])))

    by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in pair_selected:
        key = str(row["biosample"]) or f"no_biosample:{row['assembly_accession']}"
        by_sample[key].append(row)
    biosample_selected: list[dict[str, object]] = []
    for group in by_sample.values():
        biosample_selected.append(max(group, key=lambda row: float(row["selection_score"])))

    for row in flattened:
        length = int(row["total_length"])
        gc = float(row["gc_percent"])
        contigs = int(row["number_of_contigs"])
        n50 = int(row["contig_n50"])
        completeness = float(row["checkm_completeness"] or 0)
        contamination = float(row["checkm_contamination"] or 100)
        failures = []
        if not MIN_LENGTH <= length <= MAX_LENGTH:
            failures.append("genome_size")
        if not MIN_GC <= gc <= MAX_GC:
            failures.append("gc")
        if contigs > MAX_CONTIGS:
            failures.append("contigs")
        if n50 < MIN_N50:
            failures.append("n50")
        if completeness < MIN_COMPLETENESS:
            failures.append("checkm_completeness")
        if contamination > MAX_CONTAMINATION:
            failures.append("checkm_contamination")
        row["qc_pass"] = "true" if not failures else "false"
        row["qc_failure_reasons"] = ";".join(failures)

    pair_accessions = {str(row["assembly_accession"]) for row in pair_selected}
    biosample_accessions = {str(row["assembly_accession"]) for row in biosample_selected}
    for row in flattened:
        accession = str(row["assembly_accession"])
        row["selected_after_pair_dedup"] = "true" if accession in pair_accessions else "false"
        row["selected_after_biosample_dedup"] = "true" if accession in biosample_accessions else "false"
        row["selected_for_download"] = (
            "true"
            if accession in biosample_accessions and row["qc_pass"] == "true"
            else "false"
        )

    flattened.sort(key=lambda row: str(row["assembly_accession"]))
    selected = [row for row in flattened if row["selected_for_download"] == "true"]
    fields = list(flattened[0].keys())
    write_tsv(output_dir / "ncbi_mac_snapshot_all.tsv", flattened, fields)
    write_tsv(output_dir / "ncbi_mac_atlas_qc_selected.tsv", selected, fields)
    (output_dir / "ncbi_mac_atlas_download_accessions.txt").write_text(
        "\n".join(str(row["assembly_accession"]) for row in selected) + "\n"
    )

    species = Counter(str(row["reporting_species"]) for row in selected)
    fine = Counter(str(row["fine_label"]) for row in selected)
    report = [
        "# Current NCBI MAC atlas metadata",
        "",
        f"Raw current records: {len(flattened)}",
        f"After paired-accession deduplication: {len(pair_selected)}",
        f"After BioSample deduplication: {len(biosample_selected)}",
        f"QC-selected for download: {len(selected)}",
        "",
        "## Current reporting species",
        "",
    ]
    report.extend(f"- {label}: {count}" for label, count in species.most_common())
    report.extend(["", "## Fine labels retained for boundary analysis", ""])
    report.extend(f"- {label}: {count}" for label, count in fine.most_common())
    (output_dir / "ncbi_mac_atlas_metadata_report.md").write_text("\n".join(report) + "\n")
    print(
        f"raw={len(flattened)} pair_dedup={len(pair_selected)} "
        f"biosample_dedup={len(biosample_selected)} selected={len(selected)}"
    )


if __name__ == "__main__":
    main()
