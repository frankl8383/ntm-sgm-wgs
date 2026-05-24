#!/usr/bin/env python3
"""Summarize geNomad plasmid/virus/provirus outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_data_rows(path: Path) -> int:
    return len(read_rows(path))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-summary", type=Path, default=Path("results/tables/genomad_run_summary.tsv"))
    parser.add_argument("--sample-metadata", type=Path, default=Path("results/tables/table1_high_confidence_isolates.tsv"))
    parser.add_argument("--sample-output", type=Path, default=Path("results/tables/genomad_mobilome_sample_summary.tsv"))
    parser.add_argument("--contig-output", type=Path, default=Path("results/tables/genomad_mobilome_contig_calls.tsv"))
    args = parser.parse_args()

    metadata = {row["sample_id"]: row for row in read_rows(args.sample_metadata)}
    run_rows = read_rows(args.run_summary)
    sample_rows: list[dict[str, str | int]] = []
    contig_rows: list[dict[str, str]] = []

    for row in run_rows:
        sample_id = row["sample_id"]
        meta = metadata.get(sample_id, {})
        plasmid_path = Path(row["plasmid_summary"])
        virus_path = Path(row["virus_summary"])
        provirus_path = Path(row["provirus_summary"])
        plasmids = read_rows(plasmid_path)
        viruses = read_rows(virus_path)
        proviruses = read_rows(provirus_path)
        sample_rows.append(
            {
                "sample_id": sample_id,
                "public_context_clade": meta.get("public_context_clade", ""),
                "species_confidence_tier": meta.get("species_level_confidence_tier", ""),
                "genomad_status": row["status"],
                "n_conservative_plasmids": len(plasmids),
                "n_conservative_viruses": len(viruses),
                "n_find_provirus_regions": len(proviruses),
                "interpretation_note": "geNomad conservative calls are putative plasmid/virus/provirus signals; short-read draft assemblies do not prove complete plasmids or horizontal transfer.",
            }
        )
        for source, records in [
            ("plasmid_summary", plasmids),
            ("virus_summary", viruses),
            ("find_proviruses", proviruses),
        ]:
            for record in records:
                out = {
                    "sample_id": sample_id,
                    "public_context_clade": meta.get("public_context_clade", ""),
                    "call_source": source,
                    "interpretation_note": "Putative geNomad signal; not proof of complete mobile element.",
                }
                out.update(record)
                contig_rows.append(out)

    args.sample_output.parent.mkdir(parents=True, exist_ok=True)
    sample_fields = list(sample_rows[0]) if sample_rows else ["sample_id"]
    with args.sample_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_rows)

    contig_fields = sorted({key for row in contig_rows for key in row}) if contig_rows else ["sample_id", "call_source"]
    preferred = ["sample_id", "public_context_clade", "call_source", "seq_name", "source_seq", "start", "end", "length"]
    contig_fields = [field for field in preferred if field in contig_fields] + [field for field in contig_fields if field not in preferred]
    with args.contig_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=contig_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(contig_rows)


if __name__ == "__main__":
    main()
