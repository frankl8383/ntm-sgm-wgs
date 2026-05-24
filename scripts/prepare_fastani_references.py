#!/usr/bin/env python3
"""Prepare NCBI Datasets genomes as a FastANI reference panel."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd


def biosample_attribute(row: dict, name: str) -> str:
    biosample = (row.get("assemblyInfo") or {}).get("biosample") or {}
    for attr in biosample.get("attributes") or []:
        if attr.get("name") == name:
            return attr.get("value", "NA")
    return "NA"


def find_fasta(data_dir: Path, accession: str) -> Path | None:
    accession_dir = data_dir / accession
    candidates = sorted(accession_dir.glob("*_genomic.fna"))
    if candidates:
        return candidates[0]
    candidates = sorted(data_dir.glob(f"{accession}/**/*_genomic.fna"))
    return candidates[0] if candidates else None


def safe_link_name(accession: str, organism: str, strain: str) -> str:
    cleaned = "_".join([organism, strain]).replace("[", "").replace("]", "")
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in cleaned)
    while "__" in cleaned:
        cleaned = cleaned.replace("__", "_")
    return f"{accession}_{cleaned}.fna"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets-root", required=True, type=Path)
    parser.add_argument("--reference-dir", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--reference-list", required=True, type=Path)
    args = parser.parse_args()

    args.reference_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []

    for report in sorted(args.datasets_root.glob("*/ncbi_dataset/data/assembly_data_report.jsonl")):
        data_dir = report.parent
        source_dataset = report.parents[2].name
        for line in report.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            accession = row.get("accession") or row.get("currentAccession")
            if not accession:
                continue
            fasta = find_fasta(data_dir, accession)
            if fasta is None:
                continue
            organism = (row.get("organism") or {}).get("organismName", "NA")
            infraspecific = (row.get("organism") or {}).get("infraspecificNames") or {}
            strain = infraspecific.get("strain") or infraspecific.get("isolate") or biosample_attribute(row, "strain")
            strain = strain or "NA"
            assembly_info = row.get("assemblyInfo") or {}
            checkm_info = row.get("checkmInfo") or {}
            type_material = row.get("typeMaterial") or {}
            link = args.reference_dir / safe_link_name(accession, organism, strain)
            if link.exists() or link.is_symlink():
                link.unlink()
            os.symlink(os.path.relpath(fasta.resolve(), link.parent), link)
            rows.append(
                {
                    "accession": accession,
                    "organism_name": organism,
                    "strain": strain,
                    "assembly_level": assembly_info.get("assemblyLevel", "NA"),
                    "assembly_name": assembly_info.get("assemblyName", "NA"),
                    "source_database": row.get("sourceDatabase", "NA"),
                    "source_dataset": source_dataset,
                    "type_material": json.dumps(type_material, ensure_ascii=False) if type_material else "NA",
                    "checkm_completeness": checkm_info.get("completeness", "NA"),
                    "checkm_contamination": checkm_info.get("contamination", "NA"),
                    "reference_fasta": str(link),
                    "original_fasta": str(fasta),
                }
            )

    df = pd.DataFrame(rows).drop_duplicates(subset=["accession"]).sort_values(["organism_name", "accession"])
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.metadata, sep="\t", index=False)
    args.reference_list.parent.mkdir(parents=True, exist_ok=True)
    args.reference_list.write_text("\n".join(df["reference_fasta"].astype(str)) + ("\n" if not df.empty else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
