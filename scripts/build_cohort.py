#!/usr/bin/env python3
"""Build the updated 21-genome local cohort and all-candidate rescue ledger."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--direct-manifest", required=True)
    parser.add_argument("--full-rescue-evidence", required=True)
    parser.add_argument("--pilot-evidence", required=True)
    parser.add_argument("--batch-root", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    direct = read_tsv(Path(args.direct_manifest))
    rescue = read_tsv(Path(args.full_rescue_evidence))
    pilot = {row["sample_id"]: row for row in read_tsv(Path(args.pilot_evidence))}
    batch_root = Path(args.batch_root).resolve()
    output_dir = Path(args.output_dir)

    cohort: list[dict[str, object]] = []
    for row in direct:
        assembly = Path(row["fasta"]).resolve()
        if not assembly.is_file():
            raise FileNotFoundError(assembly)
        cohort.append(
            {
                "sample_id": row["sample_id"],
                "cohort_source": "directly_retained_original",
                "selected_assembly_route": "original_isolate_assembly",
                "selected_assembly_path": str(assembly),
                "selected_assembly_size_bytes": assembly.stat().st_size,
                "selected_assembly_md5": md5sum(assembly),
                "analysis_label": row["species"],
                "species_group": row["species_group"],
                "taxonomic_caution": "use_existing_tiered_species_wording",
                "checkm2_completeness": row["checkm_completeness"],
                "checkm2_contamination": row["checkm_contamination"],
                "rescue_decision": "not_applicable_directly_retained",
            }
        )

    rescued_rows = [
        row for row in rescue if row["final_rescue_decision"].startswith("rescued_interpretable")
    ]
    for row in rescued_rows:
        sample = row["sample_id"]
        assembly = batch_root / "meta" / sample / f"{sample}.metaspades_mac_bin.min500.fasta"
        if not assembly.is_file():
            raise FileNotFoundError(assembly)
        reporting_name = row["type_anchor_reporting_name"]
        cohort.append(
            {
                "sample_id": sample,
                "cohort_source": "rescued_from_mixed_or_failed_original_assembly",
                "selected_assembly_route": "metaspades_mac_bin_validated_by_strict_route",
                "selected_assembly_path": str(assembly),
                "selected_assembly_size_bytes": assembly.stat().st_size,
                "selected_assembly_md5": md5sum(assembly),
                "analysis_label": reporting_name,
                "species_group": reporting_name.replace(" ", "_"),
                "taxonomic_caution": (
                    "route_specific_type_anchor_boundary; current_reporting_name_stable"
                    if row["final_rescue_decision"].endswith("type_boundary_caution")
                    else "type_anchor_group_requires_public_context_confirmation"
                ),
                "checkm2_completeness": row["meta_checkm2_completeness"],
                "checkm2_contamination": row["meta_checkm2_contamination"],
                "rescue_decision": row["final_rescue_decision"],
            }
        )

    cohort.sort(key=lambda row: str(row["sample_id"]))
    if len(cohort) != 21 or len({str(row["sample_id"]) for row in cohort}) != 21:
        raise ValueError(f"Expected 21 unique genomes, found {len(cohort)}")

    ledger: list[dict[str, object]] = []
    for row in rescue:
        ledger.append(
            {
                "sample_id": row["sample_id"],
                "evaluated_stage": "full_read_extension",
                "strict_recruitment_rate_percent": row["strict_recruitment_rate_percent"],
                "strict_mixed_sites_per_mbp": row["strict_mixed_sites_per_mbp"],
                "meta_mixed_sites_per_mbp": row["meta_mixed_sites_per_mbp"],
                "reporting_name_or_anchor": row["type_anchor_reporting_name"],
                "final_rescue_decision": row["final_rescue_decision"],
            }
        )
    mix16 = pilot["Mix16"]
    ledger.append(
        {
            "sample_id": "Mix16",
            "evaluated_stage": "locked_three_sample_pilot",
            "strict_recruitment_rate_percent": mix16["strict_recruitment_rate_percent"],
            "strict_mixed_sites_per_mbp": mix16["strict_mixed_sites_per_mbp"],
            "meta_mixed_sites_per_mbp": mix16["meta_mixed_sites_per_mbp"],
            "reporting_name_or_anchor": mix16["best_type_anchor_label"],
            "final_rescue_decision": "unrecoverable_residual_mac_mixture",
        }
    )
    ledger.sort(key=lambda row: str(row["sample_id"]))
    if len(ledger) != 11:
        raise ValueError(f"Expected 11 rescue candidates, found {len(ledger)}")

    write_tsv(output_dir / "updated_interpretable_21_manifest.tsv", cohort)
    write_tsv(output_dir / "rescue_outcomes_all11.tsv", ledger)
    (output_dir / "updated_interpretable_21_fasta_list.txt").write_text(
        "\n".join(str(row["selected_assembly_path"]) for row in cohort) + "\n"
    )
    print(f"direct={len(direct)} rescued={len(rescued_rows)} cohort={len(cohort)} candidates={len(ledger)}")


if __name__ == "__main__":
    main()
