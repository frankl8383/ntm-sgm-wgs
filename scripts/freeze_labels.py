#!/usr/bin/env python3
"""Freeze conservative analysis labels for the updated 21-genome cohort."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def classify(row: dict[str, str]) -> dict[str, str]:
    sample = row["sample_id"]
    public_species = row["best_public_reporting_species"]
    public_fine = row["best_public_fine_label"]
    type_broad = row["best_type_broad_group"]
    type_fine = row["best_type_fine_label"]
    type_fine_margin = float(row["type_distinct_fine_label_ani_margin"])

    if public_species == "Mycobacterium colombiense" and type_broad == "Mycobacterium colombiense":
        return {
            "broad_analysis_panel": "M_colombiense",
            "provisional_genomic_lineage": "M_colombiense",
            "figure_label": "M. colombiense",
            "manuscript_wording": "Mycobacterium colombiense",
            "wording_tier": "species_supported_public_and_type_anchor",
            "boundary_note": "public and type-anchor evidence concordant",
            "phylogeny_panel": "colombiense_context",
        }

    if public_species == "Mycobacterium avium":
        subgroup = (
            "M_avium_hominissuis_adjacent_public_context"
            if "hominissuis" in public_fine
            else "M_avium_public_group"
        )
        return {
            "broad_analysis_panel": "M_avium_timonense_boundary",
            "provisional_genomic_lineage": subgroup,
            "figure_label": "M. avium context",
            "manuscript_wording": "Mycobacterium avium public-context group",
            "wording_tier": "public_context_supported_anchor_boundary",
            "boundary_note": (
                "nearest public genomes are M. avium; hominissuis-adjacent context is descriptive "
                "rather than a strict subspecies assignment, and timonense-anchor proximity is retained "
                "because the anchor panel does not represent all M. avium subspecies diversity"
            ),
            "phylogeny_panel": "avium_timonense_context",
        }

    if type_broad == "Mycobacterium paraintracellulare":
        return {
            "broad_analysis_panel": "M_intracellulare_complex",
            "provisional_genomic_lineage": "MI_MP_MIP_lineage",
            "figure_label": "MI: MP-MIP",
            "manuscript_wording": "Mycobacterium intracellulare MP-MIP genomic lineage",
            "wording_tier": "lineage_wording_required",
            "boundary_note": (
                "M. paraintracellulare/M. indicus pranii-associated genomic group; "
                "M. paraintracellulare is retained only as a historical boundary label"
            ),
            "phylogeny_panel": "intracellulare_complex_context",
        }

    if "chimaera" in public_fine.lower():
        return {
            "broad_analysis_panel": "M_intracellulare_complex",
            "provisional_genomic_lineage": "MI_chimaera_adjacent_lineage",
            "figure_label": "MI: chimaera-adjacent",
            "manuscript_wording": "Mycobacterium intracellulare chimaera-adjacent public-context lineage",
            "wording_tier": "legacy_subgroup_context_only",
            "boundary_note": "public chimaera-adjacent placement with yongonense/typical-MI anchor proximity",
            "phylogeny_panel": "intracellulare_complex_context",
        }

    if "yongonense" in type_fine.lower() and type_fine_margin < 0.15:
        return {
            "broad_analysis_panel": "M_intracellulare_complex",
            "provisional_genomic_lineage": "MI_yongonense_chimaera_boundary",
            "figure_label": "MI: yongonense/chimaera boundary",
            "manuscript_wording": "Mycobacterium intracellulare yongonense/chimaera boundary lineage",
            "wording_tier": "legacy_subgroup_boundary_wording_required",
            "boundary_note": "fine-anchor margin below 0.15 percentage points",
            "phylogeny_panel": "intracellulare_complex_context",
        }

    if type_broad == "Mycobacterium intracellulare":
        return {
            "broad_analysis_panel": "M_intracellulare_complex",
            "provisional_genomic_lineage": "MI_TMI_lineage",
            "figure_label": "MI: TMI",
            "manuscript_wording": "typical Mycobacterium intracellulare (TMI) genomic lineage",
            "wording_tier": "lineage_wording_preferred",
            "boundary_note": "typical-MI public/type context; tree confirmation required",
            "phylogeny_panel": "intracellulare_complex_context",
        }

    raise SystemExit(
        f"No analysis-label rule for {sample}: public={public_species}; type={type_broad}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ani-summary", required=True)
    parser.add_argument("--local-manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = {row["sample_id"]: row for row in read_tsv(Path(args.local_manifest))}
    rows: list[dict[str, object]] = []
    for evidence in read_tsv(Path(args.ani_summary)):
        sample = evidence["sample_id"]
        decision = classify(evidence)
        source = manifest[sample]
        rows.append(
            {
                "sample_id": sample,
                "cohort_source": source["cohort_source"],
                **decision,
                "best_current_public_species": evidence["best_public_reporting_species"],
                "best_current_public_fine_label": evidence["best_public_fine_label"],
                "best_current_public_ani": evidence["best_public_ani"],
                "public_distinct_species_ani_margin": evidence["public_distinct_species_ani_margin"],
                "best_direct_type_broad_group": evidence["best_type_broad_group"],
                "best_direct_type_fine_label": evidence["best_type_fine_label"],
                "best_direct_type_ani": evidence["best_type_ani"],
                "type_distinct_broad_group_ani_margin": evidence["type_distinct_broad_group_ani_margin"],
                "type_distinct_fine_label_ani_margin": evidence["type_distinct_fine_label_ani_margin"],
                "label_status": "provisional_until_context_tree_review",
                "included_in_updated_21": "true",
            }
        )
    rows.sort(key=lambda row: str(row["sample_id"]))
    if len(rows) != 21:
        raise SystemExit(f"Expected 21 local labels, observed {len(rows)}")
    write_tsv(output_dir / "local21_analysis_labels_pre_phylogeny.tsv", rows)

    broad_counts = Counter(str(row["broad_analysis_panel"]) for row in rows)
    lineage_counts = Counter(str(row["provisional_genomic_lineage"]) for row in rows)
    report = [
        "# Updated 21-genome analysis labels before phylogenetic confirmation",
        "",
        "These labels are analysis categories, not formal taxonomic acts. Current NCBI public labels are retained as context; type-anchor proximity is used to identify boundaries rather than to overwrite current names.",
        "",
        "## Broad panels",
        "",
    ]
    report.extend(f"- {label}: {count}" for label, count in broad_counts.most_common())
    report.extend(["", "## Provisional genomic lineages", ""])
    report.extend(f"- {label}: {count}" for label, count in lineage_counts.most_common())
    report.extend(
        [
            "",
            "## Naming policy",
            "",
            "- M. paraintracellulare is not used as an unconditional current species assignment. The MP-MIP label denotes the paraintracellulare/indicus-pranii-associated genomic lineage within the M. intracellulare complex.",
            "- M. avium labels are governed by high-ANI public context; proximity to the M. timonense anchor is reported as a boundary signal because a single M. avium subsp. avium anchor does not capture M. avium subspecies diversity.",
            "- Chimaera and yongonense wording is retained only at legacy-subgroup or boundary level pending the dedicated context tree.",
            "- M. colombiense wording is permitted for Mi24 because public and direct-anchor evidence are concordant and separated from the nearest distinct anchor.",
            "",
            "## Literature basis",
            "",
            "- Comparative WGS has separated clinical M. intracellulare into typical-MI and MP-MIP genomic groups: https://pmc.ncbi.nlm.nih.gov/articles/PMC8025370/",
            "- LPSN currently lists M. paraintracellulare as a later heterotypic synonym of M. intracellulare and not the recommended medical name: https://lpsn.dsmz.de/species/mycobacterium-paraintracellulare",
            "- M. avium subspecies share high ANI and require representative subspecies context rather than a single species anchor: https://doi.org/10.3389/fmicb.2020.01701",
        ]
    )
    (output_dir / "local21_analysis_label_policy.md").write_text("\n".join(report) + "\n")
    print(
        " ".join(f"{key}={value}" for key, value in broad_counts.most_common())
        + " | "
        + " ".join(f"{key}={value}" for key, value in lineage_counts.most_common())
    )


if __name__ == "__main__":
    main()
