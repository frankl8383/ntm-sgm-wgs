#!/usr/bin/env python3
"""Build MAC public-context genome sets for local high-confidence SGM isolates."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import pandas as pd


SPECIES_GROUPS = {
    "Mycobacterium avium": "M_avium",
    "Mycobacterium intracellulare": "M_intracellulare",
    "Mycobacterium paraintracellulare": "M_paraintracellulare",
    "Mycobacterium colombiense": "M_colombiense",
}


def biosample_attribute(row: dict[str, Any], name: str) -> str:
    biosample = (row.get("assemblyInfo") or {}).get("biosample") or {}
    for attr in biosample.get("attributes") or []:
        if attr.get("name") == name:
            return str(attr.get("value", "NA") or "NA")
    return "NA"


def organism_name(row: dict[str, Any]) -> str:
    organism = row.get("organism") or {}
    if organism.get("organismName"):
        return str(organism["organismName"])
    biosample = (row.get("assemblyInfo") or {}).get("biosample") or {}
    description = biosample.get("description") or {}
    organism = description.get("organism") or {}
    return str(organism.get("organismName", "NA") or "NA")


def strain_name(row: dict[str, Any]) -> str:
    organism = row.get("organism") or {}
    names = organism.get("infraspecificNames") or organism.get("infraspecific_names") or {}
    for key in ("strain", "isolate"):
        if names.get(key):
            return str(names[key])
    for key in ("strain", "isolate", "Sample name"):
        value = biosample_attribute(row, key)
        if value != "NA":
            return value
    return "NA"


def find_fasta(data_dir: Path, accession: str) -> Path | None:
    candidates = sorted((data_dir / accession).glob("*_genomic.fna"))
    if candidates:
        return candidates[0]
    candidates = sorted(data_dir.glob(f"{accession}/**/*_genomic.fna"))
    return candidates[0] if candidates else None


def source_dataset(report: Path) -> str:
    parts = report.parts
    if "extracted" in parts:
        index = parts.index("extracted")
        if index + 1 < len(parts):
            return parts[index + 1]
    return report.parents[2].name


def safe_name(*parts: str) -> str:
    cleaned = "_".join(str(part) for part in parts if str(part))
    cleaned = cleaned.replace("[", "").replace("]", "")
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in cleaned)
    while "__" in cleaned:
        cleaned = cleaned.replace("__", "_")
    return cleaned.strip("_")


def load_public_genomes(datasets_root: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for report in sorted(datasets_root.glob("*/ncbi_dataset/data/assembly_data_report.jsonl")):
        data_dir = report.parent
        dataset = source_dataset(report)
        for line in report.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            accession = row.get("accession") or row.get("currentAccession") or row.get("current_accession")
            if not accession:
                continue
            fasta = find_fasta(data_dir, accession)
            if fasta is None:
                continue
            name = organism_name(row)
            assembly_info = row.get("assemblyInfo") or {}
            assembly_stats = row.get("assemblyStats") or row.get("assembly_stats") or {}
            checkm_info = row.get("checkmInfo") or row.get("checkm_info") or {}
            type_material = row.get("typeMaterial") or row.get("type_material") or {}
            rows.append(
                {
                    "sample_id": accession,
                    "sample_type": "public",
                    "accession": accession,
                    "species": name,
                    "species_group": SPECIES_GROUPS.get(name, "MAC_other_or_unresolved"),
                    "strain": strain_name(row),
                    "assembly_level": assembly_info.get("assemblyLevel")
                    or assembly_info.get("assembly_level")
                    or "NA",
                    "source_database": row.get("sourceDatabase") or row.get("source_database") or "NA",
                    "source_dataset": dataset,
                    "type_material": json.dumps(type_material, ensure_ascii=False) if type_material else "NA",
                    "checkm_completeness": checkm_info.get("completeness", "NA"),
                    "checkm_contamination": checkm_info.get("contamination", "NA"),
                    "genome_size": assembly_stats.get("totalSequenceLength")
                    or assembly_stats.get("total_sequence_length")
                    or "NA",
                    "gc_percent": assembly_stats.get("gcPercent") or assembly_stats.get("gc_percent") or "NA",
                    "fasta": str(fasta),
                }
            )
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows).drop_duplicates(subset=["accession"], keep="first")
    return df


def load_local(local_report: Path) -> pd.DataFrame:
    df = pd.read_csv(local_report, sep="\t", keep_default_na=False)
    df = df[df["exclude_from_downstream_reason"].isin(["NA", ""])]
    rows: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        sample_id = row["sample_id"]
        rows.append(
            {
                "sample_id": sample_id,
                "sample_type": "local",
                "accession": sample_id,
                "species": row["final_wgs_species_call"],
                "species_group": SPECIES_GROUPS.get(row["final_wgs_species_call"], "MAC_other_or_unresolved"),
                "strain": sample_id,
                "assembly_level": "local_assembly",
                "source_database": "local",
                "source_dataset": "local_high_confidence_sgm",
                "type_material": "NA",
                "checkm_completeness": row.get("checkm2_completeness", "NA"),
                "checkm_contamination": row.get("checkm2_contamination", "NA"),
                "genome_size": row.get("assembly_size", "NA"),
                "gc_percent": row.get("gc_content", "NA"),
                "fasta": f"results/assemblies/{sample_id}/{sample_id}.assembly.fasta",
            }
        )
    return pd.DataFrame(rows)


def build_symlink_dir(metadata: pd.DataFrame, fasta_dir: Path) -> pd.DataFrame:
    fasta_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, row in metadata.iterrows():
        link_name = safe_name(row["sample_id"], row["species"], row["strain"]) + ".fna"
        link = fasta_dir / link_name
        if link.exists() or link.is_symlink():
            link.unlink()
        fasta = Path(str(row["fasta"]))
        os.symlink(os.path.relpath(fasta.resolve(), link.parent), link)
        new_row = row.to_dict()
        new_row["context_fasta"] = str(link)
        rows.append(new_row)
    return pd.DataFrame(rows)


def select_tier1(local: pd.DataFrame, fastani_hits: pd.DataFrame, public: pd.DataFrame, n: int) -> pd.DataFrame:
    selected: list[pd.DataFrame] = []
    public_by_accession = public.set_index("accession", drop=False)
    for sample_id in local["sample_id"]:
        hits = fastani_hits[fastani_hits["sample_id"] == sample_id].sort_values(
            ["ani", "alignment_fraction"], ascending=False
        )
        hits = hits[hits["accession"].isin(public_by_accession.index)].head(n)
        if hits.empty:
            continue
        subset = public_by_accession.loc[hits["accession"]].copy()
        subset["nearest_local_sample_id"] = sample_id
        subset["nearest_local_rank"] = range(1, len(subset) + 1)
        subset["nearest_local_ani"] = hits["ani"].to_numpy()
        subset["nearest_local_alignment_fraction"] = hits["alignment_fraction"].to_numpy()
        selected.append(subset)
    if not selected:
        return pd.DataFrame()
    tier1 = pd.concat(selected, ignore_index=True)
    tier1 = tier1.sort_values(["nearest_local_sample_id", "nearest_local_rank"])
    return tier1


def select_tier2(public: pd.DataFrame, local_species: list[str], max_per_species: int) -> pd.DataFrame:
    rows = []
    public = public[public["species"].isin(local_species)].copy()
    public["_type_priority"] = public["type_material"].map(lambda v: 0 if str(v) != "NA" else 1)
    public["_level_priority"] = public["assembly_level"].map(
        lambda v: {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}.get(str(v), 4)
    )
    public["_contamination"] = pd.to_numeric(public["checkm_contamination"], errors="coerce").fillna(999)
    public["_completeness"] = pd.to_numeric(public["checkm_completeness"], errors="coerce").fillna(0)
    for species, group in public.groupby("species"):
        chosen = group.sort_values(
            ["_type_priority", "_level_priority", "_contamination", "_completeness", "accession"],
            ascending=[True, True, True, False, True],
        ).head(max_per_species)
        rows.append(chosen)
    if not rows:
        return pd.DataFrame()
    tier2 = pd.concat(rows, ignore_index=True)
    return tier2.drop(columns=[c for c in tier2.columns if c.startswith("_")])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets-root", required=True, type=Path)
    parser.add_argument("--local-report", required=True, type=Path)
    parser.add_argument("--fastani-hits", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--fasta-dir", required=True, type=Path)
    parser.add_argument("--tier1-neighbors-per-local", type=int, default=20)
    parser.add_argument("--tier2-max-per-species", type=int, default=60)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    public = load_public_genomes(args.datasets_root)
    local = load_local(args.local_report)
    fastani_hits = pd.read_csv(args.fastani_hits, sep="\t")

    public = public[public["species"].isin(SPECIES_GROUPS)].copy()
    tier1 = select_tier1(local, fastani_hits, public, args.tier1_neighbors_per_local)
    tier2 = select_tier2(public, sorted(local["species"].unique()), args.tier2_max_per_species)

    combined_ids = set(local["sample_id"]) | set(tier1.get("sample_id", [])) | set(tier2.get("sample_id", []))
    combined = pd.concat([local, public[public["sample_id"].isin(combined_ids - set(local["sample_id"]))]])
    combined = combined.drop_duplicates(subset=["sample_id"], keep="first").sort_values(
        ["sample_type", "species", "sample_id"]
    )
    combined = build_symlink_dir(combined, args.fasta_dir)

    local.to_csv(args.outdir / "local_high_confidence_mac_sgm_for_public_context.tsv", sep="\t", index=False)
    public.to_csv(args.outdir / "public_mac_candidate_metadata.tsv", sep="\t", index=False)
    tier1.to_csv(args.outdir / "public_tier1_nearest_neighbors.tsv", sep="\t", index=False)
    tier2.to_csv(args.outdir / "public_tier2_species_background.tsv", sep="\t", index=False)
    combined.to_csv(args.outdir / "mac_public_context_metadata.tsv", sep="\t", index=False)
    (args.outdir / "mac_public_context_fasta_list.txt").write_text(
        "\n".join(combined["context_fasta"].astype(str)) + "\n", encoding="utf-8"
    )
    summary = pd.DataFrame(
        [
            {"table": "local_high_confidence", "n": len(local)},
            {"table": "public_candidates", "n": len(public)},
            {"table": "tier1_rows_with_duplicates_by_local", "n": len(tier1)},
            {"table": "tier1_unique_public", "n": tier1["sample_id"].nunique() if not tier1.empty else 0},
            {"table": "tier2_unique_public", "n": tier2["sample_id"].nunique() if not tier2.empty else 0},
            {"table": "combined_context_genomes", "n": len(combined)},
        ]
    )
    summary.to_csv(args.outdir / "mac_public_context_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
