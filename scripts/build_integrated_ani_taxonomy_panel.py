#!/usr/bin/env python3
"""Build a local + type-strain + public-neighbor panel for integrated ANI taxonomy."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--priority-final", required=True, type=Path)
    parser.add_argument("--species-confidence", required=True, type=Path)
    parser.add_argument("--type-metadata", required=True, type=Path)
    parser.add_argument("--public-metadata", required=True, type=Path)
    parser.add_argument("--local-vs-public-all-hits", required=True, type=Path)
    parser.add_argument("--top-public-per-local", type=int, default=5)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def clean_label(value: str, max_len: int = 80) -> str:
    text = re.sub(r"\s+", " ", str(value)).strip()
    text = re.sub(r"[^A-Za-z0-9_.=:+/ -]+", "", text)
    return text[:max_len].strip() or "NA"


def short_species(value: str) -> str:
    text = str(value)
    replacements = {
        "Mycobacterium avium": "M. avium",
        "Mycobacterium intracellulare": "M. intracellulare",
        "Mycobacterium paraintracellulare": "M. paraintracellulare",
        "Mycobacterium colombiense": "M. colombiense",
        "Mycobacterium timonense": "M. timonense",
        "Mycobacterium bouchedurhonense": "M. bouchedurhonense",
        "Mycobacterium yongonense": "M. yongonense",
        "Mycobacterium chimaera": "M. chimaera",
        "Mycobacterium marseillense": "M. marseillense",
        "Mycobacterium arosiense": "M. arosiense",
        "Mycobacterium xenopi": "M. xenopi",
        "Mycobacterium indicus pranii": "M. indicus pranii",
    }
    for key, val in replacements.items():
        if key in text:
            return val
    return text.replace("Mycobacterium ", "M. ")


def make_fasta_list(df: pd.DataFrame, path: Path) -> None:
    with path.open("w") as handle:
        for fasta in df["fasta"]:
            handle.write(f"{fasta}\n")


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    priority = read_tsv(args.priority_final)
    species_conf = read_tsv(args.species_confidence)
    type_meta = read_tsv(args.type_metadata)
    public_meta = read_tsv(args.public_metadata)
    public_hits = read_tsv(args.local_vs_public_all_hits)

    conf_lookup = species_conf.set_index("sample_id").to_dict("index")

    rows: list[dict[str, object]] = []
    for _, row in priority.iterrows():
        sample_id = str(row["sample_id"])
        conf = conf_lookup.get(sample_id, {})
        fasta = Path("results/assemblies") / sample_id / f"{sample_id}.assembly.fasta"
        rows.append(
            {
                "node_id": f"LOCAL__{sample_id}",
                "display_label": sample_id,
                "sample_id": sample_id,
                "accession": sample_id,
                "species": conf.get("analysis_group_used_for_downstream", row.get("final_wgs_species_call", "NA")),
                "species_short": short_species(conf.get("analysis_group_used_for_downstream", row.get("final_wgs_species_call", "NA"))),
                "strain": sample_id,
                "node_role": "local_downstream" if conf.get("downstream_included", False) is True else "local_qc_warning",
                "analysis_group": conf.get("analysis_group_used_for_downstream", row.get("final_wgs_species_call", "NA")),
                "species_confidence_tier": conf.get("species_level_confidence", conf.get("species_level_confidence_tier", "NA")),
                "strict_species_claim_allowed": conf.get("strict_species_claim_allowed", False),
                "source": "local",
                "fasta": str(fasta),
                "selection_reason": "priority14_local_isolate",
            }
        )

    for _, row in type_meta.iterrows():
        accession = str(row["accession"])
        species = str(row.get("requested_taxon_name", "")).strip()
        if not species:
            species = str(row["organism_name"]).split(" strain ")[0]
        rows.append(
            {
                "node_id": f"TYPE__{accession}",
                "display_label": f"{short_species(species)} {clean_label(row.get('strain', accession), 24)}",
                "sample_id": accession,
                "accession": accession,
                "species": species,
                "species_short": short_species(species),
                "strain": row.get("strain", "NA"),
                "node_role": "type_or_representative",
                "analysis_group": row.get("panel_group", "type_panel"),
                "species_confidence_tier": "reference_anchor",
                "strict_species_claim_allowed": True,
                "source": "type_panel",
                "fasta": row["reference_fasta"],
                "selection_reason": "v2_rescued_type_strain_panel",
            }
        )

    selected_public_ids: set[str] = set()
    if not public_hits.empty:
        public_hits = public_hits.sort_values(["sample_id", "ani", "alignment_fraction"], ascending=[True, False, False])
        top_hits = public_hits.groupby("sample_id").head(args.top_public_per_local)
        selected_public_ids.update(top_hits["reference_sample_id"].astype(str))
    # Ensure the best public hit stored in the species-confidence table is retained.
    if "public_final_best_accession" in species_conf.columns:
        selected_public_ids.update(species_conf["public_final_best_accession"].astype(str))

    type_accessions = {str(x) for x in type_meta["accession"]}
    public_sub = public_meta[
        public_meta["sample_type"].eq("public") & public_meta["sample_id"].astype(str).isin(selected_public_ids)
    ].copy()
    for _, row in public_sub.iterrows():
        accession = str(row["sample_id"])
        if accession in type_accessions:
            continue
        species = str(row["species"])
        rows.append(
            {
                "node_id": f"PUBLIC__{accession}",
                "display_label": f"{short_species(species)} {clean_label(row.get('strain', accession), 24)}",
                "sample_id": accession,
                "accession": accession,
                "species": species,
                "species_short": short_species(species),
                "strain": row.get("strain", "NA"),
                "node_role": "public_near_neighbor",
                "analysis_group": row.get("species_group", "public_context"),
                "species_confidence_tier": "public_context_reference",
                "strict_species_claim_allowed": True,
                "source": "public_context",
                "fasta": row["context_fasta"],
                "selection_reason": f"top_{args.top_public_per_local}_public_neighbor_or_best_hit",
            }
        )

    panel = pd.DataFrame(rows)
    panel = panel.drop_duplicates("node_id").copy()
    # Also deduplicate identical FASTA files, keeping local > type > public roles.
    role_rank = {"local_downstream": 0, "local_qc_warning": 1, "type_or_representative": 2, "public_near_neighbor": 3}
    panel["role_rank"] = panel["node_role"].map(role_rank).fillna(9).astype(int)
    panel = panel.sort_values(["role_rank", "species", "display_label"]).drop_duplicates("fasta", keep="first")
    panel = panel.drop(columns=["role_rank"]).sort_values(["node_role", "species", "display_label"]).reset_index(drop=True)
    panel["fasta_exists"] = panel["fasta"].map(lambda p: Path(p).exists())

    missing = panel[~panel["fasta_exists"]]
    if not missing.empty:
        missing.to_csv(args.outdir / "integrated_ani_panel_missing_fastas.tsv", sep="\t", index=False)
        raise SystemExit(f"Missing FASTA files in integrated panel: {args.outdir / 'integrated_ani_panel_missing_fastas.tsv'}")

    panel_path = args.outdir / "integrated_ani_panel_metadata.tsv"
    list_path = args.outdir / "integrated_ani_panel_fastani_list.txt"
    panel.to_csv(panel_path, sep="\t", index=False)
    make_fasta_list(panel, list_path)

    summary = panel.groupby(["node_role", "species_short"], dropna=False).size().reset_index(name="n")
    summary.to_csv(args.outdir / "integrated_ani_panel_summary.tsv", sep="\t", index=False)

    print(f"Wrote {panel_path} with {len(panel)} genomes")
    print(f"Wrote {list_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
