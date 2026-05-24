#!/usr/bin/env python3
"""Summarize local isolate placement in a MAC public-context tree."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import Phylo


def tree_label(path_or_name: str) -> str:
    return Path(path_or_name).name.replace(".fna", "")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--nearest-output", required=True, type=Path)
    parser.add_argument("--local-distance-output", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    tree = Phylo.read(args.tree, "newick")
    metadata = pd.read_csv(args.metadata, sep="\t", keep_default_na=False)
    metadata["tree_label"] = metadata["context_fasta"].map(tree_label)
    meta = metadata.set_index("tree_label").to_dict("index")
    terminals = {tree_label(t.name): t for t in tree.get_terminals()}

    local_labels = [label for label, row in meta.items() if row.get("sample_type") == "local" and label in terminals]
    public_labels = [label for label, row in meta.items() if row.get("sample_type") == "public" and label in terminals]

    nearest_rows = []
    for local_label in local_labels:
        local_row = meta[local_label]
        distances = []
        for public_label in public_labels:
            public_row = meta[public_label]
            distances.append(
                {
                    "local_sample_id": local_row["sample_id"],
                    "local_species": local_row["species"],
                    "public_sample_id": public_row["sample_id"],
                    "public_species": public_row["species"],
                    "public_strain": public_row["strain"],
                    "public_assembly_level": public_row["assembly_level"],
                    "public_type_material": public_row["type_material"],
                    "tree_distance": tree.distance(terminals[local_label], terminals[public_label]),
                }
            )
        for rank, row in enumerate(sorted(distances, key=lambda item: item["tree_distance"])[: args.top_n], start=1):
            row["rank"] = rank
            nearest_rows.append(row)

    local_distance_rows = []
    for i, label_i in enumerate(local_labels):
        for label_j in local_labels[i + 1 :]:
            row_i = meta[label_i]
            row_j = meta[label_j]
            local_distance_rows.append(
                {
                    "sample_id_1": row_i["sample_id"],
                    "species_1": row_i["species"],
                    "sample_id_2": row_j["sample_id"],
                    "species_2": row_j["species"],
                    "tree_distance": tree.distance(terminals[label_i], terminals[label_j]),
                }
            )

    args.nearest_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(nearest_rows).to_csv(args.nearest_output, sep="\t", index=False)
    pd.DataFrame(local_distance_rows).to_csv(args.local_distance_output, sep="\t", index=False)
    print(f"Wrote {len(nearest_rows)} nearest rows: {args.nearest_output}")
    print(f"Wrote {len(local_distance_rows)} local-distance rows: {args.local_distance_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
