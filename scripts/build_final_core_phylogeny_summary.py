#!/usr/bin/env python3
"""Create final Figure 5 phylogeny summary with expanded M. paraintracellulare."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revised-summary", required=True, type=Path)
    parser.add_argument("--expanded-metadata", required=True, type=Path)
    parser.add_argument("--expanded-site-summary", dest="expanded_site_summary", required=True, type=Path)
    parser.add_argument("--expanded-tree", required=True, type=Path)
    parser.add_argument("--expanded-iqtree", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_iqtree_report(path: Path) -> tuple[str, float | None]:
    model = "GTR+F+ASC"
    log_likelihood = None
    text = path.read_text(errors="ignore")
    for line in text.splitlines():
        if line.startswith("Model of substitution:"):
            model = line.split(":", 1)[1].strip()
        match = re.search(r"Log-likelihood of the tree:\s+(-?\d+(?:\.\d+)?)", line)
        if match:
            log_likelihood = float(match.group(1))
    return model, log_likelihood


def main() -> int:
    args = parse_args()
    summary = pd.read_csv(args.revised_summary, sep="\t")
    expanded_meta = pd.read_csv(args.expanded_metadata, sep="\t", keep_default_na=False)
    site = pd.read_csv(args.expanded_site_summary, sep="\t").iloc[0]
    model, log_likelihood = parse_iqtree_report(args.expanded_iqtree)

    expanded_row = {
        "species_group": "M_paraintracellulare",
        "n_total_genomes": int(expanded_meta.shape[0]),
        "n_local_genomes": int(expanded_meta["sample_type"].eq("local").sum()),
        "n_public_genomes": int(expanded_meta["sample_type"].eq("public").sum()),
        "alignment_sites": int(site["alignment_sites"]),
        "alignment_sequences": int(site["n_sequences"]),
        "tree_terminals": int(expanded_meta.shape[0]),
        "iqtree_model": model,
        "log_likelihood": log_likelihood,
        "alignment_path": site["alignment_path"],
        "tree_path": str(args.expanded_tree),
        "iqtree_report": str(args.expanded_iqtree),
        "constant_sites": int(site["constant_sites"]),
        "variable_sites": int(site["variable_sites"]),
        "gap_or_ambiguous_sites": int(site["gap_or_ambiguous_sites"]),
        "constant_site_fraction": float(site["constant_site_fraction"]),
        "variable_site_fraction": float(site["variable_site_fraction"]),
        "gap_or_ambiguous_site_fraction": float(site["gap_or_ambiguous_site_fraction"]),
        "phylogeny_caution": (
            "Expanded curated public background (11 public genomes). SNP-only SKA alignment; "
            "ASC model used; do not compare branch lengths across species panels."
        ),
    }

    final = summary[summary["species_group"] != "M_paraintracellulare"].copy()
    final = pd.concat([final, pd.DataFrame([expanded_row])], ignore_index=True)
    order = {
        "M_avium": 0,
        "M_colombiense": 1,
        "M_intracellulare": 2,
        "M_paraintracellulare": 3,
    }
    final["_order"] = final["species_group"].map(order).fillna(99).astype(int)
    final = final.sort_values("_order").drop(columns=["_order"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
