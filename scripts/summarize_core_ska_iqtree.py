#!/usr/bin/env python3
"""Summarize per-species SKA alignments and IQ-TREE outputs."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import AlignIO, Phylo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-dir", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def iqtree_value(path: Path, pattern: str) -> str:
    if not path.exists():
        return "NA"
    text = path.read_text(errors="ignore")
    match = re.search(pattern, text)
    return match.group(1) if match else "NA"


def main() -> int:
    args = parse_args()
    metadata = pd.read_csv(args.metadata, sep="\t", keep_default_na=False)
    rows = []
    for group_dir in sorted(args.core_dir.glob("M_*")):
        if not group_dir.is_dir():
            continue
        group = group_dir.name
        aln_path = group_dir / f"{group}.ska_snps.aln"
        tree_path = group_dir / f"{group}.iqtree_gtr.treefile"
        report_path = group_dir / f"{group}.iqtree_gtr.iqtree"
        group_meta = metadata[metadata["species_group"] == group]
        aln = AlignIO.read(aln_path, "fasta")
        tree = Phylo.read(tree_path, "newick")
        rows.append(
            {
                "species_group": group,
                "n_total_genomes": len(group_meta),
                "n_local_genomes": int((group_meta["sample_type"] == "local").sum()),
                "n_public_genomes": int((group_meta["sample_type"] == "public").sum()),
                "alignment_sites": aln.get_alignment_length(),
                "alignment_sequences": len(aln),
                "tree_terminals": len(tree.get_terminals()),
                "iqtree_model": iqtree_value(report_path, r"Model of substitution:\s+(.+)"),
                "log_likelihood": iqtree_value(report_path, r"Log-likelihood of the tree:\s+(-?[0-9.]+)"),
                "alignment_path": str(aln_path),
                "tree_path": str(tree_path),
                "iqtree_report": str(report_path),
            }
        )
    df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
