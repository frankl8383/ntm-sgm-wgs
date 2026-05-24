#!/usr/bin/env python3
"""Summarize conservative mobilome-marker signals from Bakta annotations."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


MARKER_PATTERNS: dict[str, list[re.Pattern[str]]] = {
    "integrase": [re.compile(r"\bintegrase\b|site-specific integrase", re.I)],
    "transposase": [re.compile(r"transposase|insertion sequence|IS\d+|tnp[A-Z0-9]*", re.I)],
    "recombinase": [re.compile(r"recombinase|resolvase|invertase|xer[CD]", re.I)],
    "relaxase_mobilization": [re.compile(r"relaxase|mobilization protein", re.I)],
    "conjugation_t4ss": [re.compile(r"conjugative|conjugation|type IV secretion|\bvirB\d*\b", re.I)],
    "phage_structural": [
        re.compile(
            r"phage|prophage|capsid|tail|portal|terminase|baseplate|head protein|holin|endolysin|lysin",
            re.I,
        )
    ],
    "plasmid_associated": [
        re.compile(r"plasmid|partition protein|\bpar[AB]\b|replication initiator|\brep[AB]\b", re.I)
    ],
}


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open() as raw_handle:
        lines = []
        for line in raw_handle:
            if line.startswith("#Sequence Id"):
                lines.append(line.lstrip("#"))
            elif line.startswith("#"):
                continue
            else:
                lines.append(line)
    if not lines:
        return []
    import io

    with io.StringIO("".join(lines)) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def normalize_column(row: dict[str, str], *names: str) -> str:
    lower_map = {key.lower().replace(" ", "_").replace("-", "_"): key for key in row if key is not None}
    for name in names:
        key = lower_map.get(name.lower().replace(" ", "_").replace("-", "_"))
        if key is not None:
            return row.get(key, "")
    return ""


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    if not path.exists():
        return lengths
    name = None
    seq_len = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
    if name is not None:
        lengths[name] = seq_len
    return lengths


def load_sample_metadata(path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row["sample_id"]
            metadata[sample_id] = row
    return metadata


def load_amrfinder_counts(path: Path) -> dict[str, Counter[str]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    if not path.exists():
        return counts
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row.get("sample_id", "")
            element_class = row.get("element_class", row.get("class", ""))
            if sample_id and element_class:
                counts[sample_id][element_class] += 1
    return counts


def classify_marker(text: str) -> list[str]:
    categories = []
    for category, patterns in MARKER_PATTERNS.items():
        if any(pattern.search(text) for pattern in patterns):
            categories.append(category)
    return categories


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=Path("results/tables/table1_high_confidence_isolates.tsv"))
    parser.add_argument("--bakta-summary", type=Path, default=Path("results/tables/bakta_annotation_run_summary.tsv"))
    parser.add_argument("--assemblies-dir", type=Path, default=Path("results/assemblies"))
    parser.add_argument("--amrfinder-summary", type=Path, default=Path("results/tables/amrfinderplus_amr_virulence_stress_summary.tsv"))
    parser.add_argument("--contig-output", type=Path, default=Path("results/tables/mobilome_contig_level_calls.tsv"))
    parser.add_argument("--sample-output", type=Path, default=Path("results/tables/mobilome_feature_summary.tsv"))
    args = parser.parse_args()

    metadata = load_sample_metadata(args.metadata)
    bakta_rows = read_table(args.bakta_summary)
    amr_counts = load_amrfinder_counts(args.amrfinder_summary)

    contig_calls: list[dict[str, str]] = []
    sample_category_counts: dict[str, Counter[str]] = defaultdict(Counter)
    sample_contigs_by_category: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))

    for row in bakta_rows:
        sample_id = row.get("sample_id", "")
        if not sample_id or row.get("status") not in {"ok", "already_done"}:
            continue
        tsv_path = Path(row.get("tsv", ""))
        assembly_path = args.assemblies_dir / sample_id / f"{sample_id}.assembly.fasta"
        lengths = fasta_lengths(assembly_path)
        for feature in read_table(tsv_path):
            seq_id = normalize_column(feature, "sequence_id", "sequence id", "seq_id")
            feature_type = normalize_column(feature, "type", "feature")
            locus_tag = normalize_column(feature, "locus_tag", "locus tag")
            gene = normalize_column(feature, "gene")
            product = normalize_column(feature, "product")
            dbxrefs = normalize_column(feature, "dbxrefs", "db_xrefs")
            start = normalize_column(feature, "start")
            stop = normalize_column(feature, "stop", "end")
            search_text = " ".join([feature_type, locus_tag, gene, product, dbxrefs])
            categories = classify_marker(search_text)
            for category in categories:
                sample_category_counts[sample_id][category] += 1
                sample_contigs_by_category[sample_id][category].add(seq_id)
                contig_calls.append(
                    {
                        "sample_id": sample_id,
                        "public_context_clade": metadata.get(sample_id, {}).get("public_context_clade", ""),
                        "species_confidence_tier": metadata.get(sample_id, {}).get("species_level_confidence_tier", ""),
                        "contig_id": seq_id,
                        "contig_length": lengths.get(seq_id, ""),
                        "feature_type": feature_type,
                        "start": start,
                        "stop": stop,
                        "locus_tag": locus_tag,
                        "gene": gene,
                        "product": product,
                        "dbxrefs": dbxrefs,
                        "mobilome_marker_category": category,
                        "interpretation_confidence": "keyword_annotation_screen",
                        "interpretation_note": "Conservative marker screen from Bakta annotation; not proof of complete plasmid, prophage, or horizontal transfer.",
                    }
                )

    marker_categories = list(MARKER_PATTERNS)
    summary_rows: list[dict[str, str | int]] = []
    for sample_id in metadata:
        meta = metadata[sample_id]
        out: dict[str, str | int] = {
            "sample_id": sample_id,
            "public_context_clade": meta.get("public_context_clade", ""),
            "species_level_confidence_tier": meta.get("species_level_confidence_tier", ""),
            "strict_species_claim_allowed": meta.get("strict_species_claim_allowed", ""),
        }
        total_markers = 0
        total_contigs = set()
        for category in marker_categories:
            marker_count = sample_category_counts[sample_id][category]
            contig_count = len(sample_contigs_by_category[sample_id][category])
            out[f"{category}_marker_count"] = marker_count
            out[f"{category}_contig_count"] = contig_count
            total_markers += marker_count
            total_contigs.update(sample_contigs_by_category[sample_id][category])
        out["total_mobilome_marker_count"] = total_markers
        out["total_mobilome_marker_contig_count"] = len(total_contigs)
        out["amrfinder_amr_feature_count"] = amr_counts[sample_id].get("AMR", 0)
        out["amrfinder_stress_feature_count"] = amr_counts[sample_id].get("STRESS", 0)
        out["interpretation_note"] = "Bakta keyword screen only; use as conservative mobile-element-associated feature review."
        summary_rows.append(out)

    args.contig_output.parent.mkdir(parents=True, exist_ok=True)
    contig_fields = [
        "sample_id",
        "public_context_clade",
        "species_confidence_tier",
        "contig_id",
        "contig_length",
        "feature_type",
        "start",
        "stop",
        "locus_tag",
        "gene",
        "product",
        "dbxrefs",
        "mobilome_marker_category",
        "interpretation_confidence",
        "interpretation_note",
    ]
    with args.contig_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=contig_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(contig_calls)

    sample_fields = list(summary_rows[0]) if summary_rows else ["sample_id"]
    with args.sample_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)


if __name__ == "__main__":
    main()
