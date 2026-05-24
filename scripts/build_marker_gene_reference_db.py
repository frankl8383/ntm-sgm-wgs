#!/usr/bin/env python3
"""Extract marker genes from annotated reference genomes for BLAST."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from Bio import SeqIO


MARKERS = ("rrs_16s", "rpoB", "hsp65_groEL", "sodA_candidate")


def parse_attrs(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in value.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key] = unquote(val)
    return attrs


def marker_for(feature_type: str, attrs: dict[str, str]) -> str | None:
    text = " ".join([attrs.get("gene", ""), attrs.get("Name", ""), attrs.get("product", "")]).lower()
    if feature_type.lower() == "rrna" and ("16s" in text or "rrs" in text):
        return "rrs_16s"
    if attrs.get("gene", "").lower() == "rrs" or "16s ribosomal rna" in text:
        return "rrs_16s"
    if attrs.get("gene", "").lower() == "rpob":
        return "rpoB"
    if attrs.get("gene", "").lower() in {"groel", "groel2", "hsp65"}:
        return "hsp65_groEL"
    if "chaperonin groel" in text or "molecular chaperone groel" in text:
        return "hsp65_groEL"
    if attrs.get("gene", "").lower() == "soda":
        return "sodA_candidate"
    if "mn" in text and "superoxide dismutase" in text:
        return "sodA_candidate"
    if (
        "superoxide dismutase" in text
        and "cu-zn" not in text
        and "copper" not in text
        and attrs.get("gene", "").lower() != "sodc"
    ):
        return "sodA_candidate"
    return None


def revcomp(seq: str) -> str:
    return str(seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1])


def accession_from_name(path: Path) -> str:
    match = re.match(r"(GC[AF]_\d+\.\d+)_", path.name)
    if match:
        return match.group(1)
    if re.match(r"GC[AF]_\d+\.\d+", path.parent.name):
        return path.parent.name
    raise ValueError(f"cannot parse accession from {path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets-dir", required=True, type=Path)
    parser.add_argument("--reference-metadata", required=True, type=Path)
    parser.add_argument("--fasta-output", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    args = parser.parse_args()

    ref_meta = pd.read_csv(args.reference_metadata, sep="\t").drop_duplicates(subset=["accession"]).set_index("accession")
    rows: list[dict[str, object]] = []
    records: list[tuple[str, str]] = []

    gffs = sorted(args.datasets_dir.glob("**/*_genomic.gff"))
    if not gffs:
        gffs = sorted(args.datasets_dir.glob("**/*.gff")) + sorted(args.datasets_dir.glob("**/*.gff3"))

    processed_accessions: set[str] = set()
    for gff in gffs:
        accession = accession_from_name(gff)
        if accession not in ref_meta.index or accession in processed_accessions:
            continue
        processed_accessions.add(accession)
        fasta_candidates = sorted(gff.parent.glob("*_genomic.fna"))
        if not fasta_candidates:
            continue
        genome = SeqIO.to_dict(SeqIO.parse(fasta_candidates[0], "fasta"))
        seen: dict[str, int] = {marker: 0 for marker in MARKERS}

        for line in gff.read_text(errors="replace").splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            seqid, _source, feature_type, start, end, _score, strand, _phase, attr_text = parts
            attrs = parse_attrs(attr_text)
            marker = marker_for(feature_type, attrs)
            if marker is None:
                continue
            if marker in {"rpoB", "hsp65_groEL", "sodA_candidate"} and feature_type != "CDS":
                continue
            if marker == "rrs_16s" and feature_type.lower() not in {"rrna", "gene"}:
                continue
            seq = str(genome[seqid].seq[int(start) - 1 : int(end)])
            if strand == "-":
                seq = revcomp(seq)
            seen[marker] += 1
            ref = ref_meta.loc[accession].to_dict() if accession in ref_meta.index else {}
            seq_id = f"{marker}|{accession}|{seen[marker]}"
            records.append((seq_id, seq))
            rows.append(
                {
                    "marker": marker,
                    "marker_copy": seen[marker],
                    "marker_sequence_id": seq_id,
                    "accession": accession,
                    "organism_name": ref.get("organism_name", "NA"),
                    "strain": ref.get("strain", "NA"),
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene": attrs.get("gene", "NA"),
                    "product": attrs.get("product", "NA"),
                    "length": len(seq),
                }
            )

    args.fasta_output.parent.mkdir(parents=True, exist_ok=True)
    with args.fasta_output.open("w") as out:
        for seq_id, seq in records:
            out.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")
    pd.DataFrame(rows).to_csv(args.metadata_output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
