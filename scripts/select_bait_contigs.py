#!/usr/bin/env python3
"""Select de novo contigs supported by broad MAC bait alignments."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def fasta_records(path: Path):
    name = None
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--min-identity", type=float, default=0.90)
    parser.add_argument("--min-aligned-fraction", type=float, default=0.25)
    parser.add_argument("--min-length", type=int, default=500)
    args = parser.parse_args()

    hits: dict[str, list[dict[str, float | int | str]]] = defaultdict(list)
    with Path(args.paf).open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            target = fields[5]
            nmatch = int(fields[9])
            block = int(fields[10])
            aligned_fraction = (qend - qstart) / qlen if qlen else 0
            identity = nmatch / block if block else 0
            if identity >= args.min_identity:
                hits[query].append(
                    {
                        "qlen": qlen,
                        "qstart": qstart,
                        "qend": qend,
                        "target": target,
                        "identity": identity,
                        "nmatch": nmatch,
                        "single_aligned_fraction": aligned_fraction,
                    }
                )

    output = Path(args.output)
    manifest = Path(args.manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    with output.open("w") as out:
        for name, seq in fasta_records(Path(args.assembly)):
            contig_hits = hits.get(name, [])
            intervals = sorted((int(hit["qstart"]), int(hit["qend"])) for hit in contig_hits)
            merged: list[list[int]] = []
            for start, end in intervals:
                if not merged or start > merged[-1][1]:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(merged[-1][1], end)
            union_aligned = sum(end - start for start, end in merged)
            aligned_fraction = union_aligned / len(seq) if seq else 0
            best_hit = max(contig_hits, key=lambda hit: int(hit["nmatch"]), default=None)
            weighted_identity = (
                sum(float(hit["identity"]) * int(hit["nmatch"]) for hit in contig_hits)
                / sum(int(hit["nmatch"]) for hit in contig_hits)
                if contig_hits
                else 0
            )
            keep = bool(
                len(seq) >= args.min_length
                and contig_hits
                and aligned_fraction >= args.min_aligned_fraction
            )
            rows.append(
                {
                    "contig": name,
                    "length": len(seq),
                    "best_bait_contig": best_hit["target"] if best_hit else "",
                    "weighted_identity": f"{weighted_identity:.6f}" if best_hit else "",
                    "aligned_fraction": f"{aligned_fraction:.6f}" if best_hit else "",
                    "supporting_alignments": len(contig_hits),
                    "retained": keep,
                }
            )
            if keep:
                out.write(f">{name}\n")
                for start in range(0, len(seq), 80):
                    out.write(seq[start : start + 80] + "\n")

    with manifest.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
