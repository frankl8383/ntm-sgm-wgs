#!/usr/bin/env python3
"""Filter short contigs and report assembly statistics."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def records(path: Path):
    name = None
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:]
                seq = []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def n50(lengths: list[int]) -> int:
    half = sum(lengths) / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--route", required=True)
    parser.add_argument("--min-length", type=int, default=500)
    args = parser.parse_args()

    kept = [(name, seq) for name, seq in records(Path(args.input)) if len(seq) >= args.min_length]
    if not kept:
        raise ValueError("No contigs passed the length threshold")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as out:
        for name, seq in kept:
            out.write(f">{name}\n")
            for start in range(0, len(seq), 80):
                out.write(seq[start : start + 80] + "\n")

    lengths = [len(seq) for _, seq in kept]
    total = sum(lengths)
    gc = sum(seq.upper().count("G") + seq.upper().count("C") for _, seq in kept)
    row = {
        "sample_id": args.sample,
        "route": args.route,
        "min_contig_length": args.min_length,
        "genome_size_bp": total,
        "contigs": len(lengths),
        "n50_bp": n50(lengths),
        "longest_contig_bp": max(lengths),
        "gc_percent": 100 * gc / total,
    }
    summary = Path(args.summary)
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
