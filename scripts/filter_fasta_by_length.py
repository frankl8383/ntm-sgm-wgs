#!/usr/bin/env python3
"""Filter FASTA records by minimum sequence length."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-length", type=int, default=500)
    return parser.parse_args()


def fasta_records(path: Path):
    name = None
    seq_parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts)
                name = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if name is not None:
        yield name, "".join(seq_parts)


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    with output_path.open("w", encoding="utf-8") as out:
        for name, seq in fasta_records(input_path):
            if len(seq) < args.min_length:
                continue
            kept += 1
            out.write(f"{name}\n{wrap(seq)}\n")
    if kept == 0:
        raise SystemExit(f"No FASTA records >= {args.min_length} bp in {input_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
