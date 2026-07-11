#!/usr/bin/env python3
"""Copy at most N paired records from gzipped FASTQ files."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--max-pairs", type=int, required=True)
    parser.add_argument("--output-r1", required=True)
    parser.add_argument("--output-r2", required=True)
    args = parser.parse_args()

    out1 = Path(args.output_r1)
    out2 = Path(args.output_r2)
    out1.parent.mkdir(parents=True, exist_ok=True)
    copied = 0
    with gzip.open(args.r1, "rt") as h1, gzip.open(args.r2, "rt") as h2, gzip.open(
        out1, "wt", compresslevel=5
    ) as w1, gzip.open(out2, "wt", compresslevel=5) as w2:
        while copied < args.max_pairs:
            rec1 = [h1.readline() for _ in range(4)]
            rec2 = [h2.readline() for _ in range(4)]
            if not rec1[0] and not rec2[0]:
                break
            if not rec1[0] or not rec2[0]:
                raise ValueError("FASTQ mates contain different numbers of records")
            w1.writelines(rec1)
            w2.writelines(rec2)
            copied += 1
    print(f"pairs_copied={copied}")


if __name__ == "__main__":
    main()
