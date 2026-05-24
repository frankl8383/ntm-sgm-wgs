#!/usr/bin/env python3
"""Run CheckM2 with a macOS-compatible multiprocessing start method."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--extension", default=".fasta")
    parser.add_argument("--output-directory", required=True)
    parser.add_argument("--database-path", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--lowmem", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--remove-intermediates", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    mp.set_start_method("fork", force=True)

    from checkm2.main import main as checkm2_main

    sys.argv = [
        "checkm2",
        "predict",
        "--threads",
        str(args.threads),
        "--input",
        args.input,
        "--extension",
        args.extension,
        "--output-directory",
        args.output_directory,
        "--database_path",
        args.database_path,
    ]
    if args.lowmem:
        sys.argv.append("--lowmem")
    if args.force:
        sys.argv.append("--force")
    if args.remove_intermediates:
        sys.argv.append("--remove_intermediates")
    return int(checkm2_main() or 0)


if __name__ == "__main__":
    raise SystemExit(main())
