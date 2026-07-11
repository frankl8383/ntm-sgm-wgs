#!/usr/bin/env python3
"""Run predeclared near-MAC conditions through both recovery routes.

Stage 1 estimates dilution curves against clean major-source assemblies. This
stage checks 0%, 20% and 50% calibration conditions and reciprocal 20%
validation conditions with the complete bait-recruitment and metaSPAdes-bin
routes.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


FROZEN_PRIMARY_THRESHOLD = 161.8
PAIR_DEFINITIONS = {
    "pairA": ("Mi1", "Ma20", "calibration_cross_lineage"),
    "pairB": ("Mi2", "Mi22", "calibration_cross_lineage"),
    "pairC": ("Mi23", "Mi8", "untouched_validation_within_mp_mip"),
    "pairD": ("Mi25", "Mi4", "untouched_validation_within_mp_mip"),
}
SELECTED_SEED = 202


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--pairs", default="pairA,pairB,pairC,pairD")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--threads-per-worker", type=int, default=8)
    parser.add_argument("--tool-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def load_stage1_helpers(script_dir: Path):
    path = script_dir / "near_mac_dilution.py"
    spec = importlib.util.spec_from_file_location("nearmac_stage1", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def read_fastani(path: Path) -> tuple[float, float]:
    fields = path.read_text().strip().split("\t")
    return float(fields[2]), int(fields[3]) / int(fields[4])


def alignment_rate(path: Path) -> float:
    match = re.search(r"([0-9.]+)% overall alignment rate", path.read_text())
    if not match:
        raise ValueError(f"No alignment rate in {path}")
    return float(match.group(1))


def burden_10_90(row: dict[str, str]) -> float:
    callable_positions = int(row["callable_positions_depth_ge_20"])
    sites = int(row["mixed_sites_maf_0.10_0.90"])
    return sites * 1_000_000 / callable_positions if callable_positions else 0.0


def executable(tool_dir: Path | None, name: str) -> str:
    if tool_dir is not None:
        candidate = tool_dir / name
        if candidate.exists():
            return str(candidate)
    resolved = shutil.which(name)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found: {name}")
    return resolved


def build_bait_index(bait: Path, prefix: Path, tool_env: Path | None, threads: int) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if Path(f"{prefix}.1.bt2").exists():
        return
    subprocess.run(
        [executable(tool_env, "bowtie2-build"), "--threads", str(threads), str(bait), str(prefix)],
        check=True,
    )


def run_condition(
    condition: dict[str, object],
    windows,
    stage1,
    bait: Path,
    bait_index: Path,
    outroot: Path,
    script_dir: Path,
    tool_env: Path | None,
    threads: int,
    force: bool,
) -> Path:
    pair_id = str(condition["pair_id"])
    major = str(condition["major"])
    minor = str(condition["minor"])
    fraction = int(condition["minor_percent"])
    seed = int(condition["seed"])
    condition_id = f"{pair_id}__{major}_major__{minor}_minor__p{fraction:02d}__s{seed}"
    condition_root = outroot / "conditions" / condition_id
    result = condition_root / "condition_summary.tsv"
    if result.exists() and result.stat().st_size > 0 and not force:
        return result
    if force:
        shutil.rmtree(condition_root, ignore_errors=True)

    temp = condition_root / "mixture_input"
    mix1 = temp / "mixture_R1.fastq.gz"
    mix2 = temp / "mixture_R2.fastq.gz"
    stage1.write_mixture(
        windows[(major, seed)],
        windows[(minor, seed)],
        fraction,
        mix1,
        mix2,
        major,
        minor,
    )

    condition_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    if tool_env is not None:
        env["PATH"] = f"{tool_env}:{env.get('PATH', '')}"
    subprocess.run(
        [
            str(script_dir / "run_recovery.sh"),
            condition_id,
            str(mix1),
            str(mix2),
            str(bait_index),
            str(bait),
            str(condition_root),
            str(threads),
        ],
        check=True,
        env=env,
    )

    strict_dir = condition_root / "strict"
    meta_dir = condition_root / "meta"
    mixture_dir = condition_root / "mixture"
    strict_fasta = strict_dir / f"{condition_id}.strict.min500.fasta"
    unfiltered_meta_fasta = meta_dir / f"{condition_id}.meta_mac_bin.fasta"
    meta_fasta = meta_dir / f"{condition_id}.meta_mac_bin.min500.fasta"
    subprocess.run(
        [
            sys.executable,
            str(script_dir / "filter_assembly.py"),
            "--input",
            str(unfiltered_meta_fasta),
            "--output",
            str(meta_fasta),
            "--summary",
            str(meta_dir / "assembly_stats.tsv"),
            "--sample",
            condition_id,
            "--route",
            "meta",
            "--min-length",
            "500",
        ],
        check=True,
        env=env,
    )

    meta_to_strict = condition_root / "meta_to_strict.tsv"
    strict_to_meta = condition_root / "strict_to_meta.tsv"
    shutil.copy2(condition_root / f"{condition_id}.meta_to_strict.tsv", meta_to_strict)
    subprocess.run(
        [executable(tool_env, "fastANI"), "-q", str(strict_fasta), "-r", str(meta_fasta), "-o", str(strict_to_meta)],
        check=True,
        env=env,
    )

    strict_stats = read_row(strict_dir / f"{condition_id}.assembly.tsv")
    meta_stats = read_row(meta_dir / "assembly_stats.tsv")
    strict_mix = read_row(mixture_dir / f"{condition_id}.strict.tsv")
    meta_mix = read_row(mixture_dir / f"{condition_id}.meta.tsv")
    m2s_ani, m2s_af = read_fastani(meta_to_strict)
    s2m_ani, s2m_af = read_fastani(strict_to_meta)
    strict_burden = float(strict_mix["mixed_sites_20_80_per_mbp_callable"])
    meta_burden = float(meta_mix["mixed_sites_20_80_per_mbp_callable"])
    row = {
        **condition,
        "condition_id": condition_id,
        "selected_stage": "two_route_operating_point",
        "strict_recruitment_rate_percent": alignment_rate(
            strict_dir / f"{condition_id}.bowtie2.log"
        ),
        "strict_genome_size_bp": strict_stats["genome_size_bp"],
        "meta_genome_size_bp": meta_stats["genome_size_bp"],
        "strict_contigs": strict_stats["contigs"],
        "meta_contigs": meta_stats["contigs"],
        "strict_n50_bp": strict_stats["n50_bp"],
        "meta_n50_bp": meta_stats["n50_bp"],
        "strict_gc_percent": strict_stats["gc_percent"],
        "meta_gc_percent": meta_stats["gc_percent"],
        "meta_to_strict_ani": m2s_ani,
        "strict_to_meta_ani": s2m_ani,
        "meta_to_strict_af": m2s_af,
        "strict_to_meta_af": s2m_af,
        "strict_burden_10_90_per_mbp": burden_10_90(strict_mix),
        "meta_burden_10_90_per_mbp": burden_10_90(meta_mix),
        "strict_burden_20_80_per_mbp": strict_burden,
        "meta_burden_20_80_per_mbp": meta_burden,
        "frozen_primary_threshold": FROZEN_PRIMARY_THRESHOLD,
        "either_route_primary_positive": max(strict_burden, meta_burden)
        > FROZEN_PRIMARY_THRESHOLD,
    }
    with result.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    shutil.rmtree(temp, ignore_errors=True)
    return result


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    script_dir = Path(__file__).resolve().parent
    outroot = project / "results" / "near_mac_dilution" / "stage2_two_route"
    outroot.mkdir(parents=True, exist_ok=True)
    selected_pairs = [item.strip() for item in args.pairs.split(",") if item.strip()]
    unknown = set(selected_pairs) - set(PAIR_DEFINITIONS)
    if unknown:
        raise ValueError(f"Unknown pairs: {sorted(unknown)}")

    stage1 = load_stage1_helpers(script_dir)
    samples = {
        sample for pair_id in selected_pairs for sample in PAIR_DEFINITIONS[pair_id][:2]
    }
    windows = stage1.prepare_windows(
        project,
        samples,
        [SELECTED_SEED],
        project / "results" / "near_mac_dilution" / "stage1_clean_reference" / "seeded_windows",
        min(args.workers, 2),
        args.force,
    )

    tool_env = args.tool_dir.resolve() if args.tool_dir else None
    bait = outroot / "type_anchor_bait.fasta"
    if not bait.exists():
        raise FileNotFoundError(bait)
    bait_index = outroot / "type_anchor_bait"
    build_bait_index(bait, bait_index, tool_env, args.threads_per_worker)

    conditions: list[dict[str, object]] = []
    for pair_id in selected_pairs:
        left, right, role = PAIR_DEFINITIONS[pair_id]
        if pair_id in {"pairA", "pairB"}:
            selected = [(left, right, fraction) for fraction in (0, 20, 50)]
            design_role = "calibration_route_shape"
        else:
            selected = [(left, right, 20), (right, left, 20)]
            design_role = "untouched_validation_operating_point"
        for major, minor, fraction in selected:
            conditions.append(
                {
                    "pair_id": pair_id,
                    "pair_role": role,
                    "selected_design_role": design_role,
                    "major": major,
                    "minor": minor,
                    "minor_percent": fraction,
                    "seed": SELECTED_SEED,
                }
            )
    manifest = outroot / "selected_condition_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=conditions[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(conditions)

    outputs: list[Path] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(
                run_condition,
                condition,
                windows,
                stage1,
                bait,
                bait_index,
                outroot,
                script_dir,
                tool_env,
                args.threads_per_worker,
                args.force,
            ): condition
            for condition in conditions
        }
        for future in as_completed(futures):
            condition = futures[future]
            outputs.append(future.result())
            print(
                "condition_completed\t"
                + "\t".join(f"{key}={value}" for key, value in condition.items()),
                flush=True,
            )

    rows = [read_row(path) for path in sorted(outputs)]
    summary = outroot / "selected_two_route_summary.tsv"
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"summary\t{summary}")


if __name__ == "__main__":
    main()
