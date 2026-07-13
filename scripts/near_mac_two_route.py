#!/usr/bin/env python3
"""Run the frozen 40-condition near-MAC strict/meta reconstruction benchmark."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


PRIMARY_THRESHOLD = 161.8
AUXILIARY_THRESHOLD = 952.42
SEED = 202
PAIRS = {
    "pairA": ("Mi1", "Ma20", "calibration_cross_lineage"),
    "pairB": ("Mi2", "Mi22", "calibration_cross_lineage"),
    "pairC": ("Mi23", "Mi8", "untouched_validation_within_mp_mip"),
    "pairD": ("Mi25", "Mi4", "untouched_validation_within_mp_mip"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/near_mac_dilution/expanded_two_route"),
    )
    parser.add_argument(
        "--bait-fasta",
        type=Path,
        default=Path("results/near_mac_dilution/stage2_two_route/type_anchor_bait.fasta"),
    )
    parser.add_argument("--condition-regex", default="")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--threads-per-worker", type=int, default=8)
    parser.add_argument("--tool-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--keep-intermediates", action="store_true")
    return parser.parse_args()


def load_helpers(script_dir: Path):
    path = script_dir / "near_mac_dilution.py"
    spec = importlib.util.spec_from_file_location("near_mac_dilution", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frozen_conditions() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for pair_id, (left, right, role) in PAIRS.items():
        directions = ((left, right), (right, left))
        if pair_id in {"pairA", "pairD"}:
            for major, minor in directions:
                for total_pairs in (1_000_000, 2_000_000, 4_000_000):
                    for fraction in (0, 10, 20):
                        design = "depth_sensitivity"
                        if fraction == 10 and total_pairs == 2_000_000:
                            design = "all_pair_10_percent_and_depth_sensitivity"
                        rows.append(
                            condition_row(
                                pair_id,
                                role,
                                major,
                                minor,
                                fraction,
                                total_pairs,
                                design,
                            )
                        )
        else:
            for major, minor in directions:
                rows.append(
                    condition_row(
                        pair_id,
                        role,
                        major,
                        minor,
                        10,
                        2_000_000,
                        "all_pair_10_percent",
                    )
                )
    if len(rows) != 40:
        raise AssertionError(f"Frozen design contains {len(rows)} conditions")
    return rows


def condition_row(
    pair_id: str,
    role: str,
    major: str,
    minor: str,
    fraction: int,
    total_pairs: int,
    design_role: str,
) -> dict[str, object]:
    if fraction == 0:
        expected = "combined_negative"
    elif fraction == 10:
        expected = "auxiliary_or_primary_positive"
    else:
        expected = "primary_positive"
    condition_id = (
        f"{pair_id}__{major}_major__{minor}_minor__p{fraction:02d}__"
        f"n{total_pairs}__s{SEED}"
    )
    return {
        "condition_id": condition_id,
        "pair_id": pair_id,
        "pair_role": role,
        "major": major,
        "minor": minor,
        "minor_percent": fraction,
        "total_pairs": total_pairs,
        "seed": SEED,
        "design_role": design_role,
        "expected_rule": expected,
    }


def executable(tool_dir: Path | None, name: str) -> str:
    if tool_dir is not None:
        candidate = tool_dir / name
        if candidate.exists():
            return str(candidate)
    resolved = shutil.which(name)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found: {name}")
    return resolved


def read_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def read_fastani(path: Path) -> tuple[float, float]:
    fields = path.read_text().strip().split("\t")
    if len(fields) < 5:
        raise ValueError(f"Malformed FastANI result: {path}")
    return float(fields[2]), int(fields[3]) / int(fields[4])


def alignment_rate(path: Path) -> float:
    match = re.search(r"([0-9.]+)% overall alignment rate", path.read_text())
    if not match:
        raise ValueError(f"No alignment rate in {path}")
    return float(match.group(1))


def burden(row: dict[str, str], lower: int) -> float:
    callable_positions = int(row["callable_positions_depth_ge_20"])
    key = (
        "mixed_sites_maf_0.10_0.90"
        if lower == 10
        else "mixed_sites_maf_0.20_0.80"
    )
    return (
        int(row[key]) * 1_000_000 / callable_positions
        if callable_positions
        else 0.0
    )


def satisfies(rule: str, primary: bool, auxiliary: bool) -> bool:
    if rule == "combined_negative":
        return not primary and not auxiliary
    if rule == "auxiliary_or_primary_positive":
        return primary or auxiliary
    if rule == "primary_positive":
        return primary
    raise ValueError(rule)


def build_bait_index(
    bait: Path, prefix: Path, tool_dir: Path | None, threads: int
) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if Path(f"{prefix}.1.bt2").exists():
        return
    subprocess.run(
        [
            executable(tool_dir, "bowtie2-build"),
            "--threads",
            str(threads),
            str(bait),
            str(prefix),
        ],
        check=True,
    )


def cleanup_condition(condition_root: Path, condition_id: str) -> None:
    strict = condition_root / "strict" / condition_id
    meta = condition_root / "meta" / condition_id
    mixture = condition_root / "mixture_metrics" / condition_id
    shutil.rmtree(strict / "spades", ignore_errors=True)
    shutil.rmtree(meta / "metaspades", ignore_errors=True)
    shutil.rmtree(condition_root / "mixture_input", ignore_errors=True)
    shutil.rmtree(condition_root / "minor_inputs", ignore_errors=True)
    for path in (strict / f"{condition_id}.R1.fastq.gz", strict / f"{condition_id}.R2.fastq.gz"):
        path.unlink(missing_ok=True)
    for pattern in ("*.bam", "*.bam.bai", "*.bt2"):
        for path in mixture.glob(pattern):
            path.unlink(missing_ok=True)


def run_condition(
    condition: dict[str, object],
    windows,
    helpers,
    bait: Path,
    bait_index: Path,
    output_dir: Path,
    script_dir: Path,
    tool_dir: Path | None,
    threads: int,
    force: bool,
    keep_intermediates: bool,
) -> Path:
    condition_id = str(condition["condition_id"])
    condition_root = output_dir / "conditions" / condition_id
    result = condition_root / "condition_summary.tsv"
    if result.exists() and result.stat().st_size and not force:
        return result
    if force:
        shutil.rmtree(condition_root, ignore_errors=True)
    major = str(condition["major"])
    minor = str(condition["minor"])
    seed = int(condition["seed"])
    total_pairs = int(condition["total_pairs"])
    fraction = int(condition["minor_percent"])
    mixdir = condition_root / "mixture_input"
    mix1 = mixdir / "mixture_R1.fastq.gz"
    mix2 = mixdir / "mixture_R2.fastq.gz"
    helpers.write_mixture(
        windows[(major, seed)],
        windows[(minor, seed)],
        fraction,
        mix1,
        mix2,
        major,
        minor,
        total_pairs=total_pairs,
    )
    env = os.environ.copy()
    if tool_dir is not None:
        env["PATH"] = f"{tool_dir}:{env.get('PATH', '')}"
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
            str(total_pairs),
        ],
        check=True,
        env=env,
    )
    strict = condition_root / "strict" / condition_id
    meta = condition_root / "meta" / condition_id
    mixture = condition_root / "mixture_metrics" / condition_id
    strict_stats = read_row(strict / f"{condition_id}.assembly_stats.tsv")
    meta_stats = read_row(meta / f"{condition_id}.metaspades_mac_bin_stats.tsv")
    strict_mix = read_row(mixture / f"{condition_id}.strict.tsv")
    meta_mix = read_row(mixture / f"{condition_id}.meta.tsv")
    m2s_ani, m2s_af = read_fastani(meta / f"{condition_id}.meta_to_strict.tsv")
    s2m_ani, s2m_af = read_fastani(meta / f"{condition_id}.strict_to_meta.tsv")
    strict_primary = burden(strict_mix, 20)
    meta_primary = burden(meta_mix, 20)
    strict_auxiliary = burden(strict_mix, 10)
    meta_auxiliary = burden(meta_mix, 10)
    primary_positive = max(strict_primary, meta_primary) > PRIMARY_THRESHOLD
    auxiliary_positive = max(strict_auxiliary, meta_auxiliary) > AUXILIARY_THRESHOLD
    row: dict[str, object] = {
        **condition,
        "window_policy": "nested_prefix_of_seeded_maximum_depth_window",
        "strict_recruitment_rate_percent": alignment_rate(
            strict / "logs" / f"{condition_id}.bowtie2.log"
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
        "strict_burden_10_90_per_mbp": strict_auxiliary,
        "meta_burden_10_90_per_mbp": meta_auxiliary,
        "strict_burden_20_80_per_mbp": strict_primary,
        "meta_burden_20_80_per_mbp": meta_primary,
        "frozen_auxiliary_threshold": AUXILIARY_THRESHOLD,
        "frozen_primary_threshold": PRIMARY_THRESHOLD,
        "either_route_auxiliary_positive": auxiliary_positive,
        "either_route_primary_positive": primary_positive,
        "combined_positive": primary_positive or auxiliary_positive,
        "prespecified_rule_satisfied": satisfies(
            str(condition["expected_rule"]), primary_positive, auxiliary_positive
        ),
    }
    with result.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    if not keep_intermediates:
        cleanup_condition(condition_root, condition_id)
    return result


def write_design(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    script_dir = Path(__file__).resolve().parent
    output_dir = args.output_dir if args.output_dir.is_absolute() else project / args.output_dir
    bait = args.bait_fasta if args.bait_fasta.is_absolute() else project / args.bait_fasta
    tool_dir = args.tool_dir.resolve() if args.tool_dir else None
    output_dir.mkdir(parents=True, exist_ok=True)
    all_conditions = frozen_conditions()
    write_design(output_dir / "frozen_condition_manifest.tsv", all_conditions)
    selected = all_conditions
    if args.condition_regex:
        pattern = re.compile(args.condition_regex)
        selected = [row for row in selected if pattern.search(str(row["condition_id"]))]
    if not selected:
        raise ValueError("No conditions selected")
    helpers = load_helpers(script_dir)
    windows = helpers.prepare_nested_windows(
        project,
        all_conditions,
        output_dir / "source_windows",
        args.workers,
        args.force,
    )
    if not bait.exists():
        raise FileNotFoundError(bait)
    bait_copy = output_dir / "type_anchor_bait.fasta"
    if not bait_copy.exists():
        shutil.copy2(bait, bait_copy)
    bait_index = output_dir / "type_anchor_bait"
    build_bait_index(bait_copy, bait_index, tool_dir, args.threads_per_worker)

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                run_condition,
                row,
                windows,
                helpers,
                bait_copy,
                bait_index,
                output_dir,
                script_dir,
                tool_dir,
                args.threads_per_worker,
                args.force,
                args.keep_intermediates,
            ): row
            for row in selected
        }
        for future in as_completed(futures):
            row = futures[future]
            future.result()
            print(f"condition_completed\t{row['condition_id']}", flush=True)

    completed: list[dict[str, str]] = []
    audit: list[dict[str, object]] = []
    for row in all_conditions:
        result = output_dir / "conditions" / str(row["condition_id"]) / "condition_summary.tsv"
        status = "complete" if result.exists() and result.stat().st_size else "missing"
        audit.append(
            {
                "condition_id": row["condition_id"],
                "status": status,
                "result_path": result,
            }
        )
        if status == "complete":
            completed.append(read_row(result))
    with (output_dir / "condition_completion_audit.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit)
    if completed:
        with (output_dir / "expanded_two_route_summary.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=completed[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(completed)
    print(f"completed={len(completed)}/40")


if __name__ == "__main__":
    main()
