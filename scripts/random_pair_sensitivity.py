#!/usr/bin/env python3
"""Run a random-pair sensitivity check from the cached near-MAC windows.

For each source, the three cached windows are merged by original FASTQ
coordinate, overlap is removed, and two million unique pairs are sampled
without replacement. Conditions at 0%, 10% and 20% minor source are then
mapped to the clean major-source assembly under the frozen burden rules.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import random
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


TOTAL_PAIRS = 2_000_000
PRIMARY_THRESHOLD = 161.8
POOL_SEED = 404
PAIR_DEFINITIONS = {
    "pairA": ("Mi1", "Ma20", "calibration_cross_lineage"),
    "pairB": ("Mi2", "Mi22", "calibration_cross_lineage"),
    "pairC": ("Mi23", "Mi8", "held_back_within_mp_mip"),
    "pairD": ("Mi25", "Mi4", "held_back_within_mp_mip"),
}
EXPECTED_PURE_CONTROLS_PER_PAIR = 6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--threads-per-worker", type=int, default=4)
    parser.add_argument("--tool-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--keep-pools", action="store_true")
    return parser.parse_args()


def executable(tool_dir: Path | None, name: str) -> str:
    if tool_dir is not None:
        candidate = tool_dir / name
        if candidate.exists():
            return str(candidate)
    resolved = shutil.which(name)
    if resolved is None:
        raise FileNotFoundError(f"Required executable not found: {name}")
    return resolved


def workflow_paths(project: Path) -> tuple[Path, Path]:
    results = project / "results"
    stage1 = results / "near_mac_dilution" / "stage1_clean_reference"
    outdir = results / "near_mac_random_pair_sensitivity"
    return stage1, outdir


def stable_seed(label: str) -> int:
    return int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")


def pair_key(header: str) -> str:
    key = header.split()[0]
    if key.endswith("/1") or key.endswith("/2"):
        key = key[:-2]
    return key


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def subtract_covered(
    interval: tuple[int, int], covered: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    pieces = [interval]
    for cover_start, cover_end in covered:
        next_pieces: list[tuple[int, int]] = []
        for start, end in pieces:
            if cover_end <= start or cover_start >= end:
                next_pieces.append((start, end))
                continue
            if start < cover_start:
                next_pieces.append((start, cover_start))
            if cover_end < end:
                next_pieces.append((cover_end, end))
        pieces = next_pieces
    return pieces


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_windows(cache: Path, sample: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for metadata in sorted((cache / sample).glob("seed_*/window.tsv")):
        with metadata.open(newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        seed = int(row["seed"])
        start = int(row["window_start_zero_based"])
        pairs = int(row["window_pairs"])
        rows.append(
            {
                "seed": seed,
                "start": start,
                "end": start + pairs,
                "total_source_pairs": int(row["total_source_pairs"]),
                "r1": metadata.parent / f"{sample}.R1.window2m.fastq.gz",
                "r2": metadata.parent / f"{sample}.R2.window2m.fastq.gz",
            }
        )
    if len(rows) != 3:
        raise ValueError(f"{sample}: expected three cached windows, found {len(rows)}")
    return rows


def selected_source_positions(
    union_intervals: list[tuple[int, int]], sample: str
) -> list[int]:
    union_pairs = sum(end - start for start, end in union_intervals)
    if union_pairs < TOTAL_PAIRS:
        raise ValueError(f"{sample}: union contains only {union_pairs} unique pairs")
    rng = random.Random(stable_seed(f"{sample}:union_pool:{POOL_SEED}"))
    offsets = sorted(rng.sample(range(union_pairs), TOTAL_PAIRS))
    selected: list[int] = []
    offset_index = 0
    union_cursor = 0
    for start, end in union_intervals:
        length = end - start
        limit = union_cursor + length
        while offset_index < len(offsets) and offsets[offset_index] < limit:
            selected.append(start + offsets[offset_index] - union_cursor)
            offset_index += 1
        union_cursor = limit
    if len(selected) != TOTAL_PAIRS:
        raise RuntimeError(f"{sample}: selected {len(selected)} positions")
    return selected


def copy_selected_pairs(
    r1: Path,
    r2: Path,
    selected_indices: set[int],
    out1,
    out2,
    sample: str,
) -> int:
    copied = 0
    with gzip.open(r1, "rt") as h1, gzip.open(r2, "rt") as h2:
        for index in range(TOTAL_PAIRS):
            rec1 = [h1.readline() for _ in range(4)]
            rec2 = [h2.readline() for _ in range(4)]
            if not rec1[0] or not rec2[0]:
                raise ValueError(f"{sample}: cached window ended at pair {index}")
            if pair_key(rec1[0]) != pair_key(rec2[0]):
                raise ValueError(f"{sample}: pair mismatch at cached pair {index}")
            if index in selected_indices:
                out1.writelines(rec1)
                out2.writelines(rec2)
                copied += 1
    return copied


def build_random_union_pool(
    cache: Path, pool_dir: Path, sample: str, force: bool
) -> dict[str, object]:
    out1 = pool_dir / sample / f"{sample}.R1.random_union_2m.fastq.gz"
    out2 = pool_dir / sample / f"{sample}.R2.random_union_2m.fastq.gz"
    metadata_path = pool_dir / sample / "pool_manifest.tsv"
    if out1.exists() and out2.exists() and metadata_path.exists() and not force:
        with metadata_path.open(newline="") as handle:
            return next(csv.DictReader(handle, delimiter="\t"))

    windows = load_windows(cache, sample)
    union_intervals = merge_intervals(
        [(int(row["start"]), int(row["end"])) for row in windows]
    )
    selected_positions = selected_source_positions(union_intervals, sample)

    covered: list[tuple[int, int]] = []
    ownership: list[tuple[int, int, dict[str, object]]] = []
    for window in sorted(windows, key=lambda row: int(row["seed"])):
        interval = (int(window["start"]), int(window["end"]))
        for start, end in subtract_covered(interval, covered):
            ownership.append((start, end, window))
        covered = merge_intervals([*covered, interval])

    by_seed: dict[int, set[int]] = {int(row["seed"]): set() for row in windows}
    owner_index = 0
    ownership = sorted(ownership, key=lambda item: item[0])
    for source_position in selected_positions:
        while not (
            ownership[owner_index][0]
            <= source_position
            < ownership[owner_index][1]
        ):
            owner_index += 1
        window = ownership[owner_index][2]
        by_seed[int(window["seed"])].add(source_position - int(window["start"]))

    out1.parent.mkdir(parents=True, exist_ok=True)
    copied = 0
    with gzip.open(out1, "wt", compresslevel=4) as h1, gzip.open(
        out2, "wt", compresslevel=4
    ) as h2:
        for window in sorted(windows, key=lambda row: int(row["seed"])):
            copied += copy_selected_pairs(
                Path(window["r1"]),
                Path(window["r2"]),
                by_seed[int(window["seed"])],
                h1,
                h2,
                sample,
            )
    if copied != TOTAL_PAIRS:
        raise RuntimeError(f"{sample}: wrote {copied} random-union pairs")

    row: dict[str, object] = {
        "sample_id": sample,
        "pool_seed": POOL_SEED,
        "source_total_pairs": int(windows[0]["total_source_pairs"]),
        "cached_windows": len(windows),
        "unique_union_pairs": sum(end - start for start, end in union_intervals),
        "sampled_unique_pairs": TOTAL_PAIRS,
        "sampling_without_replacement": True,
        "union_intervals_zero_based_half_open": ";".join(
            f"{start}-{end}" for start, end in union_intervals
        ),
        "source_window_starts": ";".join(
            f"{row['seed']}:{row['start']}" for row in sorted(windows, key=lambda x: int(x["seed"]))
        ),
        "r1_sha256": sha256(out1),
        "r2_sha256": sha256(out2),
    }
    with metadata_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    print(
        f"pool_completed\t{sample}\tunique_union={row['unique_union_pairs']}",
        flush=True,
    )
    return row


def prefixed_record(record: list[str], prefix: str) -> list[str]:
    record[0] = f"@{prefix}|{record[0].rstrip()[1:]}\n"
    return record


def write_random_subset(
    source1: Path,
    source2: Path,
    count: int,
    label: str,
    prefix: str,
    out1,
    out2,
) -> None:
    if count == 0:
        return
    rng = random.Random(stable_seed(label))
    if count <= TOTAL_PAIRS // 2:
        selected = set(rng.sample(range(TOTAL_PAIRS), count))
        include = lambda index: index in selected
    else:
        excluded = set(rng.sample(range(TOTAL_PAIRS), TOTAL_PAIRS - count))
        include = lambda index: index not in excluded

    copied = 0
    with gzip.open(source1, "rt") as h1, gzip.open(source2, "rt") as h2:
        for index in range(TOTAL_PAIRS):
            rec1 = [h1.readline() for _ in range(4)]
            rec2 = [h2.readline() for _ in range(4)]
            if not rec1[0] or not rec2[0]:
                raise ValueError(f"{prefix}: random pool ended at pair {index}")
            if pair_key(rec1[0]) != pair_key(rec2[0]):
                raise ValueError(f"{prefix}: pair mismatch at pool pair {index}")
            if include(index):
                out1.writelines(prefixed_record(rec1, prefix))
                out2.writelines(prefixed_record(rec2, prefix))
                copied += 1
    if copied != count:
        raise RuntimeError(f"{prefix}: wrote {copied}/{count} pairs")


def write_mixture(
    major_pool: tuple[Path, Path],
    minor_pool: tuple[Path, Path],
    minor_percent: int,
    output1: Path,
    output2: Path,
    condition_id: str,
    major: str,
    minor: str,
) -> None:
    minor_count = round(TOTAL_PAIRS * minor_percent / 100)
    major_count = TOTAL_PAIRS - minor_count
    output1.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output1, "wt", compresslevel=4) as out1, gzip.open(
        output2, "wt", compresslevel=4
    ) as out2:
        write_random_subset(
            *major_pool,
            major_count,
            f"{condition_id}:{major}:major",
            f"MAJOR_{major}",
            out1,
            out2,
        )
        write_random_subset(
            *minor_pool,
            minor_count,
            f"{condition_id}:{minor}:minor",
            f"MINOR_{minor}",
            out1,
            out2,
        )


def build_indexes(
    project: Path,
    samples: set[str],
    stage1: Path,
    tool_env: Path | None,
    threads: int,
) -> tuple[dict[str, Path], dict[str, Path]]:
    assemblies: dict[str, Path] = {}
    indexes: dict[str, Path] = {}
    for sample in sorted(samples):
        assembly = project / "results" / "assemblies" / sample / f"{sample}.assembly.fasta"
        index = stage1 / "reference_indexes" / sample / sample
        assemblies[sample] = assembly
        indexes[sample] = index
        if not Path(f"{index}.1.bt2").exists():
            index.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                [
                    executable(tool_env, "bowtie2-build"),
                    "--threads",
                    str(threads),
                    str(assembly),
                    str(index),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        if not Path(f"{assembly}.fai").exists():
            subprocess.run(
                [executable(tool_env, "samtools"), "faidx", str(assembly)], check=True
            )
    return assemblies, indexes


def run_condition(
    condition: dict[str, object],
    pools: dict[str, tuple[Path, Path]],
    assemblies: dict[str, Path],
    indexes: dict[str, Path],
    outdir: Path,
    script_dir: Path,
    tool_env: Path | None,
    threads: int,
    force: bool,
) -> Path:
    pair_id = str(condition["pair_id"])
    major = str(condition["major"])
    minor = str(condition["minor"])
    fraction = int(condition["minor_percent"])
    condition_id = (
        f"{pair_id}__{major}_major__{minor}_minor__p{fraction:02d}__random_union"
    )
    sample_out = outdir / "conditions" / condition_id
    result = sample_out / f"{condition_id}.clean_major_reference.minor_allele_burden.tsv"
    if result.exists() and result.stat().st_size > 0 and not force:
        return result

    temp = outdir / "tmp" / condition_id
    mix1 = temp / "mixture_R1.fastq.gz"
    mix2 = temp / "mixture_R2.fastq.gz"
    write_mixture(
        pools[major],
        pools[minor],
        fraction,
        mix1,
        mix2,
        condition_id,
        major,
        minor,
    )

    sample_out.mkdir(parents=True, exist_ok=True)
    logs = sample_out / "logs"
    logs.mkdir(exist_ok=True)
    bam = sample_out / f"{condition_id}.sorted.bam"
    bowtie2 = executable(tool_env, "bowtie2")
    samtools = executable(tool_env, "samtools")
    env = os.environ.copy()
    if tool_env is not None:
        env["PATH"] = f"{tool_env}:{env.get('PATH', '')}"
    with (logs / "bowtie2.log").open("w") as err:
        align = subprocess.Popen(
            [
                bowtie2,
                "--very-sensitive",
                "--no-unal",
                "-p",
                str(threads),
                "-x",
                str(indexes[major]),
                "-1",
                str(mix1),
                "-2",
                str(mix2),
            ],
            stdout=subprocess.PIPE,
            stderr=err,
            env=env,
        )
        assert align.stdout is not None
        sort = subprocess.run(
            [samtools, "sort", "-@", "2", "-o", str(bam), "-"],
            stdin=align.stdout,
            check=True,
            env=env,
        )
        align.stdout.close()
        if align.wait() != 0 or sort.returncode != 0:
            raise subprocess.CalledProcessError(align.returncode, "bowtie2")

    with (logs / "mpileup.log").open("w") as err:
        pileup = subprocess.Popen(
            [
                samtools,
                "mpileup",
                "-aa",
                "-q",
                "20",
                "-Q",
                "20",
                "-d",
                "100000",
                "-f",
                str(assemblies[major]),
                str(bam),
            ],
            stdout=subprocess.PIPE,
            stderr=err,
            env=env,
        )
        assert pileup.stdout is not None
        summarize = subprocess.run(
            [
                sys.executable,
                str(script_dir / "residual_mixture.py"),
                "--sample",
                condition_id,
                "--route",
                "clean_major_reference_random_union",
                "--output",
                str(result),
            ],
            stdin=pileup.stdout,
            check=True,
            env=env,
        )
        pileup.stdout.close()
        if pileup.wait() != 0 or summarize.returncode != 0:
            raise subprocess.CalledProcessError(pileup.returncode, "samtools mpileup")

    bam.unlink(missing_ok=True)
    Path(f"{bam}.bai").unlink(missing_ok=True)
    shutil.rmtree(temp, ignore_errors=True)
    return result


def read_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def secondary_threshold(stage1: Path) -> float:
    pure_controls: list[tuple[str, float]] = []
    for path in sorted(
        stage1.glob(
            "conditions/pair[AB]__*_major__*_minor__p00__s*"
            "/*.minor_allele_burden.tsv"
        )
    ):
        pair_id = path.parent.name.split("__", 1)[0]
        row = read_row(path)
        callable_positions = int(row["callable_positions_depth_ge_20"])
        if callable_positions == 0:
            raise ValueError(f"No callable positions in {path}")
        burden = (
            int(row["mixed_sites_maf_0.10_0.90"])
            * 1_000_000
            / callable_positions
        )
        pure_controls.append((pair_id, burden))
    counts = {
        pair_id: sum(row_pair == pair_id for row_pair, _ in pure_controls)
        for pair_id in ("pairA", "pairB")
    }
    if any(
        count != EXPECTED_PURE_CONTROLS_PER_PAIR
        for count in counts.values()
    ):
        raise ValueError(
            "Expected six pairA and six pairB pure controls to freeze the "
            f"auxiliary threshold; found {counts}"
        )
    return 5 * max(burden for _, burden in pure_controls)


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    stage1, outdir = workflow_paths(project)
    cache = stage1 / "seeded_windows"
    pool_dir = outdir / "random_union_pools"
    tool_env = args.tool_dir.resolve() if args.tool_dir else None
    script_dir = Path(__file__).resolve().parent
    outdir.mkdir(parents=True, exist_ok=True)

    samples = {
        sample for left, right, _ in PAIR_DEFINITIONS.values() for sample in (left, right)
    }
    pool_rows = [
        build_random_union_pool(cache, pool_dir, sample, args.force)
        for sample in sorted(samples)
    ]
    pool_manifest = outdir / "random_union_pool_manifest.tsv"
    with pool_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pool_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(pool_rows)
    pools = {
        sample: (
            pool_dir / sample / f"{sample}.R1.random_union_2m.fastq.gz",
            pool_dir / sample / f"{sample}.R2.random_union_2m.fastq.gz",
        )
        for sample in samples
    }
    assemblies, indexes = build_indexes(
        project, samples, stage1, tool_env, args.threads_per_worker
    )

    conditions: list[dict[str, object]] = []
    for pair_id, (left, right, role) in PAIR_DEFINITIONS.items():
        for major, minor in ((left, right), (right, left)):
            for fraction in (0, 10, 20):
                conditions.append(
                    {
                        "pair_id": pair_id,
                        "pair_role": role,
                        "major": major,
                        "minor": minor,
                        "minor_percent": fraction,
                    }
                )

    output_paths: dict[tuple[object, ...], Path] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {}
        for condition in conditions:
            key = tuple(condition.values())
            futures[
                pool.submit(
                    run_condition,
                    condition,
                    pools,
                    assemblies,
                    indexes,
                    outdir,
                    script_dir,
                    tool_env,
                    args.threads_per_worker,
                    args.force,
                )
            ] = (key, condition)
        for future in as_completed(futures):
            key, condition = futures[future]
            output_paths[key] = future.result()
            print(
                "condition_completed\t"
                + "\t".join(f"{key}={value}" for key, value in condition.items()),
                flush=True,
            )

    auxiliary_threshold = secondary_threshold(stage1)
    rows: list[dict[str, object]] = []
    for condition in conditions:
        result = read_row(output_paths[tuple(condition.values())])
        callable_positions = int(result["callable_positions_depth_ge_20"])
        mixed_10_90 = int(result["mixed_sites_maf_0.10_0.90"])
        burden_10_90 = (
            mixed_10_90 * 1_000_000 / callable_positions if callable_positions else 0.0
        )
        burden_20_80 = float(result["mixed_sites_20_80_per_mbp_callable"])
        primary_positive = burden_20_80 > PRIMARY_THRESHOLD
        auxiliary_positive = burden_10_90 > auxiliary_threshold
        rows.append(
            {
                **condition,
                "sampling_policy": "overlap_deduplicated_random_pairs_from_cached_window_union",
                "pool_seed": POOL_SEED,
                "total_pairs": TOTAL_PAIRS,
                "callable_positions_depth_ge_20": callable_positions,
                "median_effective_depth": int(result["median_effective_depth"]),
                "mixed_sites_10_90_per_mbp_callable": burden_10_90,
                "mixed_sites_20_80_per_mbp_callable": burden_20_80,
                "primary_20_80_threshold": PRIMARY_THRESHOLD,
                "auxiliary_10_90_threshold": auxiliary_threshold,
                "primary_positive": primary_positive,
                "auxiliary_positive": auxiliary_positive,
                "combined_positive": primary_positive or auxiliary_positive,
            }
        )
    metrics_path = outdir / "random_pair_condition_metrics.tsv"
    with metrics_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    aggregate_rows: list[dict[str, object]] = []
    for fraction in (0, 10, 20):
        subset = [row for row in rows if row["minor_percent"] == fraction]
        aggregate_rows.append(
            {
                "minor_percent": fraction,
                "conditions": len(subset),
                "primary_positive": sum(bool(row["primary_positive"]) for row in subset),
                "auxiliary_positive": sum(bool(row["auxiliary_positive"]) for row in subset),
                "combined_positive": sum(bool(row["combined_positive"]) for row in subset),
                "primary_burden_min": min(
                    float(row["mixed_sites_20_80_per_mbp_callable"]) for row in subset
                ),
                "primary_burden_max": max(
                    float(row["mixed_sites_20_80_per_mbp_callable"]) for row in subset
                ),
                "auxiliary_burden_min": min(
                    float(row["mixed_sites_10_90_per_mbp_callable"]) for row in subset
                ),
                "auxiliary_burden_max": max(
                    float(row["mixed_sites_10_90_per_mbp_callable"]) for row in subset
                ),
            }
        )
    aggregate_path = outdir / "random_pair_aggregate.tsv"
    with aggregate_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=aggregate_rows[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(aggregate_rows)

    report = [
        "# Near-MAC random-pair sensitivity",
        "",
        "Two million unique pairs per source were sampled without replacement from the overlap-deduplicated union of the three cached windows. The frozen primary and auxiliary thresholds were retained.",
        "",
    ]
    for row in aggregate_rows:
        report.append(
            f"- {row['minor_percent']}% minor source: primary "
            f"{row['primary_positive']}/{row['conditions']}, auxiliary "
            f"{row['auxiliary_positive']}/{row['conditions']}, combined "
            f"{row['combined_positive']}/{row['conditions']}."
        )
    (outdir / "README.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))

    if not args.keep_pools:
        shutil.rmtree(pool_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
