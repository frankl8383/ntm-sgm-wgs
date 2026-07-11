#!/usr/bin/env python3
"""Run seeded near-MAC dilution curves against clean major-source assemblies.

This is the detection-curve stage. It measures the frozen residual-mixture
statistics before selected conditions are taken through complete two-route
reconstruction.
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
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


FROZEN_THRESHOLD = 161.8
TOTAL_PAIRS = 2_000_000
PAIR_DEFINITIONS = {
    "pairA": ("Mi1", "Ma20", "calibration_cross_lineage"),
    "pairB": ("Mi2", "Mi22", "calibration_cross_lineage"),
    "pairC": ("Mi23", "Mi8", "untouched_validation_within_mp_mip"),
    "pairD": ("Mi25", "Mi4", "untouched_validation_within_mp_mip"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--pairs", default="pairA")
    parser.add_argument("--fractions", default="0,5,10,20,30,50")
    parser.add_argument("--seeds", default="101,202,303")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--threads-per-worker", type=int, default=4)
    parser.add_argument("--tool-dir", type=Path)
    parser.add_argument("--force", action="store_true")
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


def pair_key(header: str) -> str:
    key = header.split()[0]
    if key.endswith("/1") or key.endswith("/2"):
        key = key[:-2]
    return key


def stable_sample_seed(sample: str, seed: int) -> int:
    digest = hashlib.sha256(f"{sample}:{seed}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def total_sequences(project: Path, sample: str, r1: Path) -> int:
    archive = project / "results" / "qc" / "fastqc" / "clean" / f"{sample}_R1_fastqc.zip"
    if archive.exists():
        with zipfile.ZipFile(archive) as handle:
            member = next(name for name in handle.namelist() if name.endswith("fastqc_data.txt"))
            for line in handle.read(member).decode().splitlines():
                if line.startswith("Total Sequences\t"):
                    return int(line.split("\t")[1])
        raise ValueError(f"Total Sequences missing from {archive}")
    with gzip.open(r1, "rt") as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise ValueError(f"Incomplete FASTQ record in {r1}")
    return lines // 4


def copy_seeded_window(
    r1: Path,
    r2: Path,
    out1: Path,
    out2: Path,
    total: int,
    sample: str,
    seed: int,
) -> None:
    if total < TOTAL_PAIRS:
        raise ValueError(f"{sample} has only {total} pairs")
    rng = random.Random(stable_sample_seed(sample, seed))
    start = rng.randint(0, total - TOTAL_PAIRS)
    out1.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(r1, "rt") as h1, gzip.open(r2, "rt") as h2, gzip.open(
        out1, "wt", compresslevel=4
    ) as o1, gzip.open(out2, "wt", compresslevel=4) as o2:
        copied = 0
        for index in range(total):
            rec1 = [h1.readline() for _ in range(4)]
            rec2 = [h2.readline() for _ in range(4)]
            if not rec1[0] or not rec2[0]:
                raise ValueError(f"Unexpected FASTQ end for {sample}")
            if pair_key(rec1[0]) != pair_key(rec2[0]):
                raise ValueError(f"Pair mismatch for {sample} at record {index + 1}")
            if start <= index < start + TOTAL_PAIRS:
                o1.writelines(rec1)
                o2.writelines(rec2)
                copied += 1
            if copied == TOTAL_PAIRS:
                break
    if copied != TOTAL_PAIRS:
        raise ValueError(f"{sample} seed {seed}: copied {copied} pairs")
    (out1.parent / "window.tsv").write_text(
        "sample_id\tseed\ttotal_source_pairs\twindow_start_zero_based\twindow_pairs\n"
        f"{sample}\t{seed}\t{total}\t{start}\t{TOTAL_PAIRS}\n",
        encoding="utf-8",
    )


def prepare_windows(
    project: Path,
    samples: set[str],
    seeds: list[int],
    cache: Path,
    workers: int,
    force: bool,
) -> dict[tuple[str, int], tuple[Path, Path]]:
    outputs: dict[tuple[str, int], tuple[Path, Path]] = {}
    futures = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for sample in sorted(samples):
            source1 = project / "data" / "clean_fastq" / f"{sample}_R1.fastq.gz"
            source2 = project / "data" / "clean_fastq" / f"{sample}_R2.fastq.gz"
            total = total_sequences(project, sample, source1)
            for seed in seeds:
                outdir = cache / sample / f"seed_{seed}"
                out1 = outdir / f"{sample}.R1.window2m.fastq.gz"
                out2 = outdir / f"{sample}.R2.window2m.fastq.gz"
                outputs[(sample, seed)] = (out1, out2)
                if force or not out1.exists() or not out2.exists():
                    futures[
                        pool.submit(
                            copy_seeded_window,
                            source1,
                            source2,
                            out1,
                            out2,
                            total,
                            sample,
                            seed,
                        )
                    ] = (sample, seed)
        for future in as_completed(futures):
            sample, seed = futures[future]
            future.result()
            print(f"window_completed\t{sample}\tseed={seed}", flush=True)
    return outputs


def prefixed_record(record: list[str], prefix: str) -> list[str]:
    body = record[0].rstrip()[1:]
    record[0] = f"@{prefix}|{body}\n"
    return record


def append_pairs(source1: Path, source2: Path, count: int, prefix: str, out1, out2) -> None:
    if count == 0:
        return
    copied = 0
    with gzip.open(source1, "rt") as h1, gzip.open(source2, "rt") as h2:
        while copied < count:
            rec1 = [h1.readline() for _ in range(4)]
            rec2 = [h2.readline() for _ in range(4)]
            if not rec1[0] or not rec2[0]:
                raise ValueError(f"Window ended after {copied}/{count} pairs")
            if pair_key(rec1[0]) != pair_key(rec2[0]):
                raise ValueError(f"Pair mismatch while mixing {prefix}")
            out1.writelines(prefixed_record(rec1, prefix))
            out2.writelines(prefixed_record(rec2, prefix))
            copied += 1


def write_mixture(
    major_window: tuple[Path, Path],
    minor_window: tuple[Path, Path],
    minor_percent: int,
    output1: Path,
    output2: Path,
    major: str,
    minor: str,
) -> None:
    minor_count = round(TOTAL_PAIRS * minor_percent / 100)
    major_count = TOTAL_PAIRS - minor_count
    output1.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output1, "wt", compresslevel=4) as out1, gzip.open(
        output2, "wt", compresslevel=4
    ) as out2:
        append_pairs(*major_window, major_count, f"MAJOR_{major}", out1, out2)
        append_pairs(*minor_window, minor_count, f"MINOR_{minor}", out1, out2)


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def run_condition(
    condition: dict[str, object],
    windows: dict[tuple[str, int], tuple[Path, Path]],
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
    seed = int(condition["seed"])
    condition_id = f"{pair_id}__{major}_major__{minor}_minor__p{fraction:02d}__s{seed}"
    sample_out = outdir / "conditions" / condition_id
    result = sample_out / f"{condition_id}.clean_major_reference.minor_allele_burden.tsv"
    if result.exists() and result.stat().st_size > 0 and not force:
        return result

    temp = outdir / "tmp" / condition_id
    mix1 = temp / "mixture_R1.fastq.gz"
    mix2 = temp / "mixture_R2.fastq.gz"
    write_mixture(
        windows[(major, seed)],
        windows[(minor, seed)],
        fraction,
        mix1,
        mix2,
        major,
        minor,
    )

    sample_out.mkdir(parents=True, exist_ok=True)
    logs = sample_out / "logs"
    logs.mkdir(exist_ok=True)
    bam = sample_out / f"{condition_id}.sorted.bam"
    bowtie2 = executable(tool_env, "bowtie2")
    samtools = executable(tool_env, "samtools")
    python = sys.executable
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
                python,
                str(script_dir / "residual_mixture.py"),
                "--sample",
                condition_id,
                "--route",
                "clean_major_reference",
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


def build_indexes(
    project: Path,
    samples: set[str],
    outdir: Path,
    tool_env: Path | None,
    threads: int,
) -> tuple[dict[str, Path], dict[str, Path]]:
    assemblies: dict[str, Path] = {}
    indexes: dict[str, Path] = {}
    for sample in sorted(samples):
        assembly = project / "results" / "assemblies" / sample / f"{sample}.assembly.fasta"
        index = outdir / "reference_indexes" / sample / sample
        index.parent.mkdir(parents=True, exist_ok=True)
        assemblies[sample] = assembly
        indexes[sample] = index
        if not Path(f"{index}.1.bt2").exists():
            subprocess.run(
                [executable(tool_env, "bowtie2-build"), "--threads", str(threads), str(assembly), str(index)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        if not Path(f"{assembly}.fai").exists():
            subprocess.run([executable(tool_env, "samtools"), "faidx", str(assembly)], check=True)
    return assemblies, indexes


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    analysis = project
    script_dir = Path(__file__).resolve().parent
    tool_env = args.tool_dir.resolve() if args.tool_dir else None
    selected_pair_ids = [item.strip() for item in args.pairs.split(",") if item.strip()]
    if selected_pair_ids == ["all"]:
        selected_pair_ids = list(PAIR_DEFINITIONS)
    unknown = set(selected_pair_ids) - set(PAIR_DEFINITIONS)
    if unknown:
        raise ValueError(f"Unknown pairs: {sorted(unknown)}")
    fractions = [int(value) for value in args.fractions.split(",")]
    seeds = [int(value) for value in args.seeds.split(",")]
    if any(value < 0 or value > 50 for value in fractions):
        raise ValueError("Minor fractions must be between 0 and 50 percent")

    outdir = analysis / "results" / "near_mac_dilution" / "stage1_clean_reference"
    outdir.mkdir(parents=True, exist_ok=True)
    samples = {
        sample
        for pair_id in selected_pair_ids
        for sample in PAIR_DEFINITIONS[pair_id][:2]
    }
    windows = prepare_windows(
        project,
        samples,
        seeds,
        outdir / "seeded_windows",
        min(args.workers, 2),
        args.force,
    )
    assemblies, indexes = build_indexes(
        project, samples, outdir, tool_env, args.threads_per_worker
    )

    conditions: list[dict[str, object]] = []
    for pair_id in selected_pair_ids:
        left, right, role = PAIR_DEFINITIONS[pair_id]
        for major, minor in ((left, right), (right, left)):
            for fraction in fractions:
                for seed in seeds:
                    conditions.append(
                        {
                            "pair_id": pair_id,
                            "pair_role": role,
                            "major": major,
                            "minor": minor,
                            "minor_percent": fraction,
                            "seed": seed,
                        }
                    )

    outputs: dict[tuple[object, ...], Path] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {}
        for condition in conditions:
            key = tuple(condition.values())
            futures[
                pool.submit(
                    run_condition,
                    condition,
                    windows,
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
            outputs[key] = future.result()
            print(
                "condition_completed\t"
                + "\t".join(f"{name}={value}" for name, value in condition.items()),
                flush=True,
            )

    summary: list[dict[str, object]] = []
    for condition in conditions:
        key = tuple(condition.values())
        result = read_single_row(outputs[key])
        callable_positions = int(result["callable_positions_depth_ge_20"])
        mixed_10_90 = int(result["mixed_sites_maf_0.10_0.90"])
        burden_20_80 = float(result["mixed_sites_20_80_per_mbp_callable"])
        row = dict(condition)
        row.update(
            {
                "total_pairs": TOTAL_PAIRS,
                "reference_policy": "clean_major_source_assembly",
                "seed_policy": "seeded_random_start_window",
                "callable_positions_depth_ge_20": callable_positions,
                "median_effective_depth": int(result["median_effective_depth"]),
                "mixed_sites_maf_0.10_0.90": mixed_10_90,
                "mixed_sites_10_90_per_mbp_callable": (
                    mixed_10_90 * 1_000_000 / callable_positions if callable_positions else 0.0
                ),
                "mixed_sites_maf_0.20_0.80": int(result["mixed_sites_maf_0.20_0.80"]),
                "mixed_sites_20_80_per_mbp_callable": burden_20_80,
                "frozen_20_80_threshold_per_mbp": FROZEN_THRESHOLD,
                "exceeds_frozen_threshold": burden_20_80 > FROZEN_THRESHOLD,
            }
        )
        summary.append(row)
    summary_path = outdir / "nearmac_dilution_stage1_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    print(f"summary\t{summary_path}")


if __name__ == "__main__":
    main()
