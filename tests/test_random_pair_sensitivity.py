import csv
from pathlib import Path

import scripts.random_pair_sensitivity as random_pair_sensitivity
from scripts.random_pair_sensitivity import (
    executable,
    merge_intervals,
    secondary_threshold,
    subtract_covered,
    workflow_paths,
)


def test_merge_intervals_removes_overlap():
    assert merge_intervals([(10, 20), (15, 30), (40, 50)]) == [(10, 30), (40, 50)]


def test_subtract_covered_returns_unseen_segments():
    assert subtract_covered((10, 30), [(0, 15), (20, 25)]) == [(15, 20), (25, 30)]


def test_workflow_paths_match_near_mac_producer(tmp_path):
    stage1, outdir = workflow_paths(tmp_path)
    assert stage1 == tmp_path / "results/near_mac_dilution/stage1_clean_reference"
    assert outdir == tmp_path / "results/near_mac_random_pair_sensitivity"


def test_executable_uses_configured_tool_directory(tmp_path):
    tool_dir = tmp_path / "bin"
    tool_dir.mkdir()
    tool = tool_dir / "bowtie2"
    tool.write_text("#!/bin/sh\n", encoding="utf-8")
    tool.chmod(0o755)
    assert executable(tool_dir, "bowtie2") == str(tool)


def test_checked_in_pileup_summarizer_is_available():
    script_dir = Path(random_pair_sensitivity.__file__).resolve().parent
    assert (script_dir / "residual_mixture.py").is_file()


def test_secondary_threshold_uses_pair_a_b_pure_controls(tmp_path):
    stage1 = tmp_path / "stage1"
    for pair_id, mixed_sites in (("pairA", 18), ("pairB", 19)):
        for direction in ("left", "right"):
            for seed in (101, 202, 303):
                condition = (
                    stage1
                    / "conditions"
                    / f"{pair_id}__{direction}_major__minor_minor__p00__s{seed}"
                )
                condition.mkdir(parents=True)
                row = {
                    "callable_positions_depth_ge_20": 100_000,
                    "mixed_sites_maf_0.10_0.90": mixed_sites,
                }
                with (
                    condition
                    / f"{condition.name}.minor_allele_burden.tsv"
                ).open("w", newline="") as handle:
                    writer = csv.DictWriter(
                        handle, fieldnames=row, delimiter="\t"
                    )
                    writer.writeheader()
                    writer.writerow(row)
    assert secondary_threshold(stage1) == 950.0
