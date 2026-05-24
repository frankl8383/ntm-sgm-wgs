#!/usr/bin/env python3
"""Stage BGI clean FASTQs and generate a project samplesheet."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
import tarfile
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover - handled at runtime
    yaml = None


SAMPLE_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        if yaml is not None:
            return yaml.safe_load(handle) or {}
        return parse_simple_yaml(handle.read())


def parse_scalar(value: str):
    value = value.strip()
    if value in {"", "null", "Null", "NULL", "~"}:
        return None
    if value in {"true", "True", "TRUE"}:
        return True
    if value in {"false", "False", "FALSE"}:
        return False
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def parse_simple_yaml(text: str) -> dict:
    """Small fallback parser for the simple mapping shape used by config.yaml."""
    root: dict = {}
    stack: list[tuple[int, dict]] = [(-1, root)]
    for raw_line in text.splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        if raw_line.lstrip().startswith("- "):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        line = raw_line.strip()
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if value.strip() == "":
            child: dict = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = parse_scalar(value)
    return root


def resolve_project_path(project_root: Path, value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if not path.is_absolute():
        path = project_root / path
    return path


def find_member(tar: tarfile.TarFile, sample_id: str, read_number: int) -> tarfile.TarInfo | None:
    suffixes = (
        f"/{sample_id}.IS350_Clean.{read_number}.fq.gz",
        f"/{sample_id}.Clean.{read_number}.fq.gz",
        f"/{sample_id}_{read_number}.fq.gz",
        f"/{sample_id}_R{read_number}.fastq.gz",
        f"/{sample_id}_R{read_number}.fq.gz",
    )
    for member in tar.getmembers():
        if not member.isfile():
            continue
        name = member.name
        lower = name.lower()
        if not (lower.endswith(".fq.gz") or lower.endswith(".fastq.gz")):
            continue
        if any(name.endswith(suffix) for suffix in suffixes):
            return member
    return None


def extract_member(tar: tarfile.TarFile, member: tarfile.TarInfo, output: Path, overwrite: bool) -> None:
    if output.exists() and not overwrite:
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    source = tar.extractfile(member)
    if source is None:
        raise RuntimeError(f"Could not extract {member.name}")
    tmp = output.with_suffix(output.suffix + ".tmp")
    with source, tmp.open("wb") as handle:
        shutil.copyfileobj(source, handle)
    tmp.replace(output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml", help="Project config YAML.")
    parser.add_argument("--archive-dir", default=None, help="BGI Separate directory containing per-sample tar.gz files.")
    parser.add_argument("--output-dir", default=None, help="Directory for staged FASTQ files.")
    parser.add_argument("--samplesheet", default=None, help="Output samplesheet TSV.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite staged FASTQ files.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config_path = Path(args.config)
    project_root = config_path.resolve().parent.parent
    config = load_config(config_path)
    paths = config.get("paths", {})
    defaults = config.get("metadata_defaults", {})

    archive_dir = Path(args.archive_dir) if args.archive_dir else resolve_project_path(project_root, paths.get("raw_bgi_archive_dir"))
    output_dir = Path(args.output_dir) if args.output_dir else resolve_project_path(project_root, paths.get("raw_fastq_dir"))
    samplesheet = Path(args.samplesheet) if args.samplesheet else resolve_project_path(project_root, paths.get("samplesheet"))
    if archive_dir is None or output_dir is None or samplesheet is None:
        raise SystemExit("archive-dir, output-dir, and samplesheet must be configured.")
    if not archive_dir.exists():
        raise SystemExit(f"Archive directory does not exist: {archive_dir}")

    rows: list[dict[str, str]] = []
    errors: list[str] = []
    for archive in sorted(archive_dir.glob("*.tar.gz")):
        sample_id = archive.name.removesuffix(".tar.gz")
        if not SAMPLE_RE.match(sample_id):
            errors.append(f"Invalid sample_id from archive name: {sample_id}")
            continue
        r1_out = output_dir / f"{sample_id}_R1.fastq.gz"
        r2_out = output_dir / f"{sample_id}_R2.fastq.gz"
        try:
            with tarfile.open(archive, "r:gz") as tar:
                r1_member = find_member(tar, sample_id, 1)
                r2_member = find_member(tar, sample_id, 2)
                if r1_member is None or r2_member is None:
                    errors.append(f"{sample_id}: expected paired clean FASTQ files not found in {archive}")
                    continue
                extract_member(tar, r1_member, r1_out, args.overwrite)
                extract_member(tar, r2_member, r2_out, args.overwrite)
        except (tarfile.TarError, OSError, RuntimeError) as exc:
            errors.append(f"{sample_id}: {exc}")
            continue

        rows.append(
            {
                "sample_id": sample_id,
                "fastq_r1": str(r1_out.relative_to(project_root)),
                "fastq_r2": str(r2_out.relative_to(project_root)),
                "initial_species_label": str(defaults.get("initial_species_label", "unknown")),
                "initial_growth_type": str(defaults.get("initial_growth_type", "presumed_SGM")),
                "collection_date": str(defaults.get("collection_date", "NA")),
                "isolation_source": str(defaults.get("isolation_source", "NA")),
                "patient_id_anonymized": str(defaults.get("patient_id_anonymized", "NA")),
                "sequencing_platform": str(defaults.get("sequencing_platform", "DNBSEQ")),
                "read_length": str(defaults.get("read_length", "150")),
                "notes": str(defaults.get("notes", "imported_from_BGI_clean_fastq_archives")),
            }
        )

    fieldnames = [
        "sample_id",
        "fastq_r1",
        "fastq_r2",
        "initial_species_label",
        "initial_growth_type",
        "collection_date",
        "isolation_source",
        "patient_id_anonymized",
        "sequencing_platform",
        "read_length",
        "notes",
    ]
    samplesheet.parent.mkdir(parents=True, exist_ok=True)
    with samplesheet.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Staged {len(rows)} samples into {output_dir}")
    print(f"Wrote samplesheet: {samplesheet}")
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
