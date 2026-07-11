#!/usr/bin/env python3
"""Summarize the full-read rescue extension under the locked pilot gates."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


MIN_ANI = 99.0
MIN_RECIPROCAL_AF = 0.85
MIN_COMPLETENESS = 90.0
MAX_CONTAMINATION = 5.0
MIN_GENOME_SIZE = 4_000_000
MAX_GENOME_SIZE = 7_500_000
MIN_GC = 66.0
MAX_GC = 71.0
MAX_CONTIGS = 500
MIN_N50 = 20_000
MIN_TYPE_ANI = 94.0
MIN_TYPE_AF = 0.50


def read_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def read_fastani(path: Path) -> tuple[float, float]:
    fields = path.read_text().strip().split("\t")
    return float(fields[2]), int(fields[3]) / int(fields[4])


def alignment_rate(path: Path) -> float:
    match = re.search(r"([0-9.]+)% overall alignment rate", path.read_text())
    if not match:
        raise ValueError(f"No Bowtie2 alignment rate in {path}")
    return float(match.group(1))


def as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def parse_query(path: str) -> tuple[str, str]:
    stem = Path(path).stem
    sample, route = stem.rsplit("_", 1)
    if route not in {"strict", "meta"}:
        raise ValueError(path)
    return sample, route


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--batch-root", required=True)
    parser.add_argument("--pilot-root", required=True)
    parser.add_argument("--type-metadata", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-md", required=True)
    args = parser.parse_args()

    root = Path(args.batch_root)
    pilot_root = Path(args.pilot_root)
    with Path(args.config).open(newline="") as handle:
        config = list(csv.DictReader(handle, delimiter="\t"))

    pure_rates = [
        float(
            read_row(
                pilot_root
                / "mixture_metrics"
                / "cross65"
                / "cross65.metaspades_bin.minor_allele_burden.tsv"
            )["mixed_sites_20_80_per_mbp_callable"]
        ),
        float(
            read_row(
                pilot_root
                / "mixture_metrics"
                / "cross65_strict"
                / "cross65.strict.minor_allele_burden.tsv"
            )["mixed_sites_20_80_per_mbp_callable"]
        ),
    ]
    mixture_threshold = 5 * max(pure_rates)

    checkm = {
        row["Name"]: row
        for row in csv.DictReader((root / "checkm2" / "quality_report.tsv").open(), delimiter="\t")
    }
    gunc = {
        row["genome"]: row
        for row in csv.DictReader(
            (root / "gunc" / "GUNC.progenomes_2.1.maxCSS_level.tsv").open(),
            delimiter="\t",
        )
    }
    type_metadata = {
        row["accession"]: row
        for row in csv.DictReader(Path(args.type_metadata).open(), delimiter="\t")
    }
    type_hits: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    with (root / "type_ani" / "full_rescue_vs_type_anchors.fastani.tsv").open() as handle:
        for line in handle:
            query, reference, ani, matched, total = line.rstrip().split("\t")
            sample, route = parse_query(query)
            accession = Path(reference).stem
            type_hits[(sample, route)].append(
                {
                    "accession": accession,
                    "ani": float(ani),
                    "af": int(matched) / int(total),
                    "label": type_metadata[accession]["organism_name"],
                    "requested_name": type_metadata[accession]["requested_taxon_name"],
                }
            )

    rows: list[dict[str, object]] = []
    for cfg in config:
        sample = cfg["sample_id"]
        strict_stats = read_row(root / "strict" / sample / f"{sample}.assembly_stats.tsv")
        meta_stats = read_row(root / "meta" / sample / f"{sample}.metaspades_mac_bin_stats.tsv")
        strict_mix = read_row(root / "mixture_metrics" / sample / f"{sample}.strict.minor_allele_burden.tsv")
        meta_mix = read_row(root / "mixture_metrics" / sample / f"{sample}.meta.minor_allele_burden.tsv")
        meta_to_strict_ani, meta_to_strict_af = read_fastani(
            root / "meta" / sample / f"{sample}.meta_to_strict.tsv"
        )
        strict_to_meta_ani, strict_to_meta_af = read_fastani(
            root / "meta" / sample / f"{sample}.strict_to_meta.tsv"
        )
        strict_qc = checkm[f"{sample}_strict"]
        meta_qc = checkm[f"{sample}_meta"]
        strict_gunc = gunc[f"{sample}_strict"]
        meta_gunc = gunc[f"{sample}_meta"]
        strict_type = max(type_hits[(sample, "strict")], key=lambda hit: float(hit["ani"]))
        meta_type = max(type_hits[(sample, "meta")], key=lambda hit: float(hit["ani"]))

        sizes = [int(strict_stats["genome_size_bp"]), int(meta_stats["genome_size_bp"])]
        gcs = [float(strict_stats["gc_percent"]), float(meta_stats["gc_percent"])]
        contigs = [int(strict_stats["contigs"]), int(meta_stats["contigs"])]
        n50s = [int(strict_stats["n50_bp"]), int(meta_stats["n50_bp"])]
        mixed = [
            float(strict_mix["mixed_sites_20_80_per_mbp_callable"]),
            float(meta_mix["mixed_sites_20_80_per_mbp_callable"]),
        ]
        gates = {
            "assembly": min(sizes) >= MIN_GENOME_SIZE
            and max(sizes) <= MAX_GENOME_SIZE
            and min(gcs) >= MIN_GC
            and max(gcs) <= MAX_GC
            and max(contigs) <= MAX_CONTIGS
            and min(n50s) >= MIN_N50,
            "route_concordance": min(meta_to_strict_ani, strict_to_meta_ani) >= MIN_ANI
            and min(meta_to_strict_af, strict_to_meta_af) >= MIN_RECIPROCAL_AF,
            "checkm2": min(float(strict_qc["Completeness"]), float(meta_qc["Completeness"]))
            >= MIN_COMPLETENESS
            and max(float(strict_qc["Contamination"]), float(meta_qc["Contamination"]))
            <= MAX_CONTAMINATION,
            "gunc": as_bool(strict_gunc["pass.GUNC"]) and as_bool(meta_gunc["pass.GUNC"]),
            "minor_alleles": max(mixed) <= mixture_threshold,
            "type_anchor": strict_type["requested_name"] == meta_type["requested_name"]
            and min(float(strict_type["ani"]), float(meta_type["ani"])) >= MIN_TYPE_ANI
            and min(float(strict_type["af"]), float(meta_type["af"])) >= MIN_TYPE_AF,
        }
        if all(gates.values()):
            if strict_type["accession"] == meta_type["accession"]:
                decision = "rescued_interpretable"
            else:
                decision = "rescued_interpretable_type_boundary_caution"
        elif not gates["minor_alleles"]:
            decision = "unrecoverable_residual_mac_mixture"
        elif not gates["route_concordance"]:
            decision = "unrecoverable_route_discordance"
        elif not gates["assembly"] or not gates["checkm2"] or not gates["gunc"]:
            decision = "unrecoverable_genome_qc"
        else:
            decision = "unrecoverable_type_anchor_instability"

        rows.append(
            {
                "sample_id": sample,
                "batch_role": cfg["batch_role"],
                "baseline_mycobacterium_fraction": cfg["baseline_mycobacterium_fraction"],
                "baseline_checkm2_contamination": cfg["baseline_checkm2_contamination"],
                "strict_recruitment_rate_percent": alignment_rate(
                    root / "strict" / sample / "logs" / f"{sample}.bowtie2.log"
                ),
                "strict_genome_size_bp": sizes[0],
                "meta_genome_size_bp": sizes[1],
                "strict_contigs": contigs[0],
                "meta_contigs": contigs[1],
                "strict_n50_bp": n50s[0],
                "meta_n50_bp": n50s[1],
                "strict_gc_percent": gcs[0],
                "meta_gc_percent": gcs[1],
                "strict_checkm2_completeness": strict_qc["Completeness"],
                "meta_checkm2_completeness": meta_qc["Completeness"],
                "strict_checkm2_contamination": strict_qc["Contamination"],
                "meta_checkm2_contamination": meta_qc["Contamination"],
                "strict_gunc_pass": strict_gunc["pass.GUNC"],
                "meta_gunc_pass": meta_gunc["pass.GUNC"],
                "strict_gunc_contamination_portion": strict_gunc["contamination_portion"],
                "meta_gunc_contamination_portion": meta_gunc["contamination_portion"],
                "strict_mixed_sites_per_mbp": mixed[0],
                "meta_mixed_sites_per_mbp": mixed[1],
                "benchmark_derived_mixture_threshold": f"{mixture_threshold:.6f}",
                "meta_to_strict_ani": meta_to_strict_ani,
                "strict_to_meta_ani": strict_to_meta_ani,
                "meta_to_strict_af": meta_to_strict_af,
                "strict_to_meta_af": strict_to_meta_af,
                "best_type_anchor_accession": strict_type["accession"],
                "best_type_anchor_label": strict_type["label"],
                "meta_best_type_anchor_accession": meta_type["accession"],
                "meta_best_type_anchor_label": meta_type["label"],
                "type_anchor_reporting_name": strict_type["requested_name"],
                "exact_type_anchor_concordance": strict_type["accession"] == meta_type["accession"],
                "strict_type_ani": strict_type["ani"],
                "meta_type_ani": meta_type["ani"],
                "strict_type_af": strict_type["af"],
                "meta_type_af": meta_type["af"],
                "assembly_gate": gates["assembly"],
                "route_concordance_gate": gates["route_concordance"],
                "checkm2_gate": gates["checkm2"],
                "gunc_gate": gates["gunc"],
                "minor_allele_gate": gates["minor_alleles"],
                "type_anchor_gate": gates["type_anchor"],
                "final_rescue_decision": decision,
            }
        )

    output = Path(args.output_tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(str(row["final_rescue_decision"]) for row in rows)
    rescued_count = sum(
        value for key, value in counts.items() if key.startswith("rescued_interpretable")
    )
    report = [
        "# Full rescue extension",
        "",
        f"Ten samples were processed with the two locked rescue routes. Residual-mixture review used the same deterministic two-million-pair depth as the synthetic benchmark; the fixed threshold was {mixture_threshold:.1f} intermediate-frequency sites per Mb callable sequence.",
        "",
        f"- Rescued interpretable: {rescued_count}",
        f"- Unrecoverable residual MAC mixture: {counts['unrecoverable_residual_mac_mixture']}",
        f"- Other unrecoverable categories: {len(rows) - rescued_count - counts['unrecoverable_residual_mac_mixture']}",
        "",
        "## Sample decisions",
        "",
    ]
    for row in rows:
        report.append(
            f"- {row['sample_id']}: **{row['final_rescue_decision']}**; recruitment "
            f"{float(row['strict_recruitment_rate_percent']):.2f}%; strict/meta contigs "
            f"{row['strict_contigs']}/{row['meta_contigs']}; mixed-site burden "
            f"{float(row['strict_mixed_sites_per_mbp']):.1f}/{float(row['meta_mixed_sites_per_mbp']):.1f} per Mb; "
            f"type anchor {row['type_anchor_reporting_name']}."
        )
    rescued = [
        str(row["sample_id"])
        for row in rows
        if str(row["final_rescue_decision"]).startswith("rescued_interpretable")
    ]
    failed = [
        str(row["sample_id"])
        for row in rows
        if not str(row["final_rescue_decision"]).startswith("rescued_interpretable")
    ]
    report.extend(
        [
            "",
            "## Cohort consequence",
            "",
            f"The full-read extension supports rescue of {', '.join(rescued)}. {', '.join(failed)} remain excluded. Together with the 13 directly retained genomes, the provisional interpretable local cohort contains {13 + len(rescued)} genomes before final public-context placement and nomenclature review.",
        ]
    )
    Path(args.output_md).write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
