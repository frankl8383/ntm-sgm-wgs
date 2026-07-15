#!/usr/bin/env python3
"""Recalculate recovery eligibility across residual-mixture threshold multipliers."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


NON_MIXTURE_GATES = [
    "assembly_gate",
    "route_concordance_gate",
    "checkm2_gate",
    "gunc_gate",
    "type_anchor_gate",
]


def analyse(data: pd.DataFrame, multipliers: list[float]) -> tuple[pd.DataFrame, pd.DataFrame]:
    baseline = float(data["benchmark_derived_mixture_threshold"].iloc[0]) / 5.0
    working = data.copy()
    working["maximum_route_burden"] = working[
        ["strict_mixed_sites_per_mbp", "meta_mixed_sites_per_mbp"]
    ].max(axis=1)
    working["all_non_mixture_gates"] = working[NON_MIXTURE_GATES].all(axis=1)
    expected = set(
        working.loc[
            working["final_rescue_decision"].str.startswith("rescued_interpretable"),
            "sample_id",
        ]
    )

    sample_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for multiplier in multipliers:
        threshold = baseline * multiplier
        mixture_ok = working["maximum_route_burden"].le(threshold)
        eligible = working["all_non_mixture_gates"] & mixture_ok
        for row, residual_ok, final_ok in zip(
            working.itertuples(index=False), mixture_ok, eligible, strict=True
        ):
            sample_rows.append(
                {
                    "multiplier": multiplier,
                    "threshold_sites_per_callable_mbp": threshold,
                    "sample_id": row.sample_id,
                    "strict_burden": row.strict_mixed_sites_per_mbp,
                    "meta_burden": row.meta_mixed_sites_per_mbp,
                    "maximum_route_burden": row.maximum_route_burden,
                    "all_non_mixture_gates": row.all_non_mixture_gates,
                    "residual_mixture_gate": bool(residual_ok),
                    "eligible_for_recovery": bool(final_ok),
                    "frozen_decision": row.final_rescue_decision,
                }
            )

        passing = set(working.loc[eligible, "sample_id"])
        summary_rows.append(
            {
                "multiplier": multiplier,
                "threshold_sites_per_callable_mbp": threshold,
                "eligible_count": int(eligible.sum()),
                "excluded_count": int((~eligible).sum()),
                "eligible_samples": ";".join(sorted(passing)),
                "excluded_samples": ";".join(sorted(working.loc[~eligible, "sample_id"])),
                "reproduces_frozen_classification": passing == expected,
            }
        )

    return pd.DataFrame(summary_rows), pd.DataFrame(sample_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--multipliers",
        type=float,
        nargs="+",
        default=[2.0, 3.0, 4.0, 5.0, 7.5, 10.0],
    )
    args = parser.parse_args()

    summary, sample_calls = analyse(pd.read_csv(args.input, sep="\t"), args.multipliers)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "threshold_multiplier_summary.tsv", sep="\t", index=False)
    sample_calls.to_csv(
        args.output_dir / "threshold_multiplier_sample_calls.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
