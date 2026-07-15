from __future__ import annotations

import pandas as pd

from scripts.threshold_sensitivity import analyse


def test_threshold_multiplier_classification() -> None:
    data = pd.DataFrame(
        {
            "sample_id": ["pass_a", "pass_b", "fail_mixed"],
            "benchmark_derived_mixture_threshold": [50.0, 50.0, 50.0],
            "strict_mixed_sites_per_mbp": [8.0, 25.0, 120.0],
            "meta_mixed_sites_per_mbp": [9.0, 24.0, 125.0],
            "assembly_gate": [True, True, True],
            "route_concordance_gate": [True, True, True],
            "checkm2_gate": [True, True, True],
            "gunc_gate": [True, True, True],
            "type_anchor_gate": [True, True, True],
            "final_rescue_decision": [
                "rescued_interpretable",
                "rescued_interpretable",
                "excluded_residual_mixture",
            ],
        }
    )

    summary, sample_calls = analyse(data, [2.0, 3.0, 10.0])

    assert summary["reproduces_frozen_classification"].tolist() == [False, True, True]
    assert len(sample_calls) == 9
