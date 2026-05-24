#!/usr/bin/env python3
"""Write allowed/prohibited wording guardrails for Figure 6 features."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    rows = [
        {
            "feature_or_cell": "rrl/rrs hotspot review",
            "allowed_wording": "No known curated hotspot alternate was detected under the current genome-only rRNA hotspot review.",
            "prohibited_wording": "The isolates are clinically macrolide- or aminoglycoside-susceptible.",
            "required_caveat": "Clinical resistance/susceptibility requires AST/MIC and species-specific interpretation.",
        },
        {
            "feature_or_cell": "erm screen",
            "allowed_wording": "No confident Erm-like hit was detected by the current protein-level screen.",
            "prohibited_wording": "Macrolide inducible resistance is absent.",
            "required_caveat": "Erm biology and inducible resistance are species dependent and require phenotype/literature support.",
        },
        {
            "feature_or_cell": "Mi18 partial aph(3')-IIa",
            "allowed_wording": "Mi18 carried a low-confidence partial AMRFinderPlus aph(3')-IIa hit retained for manual review.",
            "prohibited_wording": "Mi18 is aminoglycoside resistant.",
            "required_caveat": "Partial hit on a short/low-coverage contig is not a confident AMR determinant.",
        },
        {
            "feature_or_cell": "gyrA/gyrB coding differences",
            "allowed_wording": "gyrA/gyrB differences are manual-review flags relative to selected reference loci.",
            "prohibited_wording": "These are fluoroquinolone resistance mutations.",
            "required_caveat": "Do not infer fluoroquinolone resistance without AST/MIC and validated NTM species-specific mutations.",
        },
        {
            "feature_or_cell": "STRESS:arsN1",
            "allowed_wording": "AMRFinderPlus identified an arsN1 stress-feature signal in 11/13 retained genomes.",
            "prohibited_wording": "arsN1 proves increased virulence, environmental adaptation, or clinical treatment relevance.",
            "required_caveat": "Report as a recurrent stress-feature signal only; biological relevance requires targeted analysis.",
        },
        {
            "feature_or_cell": "Figure 6 overall",
            "allowed_wording": "AMR/stress feature screening and curated NTM resistance-locus review.",
            "prohibited_wording": "Clinical antimicrobial resistance prediction.",
            "required_caveat": "Genome-only review; not clinical resistance prediction until AST/MIC is integrated.",
        },
    ]
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
