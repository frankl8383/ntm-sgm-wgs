from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "residual_mixture.py"


def load_module():
    spec = importlib.util.spec_from_file_location("residual_mixture", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ResidualMixtureTests(unittest.TestCase):
    def test_mpileup_tokens_and_indels_are_skipped(self):
        module = load_module()
        self.assertEqual(
            module.base_counts("G", "^].,$.,+2agA-1tT*<>"),
            {"G": 4, "A": 1, "T": 1},
        )

    def test_known_minor_fractions_and_depth_filter(self):
        pileup = "\n".join(
            [
                "ctg\t1\tG\t20\t" + "." * 15 + "A" * 5 + "\t*",
                "ctg\t2\tG\t19\t" + "." * 14 + "A" * 5 + "\t*",
                "ctg\t3\tG\t50\t" + "." * 45 + "A" * 5 + "\t*",
                "ctg\t4\tG\t20\t" + "." * 20 + "\t*",
            ]
        ) + "\n"
        with tempfile.TemporaryDirectory() as tempdir:
            output = Path(tempdir) / "summary.tsv"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--sample",
                    "synthetic",
                    "--route",
                    "test",
                    "--output",
                    str(output),
                ],
                input=pileup,
                text=True,
                check=True,
            )
            with output.open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(int(row["callable_positions_depth_ge_20"]), 3)
        self.assertEqual(int(row["mixed_sites_maf_0.10_0.90"]), 2)
        self.assertEqual(int(row["mixed_sites_maf_0.20_0.80"]), 1)


if __name__ == "__main__":
    unittest.main()
