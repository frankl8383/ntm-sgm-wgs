import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "near_mac_two_route.py"
SPEC = importlib.util.spec_from_file_location("near_mac_two_route", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class NearMacDesignTests(unittest.TestCase):
    def test_frozen_expanded_design(self):
        rows = MODULE.frozen_conditions()
        self.assertEqual(len(rows), 40)
        self.assertEqual(len({row["condition_id"] for row in rows}), 40)
        all_pair_10 = [
            row
            for row in rows
            if row["minor_percent"] == 10 and row["total_pairs"] == 2_000_000
        ]
        self.assertEqual(len(all_pair_10), 8)
        depth = [row for row in rows if row["pair_id"] in {"pairA", "pairD"}]
        self.assertEqual(len(depth), 36)
        self.assertEqual(
            {row["total_pairs"] for row in depth},
            {1_000_000, 2_000_000, 4_000_000},
        )
        self.assertTrue(
            all(row["seed"] == 202 for row in rows)
        )


if __name__ == "__main__":
    unittest.main()
