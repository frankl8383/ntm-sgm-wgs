from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_atlas.py"


def record(accession: str, paired: str, contamination: float) -> dict:
    return {
        "accession": accession,
        "current_accession": accession,
        "paired_accession": paired,
        "assembly_info": {
            "assembly_level": "Contig",
            "biosample": {"accession": "SAMN_TEST", "attributes": []},
        },
        "assembly_stats": {
            "total_sequence_length": "5200000",
            "gc_percent": 68.5,
            "number_of_contigs": 100,
            "contig_n50": 50000,
        },
        "average_nucleotide_identity": {
            "submitted_species": "Mycobacterium avium",
            "taxonomy_check_status": "OK",
        },
        "checkm_info": {"completeness": 96.0, "contamination": contamination},
        "organism": {"organism_name": "Mycobacterium avium", "tax_id": 1764},
    }


class BuildAtlasTests(unittest.TestCase):
    def test_qc_precedes_paired_accession_selection(self):
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            input_path = root / "records.jsonl"
            records = [
                record("GCF_000000001.1", "GCA_000000001.1", 10.0),
                record("GCA_000000001.1", "GCF_000000001.1", 1.0),
            ]
            input_path.write_text(
                "".join(json.dumps(item) + "\n" for item in records),
                encoding="utf-8",
            )
            output = root / "output"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--jsonl",
                    str(input_path),
                    "--output-dir",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with (output / "ncbi_mac_atlas_qc_selected.tsv").open() as handle:
                selected = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["assembly_accession"] for row in selected], ["GCA_000000001.1"])
        self.assertEqual(selected[0]["selection_decision"], "selected_qc_first")


if __name__ == "__main__":
    unittest.main()
