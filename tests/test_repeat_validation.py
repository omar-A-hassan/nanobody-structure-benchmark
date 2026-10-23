"""Boundary checks for versioned repeat archives and checkpoint logging."""
from pathlib import Path
import hashlib
import json
import tempfile
import unittest
import warnings
import zipfile

from reproduce import verify_coordinate_archive
from src.verify_archive_reruns import checkpoint_provenance
from src.summarize_completed_repeats import compare, FIELDS


class RepeatValidationTests(unittest.TestCase):
    def test_reusable_nbb_notebook_keeps_corrected_weight_logging(self):
        root = Path(__file__).resolve().parents[1]
        notebook = json.loads((root / 'kaggle/kernels/revision-archive-20261023/nbb2/repeat.ipynb').read_text())
        code = '\n'.join(''.join(c['source']) for c in notebook['cells'] if c['cell_type'] == 'code')
        self.assertIn("os.path.basename(p).startswith('nanobody_model_')", code)
        self.assertIn("assert len([p for p in R['checkpoint_sha256']", code)
        self.assertIn("]) == 4", code)

    def test_paired_changes_are_not_differences_between_medians(self):
        primary = {t: dict.fromkeys(FIELDS, x) for t, x in [('a', 1), ('b', 3), ('c', 5)]}
        repeated = {t: dict.fromkeys(FIELDS, x) for t, x in [('a', 2), ('b', 3), ('c', 100)]}
        row = compare(primary, repeated, ['a', 'b', 'c'])['chain_cdr3']
        self.assertEqual(row['primary_median'], 3)
        self.assertEqual(row['repeat_median'], 3)
        self.assertEqual(row['median_paired_change'], 1)
        self.assertEqual(row['median_absolute_change'], 1)
        self.assertEqual(row['maximum_absolute_change'], 95)
        self.assertEqual(row['exact_values'], 1)
        self.assertEqual(row['n'], 3)

    def test_v1_missing_checkpoint_hashes_is_explicitly_accepted(self):
        result = checkpoint_provenance({}, "nbb2", required=False)
        self.assertEqual(result["status"], "not_recorded")
        with self.assertRaises(AssertionError):
            checkpoint_provenance({}, "nbb2", required=True)

    def test_new_numbered_checkpoint_hashes_are_required_when_requested(self):
        hashes = {f"weights/nanobody_model_{i}": "a" * 64 for i in range(1, 5)}
        self.assertEqual(checkpoint_provenance(hashes, "nbb2", required=True)["nbb_weight_hashes_logged"], 4)
        with self.assertRaises(AssertionError):
            checkpoint_provenance({"nanobody_model_1": "invalid"}, "nbb2", required=False)

    def test_coordinate_checks_reject_tampering_duplicates_and_traversal(self):
        payload = b"ATOM coordinate fixture\n"
        entry = {"path": "target.pdb", "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / "coordinates.zip"
            with zipfile.ZipFile(archive, "w") as stream:
                stream.writestr(entry["path"], payload)
            self.assertEqual(verify_coordinate_archive(archive, [entry]), 1)
            with self.assertRaises(AssertionError):
                verify_coordinate_archive(archive, [{**entry, "sha256": "0" * 64}])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                with zipfile.ZipFile(archive, "a") as stream:
                    stream.writestr(entry["path"], payload)
            with self.assertRaises(AssertionError):
                verify_coordinate_archive(archive, [entry])
            with zipfile.ZipFile(archive, "w") as stream:
                stream.writestr("../outside.pdb", payload)
            with self.assertRaises(AssertionError):
                verify_coordinate_archive(archive, [{**entry, "path": "../outside.pdb"}])


if __name__ == "__main__":
    unittest.main()
