"""Check regional masking, rigid fits and the supplied prepared VHH score rows."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from src.score_prepared_vhh import calculate, framework_scores, original_kabsch, read_vhh


class PreparedVHHTests(unittest.TestCase):
    def test_framework_fit_does_not_fit_the_loop(self):
        reference = np.random.default_rng(17).normal(size=(20, 4, 3))
        masks = {"cdr3": np.arange(20) >= 16, "framework": np.arange(20) < 16}
        prediction = reference.copy()
        prediction[masks["cdr3"]] += [5, 0, 0]
        prediction[0, 3] = np.nan
        scores, matched, fit_atoms, _ = framework_scores(prediction, reference, masks, original_kabsch())
        self.assertAlmostEqual(scores["cdr3"], 5)
        self.assertAlmostEqual(scores["framework"], 0)
        self.assertEqual((matched, fit_atoms), (79, 63))
        rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
        moved = prediction @ rotation + [50, -30, 11]
        rotated, _, _, _ = framework_scores(moved, reference, masks, original_kabsch())
        for key in scores:
            self.assertAlmostEqual(scores[key], rotated[key])

    def test_prepared_rows_match_the_reported_precision(self):
        report = calculate()
        expected = {
            "NanoBodyBuilder2": (2.660, .724, 1.059, .928),
            "IgFold": (3.079, 1.424, 2.613, .889),
            "SimpleFold": (2.196, .813, 5.763, .960),
            "Experimental": (0., 0., 0., 0.),
        }
        for row in report["rows"]:
            actual = tuple(round(row[k], 3) for k in ("cdr1", "cdr2", "cdr3", "framework"))
            self.assertEqual(actual, expected[row["model"]])
            self.assertEqual((row["residues"], row["heavy_atoms"], row["matched_atoms"], row["framework_fit_atoms"]),
                             (126, 973, 492, 350))

    def test_parser_rejects_wrong_sequence_and_duplicate_atoms(self):
        root = Path(__file__).resolve().parents[1]
        source = root / "data/docking/7R73_experimental_DOCK.pdb"
        with self.assertRaises(AssertionError):
            read_vhh(source, "WRONG")
        with tempfile.TemporaryDirectory() as directory:
            duplicate = Path(directory) / "duplicate.pdb"
            atom = next(line for line in source.read_text().splitlines() if line.startswith("ATOM  "))
            duplicate.write_text(atom + "\n" + atom + "\n")
            with self.assertRaisesRegex(AssertionError, "Duplicate atom"):
                read_vhh(duplicate, "")


if __name__ == "__main__":
    unittest.main()
