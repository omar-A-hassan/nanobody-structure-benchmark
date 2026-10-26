"""Checks for sequence correspondence and native-contact overlap calculations."""
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("complex_audit", ROOT / "src/audit_complex_overlap.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class ComplexOverlapTests(unittest.TestCase):
    def test_global_identity_counts_gaps(self):
        engine = audit.aligner()
        self.assertEqual(audit.identity("AAAA", "AAA", engine), .75)
        self.assertEqual(audit.identity("ACDE", "ACDE", engine), 1)
        self.assertEqual(audit.residue_map("ACDE", "ACE", engine), {0: 0, 1: 1, 2: 3})

    def test_unmapped_residues_are_not_discarded(self):
        result = audit.overlap([0, 1], [0, 1, 2], {0: 0, 1: 1})
        self.assertEqual(result["jaccard"], 2/3)
        result = audit.contact_similarity([(0, 0)], [(0, 0), (1, 1)], {0: 0}, {0: 0})
        self.assertEqual(result["jaccard"], .5)
        self.assertEqual(result["query_recall"], 1)
        self.assertIsNone(audit.overlap([], [], {})["jaccard"])
        self.assertEqual(audit.overlap([0], [], {})["jaccard"], 0)
        self.assertEqual(audit.contact_similarity([(0, 0)], [], {}, {})["jaccard"], 0)

    def test_cache_does_not_merge_chain_case(self):
        self.assertNotEqual(audit.coordinate_path("pdb_test_N_+").name.casefold(),
                            audit.coordinate_path("pdb_test_n_+").name.casefold())

    def test_heavy_atom_contact_extraction(self):
        names = ["auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_atom_id", "type_symbol",
                 "label_comp_id", "occupancy", "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num"]
        payload = "data_test\nloop_\n" + "\n".join("_atom_site."+n for n in names) + "\n"
        payload += "H 1 ? CA C GLY 1 0 0 0 1\nA 5 ? CA C GLY 1 3 0 0 1\nA 6 ? CA C ALA 1 10 0 0 1\n"
        payload += "A 6 ? H H ALA 1 1 0 0 1\nA 7 ? CA C ALA 1 1 0 0 2\n"
        entry = {"key": "test", "hchain": "H", "chain": "A", "heavy_numbering": [["G", [1, " "]]],
                 "antigen_numbering": [["G", [5, " "]], ["A", [6, " "]], ["A", [7, " "]]]}
        result = audit.extract_contacts(payload.encode(), entry)
        self.assertEqual(result["contacts"], [(0, 0)])
        self.assertEqual(result["epitope"], [0])
        entry["heavy_numbering"] = [["A", [1, " "]]]
        with self.assertRaisesRegex(AssertionError, "Heavy residue mismatch"):
            audit.extract_contacts(payload.encode(), entry)
        entry["heavy_numbering"] = [["G", [1, " "]]]
        entry["hchain"] = "h"
        with self.assertRaisesRegex(ValueError, "case-sensitive"):
            audit.extract_contacts(payload.encode(), entry)
        entry["hchain"] = "H"
        with self.assertRaisesRegex(ValueError, "positive-occupancy"):
            audit.extract_contacts(payload.replace("GLY 1 0 0 0", "GLY 0 0 0 0").encode(), entry)


if __name__ == "__main__":
    unittest.main()
