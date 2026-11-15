import unittest
from analysis_core import validate_dna, validate_rna, compare_sequences, hybrid_dg37, cfd_score_from_json

class TestSequences(unittest.TestCase):
    def test_rna_rejects_t(self):
        with self.assertRaisesRegex(ValueError, "T가 아니라 U"):
            validate_rna("A"*19+"T")

    def test_dna_rejects_u(self):
        with self.assertRaises(ValueError):
            validate_dna("A"*19+"U")

    def test_t_u_are_equivalent_for_alignment(self):
        result = compare_sequences("GAGTCCGAGCAGAAGAAGAA", "GAGUCCGAGCAGAAGAAGAA")
        self.assertEqual(result["mismatch_count"], 0)
        self.assertEqual(result["identity"], 1.0)

    def test_mismatch_position(self):
        result = compare_sequences("GAGTCCGAGCAGAAGAAGAA", "GAGUCCGAGCAGAAGCAGAA")
        self.assertEqual(result["mismatch_positions"], [16])

    def test_nn_energy_finite(self):
        self.assertLess(hybrid_dg37("GAGUCCGAGCAGAAGAAGAA"), 0)

    def test_cfd_uses_guide_rna_and_offtarget_dna(self):
        guide = "A" * 20
        off = "A" * 19 + "G"
        mm = {"rA:dC,20": 0.5}
        pam = {"GG": 1.0}
        self.assertAlmostEqual(cfd_score_from_json(guide, off, "AGG", mm, pam), 0.5)

if __name__ == "__main__":
    unittest.main()
