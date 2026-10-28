"""Scientific descriptor and FASTA regression checks (stdlib unittest)."""

from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from spotone.fasta import AMINO_ACIDS, FastaRecord, read_fasta
from spotone.features import FEATURE_COLUMNS, METADATA_COLUMNS, PROPERTY_COLUMNS, WINDOW_RADII, extract_features


ROOT = Path(__file__).resolve().parents[1]


class FastaTests(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.fasta"
            path.write_text(text, encoding="utf-8")
            return read_fasta(path)

    def test_wrapped_lowercase_blank_lines_descriptions_and_bom(self):
        self.assertEqual(
            self.parse("\ufeff\n>pdb1a4y_A descriptive header\nac d\n\neF\n>second\ngh\n"),
            [FastaRecord("pdb1a4y_A", "ACDEF"), FastaRecord("second", "GH")],
        )

    def test_invalid_fasta_is_rejected(self):
        for text, message in [
            ("", "no protein"),
            ("ACD\n", "before the first"),
            ("> \nACD\n", "Empty FASTA header"),
            (">first\n>second\nACD\n", "empty sequence"),
            (">first\nACD\n>last\n", "empty sequence"),
            (">same description\nACD\n>same other\nEFG\n", "Duplicate FASTA identifier"),
        ]:
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, message):
                self.parse(text)

    def test_unknown_ambiguous_gap_stop_and_digits_are_not_removed(self):
        for residue in "XBZJUO-*1ßı":
            with self.subTest(residue=residue), self.assertRaisesRegex(ValueError, "sequence position 3"):
                self.parse(f">protein\nAC{residue}D\n")


class FeatureTests(unittest.TestCase):
    def test_schema_matches_supplied_training_table(self):
        archived = pd.read_csv(ROOT / "datasets/1_homemade.csv", nrows=0)
        metadata = {"CPX", "PDBChain", "PDBResNo", "PDBResName", "Classe", "binary_class"}
        self.assertEqual(FEATURE_COLUMNS, [column for column in archived if column not in metadata])
        self.assertEqual(len(FEATURE_COLUMNS), 173)

    def test_property_resources_are_exact_snapshots(self):
        for name in ["amino_properties.csv", "encoding.csv"]:
            self.assertEqual((ROOT / "resources" / name).read_bytes(), (ROOT / "src/spotone/data" / name).read_bytes())

    def test_one_hot_identity_properties_and_quartile_boundaries(self):
        frame = extract_features([FastaRecord("all", AMINO_ACIDS)])
        self.assertEqual(frame.columns.tolist(), METADATA_COLUMNS + FEATURE_COLUMNS)
        np.testing.assert_array_equal(frame[FEATURE_COLUMNS[:20]].to_numpy(), np.eye(20))
        props = pd.read_csv(ROOT / "resources/amino_properties.csv").set_index("res_letter")
        np.testing.assert_array_equal(frame[list(PROPERTY_COLUMNS)], props.loc[list(AMINO_ACIDS), list(PROPERTY_COLUMNS)])
        self.assertEqual(frame["rel_dist_1"].tolist(), [1] * 4 + [2] * 5 + [3] * 5 + [4] * 6)
        self.assertEqual(frame["position"].tolist(), list(range(1, 21)))

    def test_windows_match_independent_direct_means_at_both_termini_and_center(self):
        sequence = AMINO_ACIDS * 9
        frame = extract_features([FastaRecord("protein", sequence)])
        properties = pd.read_csv(ROOT / "resources/amino_properties.csv").set_index("res_letter")
        lookup = properties.loc[list(sequence), list(PROPERTY_COLUMNS)].to_numpy(dtype=float)
        for center in [0, 1, 17, 90, len(sequence) - 2, len(sequence) - 1]:
            expected = []
            for property_index in range(len(PROPERTY_COLUMNS)):
                for radius in WINDOW_RADII:
                    expected.append(lookup[max(0, center - radius):min(len(sequence), center + radius + 1), property_index].mean())
            np.testing.assert_allclose(frame.loc[center, [str(x) for x in range(3, 136)]].to_numpy(dtype=float), expected, rtol=1e-12, atol=1e-10)

    def test_single_residue_and_multiple_sequences_preserve_coordinates(self):
        frame = extract_features([FastaRecord("first", "a"), FastaRecord("pdb1a4y_A", "WY")])
        self.assertEqual(frame["sequence_id"].tolist(), ["first", "pdb1a4y_A", "pdb1a4y_A"])
        self.assertEqual(frame["position"].tolist(), [1, 1, 2])
        self.assertEqual(frame["residue"].tolist(), ["A", "W", "Y"])
        self.assertEqual(frame["PDBResName"].tolist(), ["ALA", "TRP", "TYR"])
        self.assertEqual(frame.iloc[0]["rel_dist_1"], 4)
        for property_index, property_name in enumerate(PROPERTY_COLUMNS):
            for radius_index in range(len(WINDOW_RADII)):
                self.assertAlmostEqual(frame.iloc[0][str(3 + property_index * len(WINDOW_RADII) + radius_index)], frame.iloc[0][property_name])

    def test_compatible_historical_residue_vector(self):
        # This archived row is verified to use the same sequence coordinates.
        # Other archived rows have undocumented PDB-to-sequence mappings; they
        # must not be silently remapped or treated as regeneration fixtures.
        record = next(record for record in read_fasta(ROOT / "tests/data/pdb1a4y_A.fasta") if record.identifier == "pdb1a4y_A")
        frame = extract_features([record])
        archived = pd.read_csv(ROOT / "datasets/1_homemade.csv")
        row = archived.loc[(archived["CPX"] == "1A4Y") & (archived["PDBChain"] == "A") & (archived["PDBResNo"] == 261)].iloc[0]
        self.assertEqual(frame.iloc[260]["PDBResName"], row["PDBResName"])
        np.testing.assert_allclose(frame.iloc[260][FEATURE_COLUMNS].to_numpy(dtype=float), row[FEATURE_COLUMNS].to_numpy(dtype=float), rtol=1e-12, atol=1e-10)

    def test_direct_records_receive_the_same_validation(self):
        cases = [[], [FastaRecord("empty", "")], [FastaRecord("bad", "AX")], [FastaRecord("dup", "A"), FastaRecord("dup", "C")]]
        for records in cases:
            with self.subTest(records=records), self.assertRaises(ValueError):
                extract_features(records)


if __name__ == "__main__":
    unittest.main()
