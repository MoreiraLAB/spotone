import unittest

import numpy as np

from spotone.fasta import FastaRecord
from spotone.features import FEATURE_COLUMNS, extract_features
from spotone.published import (
    PUBLISHED_MODEL_SHA256,
    adapt_probability,
    load_published_forest,
    predict_records,
    published_model_manifest,
)


class PublishedModelTests(unittest.TestCase):
    def test_manifest_and_structure(self):
        manifest = published_model_manifest()
        self.assertEqual(manifest["source_sha256"], PUBLISHED_MODEL_SHA256)
        self.assertEqual(manifest["source_sklearn_version"], "0.21.2")
        self.assertEqual(manifest["n_estimators"], 500)
        self.assertEqual(manifest["n_features"], 173)
        forest = load_published_forest()
        self.assertEqual(forest.n_estimators, 500)
        self.assertEqual(forest.n_features, len(FEATURE_COLUMNS))

    def test_deterministic_published_fixture(self):
        result = predict_records([FastaRecord("tiny", "ACDEFGHIKLMNPQRSTVWY")])
        expected_raw = [
            0.1517391304347826,
            0.148,
            0.228,
            0.142,
            0.1738472906403941,
        ]
        np.testing.assert_allclose(result.loc[:4, "raw_probability"], expected_raw, rtol=0, atol=1e-15)
        self.assertEqual(result.iloc[-1]["residue"], "Y")
        self.assertAlmostEqual(float(result.iloc[-1]["raw_probability"]), 0.5213441976442006, places=15)
        self.assertEqual(int(result["adapted_class"].sum()), 1)

    def test_historical_adaptation_boundaries_are_preserved(self):
        # Raw classification is >= 0.5, adapted classification is > 0.5.
        corrected, original, adapted = adapt_probability(0.5, "ALA")
        self.assertEqual((corrected, original, adapted), (0.5, 1, 0))
        # Historical PRO rule gates at +0.16 but adds +0.15.
        corrected, original, adapted = adapt_probability(0.34, "PRO")
        self.assertEqual((corrected, original, adapted), (0.49, 0, 0))
        # Historical ASN asymmetry: gate +0.20, add +0.15.
        corrected, original, adapted = adapt_probability(0.30, "ASN")
        self.assertEqual((corrected, original, adapted), (0.45, 0, 0))


    def test_runtime_checksum_is_enforced(self):
        import tempfile
        from pathlib import Path
        from spotone.published import load_published_forest
        load_published_forest.cache_clear()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "published_model_v1.npz"
            path.write_bytes(b"not the published model")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                load_published_forest(path)
        load_published_forest.cache_clear()

    def test_features_are_finite_and_173_wide(self):
        frame = extract_features([FastaRecord("tiny", "ACDEFGHIKLMNPQRSTVWY")])
        self.assertEqual(frame[FEATURE_COLUMNS].shape, (20, 173))
        self.assertTrue(np.isfinite(frame[FEATURE_COLUMNS].to_numpy(float)).all())


if __name__ == "__main__":
    unittest.main()
