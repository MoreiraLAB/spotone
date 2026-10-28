"""Behavioral checks for leakage, saved decision rules and schema validation."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import pickle
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from spotone.model import (
    evaluate_predictions,
    load_model,
    predict_features,
    read_dataset,
    split_dataset,
    train_model,
    tune_residue_thresholds,
)


def small_dataset() -> pd.DataFrame:
    rows = []
    for group in range(15):
        for position in range(8):
            positive = position % 3 == 0
            rows.append({
                "CPX": f"G{group:03d}", "PDBChain": "A", "PDBResNo": position + 1,
                "PDBResName": "TRP" if position % 2 else "ALA",
                "Classe": "HS" if positive else "NS", "binary_class": int(positive),
                "signal": (position % 3) + group / 40,
                "context": ((position + group) % 7) / 7,
            })
    return pd.DataFrame(rows)


class EvaluationTests(unittest.TestCase):
    def test_auc_uses_ranking_not_thresholded_predictions(self):
        metrics = evaluate_predictions([0, 0, 1, 1], [0.05, 0.1, 0.2, 0.3])
        self.assertEqual(metrics["roc_auc"], 1.0)
        self.assertEqual(metrics["recall_hs"], 0.0)
        self.assertEqual(metrics["confusion_matrix"], [[2, 0], [2, 0]])

    def test_single_class_auc_is_explicitly_unavailable(self):
        metrics = evaluate_predictions([0, 0], [0.1, 0.2])
        self.assertIsNone(metrics["roc_auc"])
        self.assertIsNone(metrics["average_precision"])

    def test_fractional_labels_are_not_silently_truncated(self):
        for actual, predicted in [([0.2, 1], None), ([0, 1], [0, 0.8])]:
            with self.subTest(actual=actual, predicted=predicted):
                with self.assertRaisesRegex(ValueError, "labels must be 0 or 1"):
                    evaluate_predictions(actual, [0.1, 0.9], predicted)
        with self.assertRaisesRegex(ValueError, "labels must be 0 or 1"):
            tune_residue_thresholds(["ALA", "ALA"], [0.2, 1], [0.1, 0.4])

    def test_thresholds_use_historical_false_negative_rule(self):
        thresholds = tune_residue_thresholds(
            ["ALA", "ALA", "TRP", "TRP", "CYS"],
            [0, 1, 0, 1, 1], [0.2, 0.4, 0.4, 0.2, 0.1],
        )
        self.assertAlmostEqual(thresholds["ALA"], 0.4)
        self.assertAlmostEqual(thresholds["TRP"], 0.5)
        self.assertAlmostEqual(thresholds["CYS"], 0.5)
        self.assertEqual(len(thresholds), 20)

    def test_probability_prediction_disagreement_cannot_increase_threshold(self):
        thresholds = tune_residue_thresholds(
            ["ALA", "ALA"], [0, 1], [0.2, 0.8], raw_predictions=[0, 0]
        )
        self.assertEqual(thresholds["ALA"], 0.5)

    def test_group_partitions_are_disjoint_reproducible_and_exhaustive(self):
        frame = small_dataset()
        first = split_dataset(frame, seed=42)
        second = split_dataset(frame, seed=42)
        groups = []
        for name in ("train", "validation", "test"):
            np.testing.assert_array_equal(first[name], second[name])
            self.assertEqual(set(frame.iloc[first[name]]["binary_class"]), {0, 1})
            groups.append(set(frame.iloc[first[name]]["CPX"]))
        self.assertFalse(groups[0] & groups[1])
        self.assertFalse(groups[0] & groups[2])
        self.assertFalse(groups[1] & groups[2])
        self.assertEqual(sorted(np.concatenate(list(first.values()))), list(range(len(frame))))


class ModelRoundtripTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.frame = small_dataset()
        self.dataset = self.root / "dataset.csv"
        self.frame.to_csv(self.dataset, index=False)

    def tearDown(self):
        self.temporary.cleanup()

    def train(self):
        return train_model(self.dataset, self.root / "model", n_estimators=16)

    def test_artifact_roundtrip_retains_thresholds_schema_and_scope(self):
        report = self.train()
        bundle = load_model(self.root / "model")
        self.assertEqual(bundle["thresholds"], report["thresholds"])
        self.assertEqual(bundle["feature_columns"], ["signal", "context"])
        self.assertEqual(bundle["feature_schema"], "custom-feature-table-v1")
        self.assertTrue({"spotone", "scipy"}.issubset(bundle["versions"]))
        self.assertIn("train split only", bundle["training_scope"])
        predicted = predict_features(self.frame, bundle)
        reordered = predict_features(self.frame.loc[:, list(reversed(self.frame.columns))], bundle)
        np.testing.assert_array_equal(predicted["probability_hotspot"], reordered["probability_hotspot"])
        np.testing.assert_array_equal(predicted["prediction"], reordered["prediction"])
        assignments = pd.read_csv(self.root / "model" / "split_assignments.csv")
        self.assertEqual(assignments.groupby("CPX")["split"].nunique().max(), 1)
        with (self.root / "model" / "report.json").open() as handle:
            persisted = json.load(handle)
        self.assertEqual(report, persisted)
        test = pd.read_csv(self.root / "model" / "test_predictions.csv")
        self.assertEqual(len(test), report["splits"]["test"]["rows"])
        self.assertEqual(set(test["CPX"]), set(report["splits"]["test"]["groups"]))
        self.assertEqual(bundle["estimator"].n_features_in_, 2)

    def test_changed_saved_threshold_changes_prediction_after_reload(self):
        self.train()
        bundle = load_model(self.root / "model")
        bundle["thresholds"] = {name: 0.0 for name in bundle["thresholds"]}
        import joblib

        joblib.dump(bundle, self.root / "changed.joblib")
        changed = predict_features(self.frame, load_model(self.root / "changed.joblib"))
        self.assertTrue(changed["prediction"].eq(1).all())

    def test_feature_resource_hash_mismatch_fails_loading(self):
        self.train()
        bundle = load_model(self.root / "model")
        from spotone.features import FEATURE_SCHEMA
        import joblib

        bundle["feature_schema"] = FEATURE_SCHEMA
        bundle["feature_resource_sha256"]["amino_properties.csv"] = "incorrect"
        joblib.dump(bundle, self.root / "changed.joblib")
        with self.assertRaisesRegex(ValueError, "hashes do not match"):
            load_model(self.root / "changed.joblib")

    def test_missing_nonfinite_or_unexpected_features_fail(self):
        self.train()
        bundle = load_model(self.root / "model")
        with self.assertRaisesRegex(ValueError, "missing"):
            predict_features(self.frame.drop(columns=["signal"]), bundle)
        frame = self.frame.copy()
        frame.loc[0, "signal"] = float("nan")
        with self.assertRaisesRegex(ValueError, "Non-finite"):
            predict_features(frame, bundle)
        with self.assertRaisesRegex(ValueError, "Unexpected"):
            predict_features(self.frame.assign(extra=0.0), bundle)

    def test_training_outputs_cannot_overwrite_the_dataset(self):
        for index, filename in enumerate(("model.joblib", "report.json", "split_assignments.csv", "test_predictions.csv")):
            with self.subTest(filename=filename):
                output = self.root / f"conflict-{index}"
                output.mkdir()
                dataset = output / filename
                original = self.dataset.read_bytes()
                dataset.write_bytes(original)
                with patch("spotone.model._estimator", side_effect=AssertionError("should fail before fitting")):
                    with self.assertRaisesRegex(ValueError, "would overwrite the input dataset"):
                        train_model(dataset, output, n_estimators=1)
                self.assertEqual(dataset.read_bytes(), original)
                self.assertEqual(list(output.iterdir()), [dataset])

    def test_training_output_symlink_cannot_overwrite_the_dataset(self):
        output = self.root / "linked-output"
        output.mkdir()
        (output / "model.joblib").symlink_to(self.dataset)
        original = self.dataset.read_bytes()
        with self.assertRaisesRegex(ValueError, "would overwrite the input dataset"):
            train_model(self.dataset, output, n_estimators=1)
        self.assertEqual(self.dataset.read_bytes(), original)

    def test_corrupt_or_incompatible_model_bytes_have_actionable_errors(self):
        corrupted = self.root / "corrupted.joblib"
        for contents in (b"", b"not a pickle\n"):
            with self.subTest(contents=contents):
                corrupted.write_bytes(contents)
                with self.assertRaisesRegex(ValueError, "Could not load model artifact.*corrupted.joblib"):
                    load_model(corrupted)
        for exception in (
            pickle.UnpicklingError("bad stream"), KeyError("bad opcode"), IndexError("truncated"),
            AttributeError("missing class"), ModuleNotFoundError("old dependency"),
        ):
            with self.subTest(exception=type(exception).__name__):
                with patch("spotone.model.joblib.load", side_effect=exception):
                    with self.assertRaisesRegex(ValueError, "Could not load model artifact"):
                        load_model(corrupted)

    def test_malformed_model_metadata_and_estimators_are_rejected_at_load(self):
        self.train()
        valid = load_model(self.root / "model")
        import joblib

        cases = [
            ("versions", [], "invalid versions mapping"),
            ("versions", {"scikit_learn": 1}, "invalid versions mapping"),
            ("versions", {}, "invalid versions mapping"),
            ("estimator", None, "does not support predict_proba"),
        ]
        for key, value, message in cases:
            with self.subTest(key=key, value=value):
                bundle = {**valid, key: value}
                joblib.dump(bundle, self.root / "malformed.joblib")
                with self.assertRaisesRegex(ValueError, message):
                    load_model(self.root / "malformed.joblib")
        for attribute, value, message in [
            ("classes_", np.array([0, 2]), "exactly the NS=0 and HS=1 classes"),
            ("n_features_in_", 99, "feature count does not match"),
            ("feature_names_in_", np.array(["context", "signal"]), "feature order does not match"),
        ]:
            with self.subTest(attribute=attribute):
                bundle = deepcopy(valid)
                setattr(bundle["estimator"], attribute, value)
                joblib.dump(bundle, self.root / "malformed.joblib")
                with self.assertRaisesRegex(ValueError, message):
                    load_model(self.root / "malformed.joblib")

    def test_duplicate_identity_is_rejected(self):
        pd.concat([self.frame, self.frame.iloc[[0]]]).to_csv(self.dataset, index=False)
        with self.assertRaisesRegex(ValueError, "Duplicate residue"):
            read_dataset(self.dataset)

    def test_disagreeing_labels_are_rejected(self):
        self.frame.loc[0, "Classe"] = "NS"
        self.frame.to_csv(self.dataset, index=False)
        with self.assertRaisesRegex(ValueError, "disagree"):
            read_dataset(self.dataset)

    def test_duplicate_headers_are_rejected_before_pandas_renames_them(self):
        self.dataset.write_text("CPX,CPX\na,b\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate column"):
            read_dataset(self.dataset)

    def test_prediction_columns_are_not_allowed_as_training_features(self):
        self.frame.assign(prediction=1).to_csv(self.dataset, index=False)
        with self.assertRaisesRegex(ValueError, "cannot be features"):
            read_dataset(self.dataset)

    def test_tuning_does_not_use_test_complexes(self):
        report = train_model(
            self.dataset, self.root / "tuned", n_estimators=4, tune=True,
        )
        self.assertEqual(report["hyperparameter_tuning"]["partition"], "train")
        self.assertEqual(report["hyperparameter_tuning"]["candidates"], 12)
        self.assertEqual(report["threshold_scope"], "validation split only")


if __name__ == "__main__":
    unittest.main()
