"""End-to-end CLI checks from outside the source checkout.

The fixture trains a small real model once from the unchanged bundled table.
Its numerical performance is not asserted: these tests verify execution,
provenance, split separation, prediction consistency and useful failure paths.
"""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class CliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="spotone cli tests with spaces ")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.work = Path(cls.temporary.name)
        cls.environment = os.environ.copy()
        cls.environment["PYTHONPATH"] = str(ROOT / "src")
        cls.environment["PYTHONDONTWRITEBYTECODE"] = "1"
        cls.fasta = cls.work / "proteins with spaces.fasta"
        cls.fasta.write_text(
            ">first description\nACDEFGHIKL\nMNPQRSTVWY\n\n>second\nWAGTMK\n",
            encoding="utf-8",
        )
        cls.model_dir = cls.work / "trained model"
        cls.model = cls.model_dir / "model.joblib"
        cls.feature_path = cls.work / "extracted features.csv"
        for arguments in (
            ("train", "--trees", "8", "--jobs", "1", "--seed", "42", "--output-dir", cls.model_dir),
            ("features", cls.fasta, "--output", cls.feature_path),
        ):
            completed = cls.command(*arguments)
            if completed.returncode:
                raise RuntimeError(
                    f"CLI fixture failed ({completed.returncode}): {arguments}\n"
                    f"stdout: {completed.stdout}\nstderr: {completed.stderr}"
                )

    @classmethod
    def python_command(cls, *arguments, cwd=None):
        return subprocess.run(
            [sys.executable, *map(str, arguments)],
            cwd=cls.work if cwd is None else cwd,
            env=cls.environment,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=60,
            check=False,
        )

    @classmethod
    def command(cls, *arguments, cwd=None):
        return cls.python_command("-m", "spotone", *arguments, cwd=cwd)

    def assert_success(self, completed):
        self.assertEqual(
            completed.returncode, 0,
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
        )

    def assert_user_error(self, completed, expected):
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("spotone: error:", completed.stderr)
        self.assertIn(expected, completed.stderr)
        self.assertNotIn("Traceback", completed.stderr)

    def test_bundled_training_writes_reviewable_artifacts_and_disjoint_groups(self):
        for name in ["model.joblib", "report.json", "split_assignments.csv", "test_predictions.csv"]:
            self.assertTrue((self.model_dir / name).is_file(), name)
        report = json.loads((self.model_dir / "report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["model_identity"], "reconstructed-spotone")
        self.assertEqual(report["dataset"]["rows"], 534)
        self.assertEqual(report["dataset"]["features"], 173)
        assignments = pd.read_csv(self.model_dir / "split_assignments.csv")
        self.assertEqual(len(assignments), 534)
        self.assertFalse(assignments["source_row"].duplicated().any())
        self.assertEqual(set(assignments["split"]), {"train", "validation", "test"})
        self.assertTrue(assignments.groupby("CPX")["split"].nunique().eq(1).all())
        held_out = pd.read_csv(self.model_dir / "test_predictions.csv")
        self.assertEqual(len(held_out), report["splits"]["test"]["rows"])
        self.assertEqual(set(held_out["CPX"]), set(report["splits"]["test"]["groups"]))

    def test_fasta_and_feature_csv_predictions_agree_with_sequence_coordinates(self):
        direct_path = self.work / "direct predictions.csv"
        table_path = self.work / "table predictions.csv"
        self.assert_success(self.command("predict", self.fasta, "--model", self.model, "-o", direct_path))
        self.assert_success(self.command("predict-features", self.feature_path, "--model", self.model, "-o", table_path))
        direct = pd.read_csv(direct_path)
        from_table = pd.read_csv(table_path)
        pd.testing.assert_frame_equal(direct, from_table)
        self.assertEqual(len(direct), 26)
        self.assertEqual(direct["sequence_id"].tolist(), ["first"] * 20 + ["second"] * 6)
        self.assertEqual(direct["position"].tolist(), list(range(1, 21)) + list(range(1, 7)))
        self.assertEqual("".join(direct["residue"]), "ACDEFGHIKLMNPQRSTVWYWAGTMK")
        self.assertTrue(direct["probability_hotspot"].between(0, 1).all())
        self.assertTrue(direct["threshold"].between(0, 1).all())
        self.assertEqual(direct["prediction"].tolist(), (direct["probability_hotspot"] >= direct["threshold"]).astype(int).tolist())
        self.assertEqual(direct["predicted_class"].tolist(), ["HS" if value else "NS" for value in direct["prediction"]])
        feature_table = pd.read_csv(self.feature_path)
        self.assertEqual(feature_table.shape, (26, 177))

    def test_malformed_fasta_fails_without_creating_predictions(self):
        invalid = self.work / "invalid.fasta"
        invalid.write_text(">bad\nAXD\n", encoding="utf-8")
        output = self.work / "invalid predictions.csv"
        completed = self.command("predict", invalid, "--model", self.model, "-o", output)
        self.assert_user_error(completed, "Invalid amino acid 'X'")
        self.assertIn("sequence position 2", completed.stderr)
        self.assertFalse(output.exists())

    def test_nonexistent_model_fails_with_an_actionable_path(self):
        missing = self.work / "missing model.joblib"
        output = self.work / "missing model predictions.csv"
        completed = self.command("predict", self.fasta, "--model", missing, "-o", output)
        self.assert_user_error(completed, "missing model.joblib")
        self.assertFalse(output.exists())

    def test_nonfinite_feature_csv_is_rejected(self):
        frame = pd.read_csv(self.feature_path)
        frame.loc[0, "MW"] = np.inf
        invalid = self.work / "nonfinite features.csv"
        frame.to_csv(invalid, index=False)
        output = self.work / "nonfinite predictions.csv"
        completed = self.command("predict-features", invalid, "--model", self.model, "-o", output)
        self.assert_user_error(completed, "Non-finite feature")
        self.assertIn("MW", completed.stderr)
        self.assertFalse(output.exists())

    def test_outputs_cannot_overwrite_input_sequences_features_or_model(self):
        originals = {path: path.read_bytes() for path in [self.fasta, self.feature_path, self.model]}
        for arguments in (
            ("features", self.fasta, "-o", self.fasta),
            ("predict", self.fasta, "--model", self.model, "-o", self.fasta),
            ("predict", self.fasta, "--model", self.model, "-o", self.model),
            ("predict", self.fasta, "--model", self.model_dir, "-o", self.model),
            ("predict-features", self.feature_path, "--model", self.model, "-o", self.feature_path),
            ("predict-features", self.feature_path, "--model", self.model_dir, "-o", self.model),
        ):
            with self.subTest(command=arguments[0], output=arguments[-1]):
                self.assert_user_error(self.command(*arguments), "output path must differ")
        for path, original in originals.items():
            self.assertEqual(path.read_bytes(), original)

    def test_help_runs_from_another_directory(self):
        clean = self.work / "help checks"
        clean.mkdir()
        self.assert_success(self.command("--help", cwd=clean))
        self.assertEqual(list(clean.iterdir()), [])

    def test_doctor_checks_dependencies_and_bundled_data_from_another_directory(self):
        completed = self.command("doctor")
        self.assert_success(completed)
        report = json.loads(completed.stdout)
        self.assertEqual(report["bundled_training_rows"], 534)
        self.assertEqual(report["features_per_residue"], 173)
        self.assertEqual(report["feature_check_rows"], 20)
        self.assertTrue({"numpy", "pandas", "scikit-learn", "joblib"}.issubset(report["dependencies"]))


if __name__ == "__main__":
    unittest.main()
