"""Reconstructed SPOTONE models with explicit, group-separated evaluation.

The historical repository does not contain fitted weights. Models produced here
are new fits to the supplied feature tables, not the published 2020 predictor.
Only load joblib artifacts you created yourself or obtained from a trusted source:
joblib, like pickle, can execute code while loading.
"""

from __future__ import annotations

import csv
import hashlib
import json
from numbers import Integral
import pickle
import platform
from importlib.resources import files
from pathlib import Path
from typing import Any
import warnings

import joblib
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.ensemble import AdaBoostClassifier, ExtraTreesClassifier
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold, GroupShuffleSplit, RandomizedSearchCV
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from . import __version__

ARTIFACT_VERSION = 1
IDENTITY_COLUMNS = ["CPX", "PDBChain", "PDBResNo", "PDBResName"]
TARGET_COLUMNS = ["Classe", "binary_class"]
METADATA_COLUMNS = IDENTITY_COLUMNS + TARGET_COLUMNS + [
    "sequence_id", "position", "residue"
]
AMINO_ACIDS = (
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
)
MODEL_NAMES = ("extra-trees", "adaboost", "svm", "mlp")


def _numeric_features(frame: pd.DataFrame, names: list[str]) -> pd.DataFrame:
    """Validate numeric values without replacing absent/invalid data by zeros."""
    if not names:
        raise ValueError("The input table contains no feature columns.")
    try:
        values = frame.loc[:, names].apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError("Feature columns must contain only numeric values.") from exc
    finite = np.isfinite(values.to_numpy(dtype=float))
    if not finite.all():
        bad = np.argwhere(~finite)[0]
        raise ValueError(
            f"Non-finite feature at row {int(bad[0]) + 1}, column {names[int(bad[1])]!r}."
        )
    return values.astype(float)


def read_dataset(path: str | Path) -> tuple[pd.DataFrame, list[str]]:
    """Read a labeled archived table, rejecting ambiguous identifiers or labels."""
    path = Path(path)
    # pandas renames duplicate CSV headers; inspect the original header first.
    with path.open(encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle), [])
    if len(header) != len(set(header)):
        raise ValueError("The dataset has duplicate column names.")
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError("The dataset is empty.")
    required = IDENTITY_COLUMNS + TARGET_COLUMNS
    missing = [name for name in required if name not in frame]
    if missing:
        raise ValueError(f"Dataset is missing required columns: {', '.join(missing)}.")
    for name in IDENTITY_COLUMNS:
        if frame[name].isna().any() or frame[name].astype(str).str.strip().eq("").any():
            raise ValueError(f"Dataset identifier {name!r} contains an empty value.")
    frame["CPX"] = frame["CPX"].astype(str).str.strip().str.upper()
    frame["PDBChain"] = frame["PDBChain"].astype(str).str.strip()
    frame["PDBResName"] = frame["PDBResName"].astype(str).str.strip().str.upper()
    if not frame["PDBResName"].isin(AMINO_ACIDS).all():
        raise ValueError("PDBResName must use the 20 canonical three-letter amino-acid codes.")
    if frame.duplicated(["CPX", "PDBChain", "PDBResNo"]).any():
        raise ValueError("Duplicate residue identities (CPX, PDBChain, PDBResNo) in dataset.")
    frame["Classe"] = frame["Classe"].astype(str).str.strip().str.upper()
    if not frame["Classe"].isin(["NS", "HS"]).all():
        raise ValueError("Classe must contain only NS and HS.")
    try:
        numeric_labels = pd.to_numeric(frame["binary_class"], errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError("binary_class must contain 0 (NS) or 1 (HS).") from exc
    if not numeric_labels.isin([0, 1]).all():
        raise ValueError("binary_class must contain 0 (NS) or 1 (HS).")
    frame["binary_class"] = numeric_labels.astype(int)
    if not frame["Classe"].map({"NS": 0, "HS": 1}).eq(frame["binary_class"]).all():
        raise ValueError("Classe and binary_class disagree; expected NS=0 and HS=1.")
    if frame["binary_class"].nunique() != 2:
        raise ValueError("Training requires both HS and NS examples.")
    feature_names = [name for name in frame.columns if name not in METADATA_COLUMNS]
    # Reserved prediction columns are a common accidental target leak when a
    # previous prediction table is supplied as training data.
    leaked = set(feature_names) & {
        "prediction", "predicted_class", "probability_hotspot", "threshold", "split"
    }
    if leaked:
        raise ValueError(f"Prediction/output columns cannot be features: {sorted(leaked)}.")
    frame.loc[:, feature_names] = _numeric_features(frame, feature_names)
    return frame, feature_names


def split_dataset(frame: pd.DataFrame, seed: int = 42) -> dict[str, np.ndarray]:
    """Make deterministic ~60/20/20 CPX-group splits with both classes in each.

    Only label presence is checked when allocating groups. No model is fitted or
    test score inspected to choose a split. Groups are complexes, not sequence
    homology clusters; related proteins can therefore still cross partitions.
    """
    if frame["CPX"].nunique() < 5:
        raise ValueError("At least five distinct CPX groups are required for train/validation/test splits.")
    groups = frame["CPX"].to_numpy()
    labels = frame["binary_class"].to_numpy()
    positions = np.arange(len(frame))
    for attempt in range(100):
        split_seed = seed + attempt
        outer = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=split_seed)
        development, test = next(outer.split(positions, labels, groups))
        inner = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=split_seed)
        train_local, validation_local = next(
            inner.split(development, labels[development], groups[development])
        )
        split = {
            "train": np.sort(development[train_local]),
            "validation": np.sort(development[validation_local]),
            "test": np.sort(test),
        }
        if all(np.unique(labels[index]).size == 2 for index in split.values()):
            return split
    raise ValueError(
        "Could not form group-separated train/validation/test partitions containing both classes. "
        "Provide more CPX groups with HS and NS examples."
    )


def _binary_labels(values: Any) -> np.ndarray:
    """Validate binary labels before casting; never truncate fractional labels."""
    try:
        raw = np.asarray(values)
        if np.iscomplexobj(raw):
            raise ValueError("Binary labels cannot be complex numbers.")
        numeric = np.asarray(values, dtype=float).ravel()
    except (TypeError, ValueError) as exc:
        raise ValueError("Evaluation labels must be 0 or 1.") from exc
    if not np.isfinite(numeric).all() or not np.isin(numeric, [0, 1]).all():
        raise ValueError("Evaluation labels must be 0 or 1.")
    return numeric.astype(int)


def evaluate_predictions(
    actual: Any, probability: Any, predictions: Any | None = None
) -> dict[str, Any]:
    """Report HS-positive and weighted metrics; ROC-AUC always uses scores."""
    y = _binary_labels(actual)
    p = np.asarray(probability, dtype=float).ravel()
    predicted = (p >= 0.5).astype(int) if predictions is None else _binary_labels(predictions)
    if not (len(y) == len(p) == len(predicted)) or not len(y):
        raise ValueError("Labels, probabilities and predictions must have equal, nonzero lengths.")
    if not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("Probabilities must be finite values between zero and one.")
    result: dict[str, Any] = {
        "n": int(len(y)),
        "accuracy": float(accuracy_score(y, predicted)),
        "balanced_accuracy": float(np.mean([
            np.mean(predicted[y == label] == label) for label in np.unique(y)
        ])),
        "roc_auc": float(roc_auc_score(y, p)) if np.unique(y).size == 2 else None,
        "average_precision": float(average_precision_score(y, p)) if np.any(y == 1) else None,
        "precision_hs": float(precision_score(y, predicted, zero_division=0)),
        "recall_hs": float(recall_score(y, predicted, zero_division=0)),
        "f1_hs": float(f1_score(y, predicted, zero_division=0)),
        "precision_weighted": float(precision_score(y, predicted, average="weighted", zero_division=0)),
        "recall_weighted": float(recall_score(y, predicted, average="weighted", zero_division=0)),
        "f1_weighted": float(f1_score(y, predicted, average="weighted", zero_division=0)),
        "mcc": float(matthews_corrcoef(y, predicted)) if np.unique(np.concatenate([y, predicted])).size > 1 else 0.0,
        "confusion_matrix_order": ["NS", "HS"],
        "confusion_matrix": confusion_matrix(y, predicted, labels=[0, 1]).tolist(),
    }
    return result


def tune_residue_thresholds(
    residues: Any, actual: Any, probability: Any, raw_predictions: Any | None = None
) -> dict[str, float]:
    """Apply the historical FN/TN rule using validation rows only.

    If a residue type has false negatives and true negatives, and the largest
    false-negative HS score exceeds the largest true-negative score, its decision
    threshold is reduced from 0.5 to that false-negative score. Types lacking
    either group retain 0.5. These are decision thresholds, not calibration.
    """
    residues = np.asarray(residues, dtype=str)
    actual = _binary_labels(actual)
    probability = np.asarray(probability, dtype=float).ravel()
    predicted = (probability >= 0.5).astype(int) if raw_predictions is None else _binary_labels(raw_predictions)
    if not (len(residues) == len(actual) == len(probability) == len(predicted)):
        raise ValueError("Residues, labels, probabilities and predictions must have equal lengths.")
    if not np.isfinite(probability).all() or ((probability < 0) | (probability > 1)).any():
        raise ValueError("Probabilities must be finite values between zero and one.")
    thresholds = {amino: 0.5 for amino in AMINO_ACIDS}
    for amino in AMINO_ACIDS:
        negative_prediction = (residues == amino) & (predicted == 0)
        fn = probability[negative_prediction & (actual == 1)]
        tn = probability[negative_prediction & (actual == 0)]
        if len(fn) and len(tn) and float(tn.max()) < float(fn.max()) < 0.5:
            thresholds[amino] = float(fn.max())
    return thresholds


def _estimator(name: str, seed: int, n_estimators: int, n_jobs: int) -> Any:
    if name == "extra-trees":
        return ExtraTreesClassifier(
            n_estimators=n_estimators, random_state=seed, n_jobs=n_jobs,
            bootstrap=True, class_weight="balanced_subsample", max_features="sqrt",
        )
    if name == "adaboost":
        return AdaBoostClassifier(n_estimators=n_estimators, random_state=seed)
    if name == "svm":
        return make_pipeline(StandardScaler(), SVC(probability=True, random_state=seed))
    if name == "mlp":
        return make_pipeline(StandardScaler(), MLPClassifier(random_state=seed, max_iter=1000))
    raise ValueError(f"Unknown model {name!r}; choose one of {', '.join(MODEL_NAMES)}.")


def _positive_probability(estimator: Any, features: pd.DataFrame) -> np.ndarray:
    classes = np.asarray(estimator.classes_)
    indices = np.flatnonzero(classes == 1)
    if len(indices) != 1:
        raise ValueError("Model does not contain the HS class (1).")
    return estimator.predict_proba(features)[:, int(indices[0])]


def _schema_name(names: list[str]) -> str:
    from .features import FEATURE_COLUMNS, FEATURE_SCHEMA

    return FEATURE_SCHEMA if len(names) == len(FEATURE_COLUMNS) and set(names) == set(FEATURE_COLUMNS) else "custom-feature-table-v1"


def _resource_hashes() -> dict[str, str]:
    return {
        name: hashlib.sha256(files("spotone").joinpath("data", name).read_bytes()).hexdigest()
        for name in ("amino_properties.csv", "encoding.csv")
    }


def train_model(
    dataset_path: str | Path,
    output_dir: str | Path,
    *,
    seed: int = 42,
    n_estimators: int = 500,
    n_jobs: int = 1,
    tune: bool = False,
    model_name: str = "extra-trees",
) -> dict[str, Any]:
    """Fit train rows, tune thresholds on validation, and evaluate unseen CPX groups."""
    if n_estimators < 1:
        raise ValueError("n_estimators must be a positive integer.")
    if n_jobs == 0:
        raise ValueError("n_jobs cannot be zero.")
    if seed < 0 or seed > 2**32 - 101:
        raise ValueError("seed must be between 0 and 4294967195.")
    if tune and model_name != "extra-trees":
        raise ValueError("Bounded hyperparameter tuning is available for extra-trees only.")
    dataset_path = Path(dataset_path)
    output_dir = Path(output_dir)
    for filename in ("model.joblib", "report.json", "split_assignments.csv", "test_predictions.csv"):
        target = output_dir / filename
        if target.resolve() == dataset_path.resolve() or (
            target.exists() and dataset_path.exists() and target.samefile(dataset_path)
        ):
            raise ValueError(
                f"Training output {target} would overwrite the input dataset. "
                "Choose a different output directory."
            )
    frame, feature_names = read_dataset(dataset_path)
    split = split_dataset(frame, seed)
    features = _numeric_features(frame, feature_names)
    labels = frame["binary_class"].to_numpy()
    estimator = _estimator(model_name, seed, n_estimators, n_jobs)
    train = split["train"]
    tuning: dict[str, Any] = {"enabled": bool(tune)}
    if tune:
        folds = min(5, frame.iloc[train]["CPX"].nunique())
        if folds < 2:
            raise ValueError("Hyperparameter tuning requires at least two training CPX groups.")
        distributions = {
            "max_features": [1.0, "sqrt", "log2"],
            "max_depth": [None, 5, 10, 20],
            "min_samples_leaf": [1, 2, 4],
            "class_weight": [None, "balanced", "balanced_subsample"],
        }
        # Candidate jobs are parallel; each individual forest uses one worker to
        # avoid multiplying process/thread counts.
        estimator.set_params(n_jobs=1)
        search = RandomizedSearchCV(
            estimator, distributions, n_iter=12, scoring="f1",
            cv=GroupKFold(n_splits=folds), random_state=seed,
            n_jobs=n_jobs, refit=True, error_score="raise",
        )
        search.fit(features.iloc[train], labels[train], groups=frame.iloc[train]["CPX"])
        estimator = search.best_estimator_
        estimator.set_params(n_jobs=n_jobs)
        tuning.update({
            "candidates": 12, "folds": int(folds), "scoring": "f1_hs",
            "group_column": "CPX", "partition": "train",
            "best_parameters": search.best_params_,
            "best_cv_score": float(search.best_score_),
        })
    else:
        estimator.fit(features.iloc[train], labels[train])
    validation = split["validation"]
    validation_probability = _positive_probability(estimator, features.iloc[validation])
    thresholds = tune_residue_thresholds(
        frame.iloc[validation]["PDBResName"], labels[validation], validation_probability,
    )
    versions = {
        "spotone": __version__, "python": platform.python_version(), "numpy": np.__version__,
        "pandas": pd.__version__, "scikit_learn": sklearn.__version__,
        "joblib": joblib.__version__, "scipy": scipy.__version__,
    }
    source = {
        "filename": dataset_path.name,
        "sha256": hashlib.sha256(dataset_path.read_bytes()).hexdigest(),
        "rows": int(len(frame)), "features": len(feature_names),
        "complexes": int(frame["CPX"].nunique()),
    }
    resource_hashes = _resource_hashes()
    bundle: dict[str, Any] = {
        "artifact_version": ARTIFACT_VERSION,
        "model_identity": "reconstructed-spotone",
        "model_name": model_name,
        "estimator": estimator,
        "feature_columns": feature_names,
        "feature_schema": _schema_name(feature_names),
        "thresholds": thresholds,
        "positive_class": 1,
        "class_labels": {0: "NS", 1: "HS"},
        "seed": seed, "versions": versions, "dataset": source,
        "feature_resource_sha256": resource_hashes,
        "training_scope": "train split only; no refit after validation/test evaluation",
        "threshold_scope": "validation split only",
    }
    report: dict[str, Any] = {
        "model_identity": bundle["model_identity"],
        "model_name": model_name,
        "feature_schema": bundle["feature_schema"],
        "seed": seed,
        "versions": versions,
        "dataset": source,
        "feature_resource_sha256": resource_hashes,
        "residue_coverage": {
            amino: {
                "all_rows": int((frame["PDBResName"] == amino).sum()),
                "all_hs": int(((frame["PDBResName"] == amino) & (frame["binary_class"] == 1)).sum()),
                "training_rows": int((frame.iloc[train]["PDBResName"] == amino).sum()),
                "training_hs": int(((frame.iloc[train]["PDBResName"] == amino) & (frame.iloc[train]["binary_class"] == 1)).sum()),
            }
            for amino in AMINO_ACIDS
        },
        "model_parameters": {
            key: value for key, value in estimator.get_params(deep=True).items()
            if isinstance(value, (str, int, float, bool, type(None)))
        },
        "split_policy": {
            "group_column": "CPX", "fractions_by_group": {"train": 0.6, "validation": 0.2, "test": 0.2},
            "algorithm": "GroupShuffleSplit, first allocation with both classes in each partition",
            "max_allocation_attempts": 100,
            "homology_separated": False,
        },
        "training_scope": bundle["training_scope"],
        "threshold_scope": bundle["threshold_scope"],
        "threshold_method": "validation-only maximum-FN versus maximum-TN rule using p(HS) >= 0.5 as the initial decision",
        "thresholds": thresholds,
        "hyperparameter_tuning": tuning,
        "splits": {},
        "limitations": [
            "Newly reconstructed fit; historical fitted weights and exact split assignments are not supplied.",
            "These results do not reproduce or validate the metrics reported in the 2020 paper.",
            "Complex separation does not guarantee separation of homologous or identical protein sequences.",
            "Archived training features are preserved; supplied FASTA sequences do not reproduce every archived row.",
            "FASTA prediction positions are 1-based sequence offsets, not PDB author residue numbers.",
            "HS probabilities are classifier scores, not calibrated confidence or predicted binding-energy changes.",
        ],
        "artifacts": {
            "model": "model.joblib", "report": "report.json",
            "split_assignments": "split_assignments.csv", "test_predictions": "test_predictions.csv",
        },
    }
    assignments = frame[IDENTITY_COLUMNS + TARGET_COLUMNS].copy()
    assignments.insert(0, "source_row", np.arange(len(frame)))
    assignments["split"] = ""
    test_predictions: pd.DataFrame | None = None
    for split_name, positions in split.items():
        partition = frame.iloc[positions]
        predictions = predict_features(partition, bundle)
        groups = sorted(partition["CPX"].unique().tolist())
        report["splits"][split_name] = {
            "rows": len(positions), "groups": groups, "n_groups": len(groups),
            "class_counts": {"NS": int((labels[positions] == 0).sum()), "HS": int((labels[positions] == 1).sum())},
            "default_threshold": evaluate_predictions(labels[positions], predictions["probability_hotspot"]),
            "residue_thresholds": evaluate_predictions(
                labels[positions], predictions["probability_hotspot"], predictions["prediction"]
            ),
        }
        assignments.loc[positions, "split"] = split_name
        if split_name == "test":
            test_predictions = predictions
            test_predictions.insert(0, "source_row", positions)
    output_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, output_dir / "model.joblib", compress=3)
    assignments.to_csv(output_dir / "split_assignments.csv", index=False)
    assert test_predictions is not None
    test_predictions.to_csv(output_dir / "test_predictions.csv", index=False)
    (output_dir / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return report


def load_model(path: str | Path) -> dict[str, Any]:
    """Load a trusted local artifact, including thresholds and feature schema."""
    path = Path(path)
    if path.is_dir():
        path = path / "model.joblib"
    try:
        bundle = joblib.load(path)
    except (EOFError, pickle.UnpicklingError, KeyError, IndexError, AttributeError, ImportError, TypeError, ValueError) as exc:
        raise ValueError(
            f"Could not load model artifact {path}: {type(exc).__name__}: {exc}. "
            "Use a complete trusted SPOTONE artifact and its recorded environment, or retrain."
        ) from exc
    if not isinstance(bundle, dict) or bundle.get("artifact_version") != ARTIFACT_VERSION:
        raise ValueError("Unsupported SPOTONE model artifact; retrain using this version.")
    required = {"estimator", "feature_columns", "feature_schema", "thresholds", "versions"}
    if not required.issubset(bundle):
        raise ValueError(f"Model artifact is incomplete: missing {sorted(required - set(bundle))}.")
    names = bundle["feature_columns"]
    if not isinstance(names, list) or not names or not all(isinstance(n, str) for n in names) or len(set(names)) != len(names):
        raise ValueError("Model artifact has an invalid feature schema.")
    if not isinstance(bundle["feature_schema"], str) or not bundle["feature_schema"]:
        raise ValueError("Model artifact has an invalid feature schema identifier.")
    versions = bundle["versions"]
    if not isinstance(versions, dict) or not versions.get("scikit_learn") or not all(
        isinstance(key, str) and isinstance(value, str) and bool(value)
        for key, value in versions.items()
    ):
        raise ValueError("Model artifact has an invalid versions mapping; scikit_learn must be recorded.")
    estimator = bundle["estimator"]
    if not callable(getattr(estimator, "predict_proba", None)):
        raise ValueError("Model artifact estimator does not support predict_proba.")
    try:
        classes = np.asarray(estimator.classes_)
        feature_count = estimator.n_features_in_
    except (AttributeError, TypeError) as exc:
        raise ValueError("Model artifact estimator is not fitted or lacks classes_/n_features_in_.") from exc
    if not (np.array_equal(classes, [0, 1]) or np.array_equal(classes, [1, 0])):
        raise ValueError("Model artifact estimator must contain exactly the NS=0 and HS=1 classes.")
    if not isinstance(feature_count, Integral) or isinstance(feature_count, bool) or feature_count != len(names):
        raise ValueError("Model artifact estimator feature count does not match its saved feature columns.")
    if hasattr(estimator, "feature_names_in_") and list(estimator.feature_names_in_) != names:
        raise ValueError("Model artifact estimator feature order does not match its saved feature columns.")
    if not isinstance(bundle["thresholds"], dict) or set(bundle["thresholds"]) != set(AMINO_ACIDS):
        raise ValueError("Model artifact does not contain all residue thresholds.")
    if not all(isinstance(v, (int, float)) and np.isfinite(v) and 0 <= v <= 1 for v in bundle["thresholds"].values()):
        raise ValueError("Model artifact contains invalid decision thresholds.")
    from .features import FEATURE_SCHEMA

    if bundle["feature_schema"] == FEATURE_SCHEMA and bundle.get("feature_resource_sha256") != _resource_hashes():
        raise ValueError(
            "Model feature-resource hashes do not match this installation. "
            "Use the original property snapshots or retrain with the intended resources."
        )
    stored_sklearn = versions["scikit_learn"]
    if stored_sklearn != sklearn.__version__:
        warnings.warn(
            f"Model was fitted with scikit-learn {stored_sklearn}; runtime is {sklearn.__version__}. "
            "Use the recorded environment or retrain for reliable serialization compatibility.",
            RuntimeWarning,
            stacklevel=2,
        )
    return bundle


def predict_features(frame: pd.DataFrame, bundle: dict[str, Any]) -> pd.DataFrame:
    """Score in saved feature order; never substitute missing features or residues."""
    if frame.empty:
        raise ValueError("Prediction input contains no residues.")
    if frame.columns.duplicated().any():
        raise ValueError("Prediction input has duplicate column names.")
    names = bundle["feature_columns"]
    missing = [name for name in names if name not in frame]
    if missing:
        raise ValueError(f"Prediction input is missing {len(missing)} model feature(s): {', '.join(missing[:10])}.")
    unexpected = [name for name in frame if name not in names and name not in METADATA_COLUMNS]
    if unexpected:
        raise ValueError(f"Unexpected columns outside the saved feature schema: {', '.join(unexpected[:10])}.")
    if "PDBResName" not in frame:
        raise ValueError("Prediction input requires PDBResName for residue-specific decision thresholds.")
    residues = frame["PDBResName"].astype(str).str.strip().str.upper()
    if not residues.isin(AMINO_ACIDS).all():
        raise ValueError("Prediction residues must use canonical three-letter amino-acid codes.")
    values = _numeric_features(frame, names)
    probability = _positive_probability(bundle["estimator"], values)
    threshold = residues.map(bundle["thresholds"]).to_numpy(dtype=float)
    if not np.isfinite(probability).all() or not np.isfinite(threshold).all():
        raise ValueError("Model produced non-finite probabilities or thresholds.")
    if ((probability < 0) | (probability > 1) | (threshold < 0) | (threshold > 1)).any():
        raise ValueError("Model probabilities and thresholds must lie between zero and one.")
    prediction = (probability >= threshold).astype(int)
    metadata = [name for name in frame.columns if name in METADATA_COLUMNS]
    result = frame.loc[:, metadata].copy().reset_index(drop=True)
    result["probability_hotspot"] = probability
    result["threshold"] = threshold
    result["predicted_class"] = np.where(prediction == 1, "HS", "NS")
    result["prediction"] = prediction
    return result
