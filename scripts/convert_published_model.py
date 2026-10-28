#!/usr/bin/env python3
"""Convert the trusted historical SPOTONE sklearn-0.21 pickle into neutral NumPy arrays.

The converter accepts only the exact published model SHA256 by default and uses a
restricted Unpickler that permits only the old sklearn container classes and the
NumPy reconstruction helpers present in this artifact. The generated NPZ contains
only numeric tree arrays needed for inference; no executable pickle is needed at
runtime.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
from pathlib import Path

import numpy as np

EXPECTED_SHA256 = "c6746cff126bb70218cd2d3ac91764ca629a9d405d415ed269e504e92712f3af"
MODEL_VERSION = "published-2020-extratrees-v1"


class _EstimatorProxy:
    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        obj._newargs = args
        return obj

    def __setstate__(self, state):
        if not isinstance(state, dict):
            raise pickle.UnpicklingError("Unexpected estimator pickle state")
        self.__dict__.update(state)


class _TreeProxy(_EstimatorProxy):
    pass


class _RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        allowed = {
            ("sklearn.ensemble.forest", "ExtraTreesClassifier"): _EstimatorProxy,
            ("sklearn.tree.tree", "ExtraTreeClassifier"): _EstimatorProxy,
            ("sklearn.tree._tree", "Tree"): _TreeProxy,
        }
        if (module, name) in allowed:
            return allowed[(module, name)]
        if (module, name) in {
            ("numpy", "dtype"),
            ("numpy", "ndarray"),
            ("numpy.core.multiarray", "_reconstruct"),
            ("numpy.core.multiarray", "scalar"),
        }:
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f"Blocked global {module}.{name}")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_legacy_pickle(path: Path, expected_sha256: str = EXPECTED_SHA256):
    digest = _sha256(path)
    if expected_sha256 and digest != expected_sha256:
        raise ValueError(f"Unexpected model SHA256: {digest}")
    with path.open("rb") as handle:
        model = _RestrictedUnpickler(handle).load()
    if getattr(model, "_sklearn_version", None) != "0.21.2":
        raise ValueError("Expected a scikit-learn 0.21.2 model")
    if getattr(model, "n_features_", None) != 173:
        raise ValueError("Expected 173 model features")
    if len(getattr(model, "estimators_", [])) != 500:
        raise ValueError("Expected 500 ExtraTrees estimators")
    if not np.array_equal(np.asarray(model.classes_), np.array([0, 1])):
        raise ValueError("Expected classes [0, 1]")
    return model, digest


def convert(source: Path, output: Path, manifest: Path) -> dict:
    model, digest = load_legacy_pickle(source)
    offsets = [0]
    left, right, feature, threshold = [], [], [], []
    value0, value1 = [], []
    max_depth = []
    random_states = []

    for estimator in model.estimators_:
        tree = estimator.tree_
        nodes = np.asarray(tree.nodes)
        values = np.asarray(tree.values, dtype=np.float64)
        if nodes.dtype.names != (
            "left_child", "right_child", "feature", "threshold", "impurity",
            "n_node_samples", "weighted_n_node_samples",
        ):
            raise ValueError(f"Unexpected tree node dtype: {nodes.dtype}")
        if values.ndim != 3 or values.shape[1:] != (1, 2):
            raise ValueError(f"Unexpected tree value shape: {values.shape}")
        left.append(nodes["left_child"].astype(np.int64, copy=False))
        right.append(nodes["right_child"].astype(np.int64, copy=False))
        feature.append(nodes["feature"].astype(np.int64, copy=False))
        threshold.append(nodes["threshold"].astype(np.float64, copy=False))
        value0.append(values[:, 0, 0])
        value1.append(values[:, 0, 1])
        offsets.append(offsets[-1] + len(nodes))
        max_depth.append(int(tree.max_depth))
        random_states.append(int(estimator.random_state))

    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        offsets=np.asarray(offsets, dtype=np.int64),
        left_child=np.concatenate(left),
        right_child=np.concatenate(right),
        feature=np.concatenate(feature),
        threshold=np.concatenate(threshold),
        value_ns=np.concatenate(value0),
        value_hs=np.concatenate(value1),
        max_depth=np.asarray(max_depth, dtype=np.int32),
        estimator_random_state=np.asarray(random_states, dtype=np.int64),
        classes=np.asarray([0, 1], dtype=np.int8),
        n_features=np.asarray([173], dtype=np.int32),
    )
    output_hash = _sha256(output)
    info = {
        "model_identity": MODEL_VERSION,
        "source_pickle": source.name,
        "source_sha256": digest,
        "source_sklearn_version": "0.21.2",
        "source_pickle_protocol": 3,
        "estimator": "ExtraTreesClassifier",
        "n_estimators": 500,
        "n_features": 173,
        "classes": [0, 1],
        "parameters": {
            "criterion": "gini",
            "max_depth": None,
            "min_samples_split": 2,
            "min_samples_leaf": 1,
            "max_features": "auto",
            "bootstrap": True,
            "oob_score": False,
            "random_state": 42,
            "class_weight": "balanced_subsample",
        },
        "runtime_format": output.name,
        "runtime_sha256": output_hash,
        "conversion": "numeric tree arrays only; no refit and no parameter changes",
    }
    manifest.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    return info


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    info = convert(args.source, args.output, args.manifest)
    print(json.dumps(info, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
