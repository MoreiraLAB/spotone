"""Inference with the original SPOTONE ExtraTrees model deployed by MoreiraLAB.

The historical scikit-learn 0.21.2 pickle is preserved for provenance, but the
runtime uses a neutral NPZ export of its exact 500 tree structures. No tree is
refit and no threshold is changed. Feature extraction uses the original 173
sequence-only descriptors.
"""
from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .fasta import FastaRecord
from .features import FEATURE_COLUMNS, extract_features

PUBLISHED_MODEL_ID = "published-2020-extratrees-v1"
PUBLISHED_MODEL_SHA256 = "c6746cff126bb70218cd2d3ac91764ca629a9d405d415ed269e504e92712f3af"

# Exact post-model rules in the production deploy_model.py. The distinction
# between the gate increment and added increment is intentional: ASN and ILE
# have asymmetric historical rules, and PRO gates on +0.16 but adds +0.15.
ADAPTATION_RULES: dict[str, tuple[float, float]] = {
    "PRO": (0.16, 0.15),
    "PHE": (0.25, 0.25),
    "ASN": (0.20, 0.15),
    "LEU": (0.06, 0.06),
    "VAL": (0.25, 0.25),
    "ILE": (0.29, 0.15),
}


@dataclass(frozen=True)
class PublishedForest:
    offsets: np.ndarray
    left_child: np.ndarray
    right_child: np.ndarray
    feature: np.ndarray
    threshold: np.ndarray
    value_ns: np.ndarray
    value_hs: np.ndarray
    n_features: int

    @property
    def n_estimators(self) -> int:
        return int(len(self.offsets) - 1)

    def predict_proba(self, features: np.ndarray | pd.DataFrame) -> np.ndarray:
        """Return [P(NS), P(HS)] with the exact archived tree ensemble.

        Old sklearn decision trees converted input to float32 before traversal;
        we do the same so threshold comparisons reproduce that behavior.
        Tree probabilities are accumulated in forest order in float64 and then
        divided by 500, matching the historical single-worker estimator.
        """
        x = np.asarray(features, dtype=np.float32)
        if x.ndim != 2 or x.shape[1] != self.n_features:
            raise ValueError(f"Expected an (n, {self.n_features}) feature matrix; got {x.shape}.")
        if not np.isfinite(x).all():
            raise ValueError("Published-model features must be finite numeric values.")
        n_samples = x.shape[0]
        out = np.zeros((n_samples, 2), dtype=np.float64)
        for tree_index in range(self.n_estimators):
            start = int(self.offsets[tree_index])
            stop = int(self.offsets[tree_index + 1])
            nodes = np.zeros(n_samples, dtype=np.int64)
            active = np.ones(n_samples, dtype=bool)
            while active.any():
                local = nodes[active]
                global_nodes = start + local
                feats = self.feature[global_nodes]
                leaf = feats < 0
                if leaf.all():
                    active[np.flatnonzero(active)] = False
                    break
                active_rows = np.flatnonzero(active)
                leaf_rows = active_rows[leaf]
                active[leaf_rows] = False
                branch_rows = active_rows[~leaf]
                branch_global = global_nodes[~leaf]
                branch_feat = feats[~leaf]
                go_left = x[branch_rows, branch_feat] <= self.threshold[branch_global]
                next_local = np.where(
                    go_left,
                    self.left_child[branch_global],
                    self.right_child[branch_global],
                )
                nodes[branch_rows] = next_local
            leaves = start + nodes
            ns = self.value_ns[leaves]
            hs = self.value_hs[leaves]
            total = ns + hs
            total = np.where(total == 0.0, 1.0, total)
            out[:, 0] += ns / total
            out[:, 1] += hs / total
        out /= float(self.n_estimators)
        return out


def _data_file(name: str):
    return resources.files("spotone").joinpath("data", name)


@lru_cache(maxsize=4)
def load_published_forest(path: str | Path | None = None) -> PublishedForest:
    """Load the neutral numeric export bundled with SPOTONE."""
    if path is None:
        resource = _data_file("published_model_v1.npz")
        with resources.as_file(resource) as local:
            return load_published_forest(local)
    path = Path(path)
    manifest = published_model_manifest()
    if path.name == manifest["runtime_format"]:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != manifest["runtime_sha256"]:
            raise ValueError("Published-model runtime checksum mismatch.")
    with np.load(path, allow_pickle=False) as data:
        offsets = np.array(data["offsets"], dtype=np.int64, copy=True)
        n_features = int(np.asarray(data["n_features"]).ravel()[0])
        forest = PublishedForest(
            offsets=offsets,
            left_child=np.array(data["left_child"], dtype=np.int64, copy=True),
            right_child=np.array(data["right_child"], dtype=np.int64, copy=True),
            feature=np.array(data["feature"], dtype=np.int64, copy=True),
            threshold=np.array(data["threshold"], dtype=np.float64, copy=True),
            value_ns=np.array(data["value_ns"], dtype=np.float64, copy=True),
            value_hs=np.array(data["value_hs"], dtype=np.float64, copy=True),
            n_features=n_features,
        )
    if forest.n_estimators != 500 or forest.n_features != 173:
        raise ValueError("Unexpected published-model structure.")
    if forest.offsets[0] != 0 or forest.offsets[-1] != len(forest.feature):
        raise ValueError("Corrupt published-model tree offsets.")
    return forest


@lru_cache(maxsize=1)
def published_model_manifest() -> dict:
    with _data_file("published_model_manifest.json").open("r", encoding="utf-8") as handle:
        return json.load(handle)


def adapt_probability(raw_probability: float, residue_three_letter: str) -> tuple[float, int, int]:
    """Return corrected probability, historical raw class, and adapted class.

    This reproduces deploy_model.py exactly, including its boundary asymmetry:
    the raw class uses ``>= 0.50`` whereas the adapted class uses ``> 0.50``.
    The returned corrected probability is rounded to four decimals as in the
    historical CSV; class assignment uses the unrounded corrected value.
    """
    raw = float(raw_probability)
    if not 0.0 <= raw <= 1.0:
        raise ValueError("Raw probability must lie between zero and one.")
    original_class = 1 if raw >= 0.50 else 0
    gate, addition = ADAPTATION_RULES.get(residue_three_letter, (0.0, 0.0))
    adapted = raw + addition if gate and raw + gate >= 0.50 else raw
    adapted = min(adapted, 1.0)
    adapted_class = 1 if adapted > 0.50 else 0
    return round(adapted, 4), original_class, adapted_class


def predict_feature_frame(frame: pd.DataFrame, forest: PublishedForest | None = None) -> pd.DataFrame:
    """Predict a feature frame produced by :func:`extract_features`."""
    if frame.empty:
        raise ValueError("Prediction input contains no residues.")
    missing = [name for name in FEATURE_COLUMNS if name not in frame]
    if missing:
        raise ValueError(f"Prediction input is missing {len(missing)} published-model feature(s).")
    forest = forest or load_published_forest()
    matrix = frame.loc[:, FEATURE_COLUMNS].to_numpy(dtype=float)
    proba = forest.predict_proba(matrix)[:, 1]
    rows = []
    for idx, raw in enumerate(proba):
        three = str(frame.iloc[idx]["PDBResName"])
        corrected, original, adapted = adapt_probability(raw, three)
        rows.append((float(raw), corrected, original, adapted))
    values = np.asarray(rows, dtype=float)
    result = frame.loc[:, ["sequence_id", "position", "residue", "PDBResName"]].copy().reset_index(drop=True)
    result.rename(columns={"PDBResName": "residue_name"}, inplace=True)
    result["raw_probability"] = values[:, 0]
    result["prediction"] = values[:, 1]
    result["original_class"] = values[:, 2].astype(int)
    result["adapted_class"] = values[:, 3].astype(int)
    result["class"] = np.where(result["adapted_class"].to_numpy() == 1, "HS", "NS")
    return result


def predict_records(records: Iterable[FastaRecord], forest: PublishedForest | None = None) -> pd.DataFrame:
    return predict_feature_frame(extract_features(records), forest=forest)
