# Security notes

The normal SPOTONE published-model runtime does not deserialize the archived scikit-learn pickle. Inference uses the checksum-protected, non-executable NumPy representation in `src/spotone/data/published_model_v1.npz`.

The original `models/original/ensemble_model.pkl` is retained for provenance. Python pickle files can execute code during deserialization and must never be loaded when obtained from an untrusted source.

`scripts/convert_published_model.py` is intended only for the exact archived SPOTONE artifact identified by the repository checksum and uses a restricted class map during the one-time conversion.

FASTA input is treated as data and validated before feature extraction. Optional reconstruction commands may load local `joblib` model bundles; use only models from trusted sources.

The optional `ifeature` command executes a separately installed iFeature checkout. Users are responsible for validating the provenance of that external software before execution.
