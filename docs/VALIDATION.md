# Validation — SPOTONE v0.4

SPOTONE v0.4 modernizes the software used to evaluate the published predictor; it does not refit the recovered published model.

## Published-model equivalence

- Original model SHA-256: `c6746cff126bb70218cd2d3ac91764ca629a9d405d415ed269e504e92712f3af`.
- Original estimator: scikit-learn 0.21.2 `ExtraTreesClassifier`.
- Number of trees: 500.
- Number of input features: 173.
- Runtime NPZ SHA-256: `d73ffbd59684ccce83fea02fb2d32a129f71b8acd470ef29fb15ed45f41507d5`.
- No tree was refitted and no tree parameter was changed.

Direct traversal of the recovered pickle and traversal of the NPZ representation produced exactly identical probabilities on deterministic test matrices and real SPOTONE feature vectors (`max abs difference = 0.0`).

The modern and recovered production feature extractors agree to floating-point summation noise (`1.1368683772161603e-13` maximum absolute difference on the canonical 20-residue check). After the historical float32 conversion used before tree traversal, predictions are identical.

The detailed machine-readable equivalence record is `docs/validation/published-model-equivalence.json`.

## Golden example

`examples/example.fasta` contains 2 sequences and 460 residues.

The canonical output is stored at:

```text
docs/validation/published-example-prediction.csv
```

SHA-256:

```text
be8703529418b90d943f38ea783ac7a1b675092d90f67ca9b267bf35fc66dc4a
```

A wheel installed outside the source checkout reproduced this file byte-for-byte.

## Self-check

Run:

```bash
spotone self-check
```

The command verifies the model checksum and structure and evaluates a deterministic 20-residue fixture. A valid v0.4 runtime reports 500 estimators and 173 features.

## Automated tests

The public test suite covers FASTA handling, feature extraction, published-model inference, CLI behaviour, optional reconstruction workflows, and external descriptor integration.

```bash
python -m pip install '.[reconstruction]'
python -m unittest discover -s tests -v
```

GitHub Actions runs the suite on the supported Python versions and independently regenerates the golden published prediction.

## Interpretation

This validation demonstrates inference compatibility with the recovered published SPOTONE model. It does not constitute a new independent measurement of the predictive-performance metrics reported in the original publication.
