# SPOTONE

**SPOTONE** predicts protein–protein interface hot spots directly from protein sequence using an Extremely Randomized Trees (ExtraTrees) classifier and sequence-derived residue descriptors.

This repository contains the scientific implementation, the recovered published model, reproducibility assets, validation fixtures, and command-line workflows for SPOTONE.

> **Reference:** Preto AJ, Moreira IS. *SPOTONE: hot SPOTs ON protein complexes with Extremely randomized trees via sequence-only features.* International Journal of Molecular Sciences. 2020;21(19):7281. https://doi.org/10.3390/ijms21197281

## What is included

- the original published ExtraTrees model recovered from the SPOTONE deployment;
- a modern, non-executable NumPy representation of the same 500 trees for current Python environments;
- the 173 sequence-only descriptors used by the published predictor;
- exact historical residue-specific probability adjustments;
- a command-line interface for prediction, feature extraction, model inspection, and validation;
- optional model reconstruction and benchmarking workflows;
- automated tests and fixed regression fixtures for reproducibility.

The published model is **not refitted** in the default prediction workflow.

## Repository structure

```text
spotone/
├── .github/
│   └── workflows/
│       └── tests.yml                  # Continuous integration
├── docs/
│   ├── ARCHITECTURE.md                # Scientific software architecture
│   ├── PUBLISHED_MODEL.md             # Model provenance and runtime conversion
│   ├── VALIDATION.md                  # Equivalence and regression validation
│   └── validation/
│       ├── published-example-prediction.csv
│       ├── published-model-equivalence.json
│       └── published-model-info.json
├── examples/
│   └── example.fasta                  # Two-sequence, 460-residue example
├── models/
│   └── original/
│       ├── ensemble_model.pkl         # Archived original fitted model
│       └── SHA256SUMS.txt
├── datasets/
│   └── 1_homemade.csv                 # Original labelled training table
├── resources/
│   ├── amino_properties.csv           # Reference physicochemical properties
│   └── encoding.csv                   # Reference residue encoding table
├── scripts/
│   └── convert_published_model.py     # Trusted one-time pickle → NPZ conversion
├── src/
│   └── spotone/
│       ├── data/
│       │   ├── amino_properties.csv
│       │   ├── encoding.csv
│       │   ├── example.fasta
│       │   ├── published_model_manifest.json
│       │   ├── published_model_v1.npz
│       │   └── training.csv
│       ├── cli.py                     # Command-line interface
│       ├── external.py                # Optional external iFeature integration
│       ├── fasta.py                   # FASTA parsing and validation
│       ├── features.py                # 173-descriptor feature extraction
│       ├── model.py                   # Optional reconstruction workflows
│       └── published.py               # Published-model inference
├── tests/                              # Unit and regression tests
├── CITATION.cff
├── pyproject.toml
├── requirements.txt
└── README.md
```

Historical public code before the 2026 modernization is preserved in the Git tag `legacy-public-pre-v0.4.0`.

## Requirements

- Python 3.10 or newer
- NumPy
- pandas

The published SPOTONE predictor does **not** require scikit-learn at runtime. Scikit-learn and joblib are only required for the optional reconstruction workflows.

## Installation

Clone the repository and create an isolated Python environment:

```bash
git clone https://github.com/MoreiraLAB/spotone.git
cd spotone

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

Verify the installation:

```bash
spotone --version
spotone self-check
spotone published-model-info
```

A successful `self-check` verifies the bundled model checksum, the expected 500-tree/173-feature structure, and a deterministic prediction fixture.

## Predict hot spots with the published model

Input sequences must be provided in FASTA format using standard one-letter amino-acid codes.

```bash
spotone predict-published examples/example.fasta -o prediction.csv
```

The output contains one row per residue:

| Column | Meaning |
| --- | --- |
| `sequence_id` | FASTA sequence identifier |
| `position` | 1-based residue position |
| `residue` | one-letter amino-acid code |
| `residue_name` | three-letter amino-acid code |
| `raw_probability` | ExtraTrees hot-spot probability before the historical residue-specific adjustment |
| `prediction` | adjusted probability reported by SPOTONE |
| `original_class` | class obtained from the raw probability |
| `adapted_class` | class after the historical residue-specific adjustment |
| `class` | `HS` (hot spot) or `NS` (non-hot spot) |

Example:

```text
sequence_id,position,residue,residue_name,raw_probability,prediction,original_class,adapted_class,class
6Q1G:H|PDBID|CHAIN|SEQUENCE,1,A,ALA,0.1715448717948718,0.1715,0,0,NS
```

## Extract the 173 SPOTONE descriptors

```bash
spotone features examples/example.fasta -o features.csv
```

For each residue, the published representation contains:

- 20 one-hot amino-acid indicators;
- 1 relative-position descriptor;
- 19 amino-acid physicochemical properties;
- 19 property means evaluated over 7 centered sequence windows with radii 2, 5, 7, 10, 25, 50, and 75.

This gives:

```text
20 + 1 + 19 + (19 × 7) = 173 descriptors
```

## Published model

The recovered model is a scikit-learn 0.21.2 `ExtraTreesClassifier` with:

- 500 trees;
- 173 input features;
- classes `[0, 1]`;
- `criterion='gini'`;
- `bootstrap=True`;
- `random_state=42`;
- `class_weight='balanced_subsample'`.

The original fitted pickle is retained unchanged at:

```text
models/original/ensemble_model.pkl
```

SHA-256:

```text
c6746cff126bb70218cd2d3ac91764ca629a9d405d415ed269e504e92712f3af
```

Because the original artifact was serialized with an old scikit-learn version, normal prediction does not unpickle it. The exact tree topology, thresholds, and leaf class weights were exported once to:

```text
src/spotone/data/published_model_v1.npz
```

SHA-256:

```text
d73ffbd59684ccce83fea02fb2d32a129f71b8acd470ef29fb15ed45f41507d5
```

No tree was refitted and no model parameter was changed during this conversion. See [`docs/PUBLISHED_MODEL.md`](docs/PUBLISHED_MODEL.md) for details.

## Historical probability adjustment

The published SPOTONE workflow applies residue-specific post-model rules for PRO, PHE, ASN, LEU, VAL, and ILE. Version 0.4 reproduces these rules exactly, including the original threshold asymmetries, because changing them would change published-model behaviour.

The raw model probability and the adjusted prediction are both included in the output so the transformation remains explicit.

## Validation and reproducibility

The modern runtime was checked against direct traversal of the recovered original model. The recorded maximum absolute probability difference is:

```text
0.0
```

The canonical example contains 2 sequences and 460 residues. Its expected prediction file is:

```text
docs/validation/published-example-prediction.csv
```

with SHA-256:

```text
be8703529418b90d943f38ea783ac7a1b675092d90f67ca9b267bf35fc66dc4a
```

Reproduce it with:

```bash
spotone predict-published examples/example.fasta -o /tmp/spotone-example.csv
shasum -a 256 /tmp/spotone-example.csv
```

Detailed validation records are available in [`docs/VALIDATION.md`](docs/VALIDATION.md) and `docs/validation/`.

## Optional reconstruction workflows

The published predictor above is the recommended workflow. Model reconstruction, classifier comparison, and training utilities are provided separately and require scikit-learn:

```bash
python -m pip install '.[reconstruction]'
```

Examples:

```bash
# Reconstruct an ExtraTrees model from the bundled training table
spotone train --output-dir results/model

# Compare supported classifiers on grouped splits
spotone benchmark --output-dir results/benchmark

# Train a reconstruction and run the bundled FASTA example
spotone demo --output-dir results/demo
```

These workflows fit new models and must not be confused with inference using the recovered published model.

## Run the tests

Install the optional reconstruction dependencies and run the complete public test suite:

```bash
python -m pip install '.[reconstruction]'
python -m unittest discover -s tests -v
```

The CI workflow also verifies the published-model self-check and reproduces the golden example prediction.

## Scientific scope and limitations

SPOTONE is a sequence-only predictor of protein–protein interface hot spots. Its score is a model probability/classification, not a binding free energy.

The original SPOTONE training data did not contain membrane proteins. Predictions for membrane proteins and other systems substantially outside the training domain should therefore be interpreted cautiously.

Validation in this repository establishes compatibility with the recovered published predictor. It does not constitute a new independent reassessment of the predictive-performance values reported in the original paper.

## Citation

If you use SPOTONE, please cite:

> Preto AJ, Moreira IS. **SPOTONE: hot SPOTs ON protein complexes with Extremely randomized trees via sequence-only features.** *International Journal of Molecular Sciences*. 2020;21(19):7281. https://doi.org/10.3390/ijms21197281

Citation metadata is also provided in [`CITATION.cff`](CITATION.cff).
