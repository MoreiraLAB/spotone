# Architecture

SPOTONE v0.4 separates the published scientific predictor from optional reconstruction workflows.

## Published inference path

`src/spotone/fasta.py` validates and parses FASTA input.

`src/spotone/features.py` computes the 173 sequence-only descriptors required by the published model.

`src/spotone/published.py` loads the checksum-protected numeric tree representation, evaluates the 500-tree ensemble, and applies the historical residue-specific probability rules.

The published inference path does not fit a model and does not require scikit-learn at runtime.

## Model artifacts

`models/original/ensemble_model.pkl` is the archived fitted scikit-learn 0.21.2 model and is retained for provenance.

`src/spotone/data/published_model_v1.npz` contains the exact numeric tree topology, thresholds, and leaf class weights used for normal inference. Its checksum is validated against `published_model_manifest.json`.

`scripts/convert_published_model.py` implements the restricted, one-time conversion from the trusted archived model to the neutral NPZ representation.

## Reconstruction path

`src/spotone/model.py` contains optional training, evaluation, and benchmarking workflows. These require the `reconstruction` dependency extra and fit new models; their outputs are distinct from the recovered published predictor.

`src/spotone/external.py` provides optional integration with a separately installed iFeature checkout for legacy descriptor workflows.

## Compatibility boundary

The scientific compatibility boundary consists of:

1. FASTA parsing and residue indexing;
2. the ordered 173-feature representation;
3. float32 input conversion before tree traversal;
4. the 500 recovered tree structures and leaf class weights;
5. historical residue-specific post-model probability rules.

Changes to these components can alter published-model predictions and therefore require explicit scientific validation and a new compatibility assessment.
