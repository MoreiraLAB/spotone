# Published model provenance

The production SPOTONE deployment loaded `support/ensemble_model.pkl` and called `predict_proba()` on 173 sequence-only descriptors. The recovered artifact has SHA-256:

`c6746cff126bb70218cd2d3ac91764ca629a9d405d415ed269e504e92712f3af`

Static/restricted inspection identifies it as a scikit-learn 0.21.2 `ExtraTreesClassifier` with 500 estimators, 173 features, classes `[0, 1]`, `random_state=42`, `bootstrap=True`, and `class_weight='balanced_subsample'`.

## Runtime migration

Modern scikit-learn cannot directly deserialize the 0.21.2 tree state because the current tree node dtype contains an additional `missing_go_to_left` field. Rather than run an end-of-life Python/scikit-learn stack in production, `scripts/convert_published_model.py` reads only the trusted known pickle through a restricted class map and exports the tree topology, thresholds and leaf class weights to a non-executable NumPy NPZ.

The production runtime loads only that NPZ. The original pickle remains unchanged as provenance.

## Equivalence

`docs/validation/published-model-equivalence.json` verifies that direct traversal of all 500 original pickle trees and traversal of the NPZ export produce exactly identical probabilities on both a deterministic random matrix and actual 173-feature vectors.

The input matrix is cast to `float32` before tree traversal because historical scikit-learn decision-tree prediction used float32 input internally. Leaf class weights and forest accumulation remain float64.
