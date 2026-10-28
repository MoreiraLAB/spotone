"""Explicit commands for extraction, training, evaluation and prediction."""

import argparse
from contextlib import ExitStack
from importlib import resources
import importlib.metadata
import json
from pathlib import Path
import sys

import pandas as pd

from . import __version__
from .fasta import read_fasta
from .features import FEATURE_COLUMNS, FEATURE_SCHEMA, extract_features


MODEL_NAMES = ("extra-trees", "adaboost", "svm", "mlp")


def _training_options(parser, *, output_default):
    parser.add_argument("--output-dir", type=Path, default=Path(output_default))
    parser.add_argument("--seed", type=int, default=42, help="Explicit reconstruction seed (default: 42).")
    parser.add_argument("--trees", type=int, default=500, help="ExtraTrees estimator count (default: 500).")
    parser.add_argument("--jobs", type=int, default=1, help="Worker count, or -1 for all CPUs (default: 1).")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="spotone",
        description="SPOTONE: sequence-only hot-spot classification with explicit model provenance.",
    )
    parser.add_argument("--version", action="version", version=f"spotone {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    feature = commands.add_parser("features", help="Extract the original 173 descriptors from protein FASTA.")
    feature.add_argument("fasta", type=Path)
    feature.add_argument("-o", "--output", type=Path, required=True)

    train = commands.add_parser("train", help="Train and evaluate a reproducible reconstructed model.")
    train.add_argument("--dataset", type=Path, help="Labelled CSV; default is the bundled, unchanged 1_homemade.csv.")
    train.add_argument("--classifier", choices=MODEL_NAMES, default="extra-trees")
    train.add_argument("--tune", action="store_true", help="Bounded ExtraTrees search on training groups only.")
    _training_options(train, output_default="results/model")

    for name, description, source in [
        ("predict", "Predict one result per FASTA residue with a saved in-house model.", "fasta"),
        ("predict-features", "Predict from an existing CSV with the saved model's feature schema.", "features"),
    ]:
        predict = commands.add_parser(name, help=description)
        predict.add_argument(source, type=Path)
        predict.add_argument("--model", type=Path, required=True, help="Trusted local SPOTONE model.joblib bundle.")
        predict.add_argument("-o", "--output", type=Path, required=True)

    benchmark = commands.add_parser("benchmark", help="Compare classifiers on explicit, identical grouped splits.")
    benchmark.add_argument("--datasets", nargs="+", type=Path, help="CSV tables; default is bundled in-house data.")
    benchmark.add_argument("--models", nargs="+", choices=MODEL_NAMES, default=list(MODEL_NAMES))
    _training_options(benchmark, output_default="results/benchmark")

    demo = commands.add_parser("demo", help="Train a reconstruction and predict the bundled example FASTA.")
    _training_options(demo, output_default="results/demo")

    commands.add_parser("doctor", help="Check dependencies, bundled data and feature extraction.")

    published = commands.add_parser("predict-published", help="Predict FASTA residues with the original published SPOTONE model.")
    published.add_argument("fasta", type=Path)
    published.add_argument("-o", "--output", type=Path, required=True)

    commands.add_parser("published-model-info", help="Show provenance of the bundled original SPOTONE model.")
    commands.add_parser("self-check", help="Verify the bundled published-model runtime and a deterministic prediction fixture.")

    external = commands.add_parser("ifeature", help="Run an optional, separately installed iFeature checkout.")
    external.add_argument("fasta", type=Path)
    external.add_argument("--ifeature-dir", type=Path, required=True)
    external.add_argument("--output-dir", type=Path, default=Path("results/ifeature"))
    external.add_argument("--features", nargs="+", help="Standard descriptor families; default: all 21 original families.")
    external.add_argument("--no-standard", action="store_true", help="Run only the requested PseKRAAC families.")
    external.add_argument("--pseudo-features", nargs="+", help="Optional PseKRAAC types, e.g. type1 type2.")
    external.add_argument("--python", default=sys.executable, help="Interpreter for the optional iFeature environment.")
    return parser


def _save_csv(frame, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    print(f"Wrote {len(frame):,} rows to {output}")


def _protect_inputs(output, *inputs):
    output = Path(output)
    for path in inputs:
        if path is None:
            continue
        source = Path(path)
        if output.resolve() == source.resolve() or (
            output.exists() and source.exists() and output.samefile(source)
        ):
            raise ValueError("The output path must differ from the input data and model paths.")


def _training_path(stack, path=None):
    if path is not None:
        return path
    return stack.enter_context(resources.as_file(resources.files("spotone").joinpath("data", "training.csv")))


def _train(path, output, args, **kwargs):
    try:
        from .model import train_model
    except ModuleNotFoundError as exc:
        raise RuntimeError("Reconstruction commands require: python -m pip install 'spotone[reconstruction]'.") from exc
    print(f"Training a reconstructed model from {path.name}; report: {output / 'report.json'}", file=sys.stderr)
    return train_model(path, output, seed=args.seed, n_estimators=args.trees, n_jobs=args.jobs, **kwargs)


def _run(args):
    if args.command == "features":
        _protect_inputs(args.output, args.fasta)
        _save_csv(extract_features(read_fasta(args.fasta)), args.output)
    elif args.command == "train":
        with ExitStack() as stack:
            dataset = _training_path(stack, args.dataset)
            _train(dataset, args.output_dir, args, tune=args.tune, model_name=args.classifier)
        print(f"Model: {args.output_dir / 'model.joblib'}")
        print(f"Evaluation: {args.output_dir / 'report.json'}")
    elif args.command in ("predict", "predict-features"):
        try:
            from .model import load_model, predict_features
        except ModuleNotFoundError as exc:
            raise RuntimeError("Reconstruction commands require: python -m pip install 'spotone[reconstruction]'.") from exc
        input_path = args.fasta if args.command == "predict" else args.features
        model_path = args.model / "model.joblib" if args.model.is_dir() else args.model
        _protect_inputs(args.output, input_path, model_path)
        bundle = load_model(args.model)
        if args.command == "predict":
            if bundle["feature_schema"] != FEATURE_SCHEMA or set(bundle["feature_columns"]) != set(FEATURE_COLUMNS):
                raise ValueError("This model uses an external feature schema. Supply matching CSV data with predict-features.")
            frame = extract_features(read_fasta(input_path))
        else:
            frame = pd.read_csv(input_path)
        _save_csv(predict_features(frame, bundle), args.output)
    elif args.command == "demo":
        try:
            from .model import load_model, predict_features
        except ModuleNotFoundError as exc:
            raise RuntimeError("Reconstruction commands require: python -m pip install 'spotone[reconstruction]'.") from exc
        args.output_dir.mkdir(parents=True, exist_ok=True)
        with ExitStack() as stack:
            dataset = _training_path(stack)
            fasta = stack.enter_context(resources.as_file(resources.files("spotone").joinpath("data", "example.fasta")))
            example = args.output_dir / "example.fasta"
            example.write_text(fasta.read_text(encoding="utf-8"), encoding="utf-8")
            model_dir = args.output_dir / "model"
            _train(dataset, model_dir, args)
            features = extract_features(read_fasta(example))
            _save_csv(predict_features(features, load_model(model_dir / "model.joblib")), args.output_dir / "predictions.csv")
        print(f"Evaluation: {model_dir / 'report.json'}")
        print("Example predictions demonstrate execution; they are not an independent validation dataset.")
    elif args.command == "benchmark":
        reports = []
        with ExitStack() as stack:
            datasets = args.datasets or [_training_path(stack)]
            for index, dataset in enumerate(datasets, 1):
                for name in dict.fromkeys(args.models):
                    output = args.output_dir / f"{index}_{dataset.stem}" / name
                    report = _train(dataset, output, args, model_name=name)
                    reports.append({"dataset": str(dataset), "classifier": name, "report": report})
        args.output_dir.mkdir(parents=True, exist_ok=True)
        target = args.output_dir / "benchmark.json"
        target.write_text(json.dumps(reports, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(f"Comparison: {target}")
    elif args.command == "doctor":
        versions = {}
        for name in ("numpy", "pandas", "scikit-learn", "joblib"):
            try:
                versions[name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                versions[name] = None
        with ExitStack() as stack:
            dataset = pd.read_csv(_training_path(stack))
        from .fasta import FastaRecord
        frame = extract_features([FastaRecord("check", "ACDEFGHIKLMNPQRSTVWY")])
        print(json.dumps({"spotone": __version__, "python": sys.version.split()[0], "dependencies": versions,
                          "bundled_training_rows": len(dataset), "features_per_residue": len(FEATURE_COLUMNS),
                          "feature_check_rows": len(frame)}, indent=2))
    elif args.command == "predict-published":
        from .published import predict_records
        _protect_inputs(args.output, args.fasta)
        _save_csv(predict_records(read_fasta(args.fasta)), args.output)
    elif args.command == "published-model-info":
        from .published import published_model_manifest
        print(json.dumps(published_model_manifest(), indent=2))
    elif args.command == "self-check":
        from .fasta import FastaRecord
        from .published import load_published_forest, predict_records, published_model_manifest
        forest = load_published_forest()
        result = predict_records([FastaRecord("self-check", "ACDEFGHIKLMNPQRSTVWY")], forest=forest)
        expected_first = [
            0.1517391304347826,
            0.148,
            0.228,
            0.142,
            0.1738472906403941,
        ]
        actual_first = [float(value) for value in result.loc[:4, "raw_probability"]]
        expected_last = 0.5213441976442006
        import numpy as np
        np.testing.assert_allclose(actual_first, expected_first, rtol=0, atol=1e-15)
        np.testing.assert_allclose(float(result.iloc[-1]["raw_probability"]), expected_last, rtol=0, atol=1e-15)
        if int(result["adapted_class"].sum()) != 1:
            raise RuntimeError("Published-model self-check hotspot count mismatch.")
        manifest = published_model_manifest()
        print(json.dumps({
            "status": "ok",
            "spotone": __version__,
            "model_identity": manifest["model_identity"],
            "runtime_sha256": manifest["runtime_sha256"],
            "n_estimators": forest.n_estimators,
            "n_features": forest.n_features,
            "fixture_rows": len(result),
        }, indent=2))
    elif args.command == "ifeature":
        from .external import run_ifeature
        if args.no_standard and args.features:
            raise ValueError("Use either --features or --no-standard, not both.")
        paths = run_ifeature(args.fasta, args.ifeature_dir, args.output_dir, features=[] if args.no_standard else args.features,
                            pseudo_features=args.pseudo_features, python_executable=args.python)
        print(f"Wrote {len(paths)} descriptor tables to {args.output_dir}")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        _run(args)
    except (ValueError, OSError, RuntimeError, pd.errors.ParserError, pd.errors.EmptyDataError) as error:
        parser.exit(2, f"spotone: error: {error}\n")
    return 0
