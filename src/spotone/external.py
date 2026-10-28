"""Optional, explicitly invoked integration with the upstream iFeature checkout.

These descriptors are not required by the 173-feature SPOTONE predictor. PCA is
intentionally not run here: feature transforms need to be fitted on training
data and saved before they can be applied consistently to validation or input.
"""

from collections.abc import Sequence
from pathlib import Path
import subprocess
import sys
import tempfile

from .fasta import read_fasta


STANDARD_FEATURES = (
    "AAC", "CKSAAP", "DPC", "DDE", "TPC", "GAAC", "CKSAAGP", "GDPC",
    "GTPC", "NMBroto", "Moran", "Geary", "CTDC", "CTDT", "CTDD", "CTriad",
    "KSCTriad", "SOCNumber", "QSOrder", "PAAC", "APAAC",
)
PSEUDO_FEATURES = (
    "type1", "type2", "type3A", "type3B", "type4", "type6A", "type6B",
    "type6C", "type7", "type8", "type9", "type10", "type11", "type12",
    "type14", "type15", "type16",
)


def _selected_features(names: Sequence[str], allowed: tuple[str, ...], kind: str) -> list[str]:
    if isinstance(names, (str, bytes)):
        raise RuntimeError(f"{kind} features must be a list of names, for example ['{allowed[0]}'].")
    selected = []
    for name in names:
        if name not in allowed:
            raise RuntimeError(
                f"Unknown {kind} feature {name!r}. Supported names: {', '.join(allowed)}."
            )
        if name not in selected:
            selected.append(name)
    return selected


def _process_message(stdout: str | None, stderr: str | None) -> str:
    message = "\n".join(part.strip() for part in (stderr, stdout) if part and part.strip())
    return f" Upstream message: {message[-2000:]}" if message else ""


def run_ifeature(
    fasta_path: str | Path,
    ifeature_dir: str | Path,
    output_dir: str | Path,
    *,
    features: Sequence[str] | None = None,
    pseudo_features: Sequence[str] | None = None,
    python_executable: str = sys.executable,
) -> list[Path]:
    """Run selected upstream iFeature families and return absolute TSV paths.

    The default runs the original 21 standard feature families. PseKRAAC
    families are opt-in; pass ``features=[]`` for a PseKRAAC-only run. PseKRAAC
    uses the historical SPOTONE settings ``--gap_lambda 5 --raactype 5``.

    Each subprocess receives a list of arguments and an explicit working
    directory. All output is staged in a fresh temporary directory; previous
    results are replaced only after every selected family succeeds and creates
    a nonempty file. A successful process exit alone is not treated as success.
    """
    standard = _selected_features(STANDARD_FEATURES if features is None else features, STANDARD_FEATURES, "standard")
    pseudo = _selected_features(() if pseudo_features is None else pseudo_features, PSEUDO_FEATURES, "PseKRAAC")
    if not standard and not pseudo:
        raise RuntimeError("Select at least one standard or PseKRAAC iFeature family.")
    if not isinstance(python_executable, (str, Path)) or not str(python_executable):
        raise RuntimeError("Specify a Python executable for the iFeature checkout.")

    fasta = Path(fasta_path).expanduser().resolve()
    checkout = Path(ifeature_dir).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    if not fasta.is_file():
        raise RuntimeError(f"The input FASTA file does not exist: {fasta}")
    try:
        read_fasta(fasta)
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Cannot use the input FASTA file: {exc}") from exc
    if not checkout.is_dir():
        raise RuntimeError(f"The iFeature checkout directory does not exist: {checkout}")

    jobs = []
    for names, script_name, extra_arguments in (
        (standard, "iFeature.py", []),
        (pseudo, "iFeaturePseKRAAC.py", ["--gap_lambda", "5", "--raactype", "5"]),
    ):
        if names and not (checkout / script_name).is_file():
            raise RuntimeError(
                f"Cannot find {script_name} in {checkout}. Supply the directory "
                "of an installed https://github.com/Superzchen/iFeature checkout."
            )
        jobs.extend((name, checkout / script_name, extra_arguments) for name in names)

    if any((output / f"{name}.tsv").resolve() == fasta for name, _, _ in jobs):
        raise RuntimeError("An iFeature output would overwrite the input FASTA. Choose a different output directory.")

    try:
        output.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".spotone-ifeature-", dir=output) as staging_dir:
            staging = Path(staging_dir)
            results = []
            for name, script, extra_arguments in jobs:
                staged_path = staging / f"{name}.tsv"
                arguments = [
                    str(python_executable), str(script), "--file", str(fasta),
                    "--type", name, "--out", str(staged_path), *extra_arguments,
                ]
                try:
                    completed = subprocess.run(
                        arguments,
                        cwd=checkout,
                        check=True,
                        capture_output=True,
                        text=True,
                        encoding="utf-8",
                        errors="replace",
                        stdin=subprocess.DEVNULL,
                    )
                except subprocess.CalledProcessError as exc:
                    raise RuntimeError(
                        f"iFeature family {name} failed with exit code {exc.returncode}."
                        + _process_message(exc.stdout, exc.stderr)
                    ) from exc
                except OSError as exc:
                    raise RuntimeError(f"Could not launch iFeature family {name}: {exc}") from exc
                if not staged_path.is_file() or staged_path.stat().st_size == 0:
                    raise RuntimeError(
                        f"iFeature family {name} exited successfully but did not create "
                        "a nonempty output file. Check the upstream dependencies and "
                        "sequence-length requirements."
                        + _process_message(completed.stdout, completed.stderr)
                    )
                results.append((staged_path, output / f"{name}.tsv"))
            for staged_path, final_path in results:
                staged_path.replace(final_path)
            return [final_path for _, final_path in results]
    except OSError as exc:
        raise RuntimeError(f"Could not write iFeature results to {output}: {exc}") from exc
