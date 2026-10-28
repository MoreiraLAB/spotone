"""The 173 in-house SPOTONE descriptors from the archived extraction code.

Window descriptors are ordered by property, then by radius. Historical CSV
column names ``3`` through ``135`` are retained to make the model schema
explicit and compatible with the supplied training table.
"""

from collections.abc import Iterable
from functools import lru_cache
from importlib import resources

import numpy as np
import pandas as pd

from .fasta import AMINO_ACIDS, FastaRecord, validate_record


WINDOW_RADII = (2, 5, 7, 10, 25, 50, 75)
PROPERTY_COLUMNS = (
    "Helix_propensity", "Sheet_propensity", "Helix_propensity_values",
    "Sheet_propensity_values", "MW", "pka_carboxylate", "pka_amine",
    "pka_side_chain", "C_number", "H_number", "N_number", "O_number",
    "S_number", "Standard_area_free", "Standard_area_protein",
    "folded_buried_area", "mean_fractional_area_loss", "residue_mass",
    "monoisotopic_mass",
)
FEATURE_COLUMNS = (
    [f"binary_{index}" for index in range(1, 21)]
    + ["rel_dist_1"]
    + list(PROPERTY_COLUMNS)
    + [str(index) for index in range(3, 136)]
)
METADATA_COLUMNS = ["sequence_id", "position", "residue", "PDBResName"]
FEATURE_SCHEMA = "spotone-inhouse-173-v1"


@lru_cache(maxsize=1)
def _tables() -> tuple[dict[str, int], np.ndarray, np.ndarray, dict[str, str]]:
    """Load versioned resource snapshots independent of the working directory."""
    root = resources.files("spotone").joinpath("data")
    with root.joinpath("amino_properties.csv").open("r", encoding="utf-8") as handle:
        properties = pd.read_csv(handle)
    with root.joinpath("encoding.csv").open("r", encoding="utf-8") as handle:
        encoding = pd.read_csv(handle)
    if list(properties.columns[3:]) != list(PROPERTY_COLUMNS):
        raise ValueError("Bundled amino acid property columns do not match the feature schema.")
    if list(encoding.columns[1:]) != list(AMINO_ACIDS):
        raise ValueError("Bundled one-hot encoding columns do not match the feature schema.")
    letters = list(AMINO_ACIDS)
    if sorted(properties["res_letter"]) != letters or sorted(encoding["res_letter"]) != letters:
        raise ValueError("Bundled resources must contain each canonical amino acid exactly once.")
    properties = properties.set_index("res_letter").loc[letters]
    encoding = encoding.set_index("res_letter").loc[letters]
    property_values = properties.loc[:, list(PROPERTY_COLUMNS)].to_numpy(dtype=float)
    encoding_values = encoding.to_numpy(dtype=float)
    if not np.isfinite(property_values).all() or not np.array_equal(encoding_values, np.eye(20)):
        raise ValueError("Bundled amino acid resources contain invalid numeric values.")
    property_values.setflags(write=False)
    encoding_values.setflags(write=False)
    return (
        {letter: index for index, letter in enumerate(letters)},
        property_values,
        encoding_values,
        properties["three_code"].to_dict(),
    )


def extract_features(records: Iterable[FastaRecord]) -> pd.DataFrame:
    """Return one row per residue with metadata and the exact 173 descriptors.

    ``position`` is a 1-based position in the input FASTA sequence. The legacy
    ``PDBResName`` column is the three-letter residue name only; it makes no
    claim about a structure or PDB numbering. Window means include the central
    residue and clip at sequence boundaries, exactly as the original code did.
    """
    letter_indices, property_table, encoding_table, three_codes = _tables()
    frames: list[pd.DataFrame] = []
    seen: set[str] = set()
    for raw_record in records:
        record = validate_record(raw_record)
        if record.identifier in seen:
            raise ValueError(f"Duplicate FASTA identifier {record.identifier!r}.")
        seen.add(record.identifier)
        sequence = record.sequence
        length = len(sequence)
        indices = np.fromiter((letter_indices[letter] for letter in sequence), dtype=np.int64)
        properties = property_table[indices]
        position = np.arange(1, length + 1, dtype=np.int64)

        # The original quartiles use strict < comparisons at 0.25, 0.5, 0.75.
        # Integer arithmetic gives the same boundary behavior without rounding.
        quartile = np.minimum((4 * position) // length + 1, 4).astype(float)

        # Prefix sums compute all clipped centered windows in O(length), rather
        # than revisiting every residue in every window for every property.
        prefix = np.vstack((np.zeros((1, len(PROPERTY_COLUMNS))), properties.cumsum(axis=0)))
        windows = np.empty((length, len(PROPERTY_COLUMNS), len(WINDOW_RADII)), dtype=float)
        centers = position - 1
        for window_index, radius in enumerate(WINDOW_RADII):
            starts = np.maximum(centers - radius, 0)
            stops = np.minimum(centers + radius + 1, length)
            windows[:, :, window_index] = (prefix[stops] - prefix[starts]) / (stops - starts)[:, None]
        values = np.column_stack((encoding_table[indices], quartile, properties, windows.reshape(length, -1)))
        frame = pd.DataFrame(values, columns=FEATURE_COLUMNS)
        frame.insert(0, "PDBResName", [three_codes[letter] for letter in sequence])
        frame.insert(0, "residue", list(sequence))
        frame.insert(0, "position", position)
        frame.insert(0, "sequence_id", record.identifier)
        frames.append(frame)
    if not frames:
        raise ValueError("At least one protein sequence is required to extract features.")
    return pd.concat(frames, ignore_index=True)
