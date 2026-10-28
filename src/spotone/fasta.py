"""Read protein FASTA without guessing residue identities or numbering.

The first whitespace-delimited word of each FASTA header is its identifier.
Positions produced by SPOTONE always refer to the supplied sequence, starting
at one; a PDB-looking identifier does not change that coordinate system.
"""

from dataclasses import dataclass
from pathlib import Path


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


@dataclass(frozen=True)
class FastaRecord:
    """One uniquely named protein sequence."""

    identifier: str
    sequence: str


def validate_record(record: FastaRecord) -> FastaRecord:
    """Validate a record and return its sequence in upper case.

    Unknown, ambiguous, gap and stop symbols are rejected rather than removed:
    removing them would silently change all subsequent residue positions.
    """
    if not isinstance(record, FastaRecord):
        raise ValueError("Expected FastaRecord objects when extracting features.")
    if not isinstance(record.identifier, str) or not record.identifier.strip():
        raise ValueError("A FASTA record has an empty identifier.")
    if any(char.isspace() for char in record.identifier):
        raise ValueError("FASTA record identifiers must not contain whitespace.")
    if not isinstance(record.sequence, str) or not record.sequence:
        raise ValueError(f"FASTA record {record.identifier!r} has an empty sequence.")
    # Validate before case conversion: Unicode characters such as sharp-s can
    # expand into canonical-looking letters when uppercased, changing positions.
    for position, residue in enumerate(record.sequence, 1):
        if residue not in AMINO_ACIDS + AMINO_ACIDS.lower():
            raise ValueError(
                f"Invalid amino acid {residue!r} in {record.identifier!r} at "
                f"sequence position {position}. Only the 20 canonical amino "
                f"acids ({AMINO_ACIDS}) are supported."
            )
    return FastaRecord(record.identifier, record.sequence.upper())


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read one or more FASTA sequences, rejecting malformed input.

    Blank lines, wrapped sequences, whitespace within sequence lines, and
    lowercase amino acids are accepted. Header descriptions are ignored.
    Duplicate first-word identifiers are an error, even if descriptions differ.
    """
    records: list[FastaRecord] = []
    seen: set[str] = set()
    identifier: str | None = None
    chunks: list[str] = []

    def finish_record() -> None:
        if identifier is not None:
            records.append(validate_record(FastaRecord(identifier, "".join(chunks))))

    with Path(path).open("r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                finish_record()
                header = stripped[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}.")
                identifier = header.split()[0]
                if identifier in seen:
                    raise ValueError(f"Duplicate FASTA identifier {identifier!r} at line {line_number}.")
                seen.add(identifier)
                chunks = []
            else:
                if identifier is None:
                    raise ValueError(
                        f"Sequence data before the first FASTA header at line {line_number}."
                    )
                chunks.append("".join(stripped.split()))
    finish_record()
    if not records:
        raise ValueError("The FASTA file contains no protein sequences.")
    return records


def parse_fasta_text(text: str) -> list[FastaRecord]:
    """Parse FASTA content already held in memory using the same validation as read_fasta."""
    import io

    if not isinstance(text, str) or not text.strip():
        raise ValueError("The FASTA input is empty.")
    records: list[FastaRecord] = []
    seen: set[str] = set()
    identifier: str | None = None
    chunks: list[str] = []

    def finish_record() -> None:
        nonlocal chunks
        if identifier is not None:
            records.append(validate_record(FastaRecord(identifier, "".join(chunks))))

    for line_number, line in enumerate(io.StringIO(text), 1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            finish_record()
            header = stripped[1:].strip()
            if not header:
                raise ValueError(f"Empty FASTA header at line {line_number}.")
            identifier = header.split()[0]
            if identifier in seen:
                raise ValueError(f"Duplicate FASTA identifier {identifier!r} at line {line_number}.")
            seen.add(identifier)
            chunks = []
        else:
            if identifier is None:
                raise ValueError(f"Sequence data before the first FASTA header at line {line_number}.")
            chunks.append("".join(stripped.split()))
    finish_record()
    if not records:
        raise ValueError("The FASTA input contains no protein sequences.")
    return records
