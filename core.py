"""Core logic: read FASTA files and list the positions where two sequences differ."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

UNAMBIGUOUS = frozenset("ACGT")
PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CT")


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record: its header line and its cleaned sequence."""

    header: str
    sequence: str


@dataclass(frozen=True)
class Variant:
    """A single-position difference between two sequences.

    ``position`` is 1-based, the convention used by genome browsers and NCBI.
    """

    position: int
    reference: str
    alternate: str

    @property
    def kind(self) -> str:
        """Classify the change as 'transition', 'transversion' or 'other'.

        Transitions swap two purines (A<->G) or two pyrimidines (C<->T);
        transversions swap a purine and a pyrimidine. Anything involving
        ambiguity codes (N, R, Y...) or gaps is reported as 'other'.
        """
        pair = {self.reference, self.alternate}
        if not pair <= UNAMBIGUOUS:
            return "other"
        if pair <= PURINES or pair <= PYRIMIDINES:
            return "transition"
        return "transversion"


@dataclass(frozen=True)
class Comparison:
    """Result of comparing a reference and a sample sequence."""

    reference_length: int
    sample_length: int
    compared_length: int
    variants: list[Variant] = field(default_factory=list)
    skipped_ambiguous: int = 0

    @property
    def evaluated_positions(self) -> int:
        """Positions that were actually judged (compared minus skipped)."""
        return self.compared_length - self.skipped_ambiguous

    @property
    def difference_rate(self) -> float:
        """Fraction of evaluated positions that differ (0.0 if none evaluated)."""
        if self.evaluated_positions == 0:
            return 0.0
        return len(self.variants) / self.evaluated_positions

    def count(self, kind: str) -> int:
        """Number of variants of a given kind ('transition', 'transversion', 'other')."""
        return sum(1 for v in self.variants if v.kind == kind)


def read_fasta(path: str | Path) -> FastaRecord:
    """Read the first record of a FASTA file.

    Whitespace and line breaks are removed from the sequence, and it is
    converted to upper case. Further records in the file are ignored.

    Raises:
        ValueError: if the file has no header line or the record has no sequence.
    """
    header: str | None = None
    chunks: list[str] = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break  # first record is complete
                header = line[1:].strip()
            elif header is None:
                raise ValueError(f"{path}: expected a FASTA header line starting with '>'")
            else:
                chunks.append("".join(line.split()))

    if header is None:
        raise ValueError(f"{path}: no FASTA record found")
    sequence = "".join(chunks).upper()
    if not sequence:
        raise ValueError(f"{path}: record '{header}' has an empty sequence")
    return FastaRecord(header=header, sequence=sequence)


def compare_sequences(
    reference: str, sample: str, *, ignore_ambiguous: bool = False
) -> Comparison:
    """Compare two sequences position by position.

    The comparison covers only the first ``min(len(reference), len(sample))``
    positions. It is *not* an alignment: an insertion or deletion shifts every
    following position and would be reported as many differences.

    Args:
        reference: reference sequence (case-insensitive).
        sample: sequence to compare against the reference (case-insensitive).
        ignore_ambiguous: if True, positions where either base is not one of
            A, C, G, T (for example N) are skipped instead of reported.
    """
    reference = reference.upper()
    sample = sample.upper()

    variants: list[Variant] = []
    skipped = 0
    for position, (ref_base, alt_base) in enumerate(zip(reference, sample), start=1):
        if ref_base == alt_base:
            continue
        if ignore_ambiguous and not {ref_base, alt_base} <= UNAMBIGUOUS:
            skipped += 1
            continue
        variants.append(Variant(position, ref_base, alt_base))

    return Comparison(
        reference_length=len(reference),
        sample_length=len(sample),
        compared_length=min(len(reference), len(sample)),
        variants=variants,
        skipped_ambiguous=skipped,
    )
