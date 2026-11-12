"""Command-line interface: python -m genome_compare REFERENCE.fasta SAMPLE.fasta."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from .core import Comparison, compare_sequences, read_fasta

TSV_HEADER = "position\treference\talternate\ttype"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genome_compare",
        description=(
            "Compare two nucleotide sequences (FASTA) position by position and "
            "list the single-base differences."
        ),
    )
    parser.add_argument("reference", type=Path, help="FASTA file with the reference sequence")
    parser.add_argument("sample", type=Path, help="FASTA file with the sequence to compare")
    parser.add_argument(
        "--ignore-ambiguous",
        action="store_true",
        help="skip positions where either base is not A, C, G or T (e.g. N)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="also write the differences to this TSV file"
    )
    return parser


def format_summary(comparison: Comparison, ref_name: str, sample_name: str) -> str:
    lines = [
        f"Reference: {ref_name} ({comparison.reference_length:,} bp)",
        f"Sample:    {sample_name} ({comparison.sample_length:,} bp)",
        f"Compared:  first {comparison.compared_length:,} positions",
    ]
    if comparison.reference_length != comparison.sample_length:
        lines.append("Note:      lengths differ; positions beyond the shorter one were not compared")
    if comparison.skipped_ambiguous:
        lines.append(f"Skipped:   {comparison.skipped_ambiguous:,} positions with ambiguous bases")
    lines.append(
        f"Differences: {len(comparison.variants):,} "
        f"({comparison.difference_rate:.4%} of {comparison.evaluated_positions:,} evaluated positions)"
    )
    lines.append(
        f"  transitions: {comparison.count('transition')}, "
        f"transversions: {comparison.count('transversion')}, "
        f"other: {comparison.count('other')}"
    )
    return "\n".join(lines)


def format_table(comparison: Comparison) -> str:
    rows = [TSV_HEADER]
    rows += [f"{v.position}\t{v.reference}\t{v.alternate}\t{v.kind}" for v in comparison.variants]
    return "\n".join(rows)


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        reference = read_fasta(args.reference)
        sample = read_fasta(args.sample)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    comparison = compare_sequences(
        reference.sequence, sample.sequence, ignore_ambiguous=args.ignore_ambiguous
    )

    print(format_summary(comparison, reference.header.split()[0], sample.header.split()[0]))
    if comparison.variants:
        print()
        print(format_table(comparison))

    if args.output:
        args.output.write_text(format_table(comparison) + "\n", encoding="utf-8")
        print(f"\nDifferences written to {args.output}")
    return 0
