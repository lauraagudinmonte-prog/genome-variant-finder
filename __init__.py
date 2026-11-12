"""genome_compare: position-by-position comparison of two nucleotide sequences."""

from .core import Comparison, FastaRecord, Variant, compare_sequences, read_fasta

__all__ = ["Comparison", "FastaRecord", "Variant", "compare_sequences", "read_fasta"]
__version__ = "0.1.0"
