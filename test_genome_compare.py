import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from genome_compare import Variant, compare_sequences, read_fasta
from genome_compare.cli import main


class VariantKindTests(unittest.TestCase):
    def test_transitions(self):
        self.assertEqual(Variant(1, "A", "G").kind, "transition")
        self.assertEqual(Variant(1, "C", "T").kind, "transition")
        self.assertEqual(Variant(1, "T", "C").kind, "transition")

    def test_transversions(self):
        self.assertEqual(Variant(1, "A", "C").kind, "transversion")
        self.assertEqual(Variant(1, "G", "T").kind, "transversion")

    def test_ambiguous_is_other(self):
        self.assertEqual(Variant(1, "A", "N").kind, "other")
        self.assertEqual(Variant(1, "-", "G").kind, "other")


class CompareSequencesTests(unittest.TestCase):
    def test_identical_sequences_have_no_variants(self):
        result = compare_sequences("ACGTACGT", "ACGTACGT")
        self.assertEqual(result.variants, [])
        self.assertEqual(result.difference_rate, 0.0)

    def test_positions_are_one_based(self):
        result = compare_sequences("ACGT", "ACCT")
        self.assertEqual(result.variants, [Variant(3, "G", "C")])

    def test_first_position_is_reported_as_one(self):
        result = compare_sequences("ACGT", "TCGT")
        self.assertEqual(result.variants[0].position, 1)

    def test_comparison_is_case_insensitive(self):
        self.assertEqual(compare_sequences("acgt", "ACGT").variants, [])

    def test_different_lengths_compare_only_the_shorter(self):
        result = compare_sequences("ACGTAC", "ACG")
        self.assertEqual(result.compared_length, 3)
        self.assertEqual(result.reference_length, 6)
        self.assertEqual(result.sample_length, 3)
        self.assertEqual(result.variants, [])

    def test_ambiguous_bases_are_reported_by_default(self):
        result = compare_sequences("ACGT", "ACNT")
        self.assertEqual(result.variants, [Variant(3, "G", "N")])
        self.assertEqual(result.skipped_ambiguous, 0)

    def test_ignore_ambiguous_skips_and_counts_them(self):
        result = compare_sequences("ACGTA", "ACNTC", ignore_ambiguous=True)
        self.assertEqual(result.variants, [Variant(5, "A", "C")])
        self.assertEqual(result.skipped_ambiguous, 1)
        self.assertEqual(result.evaluated_positions, 4)

    def test_difference_rate_and_counts(self):
        result = compare_sequences("ACGT", "GCGA")  # A->G (ts), T->A (tv)
        self.assertEqual(len(result.variants), 2)
        self.assertEqual(result.count("transition"), 1)
        self.assertEqual(result.count("transversion"), 1)
        self.assertAlmostEqual(result.difference_rate, 0.5)

    def test_empty_comparison_does_not_divide_by_zero(self):
        self.assertEqual(compare_sequences("", "ACGT").difference_rate, 0.0)


class ReadFastaTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)

    def write(self, name, text):
        path = self.dir / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_multiline_lowercase_and_whitespace(self):
        path = self.write("a.fasta", ">seq1 demo\nacgt\nAC GT \n\nTT\n")
        record = read_fasta(path)
        self.assertEqual(record.header, "seq1 demo")
        self.assertEqual(record.sequence, "ACGTACGTTT")

    def test_only_first_record_is_read(self):
        path = self.write("a.fasta", ">one\nAAA\n>two\nCCC\n")
        self.assertEqual(read_fasta(path).sequence, "AAA")

    def test_missing_header_raises(self):
        path = self.write("a.fasta", "ACGT\n")
        with self.assertRaises(ValueError):
            read_fasta(path)

    def test_empty_file_raises(self):
        path = self.write("a.fasta", "")
        with self.assertRaises(ValueError):
            read_fasta(path)

    def test_empty_sequence_raises(self):
        path = self.write("a.fasta", ">only_header\n")
        with self.assertRaises(ValueError):
            read_fasta(path)


class CommandLineTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)
        (self.dir / "ref.fasta").write_text(">ref\nACGTACGTAC\n", encoding="utf-8")
        (self.dir / "sample.fasta").write_text(">sample\nACGAACGTCC\n", encoding="utf-8")

    def run_cli(self, *args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = main([str(a) for a in args])
        return code, out.getvalue(), err.getvalue()

    def test_reports_differences(self):
        code, out, _ = self.run_cli(self.dir / "ref.fasta", self.dir / "sample.fasta")
        self.assertEqual(code, 0)
        self.assertIn("Differences: 2", out)
        self.assertIn("4\tT\tA\ttransversion", out)
        self.assertIn("9\tA\tC\ttransversion", out)

    def test_writes_tsv(self):
        tsv = self.dir / "out.tsv"
        code, _, _ = self.run_cli(
            self.dir / "ref.fasta", self.dir / "sample.fasta", "--output", tsv
        )
        self.assertEqual(code, 0)
        lines = tsv.read_text(encoding="utf-8").splitlines()
        self.assertEqual(lines[0], "position\treference\talternate\ttype")
        self.assertEqual(len(lines), 3)

    def test_missing_file_returns_error_code(self):
        code, _, err = self.run_cli(self.dir / "nope.fasta", self.dir / "sample.fasta")
        self.assertEqual(code, 1)
        self.assertIn("error:", err)


if __name__ == "__main__":
    unittest.main()
