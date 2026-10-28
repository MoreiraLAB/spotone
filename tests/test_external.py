"""Verify iFeature subprocess boundaries without installing upstream packages."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

from spotone.external import PSEUDO_FEATURES, STANDARD_FEATURES, run_ifeature


FAKE_IFEATURE = r'''
import argparse
import json
from pathlib import Path
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--file', required=True)
parser.add_argument('--type', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--gap_lambda')
parser.add_argument('--raactype')
args = parser.parse_args()
with Path('invocations.jsonl').open('a') as handle:
    handle.write(json.dumps({'file': args.file, 'type': args.type, 'out': args.out,
                             'gap_lambda': args.gap_lambda, 'raactype': args.raactype,
                             'cwd': str(Path.cwd())}) + '\n')
mode = Path('mode').read_text().strip() if Path('mode').exists() else 'success'
if mode == 'failure' or mode == 'fail-' + args.type:
    print('upstream dependency unavailable', file=sys.stderr)
    sys.exit(7)
if mode == 'no-output':
    print('sequence is too short')
    sys.exit(0)
if mode == 'empty-output':
    Path(args.out).touch()
    sys.exit(0)
assert Path(args.file).is_absolute()
assert Path(args.out).is_absolute()
assert Path(args.file).read_text().startswith('>')
Path(args.out).write_text('#\tfeature\nprotein\t1.0\n')
'''


class ExternalTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="spotone external tests ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.checkout = self.root / "iFeature checkout with spaces"
        self.checkout.mkdir()
        for name in ["iFeature.py", "iFeaturePseKRAAC.py"]:
            (self.checkout / name).write_text(FAKE_IFEATURE, encoding="utf-8")
        self.fasta = self.root / "input with spaces.fasta"
        self.fasta.write_text(">protein\nACDEFGHIKLMNPQRSTVWY\n", encoding="utf-8")
        self.output = self.root / "output with spaces"

    def run_features(self, **kwargs):
        return run_ifeature(self.fasta, self.checkout, self.output, python_executable=sys.executable, **kwargs)

    def invocations(self):
        return [json.loads(line) for line in (self.checkout / "invocations.jsonl").read_text().splitlines()]

    def test_standard_and_pseudo_calls_handle_spaces_and_leave_cwd_unchanged(self):
        original_cwd = Path.cwd()
        paths = self.run_features(features=["AAC", "DPC", "AAC"], pseudo_features=["type1"])
        self.assertEqual(Path.cwd(), original_cwd)
        self.assertEqual(paths, [self.output / "AAC.tsv", self.output / "DPC.tsv", self.output / "type1.tsv"])
        self.assertTrue(all(path.is_file() and path.stat().st_size > 0 for path in paths))
        calls = self.invocations()
        self.assertEqual([call["type"] for call in calls], ["AAC", "DPC", "type1"])
        self.assertTrue(all(call["cwd"] == str(self.checkout) for call in calls))
        self.assertEqual((calls[-1]["gap_lambda"], calls[-1]["raactype"]), ("5", "5"))
        self.assertIsNone(calls[0]["gap_lambda"])
        self.assertFalse(list(self.output.glob(".spotone-ifeature-*")))

    def test_defaults_are_all_original_standard_families_and_no_pseudo(self):
        paths = self.run_features()
        self.assertEqual(len(paths), 21)
        self.assertEqual([call["type"] for call in self.invocations()], list(STANDARD_FEATURES))
        self.assertEqual(len(PSEUDO_FEATURES), 17)

    def test_pseudo_only_does_not_require_standard_script(self):
        (self.checkout / "iFeature.py").unlink()
        paths = self.run_features(features=[], pseudo_features=["type3B"])
        self.assertEqual(paths, [self.output / "type3B.tsv"])

    def test_nonzero_exit_has_clear_error_and_preserves_existing_results(self):
        self.output.mkdir()
        old = self.output / "AAC.tsv"
        old.write_text("previous result")
        (self.checkout / "mode").write_text("failure")
        with self.assertRaisesRegex(RuntimeError, "exit code 7.*upstream dependency unavailable"):
            self.run_features(features=["AAC"])
        self.assertEqual(old.read_text(), "previous result")

    def test_exit_zero_without_output_is_failure_even_with_stale_file(self):
        self.output.mkdir()
        old = self.output / "AAC.tsv"
        old.write_text("previous result")
        (self.checkout / "mode").write_text("no-output")
        with self.assertRaisesRegex(RuntimeError, "did not create a nonempty output.*sequence is too short"):
            self.run_features(features=["AAC"])
        self.assertEqual(old.read_text(), "previous result")

    def test_empty_output_is_failure(self):
        (self.checkout / "mode").write_text("empty-output")
        with self.assertRaisesRegex(RuntimeError, "nonempty output"):
            self.run_features(features=["AAC"])
        self.assertFalse((self.output / "AAC.tsv").exists())

    def test_later_family_failure_does_not_publish_earlier_staged_results(self):
        self.output.mkdir()
        old = self.output / "AAC.tsv"
        old.write_text("previous AAC")
        (self.checkout / "mode").write_text("fail-DPC")
        with self.assertRaisesRegex(RuntimeError, "family DPC failed"):
            self.run_features(features=["AAC", "DPC"])
        self.assertEqual(old.read_text(), "previous AAC")
        self.assertFalse((self.output / "DPC.tsv").exists())

    def test_shell_characters_in_input_paths_remain_literal_arguments(self):
        self.fasta = self.root / "input; touch INJECTED; $(touch INJECTED).fasta"
        self.fasta.write_text(">protein\nACD\n")
        self.run_features(features=["AAC"])
        self.assertEqual(self.invocations()[0]["file"], str(self.fasta))
        self.assertFalse((self.checkout / "INJECTED").exists())
        self.assertFalse((self.root / "INJECTED").exists())

    def test_invalid_family_is_rejected_before_launch(self):
        for kwargs in [{"features": ["AAC;touch INJECTED"]}, {"features": []}, {"features": "AAC"}, {"pseudo_features": ["type5"]}]:
            with self.subTest(kwargs=kwargs), self.assertRaises(RuntimeError):
                self.run_features(**kwargs)
        self.assertFalse((self.checkout / "invocations.jsonl").exists())

    def test_missing_script_or_invalid_fasta_reports_actionable_error(self):
        (self.checkout / "iFeature.py").unlink()
        with self.assertRaisesRegex(RuntimeError, "Cannot find iFeature.py"):
            self.run_features(features=["AAC"])
        self.fasta.write_text(">protein\nAXD\n")
        with self.assertRaisesRegex(RuntimeError, "Invalid amino acid 'X'"):
            self.run_features(features=["AAC"])

    def test_missing_python_executable_reports_launch_error(self):
        with self.assertRaisesRegex(RuntimeError, "Could not launch iFeature family AAC"):
            run_ifeature(self.fasta, self.checkout, self.output, features=["AAC"], python_executable=str(self.root / "no-python"))

    def test_output_cannot_overwrite_input_fasta(self):
        self.fasta = self.root / "AAC.tsv"
        self.fasta.write_text(">protein\nACD\n")
        with self.assertRaisesRegex(RuntimeError, "would overwrite the input FASTA"):
            run_ifeature(self.fasta, self.checkout, self.root, features=["AAC"])
        self.assertEqual(self.fasta.read_text(), ">protein\nACD\n")
        self.assertFalse((self.checkout / "invocations.jsonl").exists())


if __name__ == "__main__":
    unittest.main()
