import io
import argparse
import os
import tempfile
from contextlib import redirect_stderr
from pathlib import Path
import sys
import unittest
from unittest import mock


from gottcha import gottcha2


class TestGottcha2Cli(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.database = Path(self.tmp.name) / "gottcha_db.species"
        self.reads = Path(self.tmp.name) / "reads.fastq"
        self.bam = Path(self.tmp.name) / "sample.gottcha_species.bam"
        for suffix in (".mmi", ".tax.tsv", ".stats", ".syldb", ".zip"):
            Path(str(self.database) + suffix).touch()
        self.reads.touch()
        self.bam.touch()

    def test_cli_dispatches_coverage_browser(self):
        with mock.patch.object(gottcha2.coverage_browser, "main") as browser_main:
            with mock.patch.object(sys, "argv", ["gottcha2", "coverage-browser", "-r", "results", "-o", "sample.html"]):
                gottcha2.cli()
        browser_main.assert_called_once_with(["-r", "results", "-o", "sample.html"])

    def test_cli_dispatches_profile(self):
        with mock.patch.object(gottcha2.profile, "main") as profile_main:
            with mock.patch.object(sys, "argv", ["gottcha2", "profile", "-i", str(self.reads), "-d", str(self.database)]):
                gottcha2.cli()
        profile_main.assert_called_once()
        args = profile_main.call_args.args[0]
        self.assertIsInstance(args, argparse.Namespace)
        self.assertEqual(args.input, [str(self.reads)])
        self.assertEqual(args.database, str(self.database))
        self.assertEqual(args.dbLevel, "species")
        self.assertEqual(args.matchIdentity, 0.95)
        self.assertFalse(args.fast)
        self.assertFalse(args.extractOnly)

    def test_cli_dispatches_fast_profile(self):
        Path(str(self.database) + ".mmi").unlink()
        with mock.patch.object(gottcha2.profile, "main") as profile_main:
            gottcha2.cli(["fast-profile", "-i", str(self.reads), "-d", str(self.database), "-np"])
        profile_main.assert_called_once()
        args = profile_main.call_args.args[0]
        self.assertTrue(args.fast)
        self.assertFalse(args.extractOnly)
        self.assertEqual(args.m2_options, "-n2 -m25 -s120 --no-long-join")

    def test_cli_dispatches_extract_without_profiling(self):
        with mock.patch.object(gottcha2.profile, "main") as profile_main:
            gottcha2.cli(["extract", "-b", str(self.bam), "-e", "562"])
        profile_main.assert_called_once()
        args = profile_main.call_args.args[0]
        self.assertTrue(args.extractOnly)
        self.assertEqual(args.extract, "562")
        self.assertEqual(args.bam, str(self.bam))
        self.assertIsNone(args.database)
        self.assertIsNone(args.matchIdentity)

    def test_invalid_input_is_rejected_before_workflow(self):
        with mock.patch.object(gottcha2.profile, "main") as profile_main:
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                gottcha2.cli(["profile", "-i", str(self.reads)])
        self.assertEqual(error.exception.code, 2)
        profile_main.assert_not_called()

    def test_profile_command_help_and_version_do_not_run_workflow(self):
        for command in ("profile", "fast-profile", "extract"):
            for option in ("--help", "--version"):
                with self.subTest(command=command, option=option):
                    output = io.StringIO()
                    with mock.patch.object(gottcha2.profile, "main") as profile_main:
                        with mock.patch.object(sys, "stdout", output), self.assertRaises(SystemExit) as error:
                            gottcha2.cli([command, option])
                    self.assertEqual(error.exception.code, 0)
                    self.assertIn(gottcha2.__version__, output.getvalue())
                    profile_main.assert_not_called()

    def test_cli_prints_version(self):
        buf = io.StringIO()
        with mock.patch.object(sys, "stdout", buf):
            with mock.patch.object(sys, "argv", ["gottcha2", "version"]):
                gottcha2.cli()
        self.assertIn(gottcha2.__version__, buf.getvalue())

    def test_cli_invalid_command_exits(self):
        buf = io.StringIO()
        with mock.patch.object(sys, "argv", ["gottcha2", "badcmd"]), mock.patch.object(sys, "stdout", buf):
            with self.assertRaises(SystemExit) as error:
                gottcha2.cli()
        self.assertEqual(error.exception.code, 1)
        self.assertIn("'badcmd' is not a valid command", buf.getvalue())


class TestProfileArguments(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.database = Path(self.tmp.name) / "gottcha_db.species.fna"
        self.reads = Path(self.tmp.name) / "reads.fastq"
        for suffix in (".mmi", ".tax.tsv", ".stats", ".syldb", ".zip"):
            Path(str(self.database) + suffix).touch()
        self.reads.touch()

    def parse(self, *options):
        return gottcha2.parse_args(
            "test", ["profile", "-i", str(self.reads), "-d", str(self.database), *options]
        )

    def test_parse_args_defaults_for_short_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = gottcha2.parse_args("test", ["profile", "-i", read_path, "-d", db_prefix])

            self.assertEqual(args.dbLevel, "species")
            self.assertEqual(args.matchIdentity, 0.95)
            self.assertEqual(args.errorRate, 0.005)
            self.assertEqual(args.prefix, "reads")
            self.assertEqual(args.input[0], os.path.abspath(read_path))
            self.assertEqual(args.secondary, "no")
            self.assertEqual(args.max_secondary, 10)
            self.assertEqual(args.secondary_ratio, 0.9)
            self.assertEqual(args.reciprocal_groups, "no")
            self.assertEqual(args.matchFraction, 0.95)
            self.assertEqual(args.matchLength, 100)
            self.assertEqual(args.presetx, "sr")
            self.assertEqual(args.m2_options, "-s120")
            self.assertEqual(args.sniScore, "0.9,0.95,0.99")

    def test_parse_args_secondary_choices(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            for secondary in ("yes", "no"):
                with self.subTest(secondary=secondary):
                    args = gottcha2.parse_args(
                        "test",
                        ["profile", "-i", read_path, "-d", db_prefix,
                         "--secondary", secondary],
                    )

                    self.assertEqual(args.secondary, secondary)

    def test_parse_args_rejects_invalid_secondary_options(self):
        cases = (
            (["--secondary", "true"], "invalid choice"),
            (["--secondary"], "expected one argument"),
            (["--no-secondary"], "unrecognized arguments: --no-secondary"),
        )
        for options, message in cases:
            with self.subTest(options=options):
                stderr = io.StringIO()
                with redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
                    gottcha2.parse_args(
                        "test", ["profile", "-i", "reads.fastq", "-d", "db.species"] + options,
                    )

                self.assertEqual(error.exception.code, 2)
                self.assertIn(message, stderr.getvalue())

    def test_parse_args_nanopore_defaults_to_direct_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "ont.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = gottcha2.parse_args("test", ["profile", "-i", read_path, "-d", db_prefix, "-np"])

            self.assertTrue(args.nanopore)
            self.assertFalse(args.ont_chunk)
            self.assertEqual(args.matchIdentity, 0.85)
            self.assertEqual(args.matchFraction, 0)
            self.assertEqual(args.matchLength, 100)
            self.assertEqual(args.errorRate, 0.01)
            self.assertEqual(args.presetx, "lr:hq")
            self.assertEqual(args.m2_options, "-n1 -m25 -s120 --no-long-join")
            self.assertEqual(args.secondary, "no")
            self.assertEqual(args.max_secondary, 10)
            self.assertEqual(args.secondary_ratio, 0.9)

    def test_parse_args_nanopore_chunk_mode_uses_chunk_defaults(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "ont.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = gottcha2.parse_args(
                "test",
                ["profile", "-i", read_path, "-d", db_prefix, "-np", "--ont-chunk"],
            )

            self.assertTrue(args.nanopore)
            self.assertTrue(args.ont_chunk)
            self.assertEqual(args.matchIdentity, 0.85)
            self.assertEqual(args.matchFraction, 0.85)
            self.assertEqual(args.matchLength, 100)
            self.assertEqual(args.errorRate, 0.03)
            self.assertEqual(args.presetx, "sr")
            self.assertEqual(args.m2_options, "-s120")

    def test_fast_nanopore_uses_reduced_reference_defaults_without_mmi(self):
        Path(str(self.database) + ".mmi").unlink()
        args = self.parse("--fast", "-np")

        self.assertTrue(args.fast)
        self.assertEqual(args.presetx, "lr:hq")
        self.assertEqual(args.m2_options, "-n2 -m25 -s120 --no-long-join")
        self.assertEqual(args.matchFraction, 0)
        self.assertEqual(args.errorRate, 0.01)

    def test_mapping_overrides_apply_to_short_reads_and_both_ont_modes(self):
        for mode in ([], ["-np"], ["-np", "--ont-chunk"]):
            with self.subTest(mode=mode):
                args = self.parse(
                    *mode, "--secondary", "yes", "--max-secondary", "25",
                    "--secondary-ratio", "0.7", "--m2-options=-n3 -s150",
                    "-xm", "map-ont", "-mi", "0.9", "-mf", "0.2",
                    "-mg", "200", "-er", "0.02",
                )
                self.assertEqual(args.secondary, "yes")
                self.assertEqual(args.max_secondary, 25)
                self.assertEqual(args.secondary_ratio, 0.7)
                self.assertEqual(args.m2_options, "-n3 -s150")
                self.assertEqual(args.presetx, "map-ont")
                self.assertEqual(args.matchIdentity, 0.9)
                self.assertEqual(args.matchFraction, 0.2)
                self.assertEqual(args.matchLength, 200)
                self.assertEqual(args.errorRate, 0.02)

    def test_secondary_threshold_boundaries(self):
        for ratio in ("0", "1"):
            with self.subTest(ratio=ratio):
                args = self.parse("--max-secondary", "0", "--secondary-ratio", ratio)
                self.assertEqual(args.max_secondary, 0)
                self.assertEqual(args.secondary_ratio, float(ratio))

    def test_rejects_invalid_mapping_thresholds(self):
        cases = (
            (["--max-secondary", "-1"], "--max-secondary must be >= 0"),
            (["--secondary-ratio", "-0.1"], "--secondary-ratio must be between 0 and 1"),
            (["--secondary-ratio", "1.1"], "--secondary-ratio must be between 0 and 1"),
            (["-mi", "1.1"], "--matchIdentity must be between 0 and 1"),
            (["-mf", "-0.1"], "--matchFraction must be between 0 and 1"),
            (["-mg", "-1"], "--matchLength must be a non-negative integer"),
        )
        for options, message in cases:
            with self.subTest(options=options):
                stderr = io.StringIO()
                with redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
                    self.parse(*options)
                self.assertEqual(error.exception.code, 2)
                self.assertIn(message, stderr.getvalue())

    def test_reciprocal_groups_are_opt_in(self):
        for choice in ("yes", "no"):
            with self.subTest(choice=choice):
                args = self.parse("--reciprocal-groups", choice)
                self.assertEqual(args.reciprocal_groups, choice)
                self.assertEqual(args.secondary, "no")

        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            self.parse("--reciprocal-groups", "true")
        self.assertEqual(error.exception.code, 2)

    def test_sni_threshold_shorthand(self):
        for supplied, expected in (
            ("0.8", "0.8,0.8,0.8"),
            ("0.8,0.96", "0.8,0.96,0.99"),
            ("0.8,0.96,0.995", "0.8,0.96,0.995"),
        ):
            with self.subTest(supplied=supplied):
                self.assertEqual(self.parse("-ss", supplied).sniScore, expected)

    def test_nocutoff_preserves_explicit_coverage_and_read_filters(self):
        args = self.parse("-nc", "-ss", "0.99", "-Mc", "0.1", "-Mr", "10")
        self.assertEqual(args.sniScore, "0,0,0")
        self.assertEqual(args.minCov, 0.1)
        self.assertEqual(args.minReads, 10)

    def test_parse_args_rejects_ont_chunk_without_nanopore_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            stderr = io.StringIO()
            with redirect_stderr(stderr), self.assertRaises(SystemExit):
                gottcha2.parse_args(
                    "test",
                    ["profile", "-i", read_path, "-d", db_prefix, "--ont-chunk"],
                )

            self.assertIn("--ont-chunk requires --nanopore", stderr.getvalue())

    def test_parse_args_extractfullref_and_nocutoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = gottcha2.parse_args(
                "test",
                ["profile", "-i", read_path, "-d", db_prefix, "-ef", "-nc"],
            )
            self.assertEqual(args.extract, "all:20:fasta")
            self.assertEqual(args.sniScore, "0,0,0")



if __name__ == "__main__":
    unittest.main()
