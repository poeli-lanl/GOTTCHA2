import io
import sys
import unittest
from unittest import mock


from gottcha import gottcha2


class TestGottcha2Cli(unittest.TestCase):
    def test_cli_dispatches_coverage_browser(self):
        with mock.patch.object(gottcha2.coverage_browser, "main") as browser_main:
            with mock.patch.object(sys, "argv", ["gottcha2", "coverage-browser", "-r", "results", "-o", "sample.html"]):
                gottcha2.cli()
        browser_main.assert_called_once_with(["-r", "results", "-o", "sample.html"])

    def test_cli_dispatches_profile(self):
        with mock.patch.object(gottcha2.profile, "main") as profile_main:
            with mock.patch.object(sys, "argv", ["gottcha2", "profile", "-i", "reads.fq"]):
                gottcha2.cli()
        profile_main.assert_called_once_with(["profile", "-i", "reads.fq"])

    def test_cli_dispatches_fast_profile(self):
        with mock.patch.object(gottcha2.profile, "main") as profile_main:
            with mock.patch.object(sys, "argv", ["gottcha2", "fast-profile", "-i", "reads.fq"]):
                gottcha2.cli()
        profile_main.assert_called_once_with(["fast-profile", "-i", "reads.fq", "--fast"])

    def test_cli_dispatches_extract_without_profiling(self):
        with mock.patch.object(gottcha2.profile, "main") as profile_main:
            with mock.patch.object(sys, "argv", ["gottcha2", "extract", "-b", "sample.bam", "-e", "562"]):
                gottcha2.cli()
        profile_main.assert_called_once_with(["extract", "-b", "sample.bam", "-e", "562", "-eo"])

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


if __name__ == "__main__":
    unittest.main()
