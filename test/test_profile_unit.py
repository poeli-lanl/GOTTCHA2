import io
import os
import tempfile
import unittest

from gottcha.utils import profile


class TestProfileUtils(unittest.TestCase):
    """Tests for helper functions that remain in gottcha.utils.profile.

    Argument parsing (parse_args) has moved to gottcha.gottcha2; see
    test/test_gottcha2_cli_unit.py::TestProfileArguments for those tests.
    """

    def test_load_database_stats_with_header_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats_path = os.path.join(tmp, "db.stats")
            with open(stats_path, "w") as f:
                f.write("Rank\tName\tTaxid\tSK\tNum\tMax\tMin\tTotalLength\tGenomeSize\tNote\n")
                f.write("species\tX\t123\tB\t1\t0\t0\t1000\t1200\tok\n")

            df = profile.load_database_stats(stats_path)
            self.assertIn("123", df.index)
            self.assertEqual(int(df.loc["123", "TotalLength"]), 1000)
            self.assertEqual(int(df.loc["123", "GenomeSize"]), 1200)

    def test_load_acc_list_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "acc.txt")
            open(p, "w").close()
            self.assertEqual(profile.load_acc_list(p), set())


if __name__ == "__main__":
    unittest.main()
