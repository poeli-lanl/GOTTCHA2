import os
import sys
import tempfile
import types
import unittest

import pandas as pd


if "pysam" not in sys.modules:
    sys.modules["pysam"] = types.ModuleType("pysam")

from gottcha.utils import extract_reads


class TestExtractReadsUtils(unittest.TestCase):
    def test_parse_taxids_filters_note_and_selection(self):
        res_df = pd.DataFrame(
            [
                {"LEVEL": "species", "NAME": "Alpha taxa", "TAXID": "11", "NOTE": ""},
                {"LEVEL": "species", "NAME": "Beta taxa", "TAXID": "22", "NOTE": "Filtered out"},
            ]
        )
        taxa_dict, qualified_taxids = extract_reads.parse_taxids("11", res_df, "unused.tsv")

        self.assertEqual(qualified_taxids, ["11"])
        self.assertIn("11", taxa_dict)
        self.assertEqual(taxa_dict["11"]["name"], "Alpha_taxa")

    def test_parse_taxids_from_file(self):
        res_df = pd.DataFrame(
            [{"LEVEL": "species", "NAME": "Taxon Name", "TAXID": "22", "NOTE": ""}]
        )
        with tempfile.TemporaryDirectory() as tmp:
            taxid_file = os.path.join(tmp, "taxids.txt")
            with open(taxid_file, "w") as f:
                f.write("# comment\n22\n")

            taxa_dict, qualified_taxids = extract_reads.parse_taxids(f"@{taxid_file}", res_df, "unused.tsv")

        self.assertEqual(qualified_taxids, ["22"])
        self.assertIn("22", taxa_dict)

    def test_iter_tasks_aoi_filters_and_lineage(self):
        extract_reads.lineage_cache = {"111": ["11"], "222": []}
        refs = [
            "ACC1|1|100|111|A",
            "ACC2|1|100|111|A",
            "ACC3|1|100|222|A",
            "BADREF",
        ]

        tasks_out = list(extract_reads._iter_tasks(refs, 5, {"ACC2"}, "filter_out"))
        tasks_in = list(extract_reads._iter_tasks(refs, 5, {"ACC2"}, "filter_in"))

        self.assertEqual(tasks_out, [("ACC1|1|100|111|A", 5, False)])
        self.assertEqual(tasks_in, [("ACC2|1|100|111|A", 5, True)])


if __name__ == "__main__":
    unittest.main()
