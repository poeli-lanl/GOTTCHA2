import io
import unittest

import pandas as pd

from gottcha.utils import report


SUMMARY_FIELDS = [
    "LEVEL", "NAME", "TAXID", "READ_COUNT", "TOTAL_BP_MAPPED", "SNI_SCORE",
    "COVERED_SIG_LEN", "SIG_COV", "DEPTH", "REL_ABUNDANCE_GC", "REL_ABUNDANCE",
]
FULL_FIELDS = SUMMARY_FIELDS + [
    "PARENT_NAME", "PARENT_TAXID", "AOI_READ_COUNT", "TOTAL_READ_LEN",
    "TOTAL_BP_MISMATCH", "TOTAL_BP_INDEL", "ALN_IDENTITY",
    "CONSENSUS_SEQ_IDENTITY", "SNI_CI95_LH", "SIG_COV_RAW", "MAPPED_SIG_LEN",
    "TOTAL_SIG_LEN", "COVERED_SIG_DEPTH", "COVERED_MAPPED_SIG_COV", "ZSCORE",
    "GENOMIC_CONTENT_EST", "ABUNDANCE", "REL_ABUNDANCE_DEPTH", "SIG_LEVEL",
    "GENOME_COUNT", "GENOME_SIZE", "NOTE",
]


class TestReport(unittest.TestCase):
    def test_summary_filters_notes_and_full_report_preserves_all_rows(self):
        for fmt, separator in (("tsv", "\t"), ("csv", ",")):
            with self.subTest(fmt=fmt):
                rows = []
                for taxid, level, note in (
                    ("1", "genus", ""),
                    ("11", "species", "Database note, with comma"),
                    ("12", "species", "Filtered out (minReads threshold 10 > 2); "),
                    ("111", "strain", "Not shown (strain-result could be biased); "),
                ):
                    row = dict.fromkeys(FULL_FIELDS, 0)
                    row.update(
                        TAXID=taxid, LEVEL=level, NAME=f"Taxon {taxid}", NOTE=note,
                        SIG_LEVEL=7, SIG_COV=0.5, SIG_COV_RAW=0.25,
                        ALN_IDENTITY=0.98, CONSENSUS_SEQ_IDENTITY=0.99,
                        SNI_CI95_LH="[0.95-0.99]",
                    )
                    rows.append(row)
                summary_out, full_out = io.StringIO(), io.StringIO()

                self.assertTrue(report.generate_report_file(
                    pd.DataFrame(rows), summary_out, full_out, fmt
                ))

                summary = pd.read_csv(io.StringIO(summary_out.getvalue()), sep=separator, dtype={"TAXID": str})
                full = pd.read_csv(io.StringIO(full_out.getvalue()), sep=separator, dtype={"TAXID": str})
                self.assertEqual(summary.columns.tolist(), SUMMARY_FIELDS)
                self.assertEqual(full.columns.tolist(), FULL_FIELDS)
                self.assertEqual(summary.TAXID.tolist(), ["1", "11"])
                self.assertEqual(full.TAXID.tolist(), ["1", "11", "12", "111"])
                self.assertEqual(set(full.SIG_LEVEL), {"species"})
                self.assertEqual(full.loc[1, "NOTE"], "Database note, with comma")
                self.assertEqual(full.loc[0, "SIG_COV"], 0.5)
                self.assertEqual(full.loc[0, "SIG_COV_RAW"], 0.25)
                self.assertEqual(full.loc[0, "ALN_IDENTITY"], 0.98)
                self.assertEqual(full.loc[0, "CONSENSUS_SEQ_IDENTITY"], 0.99)


if __name__ == "__main__":
    unittest.main()
