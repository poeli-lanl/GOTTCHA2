import gzip
import io
import json
import shutil
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

import pysam
import pandas as pd

from gottcha.utils import coverage_browser as browser
from gottcha.utils import coverage_browser_workflow as workflow


class TestBrowserProfileData(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.full = Path(tmp.name) / 'sample.full.tsv'
        self.row = {
            'NAME': 'Example genome', 'TAXID': '00111', 'PARENT_NAME': 'Example species',
            'SIG_LEVEL': 'strain', 'SIG_COV': 0.75, 'TOTAL_SIG_LEN': 1000, 'GENOME_SIZE': 5000,
        }

    def parse(self, **metrics):
        pd.DataFrame([{**self.row, **metrics}]).to_csv(self.full, sep='\t', index=False)
        return browser.parse_full_file(self.full, ['00111'])[0]

    def test_preserves_profiling_metrics_and_taxonomy_identifiers(self):
        genome = self.parse(
            LEVEL='strain', PARENT_TAXID='00011', READ_COUNT=1234,
            SNI_SCORE=0.998765, SNI_CI95_LH='[0.995-0.999]', REL_ABUNDANCE=0.01234,
            DEPTH=12.5, SIG_COV_RAW=0.625, COVERED_SIG_LEN=750,
            TOTAL_BP_MISMATCH=12, ALN_IDENTITY=0.975, CONSENSUS_SEQ_IDENTITY=0.99,
            GENOME_COUNT=1, NOTE='Filtered out (minReads threshold); database note',
        )
        profile = genome['profile']
        self.assertEqual(profile['TAXID'], '00111')
        self.assertEqual(profile['PARENT_TAXID'], '00011')
        self.assertEqual(profile['READ_COUNT'], 1234)
        self.assertEqual(profile['SNI_SCORE'], 0.998765)
        self.assertEqual(profile['SNI_CI95_LH'], '[0.995-0.999]')
        self.assertEqual(profile['REL_ABUNDANCE'], 0.01234)
        self.assertEqual(profile['DEPTH'], 12.5)
        self.assertEqual(profile['SIG_COV_RAW'], 0.625)
        self.assertEqual(profile['ALN_IDENTITY'], 0.975)
        self.assertEqual(profile['CONSENSUS_SEQ_IDENTITY'], 0.99)
        self.assertEqual(profile['NOTE'], 'Filtered out (minReads threshold); database note')
        self.assertEqual(genome['sigCov'], profile['SIG_COV'])
        self.assertEqual(genome['totalLength'], profile['TOTAL_SIG_LEN'])

    def test_older_reports_keep_optional_fields_absent(self):
        genome = self.parse()
        self.assertEqual(genome['name'], 'Example genome')
        self.assertEqual(genome['genomeSize'], 5000)
        self.assertNotIn('READ_COUNT', genome['profile'])
        self.assertNotIn('NOTE', genome['profile'])
        self.assertNotIn('SNI_SCORE', genome['profile'])

    def test_missing_and_nonfinite_metrics_are_distinct_from_zero(self):
        genome = self.parse(
            READ_COUNT=0, AOI_READ_COUNT=0, REL_ABUNDANCE=0, SNI_SCORE='',
            DEPTH='inf', ALN_IDENTITY='not measured', TOTAL_BP_INDEL='NaN', NOTE='',
        )
        profile = genome['profile']
        for field in ('READ_COUNT', 'AOI_READ_COUNT', 'REL_ABUNDANCE'):
            self.assertEqual(profile[field], 0)
        for field in ('SNI_SCORE', 'DEPTH', 'ALN_IDENTITY', 'TOTAL_BP_INDEL'):
            self.assertIsNone(profile[field])
        self.assertEqual(profile['NOTE'], '')
        json.dumps(genome, allow_nan=False)

    def test_notes_remain_text_when_embedded_in_html(self):
        note = 'Review </script><script>alert("note")</script> & keep this text'
        genome = self.parse(NOTE=note)
        embedded = browser._compact_json([genome])
        self.assertNotIn('<script>', embedded)
        self.assertNotIn('</script>', embedded)
        self.assertEqual(json.loads(embedded)[0]['profile']['NOTE'], note)

    def test_profile_data_only_contains_selected_genomes(self):
        pd.DataFrame([
            {**self.row, 'READ_COUNT': 12},
            {**self.row, 'TAXID': '222', 'READ_COUNT': 99},
        ]).to_csv(self.full, sep='\t', index=False)
        genomes = browser.parse_full_file(self.full, ['222'])
        self.assertEqual(len(genomes), 1)
        self.assertEqual(genomes[0]['profile']['READ_COUNT'], 99)


class TestBrowserDiscovery(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="coverage browser ")
        self.addCleanup(tmp.cleanup)
        self.directory = Path(tmp.name)
        self.bam = self.touch("sample.gottcha_species.bam")
        self.full = self.touch("sample.full.tsv")

    def touch(self, name):
        path = self.directory / name
        path.touch()
        return path

    def test_discovers_matching_sample_and_reference(self):
        reference = self.touch("sample.sylph_extracted.fa.gz")
        inputs = workflow.resolve_inputs(results=self.directory)
        self.assertEqual(inputs.bam, self.bam.resolve())
        self.assertEqual(inputs.full, self.full.resolve())
        self.assertEqual(inputs.reference, reference.resolve())
        self.assertEqual(inputs.output, self.directory.resolve() / "sample.coverage.html")

    def test_multiple_samples_require_prefix_and_never_mix_files(self):
        self.touch("other.gottcha_species.bam")
        self.touch("other.full.tsv")
        self.touch("other.sylph_extracted.fa.gz")
        with self.assertRaisesRegex(ValueError, "found 2"):
            workflow.resolve_inputs(results=self.directory)
        inputs = workflow.resolve_inputs(results=self.directory, prefix="sample")
        self.assertEqual(inputs.bam, self.bam.resolve())
        self.assertIsNone(inputs.reference)

    def test_explicit_bam_and_reference_override_directory_discovery(self):
        bam = self.touch("external.bam")
        reference = self.touch("reference.fa")
        inputs = workflow.resolve_inputs(
            results=self.directory, prefix="sample", bam=bam, reference=reference,
            output=self.directory / "chosen.html",
        )
        self.assertEqual(inputs.bam, bam.resolve())
        self.assertEqual(inputs.reference, reference.resolve())
        self.assertEqual(inputs.full, self.full.resolve())
        self.assertEqual(inputs.output.name, "chosen.html")

    def test_multiple_database_levels_require_explicit_bam(self):
        self.touch("sample.gottcha_strain.bam")
        with self.assertRaisesRegex(ValueError, "Specify --bam"):
            workflow.resolve_inputs(results=self.directory, prefix="sample")

    def test_missing_reference_allows_coverage_only(self):
        inputs = workflow.resolve_inputs(results=self.directory)
        self.assertIsNone(inputs.reference)
        self.assertEqual(inputs.vcfs, [])

    def test_existing_vcf_is_used_without_recalling_variants(self):
        vcf = self.touch("sample.gottcha_species.vcf.gz")
        inputs = workflow.resolve_inputs(results=self.directory)
        self.assertEqual(inputs.vcfs, [vcf.resolve()])
        skipped = workflow.resolve_inputs(results=self.directory, no_variants=True)
        self.assertEqual(skipped.vcfs, [])

    def test_missing_input_and_output_collisions_are_rejected(self):
        cases = (
            ({}, "Provide --results"),
            ({"bam": self.bam}, "--full is required"),
            ({"results": self.directory, "prefix": "absent"}, "found 0"),
            ({"results": self.directory, "reference": self.directory / "absent.fa"}, "Reference file not found"),
            ({"results": self.directory, "output": self.full}, "must not overwrite"),
            ({"results": self.directory, "output": self.directory}, "not a directory"),
        )
        for options, message in cases:
            with self.subTest(options=options), self.assertRaisesRegex(ValueError, message):
                workflow.resolve_inputs(**options)


class TestBrowserWorkflow(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="coverage browser ")
        self.addCleanup(tmp.cleanup)
        self.directory = Path(tmp.name)
        self.bam = self.directory / "sample.gottcha_species.bam"
        self.reference = self.directory / "sample.sylph_extracted.fa.gz"
        self.full = self.directory / "sample.full.tsv"
        self.output = self.directory / "sample.coverage.html"
        contig = "chrA|101|220|111|"
        sequence = "ACGT" * 30
        with gzip.open(self.reference, "wt") as handle:
            handle.write(f">{contig}\n{sequence}\n")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": contig, "LN": len(sequence)}],
                  "RG": [{"ID": "sample", "SM": "sample"}]}
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as bam:
            for index in range(12):
                alignment = pysam.AlignedSegment(bam.header)
                alignment.query_name = f"read{index}"
                alignment.query_sequence = sequence[:50] + "T" + sequence[51:]
                alignment.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
                alignment.reference_id = 0
                alignment.reference_start = 0
                alignment.mapping_quality = 60
                alignment.cigar = [(0, len(sequence))]
                alignment.set_tags([("RG", "sample"), ("NM", 1), ("MD", "50G69")])
                bam.write(alignment)
        self.full.write_text(
            "LEVEL\tTAXID\tPARENT_NAME\tSIG_COV\tSIG_LEVEL\tNAME\tTOTAL_SIG_LEN\tGENOME_SIZE\n"
            "strain\t111\tExample species\t1.0\tstrain\tExample genome\t120\t5000\n"
        )

    def run_main(self, *args):
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            browser.main(list(map(str, args)))

    def test_real_bam_coverage_and_haploid_variant_calling(self):
        original_reference = self.reference.read_bytes()
        inputs = workflow.resolve_inputs(results=self.directory)
        with tempfile.TemporaryDirectory(dir=self.directory) as tmp:
            coverage, vcfs = workflow.prepare_inputs(inputs, tmp, threads=2)
            rows, taxids = browser.parse_coverage_file(coverage)
            variants, sources = browser.parse_vcf_files([str(vcf) for vcf in vcfs], rows)
            self.assertEqual(taxids, ["111"])
            self.assertEqual(rows[0]["numreads"], 12)
            self.assertEqual(rows[0]["coverage"], 100)
            self.assertEqual(len(variants), 1)
            self.assertEqual(variants[0]["pos"], 151)
            self.assertEqual(variants[0]["ref"], "G")
            self.assertEqual(variants[0]["alt"], "T")
            self.assertEqual(sources[0]["name"], "sample.gottcha_species.vcf.gz")
            self.assertTrue(Path(str(vcfs[0]) + ".tbi").is_file())
            filtered, _ = browser.parse_vcf_files([str(vcf) for vcf in vcfs], rows, min_depth=100)
            self.assertEqual(filtered, [])
        self.assertEqual(self.reference.read_bytes(), original_reference)
        self.assertFalse(Path(str(self.reference) + ".fai").exists())

    def test_directory_mode_generates_self_contained_html_and_cleans_intermediates(self):
        self.run_main("-r", self.directory, "-t", "2")
        html = self.output.read_text()
        self.assertTrue("Example genome" in html)
        self.assertTrue("sample.gottcha_species.vcf.gz" in html)
        self.assertTrue("browserResourcesReady" in html)
        self.assertTrue("/publicdata/" not in html)
        self.assertTrue("_PLACEHOLDER" not in html)
        self.assertEqual(list(self.directory.glob(".gottcha2-coverage-*")), [])

    def test_explicit_inputs_accept_plain_fasta_and_custom_output(self):
        reference = self.directory / "reference.fa"
        with gzip.open(self.reference, "rb") as source:
            reference.write_bytes(source.read())
        output = self.directory / "reports" / "custom.html"
        self.run_main("-b", self.bam, "--reference", reference, "-f", self.full, "-o", output)
        self.assertTrue("sample.gottcha_species.vcf.gz" in output.read_text())

    def test_coverage_only_skips_variant_calling_even_with_reference(self):
        with mock.patch.object(workflow.bcftools, "mpileup", side_effect=AssertionError("must not call")):
            self.run_main("-r", self.directory, "--no-variants")
        self.assertTrue("Example genome" in self.output.read_text())
        self.assertTrue("sample.gottcha_species.vcf.gz" not in self.output.read_text())

    def test_directory_without_reference_generates_coverage_only(self):
        self.reference.unlink()
        self.run_main("-r", self.directory)
        self.assertTrue("Example genome" in self.output.read_text())
        self.assertTrue("sample.gottcha_species.vcf.gz" not in self.output.read_text())

    def test_directory_reuses_existing_indexed_vcf(self):
        inputs = workflow.resolve_inputs(results=self.directory)
        with tempfile.TemporaryDirectory() as tmp:
            _, vcfs = workflow.prepare_inputs(inputs, tmp)
            for suffix in ("", ".tbi"):
                shutil.copyfile(str(vcfs[0]) + suffix, self.directory / (vcfs[0].name + suffix))
        with mock.patch.object(workflow.bcftools, "mpileup", side_effect=AssertionError("must not recall")):
            self.run_main("-r", self.directory)
        self.assertTrue("sample.gottcha_species.vcf.gz" in self.output.read_text())

    def test_pileup_failure_stops_variant_calling(self):
        inputs = workflow.resolve_inputs(results=self.directory)
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.object(workflow.bcftools, "mpileup", side_effect=pysam.SamtoolsError("bad pileup")):
                with mock.patch.object(workflow.bcftools, "call") as call:
                    with self.assertRaises(pysam.SamtoolsError):
                        workflow.prepare_inputs(inputs, tmp)
        call.assert_not_called()

    def test_existing_vcf_requires_an_index(self):
        (self.directory / "sample.gottcha_species.vcf.gz").touch()
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
            browser.main(["-r", str(self.directory)])
        self.assertEqual(error.exception.code, 1)
        self.assertIn("Missing tabix index", stderr.getvalue())
        self.assertFalse(self.output.exists())

    def test_legacy_coverage_and_full_inputs_still_work(self):
        coverage = self.directory / "coverage.tsv"
        pysam.coverage("-o", str(coverage), str(self.bam), catch_stdout=False)
        self.run_main("-c", coverage, "-f", self.full, "-o", self.output)
        self.assertTrue("Example genome" in self.output.read_text())

    def test_preparation_failure_preserves_existing_html_and_cleans_temporary_files(self):
        self.output.write_text("previous report")
        with mock.patch.object(workflow.pysam, "coverage", side_effect=pysam.SamtoolsError("bad BAM")):
            with self.assertRaises(SystemExit) as error:
                self.run_main("-r", self.directory)
        self.assertEqual(error.exception.code, 1)
        self.assertEqual(self.output.read_text(), "previous report")
        self.assertEqual(list(self.directory.glob(".gottcha2-coverage-*")), [])

    def test_reference_contigs_must_match_bam(self):
        with gzip.open(self.reference, "wt") as handle:
            handle.write(">other\nACGT\n")
        inputs = workflow.resolve_inputs(results=self.directory)
        with tempfile.TemporaryDirectory(dir=self.directory) as tmp:
            with self.assertRaisesRegex(ValueError, "missing or mismatched"):
                workflow.prepare_inputs(inputs, tmp)

    def test_invalid_cli_options_fail_before_preparing_data(self):
        for args in (("--threads", "0"), ("--min-depth", "-1"), ("--no-variants", "--vcf", "x.vcf.gz")):
            with self.subTest(args=args), self.assertRaises(SystemExit) as error:
                self.run_main("-r", self.directory, *args)
            self.assertEqual(error.exception.code, 2)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
