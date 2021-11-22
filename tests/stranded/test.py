import os.path
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class StrandedBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir, 'stranded',
                                      self._sub_test_name())
        self._generated_input_dir = os.path.join(self._test_dir,
                                                 'generated_input')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._tmp_dir = os.path.join(self._test_dir, 'tmp')
        tests.util.recreate_dirs([
            self._generated_input_dir, self._out_dir, self._tmp_dir,
            self._command_output_dir()
        ])

        self._read_length = 50
        self._chromosome_length = 2000
        self._task = 'both'

        self._sample_1_bams_path = os.path.join(self._generated_input_dir,
                                                'b1.txt')
        sample_1_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_1_rep_{}.bam')
        self._sample_1_bams = self._create_sample_1_bams(
            self._sample_1_bams_path, sample_1_bam_replicate_template)

        self._gtf_path = os.path.join(self._generated_input_dir, 'test.gtf')
        transcript_exons, transcript_genes, transcript_strands = (
            self._exons_genes_and_strands_by_transcript())
        self._gtf = self._create_gtf_from_transcripts(
            self._gtf_path,
            transcript_exons,
            genes=transcript_genes,
            strands=transcript_strands)

    def _command_output_dir(self):
        return os.path.join(self._test_dir, 'command_output')

    def _rmats_arguments(self):
        return [
            '--b1',
            self._sample_1_bams_path,
            '--gtf',
            self._gtf_path,
            '--od',
            self._out_dir,
            '-t',
            self._read_type,
            '--readLength',
            str(self._read_length),
            '--tmp',
            self._tmp_dir,
            '--task',
            self._task,
            '--statoff',
        ]

    def _exons_genes_and_strands_by_transcript(self):
        exons_genes_and_strands = [
            # SE gene 1 + strand
            ([(1, 100), (201, 300), (401, 500)], [1], ['+']),
            ([(1, 100), (401, 500)], [1], ['+']),
            # SE gene 2 - strand
            ([(1001, 1100), (1201, 1300), (1401, 1500)], [2], ['-']),
            ([(1001, 1100), (1401, 1500)], [2], ['-']),
        ]
        exons_by_transcript = list()
        genes_by_transcript = list()
        strands_by_transcript = list()
        for exons, genes, strands in exons_genes_and_strands:
            for gene in genes:
                for strand in strands:
                    exons_by_transcript.append(exons)
                    genes_by_transcript.append(gene)
                    strands_by_transcript.append(strand)

        return exons_by_transcript, genes_by_transcript, strands_by_transcript

    def _check_results_first_strand(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 2)
        for row in se_mats_jc_rows:
            self.assertIn(row['exonStart_0base'], ['200', '1200'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,1')
                self.assertEqual(row['strand'], '+')
            elif row['exonStart_0base'] == '1200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,0')
                self.assertEqual(row['strand'], '-')

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(se_mats_jcec_path))
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 2)
        for row in se_mats_jcec_rows:
            self.assertIn(row['exonStart_0base'], ['200', '1200'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,1')
                self.assertEqual(row['strand'], '+')
            elif row['exonStart_0base'] == '1200':
                self.assertEqual(row['IJC_SAMPLE_1'], '2,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,0')
                self.assertEqual(row['strand'], '-')

    def _check_results_second_strand(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 2)
        for row in se_mats_jc_rows:
            self.assertIn(row['exonStart_0base'], ['200', '1200'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,0')
                self.assertEqual(row['strand'], '+')
            elif row['exonStart_0base'] == '1200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,1')
                self.assertEqual(row['strand'], '-')

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(se_mats_jcec_path))
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 2)
        for row in se_mats_jcec_rows:
            self.assertIn(row['exonStart_0base'], ['200', '1200'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '2,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,0')
                self.assertEqual(row['strand'], '+')
            elif row['exonStart_0base'] == '1200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,1')
                self.assertEqual(row['strand'], '-')

    def _check_results_unstranded(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 2)
        for row in se_mats_jc_rows:
            self.assertIn(row['exonStart_0base'], ['200', '1200'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['strand'], '+')
            elif row['exonStart_0base'] == '1200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['strand'], '-')

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(se_mats_jcec_path))
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 2)
        for row in se_mats_jcec_rows:
            self.assertIn(row['exonStart_0base'], ['200', '1200'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['strand'], '+')
            elif row['exonStart_0base'] == '1200':
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['strand'], '-')


class PairedStrandedBaseTest(StrandedBaseTest):
    def setUp(self):
        super().setUp()
        self._read_type = 'paired'

    def _create_sample_1_bams(self, sample_1_bams_path,
                              sample_1_replicate_template):
        rep_1_bam_path = sample_1_replicate_template.format(1)
        rep_1_bam = self._create_bam_from_paired_read_coords(
            rep_1_bam_path,
            self._chromosome_length,
            self._read_length,
            self._paired_read_coords_1_1(),
            is_reversed_1=False,
            is_reversed_2=True)

        rep_2_bam_path = sample_1_replicate_template.format(2)
        rep_2_bam = self._create_bam_from_paired_read_coords(
            rep_2_bam_path,
            self._chromosome_length,
            self._read_length,
            self._paired_read_coords_1_2(),
            is_reversed_1=True,
            is_reversed_2=False)

        sample_1_bams = [rep_1_bam, rep_2_bam]
        self._write_bams(sample_1_bams, sample_1_bams_path)
        return sample_1_bams

    def _include_read_1(self):
        return ([[81, 100], [201, 300]], [[201, 300]])

    def _skip_read_1(self):
        return ([[81, 100], [401, 500]], [[401, 500]])

    def _include_read_2(self):
        return ([[1081, 1100], [1201, 1300]], [[1201, 1300]])

    def _skip_read_2(self):
        return ([[1081, 1100], [1401, 1500]], [[1401, 1500]])

    def _paired_read_coords_1_1(self):
        return [
            self._include_read_1(),
            self._skip_read_1(),
            self._include_read_2(),
            self._skip_read_2(),
        ]

    def _paired_read_coords_1_2(self):
        return [
            self._include_read_1(),
            self._skip_read_1(),
            self._include_read_2(),
            self._skip_read_2(),
        ]


class SingleEndStrandedBaseTest(StrandedBaseTest):
    def setUp(self):
        super().setUp()
        self._read_type = 'single'

    def _create_sample_1_bams(self, sample_1_bams_path,
                              sample_1_replicate_template):
        rep_1_bam_path = sample_1_replicate_template.format(1)
        rep_1_bam = self._create_bam_from_single_end_read_coords(
            rep_1_bam_path,
            self._chromosome_length,
            self._read_length,
            self._single_end_read_coords_1_1(),
            is_reversed=False)

        rep_2_bam_path = sample_1_replicate_template.format(2)
        rep_2_bam = self._create_bam_from_single_end_read_coords(
            rep_2_bam_path,
            self._chromosome_length,
            self._read_length,
            self._single_end_read_coords_1_2(),
            is_reversed=True)

        sample_1_bams = [rep_1_bam, rep_2_bam]
        self._write_bams(sample_1_bams, sample_1_bams_path)
        return sample_1_bams

    def _include_read_1(self):
        return [[81, 100], [201, 300]]

    def _include_exon_read_1(self):
        return [[201, 300]]

    def _skip_read_1(self):
        return [[81, 100], [401, 500]]

    def _skip_exon_read_1(self):
        return [[401, 500]]

    def _include_read_2(self):
        return [[1081, 1100], [1201, 1300]]

    def _include_exon_read_2(self):
        return [[1201, 1300]]

    def _skip_read_2(self):
        return [[1081, 1100], [1401, 1500]]

    def _skip_exon_read_2(self):
        return [[1401, 1500]]

    def _single_end_read_coords_1_1(self):
        return [
            self._include_read_1(),
            self._include_exon_read_1(),
            self._skip_read_1(),
            self._skip_exon_read_1(),
            self._include_read_2(),
            self._include_exon_read_2(),
            self._skip_read_2(),
            self._skip_exon_read_2(),
        ]

    def _single_end_read_coords_1_2(self):
        return [
            self._include_read_1(),
            self._include_exon_read_1(),
            self._skip_read_1(),
            self._skip_exon_read_1(),
            self._include_read_2(),
            self._include_exon_read_2(),
            self._skip_read_2(),
            self._skip_exon_read_2(),
        ]


class PairedFirstStrandTest(PairedStrandedBaseTest):
    def _sub_test_name(self):
        return 'paired_first_strand'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend(['--libType', 'fr-firststrand'])
        return arguments

    def _check_results(self):
        self._check_results_first_strand()


class PairedSecondStrandTest(PairedStrandedBaseTest):
    def _sub_test_name(self):
        return 'paired_second_strand'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend(['--libType', 'fr-secondstrand'])
        return arguments

    def _check_results(self):
        self._check_results_second_strand()


class PairedUnstrandedTest(PairedStrandedBaseTest):
    def _sub_test_name(self):
        return 'paired_unstranded'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend(['--libType', 'fr-unstranded'])
        return arguments

    def _check_results(self):
        self._check_results_unstranded()


class SingleEndFirstStrandTest(SingleEndStrandedBaseTest):
    def _sub_test_name(self):
        return 'single_end_first_strand'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend(['--libType', 'fr-firststrand'])
        return arguments

    def _check_results(self):
        self._check_results_first_strand()


class SingleEndSecondStrandTest(SingleEndStrandedBaseTest):
    def _sub_test_name(self):
        return 'single_end_second_strand'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend(['--libType', 'fr-secondstrand'])
        return arguments

    def _check_results(self):
        self._check_results_second_strand()


class SingleEndUnstrandedTest(SingleEndStrandedBaseTest):
    def _sub_test_name(self):
        return 'single_end_unstranded'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend(['--libType', 'fr-unstranded'])
        return arguments

    def _check_results(self):
        self._check_results_unstranded()


if __name__ == '__main__':
    unittest.main(verbosity=2)
