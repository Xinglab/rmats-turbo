import math
import os.path
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class DartsModelBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir, 'darts_model',
                                      self._sub_test_name())
        self._generated_input_dir = os.path.join(self._test_dir,
                                                 'generated_input')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._tmp_dir = os.path.join(self._test_dir, 'tmp')
        tests.util.recreate_dirs([
            self._generated_input_dir, self._out_dir, self._tmp_dir,
            self._command_output_dir()
        ])

        self._read_type = 'paired'
        self._read_length = 50
        self._chromosome_length = 2000
        self._task = 'both'

        self._sample_1_bams_path = os.path.join(self._generated_input_dir,
                                                'b1.txt')
        sample_1_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_1_rep_{}.bam')
        self._sample_1_bams = self._create_sample_1_bams(
            self._sample_1_bams_path, sample_1_bam_replicate_template)

        self._sample_2_bams_path = os.path.join(self._generated_input_dir,
                                                'b2.txt')
        sample_2_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_2_rep_{}.bam')
        self._sample_2_bams = self._create_sample_2_bams(
            self._sample_2_bams_path, sample_2_bam_replicate_template)

        self._gtf_path = os.path.join(self._generated_input_dir, 'test.gtf')
        self._gtf = self._create_gtf_from_transcripts(
            self._gtf_path, self._exons_by_transcript())

    def _command_output_dir(self):
        return os.path.join(self._test_dir, 'command_output')

    def _rmats_arguments(self):
        return [
            '--b1',
            self._sample_1_bams_path,
            '--b2',
            self._sample_2_bams_path,
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
            '--darts-model',
        ]

    def _create_sample_1_bams(self, sample_1_bams_path,
                              sample_1_replicate_template):
        rep_1_bam_path = sample_1_replicate_template.format(1)
        rep_1_bam = self._create_bam_from_paired_read_coords(
            rep_1_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords_1_1())

        rep_2_coords = self._paired_read_coords_1_2()
        if rep_2_coords:
            rep_2_bam_path = sample_1_replicate_template.format(2)
            rep_2_bam = self._create_bam_from_paired_read_coords(
                rep_2_bam_path, self._chromosome_length, self._read_length,
                rep_2_coords)
            sample_1_bams = [rep_1_bam, rep_2_bam]
        else:
            sample_1_bams = [rep_1_bam]

        self._write_bams(sample_1_bams, sample_1_bams_path)
        return sample_1_bams

    def _create_sample_2_bams(self, sample_2_bams_path,
                              sample_2_replicate_template):
        rep_1_bam_path = sample_2_replicate_template.format(1)
        rep_1_bam = self._create_bam_from_paired_read_coords(
            rep_1_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords_2_1())

        rep_2_coords = self._paired_read_coords_2_2()
        if rep_2_coords:
            rep_2_bam_path = sample_2_replicate_template.format(2)
            rep_2_bam = self._create_bam_from_paired_read_coords(
                rep_2_bam_path, self._chromosome_length, self._read_length,
                rep_2_coords)
            sample_2_bams = [rep_1_bam, rep_2_bam]
        else:
            sample_2_bams = [rep_1_bam]

        self._write_bams(sample_2_bams, sample_2_bams_path)
        return sample_2_bams

    def _read_possibly_na_float(self, float_str):
        if float_str == 'NA':
            return math.nan

        return float(float_str)

    def _exons_by_transcript(self):
        return [
            [(1, 100), (201, 300), (401, 500)],
            [(1, 100), (401, 500)],
        ]

    def _include_read(self):
        return ([[81, 100], [201, 300]], [[201, 300]])

    def _skip_read(self):
        return ([[81, 100], [401, 500]], [[401, 500]])

    def _check_row(self, row, ijc_sample_1, sjc_sample_1, ijc_sample_2,
                   sjc_sample_2, num_replicates):
        self.assertEqual(row['IJC_SAMPLE_1'], ijc_sample_1)
        self.assertEqual(row['SJC_SAMPLE_1'], sjc_sample_1)
        self.assertEqual(row['IJC_SAMPLE_2'], ijc_sample_2)
        self.assertEqual(row['SJC_SAMPLE_2'], sjc_sample_2)
        pvalue = self._read_possibly_na_float(row['PValue'])
        fdr = self._read_possibly_na_float(row['FDR'])
        tests.util.assert_within_bounds(self, pvalue, 0, 1)
        self.assertTrue(math.isnan(fdr))

        inc_level_1s = row['IncLevel1'].split(',')
        self.assertEqual(len(inc_level_1s), num_replicates)
        for inc_level_str in inc_level_1s:
            inc_level = self._read_possibly_na_float(inc_level_str)
            if not math.isnan(inc_level):
                tests.util.assert_within_bounds(self, inc_level, 0, 1)

        inc_level_2s = row['IncLevel2'].split(',')
        self.assertEqual(len(inc_level_2s), num_replicates)
        for inc_level_str in inc_level_2s:
            inc_level = self._read_possibly_na_float(inc_level_str)
            if not math.isnan(inc_level):
                tests.util.assert_within_bounds(self, inc_level, 0, 1)

        inc_level_diff = float(row['IncLevelDifference'])
        tests.util.assert_within_bounds(self, inc_level_diff, 0, 1)


class ReplicateTest(DartsModelBaseTest):
    def _sub_test_name(self):
        return 'replicate_test'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend(['--darts-cutoff', '0.1'])
        return arguments

    def _paired_read_coords_1_1(self):
        return [
            self._include_read(),
            self._skip_read(),
        ]

    def _paired_read_coords_1_2(self):
        return [
            self._include_read(),
            self._include_read(),
            self._include_read(),
            self._skip_read(),
        ]

    def _paired_read_coords_2_1(self):
        return [
            self._skip_read(),
            self._skip_read(),
        ]

    def _paired_read_coords_2_2(self):
        return [
            self._include_read(),
            self._skip_read(),
            self._skip_read(),
            self._skip_read(),
        ]

    def _check_results(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 1)
        se_mats_jc_row_0 = se_mats_jc_rows[0]
        self._check_row(se_mats_jc_row_0, '1,3', '1,1', '0,1', '2,3', 2)

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            se_mats_jcec_path)
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 1)
        se_mats_jcec_row_0 = se_mats_jcec_rows[0]
        self._check_row(se_mats_jcec_row_0, '2,6', '1,1', '0,2', '2,3', 2)


class NoReplicatesTest(DartsModelBaseTest):
    def _sub_test_name(self):
        return 'no_replicates_test'

    def test(self):
        self._run_test()

    def _paired_read_coords_1_1(self):
        return [
            self._include_read(),
            self._include_read(),
            self._skip_read(),
        ]

    def _paired_read_coords_1_2(self):
        return list()

    def _paired_read_coords_2_1(self):
        return [
            self._include_read(),
            self._skip_read(),
            self._skip_read(),
            self._skip_read(),
        ]

    def _paired_read_coords_2_2(self):
        return list()

    def _check_results(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 1)
        se_mats_jc_row_0 = se_mats_jc_rows[0]
        self._check_row(se_mats_jc_row_0, '2', '1', '1', '3', 1)

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            se_mats_jcec_path)
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 1)
        se_mats_jcec_row_0 = se_mats_jcec_rows[0]
        self._check_row(se_mats_jcec_row_0, '4', '1', '2', '3', 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
