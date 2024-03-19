import collections
import os.path
import shutil
import subprocess
import sys
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class Test(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir, 'task_stat')
        self._generated_input_dir = os.path.join(self._test_dir,
                                                 'generated_input')
        self._out_dir_all = os.path.join(self._test_dir, 'out_all')
        self._tmp_dir_all = os.path.join(self._test_dir, 'tmp_all')
        self._out_dir_all_junctions = os.path.join(self._test_dir,
                                                   'out_all_junctions')
        self._tmp_dir_all_junctions = os.path.join(self._test_dir,
                                                   'tmp_all_junctions')
        self._out_dir_select = os.path.join(self._test_dir, 'out_select')
        self._tmp_dir_select = os.path.join(self._test_dir, 'tmp_select')
        self._out_dir_select_junctions = os.path.join(self._test_dir,
                                                      'out_select_junctions')
        self._tmp_dir_select_junctions = os.path.join(self._test_dir,
                                                      'tmp_select_junctions')
        self._out_dir_just_se = os.path.join(self._test_dir, 'out_just_se')
        self._tmp_dir_just_se = os.path.join(self._test_dir, 'tmp_just_se')

        tests.util.recreate_dirs([
            self._generated_input_dir, self._out_dir_all, self._tmp_dir_all,
            self._out_dir_all_junctions, self._tmp_dir_all_junctions,
            self._out_dir_select, self._tmp_dir_select,
            self._out_dir_select_junctions, self._tmp_dir_select_junctions,
            self._out_dir_just_se, self._tmp_dir_just_se,
            self._command_output_dir()
        ])

        self._read_type = 'paired'
        self._read_length = 50
        self._chromosome_length = 4000

        self._sample_1_bams_path = os.path.join(self._generated_input_dir,
                                                'b1.txt')
        self._sample_2_bams_path = os.path.join(self._generated_input_dir,
                                                'b2.txt')
        sample_1_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_1_rep_{}.bam')
        sample_2_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_2_rep_{}.bam')
        self._sample_1_bams = self._create_sample_1_bams(
            self._sample_1_bams_path, sample_1_bam_replicate_template)
        self._sample_2_bams = self._create_sample_2_bams(
            self._sample_2_bams_path, sample_2_bam_replicate_template)

        self._gtf_path = os.path.join(self._generated_input_dir, 'test.gtf')
        self._gtf = self._create_gtf_from_transcripts(
            self._gtf_path, self._exons_by_transcript())
        self._sub_steps = [
            'statoff',
            'statoff_junctions',
            'selected_stat',
            'selected_stat_junctions',
            'just_se',
            'deferred_stat',
        ]
        self._sub_step = None

    def test(self):
        for sub_step in self._sub_steps:
            self._sub_step = sub_step
            self._setup_sub_step()
            self._run_test()

    def _command_output_dir(self):
        return os.path.join(self._test_dir, 'command_output')

    def _rmats_arguments(self):
        if self._sub_step == 'statoff':
            return [
                '--gtf', self._gtf_path, '-t', self._read_type, '--readLength',
                str(self._read_length), '--od', self._out_dir_all, '--tmp',
                self._tmp_dir_all, '--b1', self._sample_1_bams_path, '--b2',
                self._sample_2_bams_path, '--task', 'both', '--statoff'
            ]
        if self._sub_step == 'statoff_junctions':
            return [
                '--gtf', self._gtf_path, '-t', self._read_type, '--readLength',
                str(self._read_length), '--od', self._out_dir_all_junctions,
                '--tmp', self._tmp_dir_all_junctions, '--b1',
                self._sample_1_bams_path, '--b2', self._sample_2_bams_path,
                '--task', 'both', '--statoff', '--individual-counts'
            ]
        if self._sub_step == 'selected_stat':
            return [
                '--od', self._out_dir_select, '--tmp', self._tmp_dir_select,
                '--task', 'stat'
            ]
        if self._sub_step == 'selected_stat_junctions':
            return [
                '--od', self._out_dir_select_junctions, '--tmp',
                self._tmp_dir_select_junctions, '--task', 'stat'
            ]
        if self._sub_step == 'just_se':
            return [
                '--od', self._out_dir_just_se, '--tmp', self._tmp_dir_just_se,
                '--task', 'stat'
            ]
        if self._sub_step == 'deferred_stat':
            return [
                '--od', self._out_dir_all, '--tmp', self._tmp_dir_all,
                '--task', 'stat'
            ]

        return None

    def _setup_sub_step(self):
        if self._sub_step == 'selected_stat':
            self._setup_selected_stat()
        if self._sub_step == 'selected_stat_junctions':
            self._setup_selected_stat_junctions()
        if self._sub_step == 'just_se':
            self._setup_just_se()

    def _setup_selected_stat(self):
        self._prepare_stat_inputs(self._out_dir_select, self._out_dir_all, [1],
                                  [0, 3])

    def _setup_selected_stat_junctions(self):
        self._prepare_stat_inputs(self._out_dir_select_junctions,
                                  self._out_dir_all_junctions, [1], [0, 3])

    def _setup_just_se(self):
        orig_from_gtf = os.path.join(self._out_dir_all, 'fromGTF.SE.txt')
        new_from_gtf = os.path.join(self._out_dir_just_se, 'fromGTF.SE.txt')
        shutil.copy(orig_from_gtf, new_from_gtf)
        orig_raw = os.path.join(self._out_dir_all, 'JC.raw.input.SE.txt')
        new_raw = os.path.join(self._out_dir_just_se, 'JC.raw.input.SE.txt')
        shutil.copy(orig_raw, new_raw)

    def _create_sample_1_bams(self, sample_1_bams_path,
                              sample_1_replicate_template):
        rep_1_bam_path = sample_1_replicate_template.format(1)
        rep_1_bam = self._create_bam_from_paired_read_coords(
            rep_1_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords_1_1())

        rep_2_bam_path = sample_1_replicate_template.format(2)
        rep_2_bam = self._create_bam_from_paired_read_coords(
            rep_2_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords_1_2())

        sample_1_bams = [rep_1_bam, rep_2_bam]
        self._write_bams(sample_1_bams, sample_1_bams_path)
        return sample_1_bams

    def _create_sample_2_bams(self, sample_2_bams_path,
                              sample_2_replicate_template):
        rep_1_bam_path = sample_2_replicate_template.format(1)
        rep_1_bam = self._create_bam_from_paired_read_coords(
            rep_1_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords_2_1())

        rep_2_bam_path = sample_2_replicate_template.format(2)
        rep_2_bam = self._create_bam_from_paired_read_coords(
            rep_2_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords_2_2())

        sample_2_bams = [rep_1_bam, rep_2_bam]
        self._write_bams(sample_2_bams, sample_2_bams_path)
        return sample_2_bams

    def _exons_by_transcript(self):
        return [
            [(1, 100), (201, 300), (401, 500)],  # SE 1
            [(601, 700), (801, 900), (1001, 1100)],  # SE 2
            [(1201, 1300), (1401, 1500), (1801, 1900)],  # MXE
            [(1201, 1300), (1601, 1700), (1801, 1900)],  # MXE
            [(2001, 2100), (2301, 2400)],  # A5SS
            [(2001, 2200), (2301, 2400)],  # A5SS
            [(2501, 2600), (2701, 2900)],  # A3SS
            [(2501, 2600), (2801, 2900)],  # A3SS
            [(3001, 3100), (3201, 3300)],  # RI
            [(3001, 3300)],  # RI
        ]

    def _include_read_SE_1(self):
        return ([[81, 100], [201, 300]], [[201, 300]])

    def _skip_read_SE_1(self):
        return ([[81, 100], [401, 500]], [[401, 500]])

    def _include_read_SE_2(self):
        return ([[681, 700], [801, 900]], [[801, 900]])

    def _skip_read_SE_2(self):
        return ([[681, 700], [1001, 1100]], [[1001, 1100]])

    def _se_reads_from_counts(self, i1, s1, i2, s2):
        return i1 * [self._include_read_SE_1()] + s1 * [
            self._skip_read_SE_1()
        ] + i2 * [self._include_read_SE_2()] + s2 * [self._skip_read_SE_2()]

    def _mxe_read_1(self):
        return ([[1281, 1300], [1401, 1500]], [[1401, 1500]])

    def _mxe_read_2(self):
        return ([[1281, 1300], [1601, 1700]], [[1601, 1700]])

    def _a5ss_read_1(self):
        return ([[2081, 2100], [2301, 2400]], [[2301, 2400]])

    def _a5ss_read_2(self):
        return ([[2181, 2200], [2301, 2400]], [[2301, 2400]])

    def _a3ss_read_1(self):
        return ([[2581, 2600], [2701, 2900]], [[2701, 2900]])

    def _a3ss_read_2(self):
        return ([[2581, 2600], [2801, 2900]], [[2801, 2900]])

    def _ri_read_1(self):
        return ([[3081, 3100], [3201, 3300]], [[3201, 3300]])

    def _ri_read_2(self):
        return ([[3081, 3300]], [[3001, 3220]])

    def _other_reads(self):
        return [
            self._mxe_read_1(),
            self._mxe_read_2(),
            self._a5ss_read_1(),
            self._a5ss_read_2(),
            self._a3ss_read_1(),
            self._a3ss_read_2(),
            self._ri_read_1(),
            self._ri_read_2()
        ]

    def _paired_read_coords_1_1(self):
        return self._se_reads_from_counts(10, 10, 10, 0) + self._other_reads()

    def _paired_read_coords_1_2(self):
        return self._se_reads_from_counts(15, 5, 10, 0) + self._other_reads()

    def _paired_read_coords_2_1(self):
        return self._se_reads_from_counts(10, 0, 10, 10) + self._other_reads()

    def _paired_read_coords_2_2(self):
        return self._se_reads_from_counts(10, 0, 15, 5) + self._other_reads()

    def _prepare_stat_inputs(self, new_out_dir, old_out_dir, group_1_indices,
                             group_2_indices):
        command = [
            sys.executable, tests.test_config.PREPARE_STAT_INPUTS,
            '--new-output-dir', new_out_dir, '--old-output-dir', old_out_dir,
            '--group-1-indices', ','.join([str(x) for x in group_1_indices]),
            '--group-2-indices', ','.join([str(x) for x in group_2_indices])
        ]

        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)

    def _check_results(self):
        if self._sub_step == 'statoff':
            self._check_results_statoff(junctions=False)
        elif self._sub_step == 'statoff_junctions':
            self._check_results_statoff(junctions=True)
        elif self._sub_step == 'selected_stat':
            self._check_results_selected_stat(junctions=False)
        elif self._sub_step == 'selected_stat_junctions':
            self._check_results_selected_stat(junctions=True)
        elif self._sub_step == 'just_se':
            self._check_results_just_se()
        elif self._sub_step == 'deferred_stat':
            self._check_results_deferred_stat()
        else:
            self.fail('unexpected sub_step: {}'.format(self._sub_step))

    def _read_floats(self, floats_str):
        float_strs = floats_str.split(',')
        floats = list()
        for float_str in float_strs:
            try:
                floats.append(float(float_str))
            except ValueError as e:
                self.fail('could not parse {} as float from {}: {}'.format(
                    float_str, floats_str, e))

        return floats

    def _check_results_statoff(self, junctions=False):
        self._check_no_error_results()
        out_dir = self._out_dir_all
        if junctions:
            out_dir = self._out_dir_all_junctions

        se_mats_jc_path = os.path.join(out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header,
                                      has_individual=junctions)
        self.assertEqual(len(se_mats_jc_rows), 2)
        for row in se_mats_jc_rows:
            self.assertIn(row['exonStart_0base'], ['200', '800'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '10,15')
                self.assertEqual(row['SJC_SAMPLE_1'], '10,5')
                self.assertEqual(row['IJC_SAMPLE_2'], '10,10')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
                self.assertEqual(row['PValue'], 'NA')
                self.assertEqual(row['FDR'], 'NA')
                self.assertEqual(self._read_floats(row['IncLevel1']),
                                 [0.333, 0.6])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [1.0, 1.0])
                self.assertEqual(float(row['IncLevelDifference']), -0.533)
                if junctions:
                    self.assertEqual(row['upstream_to_target_count'],
                                     '10,15,10,10')
                    self.assertEqual(row['target_to_downstream_count'],
                                     '0,0,0,0')
                    self.assertEqual(row['target_count'], '10,15,10,10')
                    self.assertEqual(row['upstream_to_downstream_count'],
                                     '10,5,0,0')
                else:
                    self.assertNotIn('upstream_to_target_count', row)
            elif row['exonStart_0base'] == '800':
                self.assertEqual(row['IJC_SAMPLE_1'], '10,10')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '10,15')
                self.assertEqual(row['SJC_SAMPLE_2'], '10,5')
                self.assertEqual(row['PValue'], 'NA')
                self.assertEqual(row['FDR'], 'NA')
                self.assertEqual(self._read_floats(row['IncLevel1']),
                                 [1.0, 1.0])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [0.333, 0.6])
                self.assertEqual(float(row['IncLevelDifference']), 0.533)
                if junctions:
                    self.assertEqual(row['upstream_to_target_count'],
                                     '10,10,10,15')
                    self.assertEqual(row['target_to_downstream_count'],
                                     '0,0,0,0')
                    self.assertEqual(row['target_count'], '10,10,10,15')
                    self.assertEqual(row['upstream_to_downstream_count'],
                                     '0,0,10,5')
                else:
                    self.assertNotIn('upstream_to_target_count', row)

        se_mats_jcec_path = os.path.join(out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(se_mats_jcec_path))
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header,
                                        has_individual=junctions)
        self.assertEqual(len(se_mats_jcec_rows), 2)
        for row in se_mats_jcec_rows:
            self.assertIn(row['exonStart_0base'], ['200', '800'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '20,30')
                self.assertEqual(row['SJC_SAMPLE_1'], '10,5')
                self.assertEqual(row['IJC_SAMPLE_2'], '20,20')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
                self.assertEqual(row['PValue'], 'NA')
                self.assertEqual(row['FDR'], 'NA')
                self.assertEqual(self._read_floats(row['IncLevel1']),
                                 [0.397, 0.664])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [1.0, 1.0])
                self.assertEqual(float(row['IncLevelDifference']), -0.47)
                if junctions:
                    self.assertEqual(row['upstream_to_target_count'],
                                     '10,15,10,10')
                    self.assertEqual(row['target_to_downstream_count'],
                                     '0,0,0,0')
                    self.assertEqual(row['target_count'], '10,15,10,10')
                    self.assertEqual(row['upstream_to_downstream_count'],
                                     '10,5,0,0')
                else:
                    self.assertNotIn('upstream_to_target_count', row)

        mxe_mats_jc_path = os.path.join(out_dir, 'MXE.MATS.JC.txt')
        mxe_mats_jc_header, mxe_mats_jc_rows, error = (
            output_parser.parse_mats_jc(mxe_mats_jc_path))
        self.assertFalse(error)
        self._check_mxe_mats_jc_header(mxe_mats_jc_header,
                                       has_individual=junctions)
        self.assertEqual(len(mxe_mats_jc_rows), 1)
        self.assertEqual(mxe_mats_jc_rows[0]['FDR'], 'NA')

        mxe_mats_jcec_path = os.path.join(out_dir, 'MXE.MATS.JCEC.txt')
        mxe_mats_jcec_header, mxe_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(mxe_mats_jcec_path))
        self.assertFalse(error)
        self._check_mxe_mats_jcec_header(mxe_mats_jcec_header,
                                         has_individual=junctions)
        self.assertEqual(len(mxe_mats_jcec_rows), 1)
        self.assertEqual(mxe_mats_jcec_rows[0]['FDR'], 'NA')

        a5ss_mats_jc_path = os.path.join(out_dir, 'A5SS.MATS.JC.txt')
        a5ss_mats_jc_header, a5ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a5ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a5ss_mats_jc_header,
                                         has_individual=junctions)
        self.assertEqual(len(a5ss_mats_jc_rows), 1)
        self.assertEqual(a5ss_mats_jc_rows[0]['FDR'], 'NA')

        a5ss_mats_jcec_path = os.path.join(out_dir, 'A5SS.MATS.JCEC.txt')
        a5ss_mats_jcec_header, a5ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a5ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a5ss_mats_jcec_header,
                                           has_individual=junctions)
        self.assertEqual(len(a5ss_mats_jcec_rows), 1)
        self.assertEqual(a5ss_mats_jcec_rows[0]['FDR'], 'NA')

        a3ss_mats_jc_path = os.path.join(out_dir, 'A3SS.MATS.JC.txt')
        a3ss_mats_jc_header, a3ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a3ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a3ss_mats_jc_header,
                                         has_individual=junctions)
        self.assertEqual(len(a3ss_mats_jc_rows), 1)
        self.assertEqual(a3ss_mats_jc_rows[0]['FDR'], 'NA')

        a3ss_mats_jcec_path = os.path.join(out_dir, 'A3SS.MATS.JCEC.txt')
        a3ss_mats_jcec_header, a3ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a3ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a3ss_mats_jcec_header,
                                           has_individual=junctions)
        self.assertEqual(len(a3ss_mats_jcec_rows), 1)
        self.assertEqual(a3ss_mats_jcec_rows[0]['FDR'], 'NA')

        ri_mats_jc_path = os.path.join(out_dir, 'RI.MATS.JC.txt')
        ri_mats_jc_header, ri_mats_jc_rows, error = output_parser.parse_mats_jc(
            ri_mats_jc_path)
        self.assertFalse(error)
        self._check_ri_mats_jc_header(ri_mats_jc_header,
                                      has_individual=junctions)
        self.assertEqual(len(ri_mats_jc_rows), 1)
        self.assertEqual(ri_mats_jc_rows[0]['FDR'], 'NA')

        ri_mats_jcec_path = os.path.join(out_dir, 'RI.MATS.JCEC.txt')
        ri_mats_jcec_header, ri_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(ri_mats_jcec_path))
        self.assertFalse(error)
        self._check_ri_mats_jcec_header(ri_mats_jcec_header,
                                        has_individual=junctions)
        self.assertEqual(len(ri_mats_jcec_rows), 1)
        self.assertEqual(ri_mats_jcec_rows[0]['FDR'], 'NA')

    def _check_results_selected_stat(self, junctions=False):
        self._check_no_error_results()
        out_dir = self._out_dir_select
        if junctions:
            out_dir = self._out_dir_select_junctions

        se_mats_jc_path = os.path.join(out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header,
                                      has_individual=junctions)
        self.assertEqual(len(se_mats_jc_rows), 2)
        for row in se_mats_jc_rows:
            self.assertIn(row['exonStart_0base'], ['200', '800'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '15')
                self.assertEqual(row['SJC_SAMPLE_1'], '5')
                self.assertEqual(row['IJC_SAMPLE_2'], '10,10')
                self.assertEqual(row['SJC_SAMPLE_2'], '10,0')
                tests.util.assert_within_bounds(self, float(row['PValue']), 0,
                                                1)
                tests.util.assert_within_bounds(self, float(row['FDR']), 0, 1)
                self.assertEqual(self._read_floats(row['IncLevel1']), [0.6])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [0.333, 1.0])
                self.assertEqual(float(row['IncLevelDifference']), -0.067)
                if junctions:
                    self.assertEqual(row['upstream_to_target_count'],
                                     '15,10,10')
                    self.assertEqual(row['target_to_downstream_count'],
                                     '0,0,0')
                    self.assertEqual(row['target_count'], '15,10,10')
                    self.assertEqual(row['upstream_to_downstream_count'],
                                     '5,10,0')
                else:
                    self.assertNotIn('upstream_to_target_count', row)
            elif row['exonStart_0base'] == '800':
                self.assertEqual(row['IJC_SAMPLE_1'], '10')
                self.assertEqual(row['SJC_SAMPLE_1'], '0')
                self.assertEqual(row['IJC_SAMPLE_2'], '10,15')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,5')
                tests.util.assert_within_bounds(self, float(row['PValue']), 0,
                                                1)
                tests.util.assert_within_bounds(self, float(row['FDR']), 0, 1)
                self.assertEqual(self._read_floats(row['IncLevel1']), [1.0])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [1.0, 0.6])
                self.assertEqual(float(row['IncLevelDifference']), 0.2)
                if junctions:
                    self.assertEqual(row['upstream_to_target_count'],
                                     '10,10,15')
                    self.assertEqual(row['target_to_downstream_count'],
                                     '0,0,0')
                    self.assertEqual(row['target_count'], '10,10,15')
                    self.assertEqual(row['upstream_to_downstream_count'],
                                     '0,0,5')
                else:
                    self.assertNotIn('upstream_to_target_count', row)

        se_mats_jcec_path = os.path.join(out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(se_mats_jcec_path))
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header,
                                        has_individual=junctions)
        self.assertEqual(len(se_mats_jcec_rows), 2)
        for row in se_mats_jcec_rows:
            self.assertIn(row['exonStart_0base'], ['200', '800'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '30')
                self.assertEqual(row['SJC_SAMPLE_1'], '5')
                self.assertEqual(row['IJC_SAMPLE_2'], '20,20')
                self.assertEqual(row['SJC_SAMPLE_2'], '10,0')
                tests.util.assert_within_bounds(self, float(row['PValue']), 0,
                                                1)
                tests.util.assert_within_bounds(self, float(row['FDR']), 0, 1)
                self.assertEqual(self._read_floats(row['IncLevel1']), [0.664])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [0.397, 1.0])
                self.assertEqual(float(row['IncLevelDifference']), -0.034)
                if junctions:
                    self.assertEqual(row['upstream_to_target_count'],
                                     '15,10,10')
                    self.assertEqual(row['target_to_downstream_count'],
                                     '0,0,0')
                    self.assertEqual(row['target_count'], '15,10,10')
                    self.assertEqual(row['upstream_to_downstream_count'],
                                     '5,10,0')
                else:
                    self.assertNotIn('upstream_to_target_count', row)

        mxe_mats_jc_path = os.path.join(out_dir, 'MXE.MATS.JC.txt')
        mxe_mats_jc_header, mxe_mats_jc_rows, error = (
            output_parser.parse_mats_jc(mxe_mats_jc_path))
        self.assertFalse(error)
        self._check_mxe_mats_jc_header(mxe_mats_jc_header,
                                       has_individual=junctions)
        self.assertEqual(len(mxe_mats_jc_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(mxe_mats_jc_rows[0]['FDR']), 0,
                                        1)

        mxe_mats_jcec_path = os.path.join(out_dir, 'MXE.MATS.JCEC.txt')
        mxe_mats_jcec_header, mxe_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(mxe_mats_jcec_path))
        self.assertFalse(error)
        self._check_mxe_mats_jcec_header(mxe_mats_jcec_header,
                                         has_individual=junctions)
        self.assertEqual(len(mxe_mats_jcec_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(mxe_mats_jcec_rows[0]['FDR']), 0,
                                        1)

        a5ss_mats_jc_path = os.path.join(out_dir, 'A5SS.MATS.JC.txt')
        a5ss_mats_jc_header, a5ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a5ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a5ss_mats_jc_header,
                                         has_individual=junctions)
        self.assertEqual(len(a5ss_mats_jc_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(a5ss_mats_jc_rows[0]['FDR']), 0,
                                        1)

        a5ss_mats_jcec_path = os.path.join(out_dir, 'A5SS.MATS.JCEC.txt')
        a5ss_mats_jcec_header, a5ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a5ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a5ss_mats_jcec_header,
                                           has_individual=junctions)
        self.assertEqual(len(a5ss_mats_jcec_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(a5ss_mats_jcec_rows[0]['FDR']),
                                        0, 1)

        a3ss_mats_jc_path = os.path.join(out_dir, 'A3SS.MATS.JC.txt')
        a3ss_mats_jc_header, a3ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a3ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a3ss_mats_jc_header,
                                         has_individual=junctions)
        self.assertEqual(len(a3ss_mats_jc_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(a3ss_mats_jc_rows[0]['FDR']), 0,
                                        1)

        a3ss_mats_jcec_path = os.path.join(out_dir, 'A3SS.MATS.JCEC.txt')
        a3ss_mats_jcec_header, a3ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a3ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a3ss_mats_jcec_header,
                                           has_individual=junctions)
        self.assertEqual(len(a3ss_mats_jcec_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(a3ss_mats_jcec_rows[0]['FDR']),
                                        0, 1)

        ri_mats_jc_path = os.path.join(out_dir, 'RI.MATS.JC.txt')
        ri_mats_jc_header, ri_mats_jc_rows, error = output_parser.parse_mats_jc(
            ri_mats_jc_path)
        self.assertFalse(error)
        self._check_ri_mats_jc_header(ri_mats_jc_header,
                                      has_individual=junctions)
        self.assertEqual(len(ri_mats_jc_rows), 1)
        tests.util.assert_within_bounds(self, float(ri_mats_jc_rows[0]['FDR']),
                                        0, 1)

        ri_mats_jcec_path = os.path.join(out_dir, 'RI.MATS.JCEC.txt')
        ri_mats_jcec_header, ri_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(ri_mats_jcec_path))
        self.assertFalse(error)
        self._check_ri_mats_jcec_header(ri_mats_jcec_header,
                                        has_individual=junctions)
        self.assertEqual(len(ri_mats_jcec_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(ri_mats_jcec_rows[0]['FDR']), 0,
                                        1)

    def _check_results_just_se(self):
        self.assertEqual(self._rmats_return_code, 0)
        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        self.assertEqual(len(err_lines), 9)
        unable_to_produce_lines = collections.defaultdict(int)
        for err_line in err_lines:
            self.assertIn('Unable to produce final output files for', err_line)
            if ' SE ' in err_line:
                unable_to_produce_lines['SE'] += 1
            if ' MXE ' in err_line:
                unable_to_produce_lines['MXE'] += 1
            if ' A5SS ' in err_line:
                unable_to_produce_lines['A5SS'] += 1
            if ' A3SS ' in err_line:
                unable_to_produce_lines['A3SS'] += 1
            if ' RI ' in err_line:
                unable_to_produce_lines['RI'] += 1

        self.assertEqual(unable_to_produce_lines, {
            'SE': 1,
            'MXE': 2,
            'A5SS': 2,
            'A3SS': 2,
            'RI': 2
        })

        se_mats_jc_path = os.path.join(self._out_dir_just_se, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 2)
        for row in se_mats_jc_rows:
            self.assertIn(row['exonStart_0base'], ['200', '800'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '10,15')
                self.assertEqual(row['SJC_SAMPLE_1'], '10,5')
                self.assertEqual(row['IJC_SAMPLE_2'], '10,10')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
                tests.util.assert_within_bounds(self, float(row['PValue']), 0,
                                                1)
                tests.util.assert_within_bounds(self, float(row['FDR']), 0, 1)
                self.assertEqual(self._read_floats(row['IncLevel1']),
                                 [0.333, 0.6])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [1.0, 1.0])
                self.assertEqual(float(row['IncLevelDifference']), -0.533)

        se_mats_jcec_path = os.path.join(self._out_dir_just_se,
                                         'SE.MATS.JCEC.txt')
        self.assertFalse(os.path.exists(se_mats_jcec_path))
        mxe_mats_jc_path = os.path.join(self._out_dir_just_se,
                                        'MXE.MATS.JC.txt')
        self.assertFalse(os.path.exists(mxe_mats_jc_path))
        mxe_mats_jcec_path = os.path.join(self._out_dir_just_se,
                                          'MXE.MATS.JCEC.txt')
        self.assertFalse(os.path.exists(mxe_mats_jcec_path))
        a5ss_mats_jc_path = os.path.join(self._out_dir_just_se,
                                         'A5SS.MATS.JC.txt')
        self.assertFalse(os.path.exists(a5ss_mats_jc_path))
        a5ss_mats_jcec_path = os.path.join(self._out_dir_just_se,
                                           'A5SS.MATS.JCEC.txt')
        self.assertFalse(os.path.exists(a5ss_mats_jcec_path))
        a3ss_mats_jc_path = os.path.join(self._out_dir_just_se,
                                         'A3SS.MATS.JC.txt')
        self.assertFalse(os.path.exists(a3ss_mats_jc_path))
        a3ss_mats_jcec_path = os.path.join(self._out_dir_just_se,
                                           'A3SS.MATS.JCEC.txt')
        self.assertFalse(os.path.exists(a3ss_mats_jcec_path))
        ri_mats_jc_path = os.path.join(self._out_dir_just_se, 'RI.MATS.JC.txt')
        self.assertFalse(os.path.exists(ri_mats_jc_path))
        ri_mats_jcec_path = os.path.join(self._out_dir_just_se,
                                         'RI.MATS.JCEC.txt')
        self.assertFalse(os.path.exists(ri_mats_jcec_path))

    def _check_results_deferred_stat(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir_all, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 2)
        for row in se_mats_jc_rows:
            self.assertIn(row['exonStart_0base'], ['200', '800'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '10,15')
                self.assertEqual(row['SJC_SAMPLE_1'], '10,5')
                self.assertEqual(row['IJC_SAMPLE_2'], '10,10')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
                tests.util.assert_within_bounds(self, float(row['PValue']), 0,
                                                1)
                tests.util.assert_within_bounds(self, float(row['FDR']), 0, 1)
                self.assertEqual(self._read_floats(row['IncLevel1']),
                                 [0.333, 0.6])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [1.0, 1.0])
                self.assertEqual(float(row['IncLevelDifference']), -0.533)
            elif row['exonStart_0base'] == '800':
                self.assertEqual(row['IJC_SAMPLE_1'], '10,10')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '10,15')
                self.assertEqual(row['SJC_SAMPLE_2'], '10,5')
                tests.util.assert_within_bounds(self, float(row['PValue']), 0,
                                                1)
                tests.util.assert_within_bounds(self, float(row['FDR']), 0, 1)
                self.assertEqual(self._read_floats(row['IncLevel1']),
                                 [1.0, 1.0])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [0.333, 0.6])
                self.assertEqual(float(row['IncLevelDifference']), 0.533)

        se_mats_jcec_path = os.path.join(self._out_dir_all, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(se_mats_jcec_path))
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 2)
        for row in se_mats_jcec_rows:
            self.assertIn(row['exonStart_0base'], ['200', '800'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '20,30')
                self.assertEqual(row['SJC_SAMPLE_1'], '10,5')
                self.assertEqual(row['IJC_SAMPLE_2'], '20,20')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
                tests.util.assert_within_bounds(self, float(row['PValue']), 0,
                                                1)
                tests.util.assert_within_bounds(self, float(row['FDR']), 0, 1)
                self.assertEqual(self._read_floats(row['IncLevel1']),
                                 [0.397, 0.664])
                self.assertEqual(self._read_floats(row['IncLevel2']),
                                 [1.0, 1.0])
                self.assertEqual(float(row['IncLevelDifference']), -0.47)

        mxe_mats_jc_path = os.path.join(self._out_dir_all, 'MXE.MATS.JC.txt')
        mxe_mats_jc_header, mxe_mats_jc_rows, error = (
            output_parser.parse_mats_jc(mxe_mats_jc_path))
        self.assertFalse(error)
        self._check_mxe_mats_jc_header(mxe_mats_jc_header)
        self.assertEqual(len(mxe_mats_jc_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(mxe_mats_jc_rows[0]['FDR']), 0,
                                        1)

        mxe_mats_jcec_path = os.path.join(self._out_dir_all,
                                          'MXE.MATS.JCEC.txt')
        mxe_mats_jcec_header, mxe_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(mxe_mats_jcec_path))
        self.assertFalse(error)
        self._check_mxe_mats_jcec_header(mxe_mats_jcec_header)
        self.assertEqual(len(mxe_mats_jcec_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(mxe_mats_jcec_rows[0]['FDR']), 0,
                                        1)

        a5ss_mats_jc_path = os.path.join(self._out_dir_all, 'A5SS.MATS.JC.txt')
        a5ss_mats_jc_header, a5ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a5ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a5ss_mats_jc_header)
        self.assertEqual(len(a5ss_mats_jc_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(a5ss_mats_jc_rows[0]['FDR']), 0,
                                        1)

        a5ss_mats_jcec_path = os.path.join(self._out_dir_all,
                                           'A5SS.MATS.JCEC.txt')
        a5ss_mats_jcec_header, a5ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a5ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a5ss_mats_jcec_header)
        self.assertEqual(len(a5ss_mats_jcec_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(a5ss_mats_jcec_rows[0]['FDR']),
                                        0, 1)

        a3ss_mats_jc_path = os.path.join(self._out_dir_all, 'A3SS.MATS.JC.txt')
        a3ss_mats_jc_header, a3ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a3ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a3ss_mats_jc_header)
        self.assertEqual(len(a3ss_mats_jc_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(a3ss_mats_jc_rows[0]['FDR']), 0,
                                        1)

        a3ss_mats_jcec_path = os.path.join(self._out_dir_all,
                                           'A3SS.MATS.JCEC.txt')
        a3ss_mats_jcec_header, a3ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a3ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a3ss_mats_jcec_header)
        self.assertEqual(len(a3ss_mats_jcec_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(a3ss_mats_jcec_rows[0]['FDR']),
                                        0, 1)

        ri_mats_jc_path = os.path.join(self._out_dir_all, 'RI.MATS.JC.txt')
        ri_mats_jc_header, ri_mats_jc_rows, error = output_parser.parse_mats_jc(
            ri_mats_jc_path)
        self.assertFalse(error)
        self._check_ri_mats_jc_header(ri_mats_jc_header)
        self.assertEqual(len(ri_mats_jc_rows), 1)
        tests.util.assert_within_bounds(self, float(ri_mats_jc_rows[0]['FDR']),
                                        0, 1)

        ri_mats_jcec_path = os.path.join(self._out_dir_all, 'RI.MATS.JCEC.txt')
        ri_mats_jcec_header, ri_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(ri_mats_jcec_path))
        self.assertFalse(error)
        self._check_ri_mats_jcec_header(ri_mats_jcec_header)
        self.assertEqual(len(ri_mats_jcec_rows), 1)
        tests.util.assert_within_bounds(self,
                                        float(ri_mats_jcec_rows[0]['FDR']), 0,
                                        1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
