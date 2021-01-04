import os.path
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class IndividualCountsBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir, 'individual_counts',
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
        self._chromosome_length = 10000
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
        ]

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
            # SE
            [(1, 100), (201, 300), (401, 500)],
            [(1, 100), (401, 500)],
            [(1001, 1100), (1201, 1300), (1401, 1500)],
            [(1001, 1100), (1401, 1500)],
            # A3SS
            [(2001, 2100), (2201, 2400)],
            [(2001, 2100), (2301, 2400)],
            [(3001, 3100), (3201, 3400)],
            [(3001, 3100), (3301, 3400)],
            # A5SS
            [(4001, 4100), (4301, 4400)],
            [(4001, 4200), (4301, 4400)],
            [(5001, 5100), (5301, 5400)],
            [(5001, 5200), (5301, 5400)],
            # MXE
            [(6001, 6100), (6201, 6300), (6601, 6700)],
            [(6001, 6100), (6401, 6500), (6601, 6700)],
            [(7001, 7100), (7201, 7300), (7601, 7700)],
            [(7001, 7100), (7401, 7500), (7601, 7700)],
            # RI
            [(8001, 8100), (8201, 8300)],
            [(8001, 8300)],
            [(9001, 9100), (9201, 9300)],
            [(9001, 9300)],
        ]

    def _se_1_upstream_to_target_read(self):
        return ([[81, 100], [201, 300]], [[1, 100], [201, 220]])

    def _se_2_upstream_to_target_read(self):
        return ([[1081, 1100], [1201, 1300]], [[1001, 1100], [1201, 1220]])

    def _se_1_target_to_downstream_read(self):
        return ([[281, 300], [401, 500]], [[201, 300], [401, 420]])

    def _se_2_target_to_downstream_read(self):
        return ([[1281, 1300], [1401, 1500]], [[1201, 1300], [1401, 1420]])

    def _se_1_target_read(self):
        return ([[201, 300]], [[201, 300]])

    def _se_2_target_read(self):
        return ([[1201, 1300]], [[1201, 1300]])

    def _se_1_upstream_to_downstream_read(self):
        return ([[81, 100], [401, 500]], [[1, 100], [401, 420]])

    def _se_2_upstream_to_downstream_read(self):
        return ([[1081, 1100], [1401, 1500]], [[1001, 1100], [1401, 1420]])

    def _se_reads(self, u_t_1, t_d_1, t_1, u_d_1, u_t_2, t_d_2, t_2, u_d_2):
        return (([self._se_1_upstream_to_target_read()] * u_t_1) +
                ([self._se_1_target_to_downstream_read()] * t_d_1) +
                ([self._se_1_target_read()] * t_1) +
                ([self._se_1_upstream_to_downstream_read()] * u_d_1) +
                ([self._se_2_upstream_to_target_read()] * u_t_2) +
                ([self._se_2_target_to_downstream_read()] * t_d_2) +
                ([self._se_2_target_read()] * t_2) +
                ([self._se_2_upstream_to_downstream_read()] * u_d_2))

    def _a3ss_1_across_short_boundary_read(self):
        return ([[2281, 2400]], [[2201, 2320]])

    def _a3ss_2_across_short_boundary_read(self):
        return ([[3281, 3400]], [[3201, 3320]])

    def _a3ss_1_long_to_flanking_read(self):
        return ([[2081, 2100], [2201, 2400]], [[2001, 2100], [2201, 2220]])

    def _a3ss_2_long_to_flanking_read(self):
        return ([[3081, 3100], [3201, 3400]], [[3001, 3100], [3201, 3220]])

    def _a3ss_1_exclusive_to_long_read(self):
        return ([[2201, 2300]], [[2201, 2300]])

    def _a3ss_2_exclusive_to_long_read(self):
        return ([[3201, 3300]], [[3201, 3300]])

    def _a3ss_1_short_to_flanking_read(self):
        return ([[2081, 2100], [2301, 2400]], [[2001, 2100], [2301, 2320]])

    def _a3ss_2_short_to_flanking_read(self):
        return ([[3081, 3100], [3301, 3400]], [[3001, 3100], [3301, 3320]])

    def _a3ss_reads(self, s_b_1, l_f_1, l_1, s_f_1, s_b_2, l_f_2, l_2, s_f_2):
        return (([self._a3ss_1_across_short_boundary_read()] * s_b_1) +
                ([self._a3ss_1_long_to_flanking_read()] * l_f_1) +
                ([self._a3ss_1_exclusive_to_long_read()] * l_1) +
                ([self._a3ss_1_short_to_flanking_read()] * s_f_1) +
                ([self._a3ss_2_across_short_boundary_read()] * s_b_2) +
                ([self._a3ss_2_long_to_flanking_read()] * l_f_2) +
                ([self._a3ss_2_exclusive_to_long_read()] * l_2) +
                ([self._a3ss_2_short_to_flanking_read()] * s_f_2))

    def _a5ss_1_across_short_boundary_read(self):
        return ([[4081, 4200]], [[4001, 4120]])

    def _a5ss_2_across_short_boundary_read(self):
        return ([[5081, 5200]], [[5001, 5120]])

    def _a5ss_1_long_to_flanking_read(self):
        return ([[4181, 4200], [4301, 4400]], [[4001, 4200], [4301, 4320]])

    def _a5ss_2_long_to_flanking_read(self):
        return ([[5181, 5200], [5301, 5400]], [[5001, 5200], [5301, 5320]])

    def _a5ss_1_exclusive_to_long_read(self):
        return ([[4101, 4200]], [[4101, 4200]])

    def _a5ss_2_exclusive_to_long_read(self):
        return ([[5101, 5200]], [[5101, 5200]])

    def _a5ss_1_short_to_flanking_read(self):
        return ([[4081, 4100], [4301, 4400]], [[4001, 4100], [4301, 4320]])

    def _a5ss_2_short_to_flanking_read(self):
        return ([[5081, 5100], [5301, 5400]], [[5001, 5100], [5301, 5320]])

    def _a5ss_reads(self, s_b_1, l_f_1, l_1, s_f_1, s_b_2, l_f_2, l_2, s_f_2):
        return (([self._a5ss_1_across_short_boundary_read()] * s_b_1) +
                ([self._a5ss_1_long_to_flanking_read()] * l_f_1) +
                ([self._a5ss_1_exclusive_to_long_read()] * l_1) +
                ([self._a5ss_1_short_to_flanking_read()] * s_f_1) +
                ([self._a5ss_2_across_short_boundary_read()] * s_b_2) +
                ([self._a5ss_2_long_to_flanking_read()] * l_f_2) +
                ([self._a5ss_2_exclusive_to_long_read()] * l_2) +
                ([self._a5ss_2_short_to_flanking_read()] * s_f_2))

    def _mxe_1_upstream_to_first_read(self):
        return ([[6081, 6100], [6201, 6300]], [[6001, 6100], [6201, 6220]])

    def _mxe_2_upstream_to_first_read(self):
        return ([[7081, 7100], [7201, 7300]], [[7001, 7100], [7201, 7220]])

    def _mxe_1_first_to_downstream_read(self):
        return ([[6281, 6300], [6601, 6700]], [[6201, 6300], [6601, 6620]])

    def _mxe_2_first_to_downstream_read(self):
        return ([[7281, 7300], [7601, 7700]], [[7201, 7300], [7601, 7620]])

    def _mxe_1_first_read(self):
        return ([[6201, 6300]], [[6201, 6300]])

    def _mxe_2_first_read(self):
        return ([[7201, 7300]], [[7201, 7300]])

    def _mxe_1_upstream_to_second_read(self):
        return ([[6081, 6100], [6401, 6500]], [[6001, 6100], [6401, 6420]])

    def _mxe_2_upstream_to_second_read(self):
        return ([[7081, 7100], [7401, 7500]], [[7001, 7100], [7401, 7420]])

    def _mxe_1_second_to_downstream_read(self):
        return ([[6481, 6500], [6601, 6700]], [[6401, 6500], [6601, 6620]])

    def _mxe_2_second_to_downstream_read(self):
        return ([[7481, 7500], [7601, 7700]], [[7401, 7500], [7601, 7620]])

    def _mxe_1_second_read(self):
        return ([[6401, 6500]], [[6401, 6500]])

    def _mxe_2_second_read(self):
        return ([[7401, 7500]], [[7401, 7500]])

    def _mxe_reads(self, u_f_1, f_d_1, f_1, u_s_1, s_d_1, s_1, u_f_2, f_d_2,
                   f_2, u_s_2, s_d_2, s_2):
        return (([self._mxe_1_upstream_to_first_read()] * u_f_1) +
                ([self._mxe_1_first_to_downstream_read()] * f_d_1) +
                ([self._mxe_1_first_read()] * f_1) +
                ([self._mxe_1_upstream_to_second_read()] * u_s_1) +
                ([self._mxe_1_second_to_downstream_read()] * s_d_1) +
                ([self._mxe_1_second_read()] * s_1) +
                ([self._mxe_2_upstream_to_first_read()] * u_f_2) +
                ([self._mxe_2_first_to_downstream_read()] * f_d_2) +
                ([self._mxe_2_first_read()] * f_2) +
                ([self._mxe_2_upstream_to_second_read()] * u_s_2) +
                ([self._mxe_2_second_to_downstream_read()] * s_d_2) +
                ([self._mxe_2_second_read()] * s_2))

    def _ri_1_upstream_to_intron_read(self):
        return ([[8081, 8200]], [[8001, 8120]])

    def _ri_2_upstream_to_intron_read(self):
        return ([[9081, 9200]], [[9001, 9120]])

    def _ri_1_intron_to_downstream_read(self):
        return ([[8181, 8300]], [[8101, 8220]])

    def _ri_2_intron_to_downstream_read(self):
        return ([[9181, 9300]], [[9101, 9220]])

    def _ri_1_intron_read(self):
        return ([[8101, 8200]], [[8101, 8200]])

    def _ri_2_intron_read(self):
        return ([[9101, 9200]], [[9101, 9200]])

    def _ri_1_upstream_to_downstream_read(self):
        return ([[8081, 8100], [8201, 8300]], [[8001, 8100], [8201, 8220]])

    def _ri_2_upstream_to_downstream_read(self):
        return ([[9081, 9100], [9201, 9300]], [[9001, 9100], [9201, 9220]])

    def _ri_reads(self, u_i_1, i_d_1, i_1, u_d_1, u_i_2, i_d_2, i_2, u_d_2):
        return (([self._ri_1_upstream_to_intron_read()] * u_i_1) +
                ([self._ri_1_intron_to_downstream_read()] * i_d_1) +
                ([self._ri_1_intron_read()] * i_1) +
                ([self._ri_1_upstream_to_downstream_read()] * u_d_1) +
                ([self._ri_2_upstream_to_intron_read()] * u_i_2) +
                ([self._ri_2_intron_to_downstream_read()] * i_d_2) +
                ([self._ri_2_intron_read()] * i_2) +
                ([self._ri_2_upstream_to_downstream_read()] * u_d_2))

    def _paired_read_coords_1_1(self):
        return (
            self._se_reads(10, 11, 12, 13, 20, 1, 22, 23) +
            self._a3ss_reads(10, 11, 12, 13, 20, 1, 22, 23) +
            self._a5ss_reads(10, 11, 12, 13, 20, 1, 22, 23) +
            self._mxe_reads(10, 11, 12, 13, 14, 15, 20, 1, 22, 23, 24, 25) +
            self._ri_reads(10, 11, 12, 13, 20, 1, 22, 23))

    def _paired_read_coords_1_2(self):
        return (
            self._se_reads(15, 16, 17, 18, 25, 2, 27, 28) +
            self._a3ss_reads(15, 16, 17, 18, 25, 2, 27, 28) +
            self._a5ss_reads(15, 16, 17, 18, 25, 2, 27, 28) +
            self._mxe_reads(15, 16, 17, 18, 19, 20, 25, 2, 27, 28, 29, 30) +
            self._ri_reads(15, 16, 17, 18, 25, 2, 27, 28))

    def _paired_read_coords_2_1(self):
        return (
            self._se_reads(20, 21, 22, 23, 30, 3, 32, 33) +
            self._a3ss_reads(20, 21, 22, 23, 30, 3, 32, 33) +
            self._a5ss_reads(20, 21, 22, 23, 30, 3, 32, 33) +
            self._mxe_reads(20, 21, 22, 23, 24, 25, 30, 3, 32, 33, 34, 35) +
            self._ri_reads(20, 21, 22, 23, 30, 3, 32, 33))

    def _paired_read_coords_2_2(self):
        return (
            self._se_reads(25, 26, 27, 28, 35, 4, 37, 38) +
            self._a3ss_reads(25, 26, 27, 28, 35, 4, 37, 38) +
            self._a5ss_reads(25, 26, 27, 28, 35, 4, 37, 38) +
            self._mxe_reads(25, 26, 27, 28, 29, 30, 35, 4, 37, 38, 39, 40) +
            self._ri_reads(25, 26, 27, 28, 35, 4, 37, 38))

    def _check_se_1_jc(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '42,62')
        self.assertEqual(row['IJC_SAMPLE_2'], '82,102')
        self._check_se_1_shared(row)

    def _check_se_2_jc(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '42,54')
        self.assertEqual(row['IJC_SAMPLE_2'], '66,78')
        self._check_se_2_shared(row)

    def _check_se_1_jcec(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '66,96')
        self.assertEqual(row['IJC_SAMPLE_2'], '126,156')
        self._check_se_1_shared(row)

    def _check_se_2_jcec(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '86,108')
        self.assertEqual(row['IJC_SAMPLE_2'], '130,152')
        self._check_se_2_shared(row)

    def _check_se_1_shared(self, row):
        self.assertEqual(row['SJC_SAMPLE_1'], '26,36')
        self.assertEqual(row['SJC_SAMPLE_2'], '46,56')
        self.assertEqual(row['upstream_to_target_count'], '20,30,40,50')
        self.assertEqual(row['target_to_downstream_count'], '22,32,42,52')
        self.assertEqual(row['target_count'], '24,34,44,54')
        self.assertEqual(row['upstream_to_downstream_count'], '26,36,46,56')

    def _check_se_2_shared(self, row):
        self.assertEqual(row['SJC_SAMPLE_1'], '46,56')
        self.assertEqual(row['SJC_SAMPLE_2'], '66,76')
        self.assertEqual(row['upstream_to_target_count'], '40,50,60,70')
        self.assertEqual(row['target_to_downstream_count'], '2,4,6,8')
        self.assertEqual(row['target_count'], '44,54,64,74')
        self.assertEqual(row['upstream_to_downstream_count'], '46,56,66,76')

    def _check_alt35_1_jc(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '42,62')
        self.assertEqual(row['IJC_SAMPLE_2'], '82,102')
        self._check_alt35_1_shared(row)

    def _check_alt35_2_jc(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '42,54')
        self.assertEqual(row['IJC_SAMPLE_2'], '66,78')
        self._check_alt35_2_shared(row)

    def _check_alt35_1_jcec(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '66,96')
        self.assertEqual(row['IJC_SAMPLE_2'], '126,156')
        self._check_alt35_1_shared(row)

    def _check_alt35_2_jcec(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '86,108')
        self.assertEqual(row['IJC_SAMPLE_2'], '130,152')
        self._check_alt35_2_shared(row)

    def _check_alt35_1_shared(self, row):
        self.assertEqual(row['SJC_SAMPLE_1'], '26,36')
        self.assertEqual(row['SJC_SAMPLE_2'], '46,56')
        self.assertEqual(row['across_short_boundary_count'], '20,30,40,50')
        self.assertEqual(row['long_to_flanking_count'], '22,32,42,52')
        self.assertEqual(row['exclusive_to_long_count'], '24,34,44,54')
        self.assertEqual(row['short_to_flanking_count'], '26,36,46,56')

    def _check_alt35_2_shared(self, row):
        self.assertEqual(row['SJC_SAMPLE_1'], '46,56')
        self.assertEqual(row['SJC_SAMPLE_2'], '66,76')
        self.assertEqual(row['across_short_boundary_count'], '40,50,60,70')
        self.assertEqual(row['long_to_flanking_count'], '2,4,6,8')
        self.assertEqual(row['exclusive_to_long_count'], '44,54,64,74')
        self.assertEqual(row['short_to_flanking_count'], '46,56,66,76')

    def _check_mxe_1_jc(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '42,62')
        self.assertEqual(row['IJC_SAMPLE_2'], '82,102')
        self.assertEqual(row['SJC_SAMPLE_1'], '54,74')
        self.assertEqual(row['SJC_SAMPLE_2'], '94,114')
        self._check_mxe_1_shared(row)

    def _check_mxe_2_jc(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '42,54')
        self.assertEqual(row['IJC_SAMPLE_2'], '66,78')
        self.assertEqual(row['SJC_SAMPLE_1'], '94,114')
        self.assertEqual(row['SJC_SAMPLE_2'], '134,154')
        self._check_mxe_2_shared(row)

    def _check_mxe_1_jcec(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '66,96')
        self.assertEqual(row['IJC_SAMPLE_2'], '126,156')
        self.assertEqual(row['SJC_SAMPLE_1'], '84,114')
        self.assertEqual(row['SJC_SAMPLE_2'], '144,174')
        self._check_mxe_1_shared(row)

    def _check_mxe_2_jcec(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '86,108')
        self.assertEqual(row['IJC_SAMPLE_2'], '130,152')
        self.assertEqual(row['SJC_SAMPLE_1'], '144,174')
        self.assertEqual(row['SJC_SAMPLE_2'], '204,234')
        self._check_mxe_2_shared(row)

    def _check_mxe_1_shared(self, row):
        self.assertEqual(row['upstream_to_first_count'], '20,30,40,50')
        self.assertEqual(row['first_to_downstream_count'], '22,32,42,52')
        self.assertEqual(row['first_count'], '24,34,44,54')
        self.assertEqual(row['upstream_to_second_count'], '26,36,46,56')
        self.assertEqual(row['second_to_downstream_count'], '28,38,48,58')
        self.assertEqual(row['second_count'], '30,40,50,60')

    def _check_mxe_2_shared(self, row):
        self.assertEqual(row['upstream_to_first_count'], '40,50,60,70')
        self.assertEqual(row['first_to_downstream_count'], '2,4,6,8')
        self.assertEqual(row['first_count'], '44,54,64,74')
        self.assertEqual(row['upstream_to_second_count'], '46,56,66,76')
        self.assertEqual(row['second_to_downstream_count'], '48,58,68,78')
        self.assertEqual(row['second_count'], '50,60,70,80')

    def _check_ri_1_jc(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '42,62')
        self.assertEqual(row['IJC_SAMPLE_2'], '82,102')
        self._check_ri_1_shared(row)

    def _check_ri_2_jc(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '42,54')
        self.assertEqual(row['IJC_SAMPLE_2'], '66,78')
        self._check_ri_2_shared(row)

    def _check_ri_1_jcec(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '66,96')
        self.assertEqual(row['IJC_SAMPLE_2'], '126,156')
        self._check_ri_1_shared(row)

    def _check_ri_2_jcec(self, row):
        self.assertEqual(row['IJC_SAMPLE_1'], '86,108')
        self.assertEqual(row['IJC_SAMPLE_2'], '130,152')
        self._check_ri_2_shared(row)

    def _check_ri_1_shared(self, row):
        self.assertEqual(row['SJC_SAMPLE_1'], '26,36')
        self.assertEqual(row['SJC_SAMPLE_2'], '46,56')
        self.assertEqual(row['upstream_to_intron_count'], '20,30,40,50')
        self.assertEqual(row['intron_to_downstream_count'], '22,32,42,52')
        self.assertEqual(row['intron_count'], '24,34,44,54')
        self.assertEqual(row['upstream_to_downstream_count'], '26,36,46,56')

    def _check_ri_2_shared(self, row):
        self.assertEqual(row['SJC_SAMPLE_1'], '46,56')
        self.assertEqual(row['SJC_SAMPLE_2'], '66,76')
        self.assertEqual(row['upstream_to_intron_count'], '40,50,60,70')
        self.assertEqual(row['intron_to_downstream_count'], '2,4,6,8')
        self.assertEqual(row['intron_count'], '44,54,64,74')
        self.assertEqual(row['upstream_to_downstream_count'], '46,56,66,76')


class FilteredTest(IndividualCountsBaseTest):
    def _sub_test_name(self):
        return 'filtered'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend([
            '--imbalance-ratio',
            '5',
        ])
        return arguments

    def _check_results(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 1)
        self.assertEqual(se_mats_jc_rows[0]['upstreamES'], '0')
        self._check_se_1_jc(se_mats_jc_rows[0])

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            se_mats_jcec_path)
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 1)
        self.assertEqual(se_mats_jcec_rows[0]['upstreamES'], '0')
        self._check_se_1_jcec(se_mats_jcec_rows[0])

        a3ss_mats_jc_path = os.path.join(self._out_dir, 'A3SS.MATS.JC.txt')
        a3ss_mats_jc_header, a3ss_mats_jc_rows, error = output_parser.parse_mats_jc(
            a3ss_mats_jc_path)
        self.assertFalse(error)
        self._check_alt35_mats_jc_header(a3ss_mats_jc_header)
        self.assertEqual(len(a3ss_mats_jc_rows), 1)
        self.assertEqual(a3ss_mats_jc_rows[0]['flankingES'], '2000')
        self._check_alt35_1_jc(a3ss_mats_jc_rows[0])

        a3ss_mats_jcec_path = os.path.join(self._out_dir, 'A3SS.MATS.JCEC.txt')
        a3ss_mats_jcec_header, a3ss_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            a3ss_mats_jcec_path)
        self.assertFalse(error)
        self._check_alt35_mats_jcec_header(a3ss_mats_jcec_header)
        self.assertEqual(len(a3ss_mats_jcec_rows), 1)
        self.assertEqual(a3ss_mats_jcec_rows[0]['flankingES'], '2000')
        self._check_alt35_1_jcec(a3ss_mats_jcec_rows[0])

        a5ss_mats_jc_path = os.path.join(self._out_dir, 'A5SS.MATS.JC.txt')
        a5ss_mats_jc_header, a5ss_mats_jc_rows, error = output_parser.parse_mats_jc(
            a5ss_mats_jc_path)
        self.assertFalse(error)
        self._check_alt35_mats_jc_header(a5ss_mats_jc_header)
        self.assertEqual(len(a5ss_mats_jc_rows), 1)
        self.assertEqual(a5ss_mats_jc_rows[0]['shortES'], '4000')
        self._check_alt35_1_jc(a5ss_mats_jc_rows[0])

        a5ss_mats_jcec_path = os.path.join(self._out_dir, 'A5SS.MATS.JCEC.txt')
        a5ss_mats_jcec_header, a5ss_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            a5ss_mats_jcec_path)
        self.assertFalse(error)
        self._check_alt35_mats_jcec_header(a5ss_mats_jcec_header)
        self.assertEqual(len(a5ss_mats_jcec_rows), 1)
        self.assertEqual(a5ss_mats_jcec_rows[0]['shortES'], '4000')
        self._check_alt35_1_jcec(a5ss_mats_jcec_rows[0])

        mxe_mats_jc_path = os.path.join(self._out_dir, 'MXE.MATS.JC.txt')
        mxe_mats_jc_header, mxe_mats_jc_rows, error = output_parser.parse_mats_jc(
            mxe_mats_jc_path)
        self.assertFalse(error)
        self._check_mxe_mats_jc_header(mxe_mats_jc_header)
        self.assertEqual(len(mxe_mats_jc_rows), 1)
        self.assertEqual(mxe_mats_jc_rows[0]['upstreamES'], '6000')
        self._check_mxe_1_jc(mxe_mats_jc_rows[0])

        mxe_mats_jcec_path = os.path.join(self._out_dir, 'MXE.MATS.JCEC.txt')
        mxe_mats_jcec_header, mxe_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            mxe_mats_jcec_path)
        self.assertFalse(error)
        self._check_mxe_mats_jcec_header(mxe_mats_jcec_header)
        self.assertEqual(len(mxe_mats_jcec_rows), 1)
        self.assertEqual(mxe_mats_jcec_rows[0]['upstreamES'], '6000')
        self._check_mxe_1_jcec(mxe_mats_jcec_rows[0])

        ri_mats_jc_path = os.path.join(self._out_dir, 'RI.MATS.JC.txt')
        ri_mats_jc_header, ri_mats_jc_rows, error = output_parser.parse_mats_jc(
            ri_mats_jc_path)
        self.assertFalse(error)
        self._check_ri_mats_jc_header(ri_mats_jc_header)
        self.assertEqual(len(ri_mats_jc_rows), 1)
        self.assertEqual(ri_mats_jc_rows[0]['upstreamES'], '8000')
        self._check_ri_1_jc(ri_mats_jc_rows[0])

        ri_mats_jcec_path = os.path.join(self._out_dir, 'RI.MATS.JCEC.txt')
        ri_mats_jcec_header, ri_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            ri_mats_jcec_path)
        self.assertFalse(error)
        self._check_ri_mats_jcec_header(ri_mats_jcec_header)
        self.assertEqual(len(ri_mats_jcec_rows), 1)
        self.assertEqual(ri_mats_jcec_rows[0]['upstreamES'], '8000')
        self._check_ri_1_jcec(ri_mats_jcec_rows[0])


class NotFilteredTest(IndividualCountsBaseTest):
    def _sub_test_name(self):
        return 'not_filtered'

    def test(self):
        self._run_test()

    def _check_results(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 2)
        se_mats_jc_row_0_start = se_mats_jc_rows[0]['upstreamES']
        se_mats_jc_row_1_start = se_mats_jc_rows[1]['upstreamES']
        self.assertIn(se_mats_jc_row_0_start, ['0', '1000'])
        self.assertIn(se_mats_jc_row_1_start, ['0', '1000'])
        self.assertNotEqual(se_mats_jc_row_0_start, se_mats_jc_row_1_start)
        if se_mats_jc_row_0_start == '0':
            self._check_se_1_jc(se_mats_jc_rows[0])
            self._check_se_2_jc(se_mats_jc_rows[1])
        else:
            self._check_se_1_jc(se_mats_jc_rows[1])
            self._check_se_2_jc(se_mats_jc_rows[0])

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            se_mats_jcec_path)
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 2)
        se_mats_jcec_row_0_start = se_mats_jcec_rows[0]['upstreamES']
        se_mats_jcec_row_1_start = se_mats_jcec_rows[1]['upstreamES']
        self.assertIn(se_mats_jcec_row_0_start, ['0', '1000'])
        self.assertIn(se_mats_jcec_row_1_start, ['0', '1000'])
        self.assertNotEqual(se_mats_jcec_row_0_start, se_mats_jcec_row_1_start)
        if se_mats_jcec_row_0_start == '0':
            self._check_se_1_jcec(se_mats_jcec_rows[0])
            self._check_se_2_jcec(se_mats_jcec_rows[1])
        else:
            self._check_se_1_jcec(se_mats_jcec_rows[1])
            self._check_se_2_jcec(se_mats_jcec_rows[0])

        a3ss_mats_jc_path = os.path.join(self._out_dir, 'A3SS.MATS.JC.txt')
        a3ss_mats_jc_header, a3ss_mats_jc_rows, error = output_parser.parse_mats_jc(
            a3ss_mats_jc_path)
        self.assertFalse(error)
        self._check_alt35_mats_jc_header(a3ss_mats_jc_header)
        self.assertEqual(len(a3ss_mats_jc_rows), 2)
        a3ss_mats_jc_row_0_start = a3ss_mats_jc_rows[0]['flankingES']
        a3ss_mats_jc_row_1_start = a3ss_mats_jc_rows[1]['flankingES']
        self.assertIn(a3ss_mats_jc_row_0_start, ['2000', '3000'])
        self.assertIn(a3ss_mats_jc_row_1_start, ['2000', '3000'])
        self.assertNotEqual(a3ss_mats_jc_row_0_start, a3ss_mats_jc_row_1_start)
        if a3ss_mats_jc_row_0_start == '2000':
            self._check_alt35_1_jc(a3ss_mats_jc_rows[0])
            self._check_alt35_2_jc(a3ss_mats_jc_rows[1])
        else:
            self._check_alt35_1_jc(a3ss_mats_jc_rows[1])
            self._check_alt35_2_jc(a3ss_mats_jc_rows[0])

        a3ss_mats_jcec_path = os.path.join(self._out_dir, 'A3SS.MATS.JCEC.txt')
        a3ss_mats_jcec_header, a3ss_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            a3ss_mats_jcec_path)
        self.assertFalse(error)
        self._check_alt35_mats_jcec_header(a3ss_mats_jcec_header)
        self.assertEqual(len(a3ss_mats_jcec_rows), 2)
        a3ss_mats_jcec_row_0_start = a3ss_mats_jcec_rows[0]['flankingES']
        a3ss_mats_jcec_row_1_start = a3ss_mats_jcec_rows[1]['flankingES']
        self.assertIn(a3ss_mats_jcec_row_0_start, ['2000', '3000'])
        self.assertIn(a3ss_mats_jcec_row_1_start, ['2000', '3000'])
        self.assertNotEqual(a3ss_mats_jcec_row_0_start,
                            a3ss_mats_jcec_row_1_start)
        if a3ss_mats_jcec_row_0_start == '2000':
            self._check_alt35_1_jcec(a3ss_mats_jcec_rows[0])
            self._check_alt35_2_jcec(a3ss_mats_jcec_rows[1])
        else:
            self._check_alt35_1_jcec(a3ss_mats_jcec_rows[1])
            self._check_alt35_2_jcec(a3ss_mats_jcec_rows[0])

        a5ss_mats_jc_path = os.path.join(self._out_dir, 'A5SS.MATS.JC.txt')
        a5ss_mats_jc_header, a5ss_mats_jc_rows, error = output_parser.parse_mats_jc(
            a5ss_mats_jc_path)
        self.assertFalse(error)
        self._check_alt35_mats_jc_header(a5ss_mats_jc_header)
        self.assertEqual(len(a5ss_mats_jc_rows), 2)
        a5ss_mats_jc_row_0_start = a5ss_mats_jc_rows[0]['shortES']
        a5ss_mats_jc_row_1_start = a5ss_mats_jc_rows[1]['shortES']
        self.assertIn(a5ss_mats_jc_row_0_start, ['4000', '5000'])
        self.assertIn(a5ss_mats_jc_row_1_start, ['4000', '5000'])
        self.assertNotEqual(a5ss_mats_jc_row_0_start, a5ss_mats_jc_row_1_start)
        if a5ss_mats_jc_row_0_start == '4000':
            self._check_alt35_1_jc(a5ss_mats_jc_rows[0])
            self._check_alt35_2_jc(a5ss_mats_jc_rows[1])
        else:
            self._check_alt35_1_jc(a5ss_mats_jc_rows[1])
            self._check_alt35_2_jc(a5ss_mats_jc_rows[0])

        a5ss_mats_jcec_path = os.path.join(self._out_dir, 'A5SS.MATS.JCEC.txt')
        a5ss_mats_jcec_header, a5ss_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            a5ss_mats_jcec_path)
        self.assertFalse(error)
        self._check_alt35_mats_jcec_header(a5ss_mats_jcec_header)
        self.assertEqual(len(a5ss_mats_jcec_rows), 2)
        a5ss_mats_jcec_row_0_start = a5ss_mats_jcec_rows[0]['shortES']
        a5ss_mats_jcec_row_1_start = a5ss_mats_jcec_rows[1]['shortES']
        self.assertIn(a5ss_mats_jcec_row_0_start, ['4000', '5000'])
        self.assertIn(a5ss_mats_jcec_row_1_start, ['4000', '5000'])
        self.assertNotEqual(a5ss_mats_jcec_row_0_start,
                            a5ss_mats_jcec_row_1_start)
        if a5ss_mats_jcec_row_0_start == '4000':
            self._check_alt35_1_jcec(a5ss_mats_jcec_rows[0])
            self._check_alt35_2_jcec(a5ss_mats_jcec_rows[1])
        else:
            self._check_alt35_1_jcec(a5ss_mats_jcec_rows[1])
            self._check_alt35_2_jcec(a5ss_mats_jcec_rows[0])

        mxe_mats_jc_path = os.path.join(self._out_dir, 'MXE.MATS.JC.txt')
        mxe_mats_jc_header, mxe_mats_jc_rows, error = output_parser.parse_mats_jc(
            mxe_mats_jc_path)
        self.assertFalse(error)
        self._check_mxe_mats_jc_header(mxe_mats_jc_header)
        self.assertEqual(len(mxe_mats_jc_rows), 2)
        mxe_mats_jc_row_0_start = mxe_mats_jc_rows[0]['upstreamES']
        mxe_mats_jc_row_1_start = mxe_mats_jc_rows[1]['upstreamES']
        self.assertIn(mxe_mats_jc_row_0_start, ['6000', '7000'])
        self.assertIn(mxe_mats_jc_row_1_start, ['6000', '7000'])
        self.assertNotEqual(mxe_mats_jc_row_0_start, mxe_mats_jc_row_1_start)
        if mxe_mats_jc_row_0_start == '6000':
            self._check_mxe_1_jc(mxe_mats_jc_rows[0])
            self._check_mxe_2_jc(mxe_mats_jc_rows[1])
        else:
            self._check_mxe_1_jc(mxe_mats_jc_rows[1])
            self._check_mxe_2_jc(mxe_mats_jc_rows[0])

        mxe_mats_jcec_path = os.path.join(self._out_dir, 'MXE.MATS.JCEC.txt')
        mxe_mats_jcec_header, mxe_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            mxe_mats_jcec_path)
        self.assertFalse(error)
        self._check_mxe_mats_jcec_header(mxe_mats_jcec_header)
        self.assertEqual(len(mxe_mats_jcec_rows), 2)
        mxe_mats_jcec_row_0_start = mxe_mats_jcec_rows[0]['upstreamES']
        mxe_mats_jcec_row_1_start = mxe_mats_jcec_rows[1]['upstreamES']
        self.assertIn(mxe_mats_jcec_row_0_start, ['6000', '7000'])
        self.assertIn(mxe_mats_jcec_row_1_start, ['6000', '7000'])
        self.assertNotEqual(mxe_mats_jcec_row_0_start,
                            mxe_mats_jcec_row_1_start)
        if mxe_mats_jcec_row_0_start == '6000':
            self._check_mxe_1_jcec(mxe_mats_jcec_rows[0])
            self._check_mxe_2_jcec(mxe_mats_jcec_rows[1])
        else:
            self._check_mxe_1_jcec(mxe_mats_jcec_rows[1])
            self._check_mxe_2_jcec(mxe_mats_jcec_rows[0])

        ri_mats_jc_path = os.path.join(self._out_dir, 'RI.MATS.JC.txt')
        ri_mats_jc_header, ri_mats_jc_rows, error = output_parser.parse_mats_jc(
            ri_mats_jc_path)
        self.assertFalse(error)
        self._check_ri_mats_jc_header(ri_mats_jc_header)
        self.assertEqual(len(ri_mats_jc_rows), 2)
        ri_mats_jc_row_0_start = ri_mats_jc_rows[0]['upstreamES']
        ri_mats_jc_row_1_start = ri_mats_jc_rows[1]['upstreamES']
        self.assertIn(ri_mats_jc_row_0_start, ['8000', '9000'])
        self.assertIn(ri_mats_jc_row_1_start, ['8000', '9000'])
        self.assertNotEqual(ri_mats_jc_row_0_start, ri_mats_jc_row_1_start)
        if ri_mats_jc_row_0_start == '8000':
            self._check_ri_1_jc(ri_mats_jc_rows[0])
            self._check_ri_2_jc(ri_mats_jc_rows[1])
        else:
            self._check_ri_1_jc(ri_mats_jc_rows[1])
            self._check_ri_2_jc(ri_mats_jc_rows[0])

        ri_mats_jcec_path = os.path.join(self._out_dir, 'RI.MATS.JCEC.txt')
        ri_mats_jcec_header, ri_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            ri_mats_jcec_path)
        self.assertFalse(error)
        self._check_ri_mats_jcec_header(ri_mats_jcec_header)
        self.assertEqual(len(ri_mats_jcec_rows), 2)
        ri_mats_jcec_row_0_start = ri_mats_jcec_rows[0]['upstreamES']
        ri_mats_jcec_row_1_start = ri_mats_jcec_rows[1]['upstreamES']
        self.assertIn(ri_mats_jcec_row_0_start, ['8000', '9000'])
        self.assertIn(ri_mats_jcec_row_1_start, ['8000', '9000'])
        self.assertNotEqual(ri_mats_jcec_row_0_start, ri_mats_jcec_row_1_start)
        if ri_mats_jcec_row_0_start == '8000':
            self._check_ri_1_jcec(ri_mats_jcec_rows[0])
            self._check_ri_2_jcec(ri_mats_jcec_rows[1])
        else:
            self._check_ri_1_jcec(ri_mats_jcec_rows[1])
            self._check_ri_2_jcec(ri_mats_jcec_rows[0])


if __name__ == '__main__':
    unittest.main(verbosity=2)
