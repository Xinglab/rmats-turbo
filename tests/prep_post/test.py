import glob
import os.path
import subprocess
import sys
import unittest

import tests.bam
import tests.base_test
import tests.gtf
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class Test(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir, 'prep_post')
        self._generated_input_dir = os.path.join(self._test_dir,
                                                 'generated_input')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._prep_1_tmp_dir = os.path.join(self._test_dir, 'tmp_prep_1')
        self._prep_2_tmp_dir = os.path.join(self._test_dir, 'tmp_prep_2')
        self._post_tmp_dir = os.path.join(self._test_dir, 'tmp_post')

        self._dup_input_bam_tmp_dir = os.path.join(self._test_dir,
                                                   'tmp_dup_input_bam')
        self._dup_prep_bam_tmp_dir = os.path.join(self._test_dir,
                                                  'tmp_dup_prep_bam')
        self._miss_input_bam_tmp_dir = os.path.join(self._test_dir,
                                                    'tmp_miss_input_bam')
        self._miss_prep_bam_tmp_dir = os.path.join(self._test_dir,
                                                   'tmp_miss_prep_bam')
        self._extra_empty_dot_rmats_tmp_dir = os.path.join(
            self._test_dir, 'tmp_extra_dot_rmats')
        self._extra_empty_dot_rmats_out_dir = os.path.join(
            self._test_dir, 'out_extra_dot_rmats')

        tests.util.recreate_dirs([
            self._generated_input_dir, self._out_dir, self._prep_1_tmp_dir,
            self._prep_2_tmp_dir, self._post_tmp_dir,
            self._dup_input_bam_tmp_dir, self._dup_prep_bam_tmp_dir,
            self._miss_input_bam_tmp_dir, self._miss_prep_bam_tmp_dir,
            self._extra_empty_dot_rmats_out_dir,
            self._extra_empty_dot_rmats_tmp_dir,
            self._command_output_dir()
        ])

        self._read_type = 'paired'
        self._read_length = 50

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
        self._gtf = self._create_gtf(self._gtf_path)
        self._sub_steps = [
            'prep_1',
            'inte_1_fail',
            'inte_1_pass',
            'prep_2',
            'inte_2_fail',
            'inte_2_pass',
            'post',
            'duplicate_input_bam',
            'duplicate_prep_bam',
            'missing_input_bam',
            'missing_prep_bam',
            'extra_empty_dot_rmats',
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
        arguments = [
            '--gtf',
            self._gtf_path,
            '-t',
            self._read_type,
            '--readLength',
            str(self._read_length),
        ]

        if self._sub_step == 'prep_1':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._prep_1_tmp_dir,
                '--b1',
                self._sample_1_bams_path,
                '--task',
                'prep',
            ])
        elif self._sub_step == 'inte_1_fail':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._post_tmp_dir,
                '--b1',
                self._sample_1_bams_path,
                '--b2',
                self._sample_2_bams_path,
                '--task',
                'inte',
            ])
        elif self._sub_step == 'inte_1_pass':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._post_tmp_dir,
                '--b1',
                self._sample_1_bams_path,
                '--task',
                'inte',
                '--statoff',
            ])
        elif self._sub_step == 'prep_2':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._prep_2_tmp_dir,
                '--b1',
                self._sample_2_bams_path,
                '--task',
                'prep',
            ])
        elif self._sub_step == 'inte_2_fail':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._post_tmp_dir,
                '--b1',
                self._sample_2_bams_path,
                '--task',
                'inte',
                '--statoff',
            ])
        elif self._sub_step == 'inte_2_pass':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._post_tmp_dir,
                '--b1',
                self._sample_1_bams_path,
                '--b2',
                self._sample_2_bams_path,
                '--task',
                'inte',
            ])
        elif self._sub_step == 'post':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._post_tmp_dir,
                '--b1',
                self._sample_1_bams_path,
                '--b2',
                self._sample_2_bams_path,
                '--task',
                'post',
            ])
        elif self._sub_step == 'duplicate_input_bam':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._dup_input_bam_tmp_dir,
                '--b1',
                self._dup_input_bam_path,
                '--task',
                'post',
                '--statoff',
            ])
        elif self._sub_step == 'duplicate_prep_bam':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._dup_prep_bam_tmp_dir,
                '--b1',
                self._dup_prep_bam_path,
                '--task',
                'post',
                '--statoff',
            ])
        elif self._sub_step == 'missing_input_bam':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._miss_input_bam_tmp_dir,
                '--b1',
                self._miss_input_bam_path,
                '--task',
                'post',
                '--statoff',
            ])
        elif self._sub_step == 'missing_prep_bam':
            arguments.extend([
                '--od',
                self._out_dir,
                '--tmp',
                self._miss_prep_bam_tmp_dir,
                '--b1',
                self._miss_prep_bam_path,
                '--task',
                'post',
                '--statoff',
            ])
        elif self._sub_step == 'extra_empty_dot_rmats':
            arguments.extend([
                '--od',
                self._extra_empty_dot_rmats_out_dir,
                '--tmp',
                self._extra_empty_dot_rmats_tmp_dir,
                '--b1',
                self._sample_1_bams_path,
                '--b2',
                self._sample_2_bams_path,
                '--task',
                'post',
            ])

        return arguments

    def _setup_sub_step(self):
        if self._sub_step == 'duplicate_input_bam':
            self._setup_dup_input_bam()
        elif self._sub_step == 'duplicate_prep_bam':
            self._setup_dup_prep_bam()
        elif self._sub_step == 'missing_input_bam':
            self._setup_miss_input_bam()
        elif self._sub_step == 'missing_prep_bam':
            self._setup_miss_prep_bam()
        elif self._sub_step == 'extra_empty_dot_rmats':
            self._setup_extra_empty_dot_rmats()

    def _setup_dup_input_bam(self):
        self._dup_input_bam_path = os.path.join(self._generated_input_dir,
                                                'dup_input.txt')
        bams = self._sample_1_bams + [self._sample_1_bams[0]]
        self._write_bams(bams, self._dup_input_bam_path)
        self._cp_with_prefix('prep_1', self._prep_1_tmp_dir,
                             self._dup_input_bam_tmp_dir)

    def _setup_dup_prep_bam(self):
        self._dup_prep_bam_path = os.path.join(self._generated_input_dir,
                                               'dup_prep.txt')
        bams = self._sample_1_bams
        self._write_bams(bams, self._dup_prep_bam_path)
        self._cp_with_prefix('prep_1', self._prep_1_tmp_dir,
                             self._dup_prep_bam_tmp_dir)
        self._cp_with_prefix('prep_1_again', self._prep_1_tmp_dir,
                             self._dup_prep_bam_tmp_dir)

    def _setup_miss_input_bam(self):
        self._miss_input_bam_path = os.path.join(self._generated_input_dir,
                                                 'miss_input.txt')
        bams = [self._sample_1_bams[0]]
        self._write_bams(bams, self._miss_input_bam_path)
        self._cp_with_prefix('prep_1', self._prep_1_tmp_dir,
                             self._miss_input_bam_tmp_dir)

    def _setup_miss_prep_bam(self):
        self._miss_prep_bam_path = os.path.join(self._generated_input_dir,
                                                'miss_prep.txt')
        bams = self._sample_1_bams + self._sample_2_bams
        self._write_bams(bams, self._miss_prep_bam_path)
        self._cp_with_prefix('prep_1', self._prep_1_tmp_dir,
                             self._miss_prep_bam_tmp_dir)

    def _setup_extra_empty_dot_rmats(self):
        self._cp_with_prefix('prep_1', self._prep_1_tmp_dir,
                             self._extra_empty_dot_rmats_tmp_dir)
        self._cp_with_prefix('prep_2', self._prep_2_tmp_dir,
                             self._extra_empty_dot_rmats_tmp_dir)
        extra_dot_rmats_f_name = os.path.join(
            self._extra_empty_dot_rmats_tmp_dir, 'extra_empty.rmats')
        with open(extra_dot_rmats_f_name, 'wt'):
            pass  # create empty file

    def _create_gtf(self, gtf_path):
        gtf = tests.gtf.GTF()
        gtf.path = gtf_path
        transcript_1 = tests.gtf.Transcript()
        transcript_1.chromosome = '1'
        transcript_1.strand = '+'
        transcript_1.gene_id = tests.util.gene_id_str(1)
        transcript_1.gene_name = tests.util.gene_name_str(1)
        transcript_1.transcript_id = tests.util.transcript_id_str(1)
        transcript_1.exons = [(1, 100), (201, 300), (401, 500)]
        gtf.transcripts = [transcript_1]
        error = gtf.write()
        self.assertFalse(error)
        return gtf

    def _create_sample_1_bams(self, sample_1_bams_path,
                              sample_1_replicate_template):
        rep_1_bam = tests.bam.BAM()
        rep_1_bam.path = sample_1_replicate_template.format(1)
        rep_2_bam = tests.bam.BAM()
        rep_2_bam.path = sample_1_replicate_template.format(2)
        sample_1_bams = [rep_1_bam, rep_2_bam]

        rep_1_read_1 = tests.bam.Read()
        rep_1_read_1.ref_seq_name = '1'  # chromosome
        rep_1_read_1.ref_seq_len = 1000  # chromosome length
        rep_1_read_1.template_name = tests.util.template_name_str([1, 1])
        rep_1_read_2 = tests.bam.Read()
        error = tests.bam.set_read_pair_from_intervals(rep_1_read_1,
                                                       rep_1_read_2,
                                                       [[76, 100], [201, 300]],
                                                       [[401, 475]],
                                                       self._read_length)
        self.assertFalse(error)
        rep_1_bam.reads = [rep_1_read_1, rep_1_read_2]

        rep_2_read_1 = tests.bam.Read()
        rep_2_read_1.ref_seq_name = '1'  # chromosome
        rep_2_read_1.ref_seq_len = 1000  # chromosome length
        rep_2_read_1.template_name = tests.util.template_name_str([1, 2])
        rep_2_read_2 = tests.bam.Read()
        error = tests.bam.set_read_pair_from_intervals(
            rep_2_read_1, rep_2_read_2, [[26, 100]], [[201, 300], [401, 425]],
            self._read_length)
        self.assertFalse(error)
        rep_2_bam.reads = [rep_2_read_1, rep_2_read_2]

        self._write_bams(sample_1_bams, sample_1_bams_path)
        return sample_1_bams

    def _create_sample_2_bams(self, sample_2_bams_path,
                              sample_2_replicate_template):
        rep_1_bam = tests.bam.BAM()
        rep_1_bam.path = sample_2_replicate_template.format(1)
        rep_2_bam = tests.bam.BAM()
        rep_2_bam.path = sample_2_replicate_template.format(2)
        sample_2_bams = [rep_1_bam, rep_2_bam]

        rep_1_read_1 = tests.bam.Read()
        rep_1_read_1.ref_seq_name = '1'  # chromosome
        rep_1_read_1.ref_seq_len = 1000  # chromosome length
        rep_1_read_1.template_name = tests.util.template_name_str([2, 1])
        rep_1_read_2 = tests.bam.Read()
        error = tests.bam.set_read_pair_from_intervals(rep_1_read_1,
                                                       rep_1_read_2,
                                                       [[76, 100], [401, 500]],
                                                       [[401, 475]],
                                                       self._read_length)
        self.assertFalse(error)
        rep_1_bam.reads = [rep_1_read_1, rep_1_read_2]

        rep_2_read_1 = tests.bam.Read()
        rep_2_read_1.ref_seq_name = '1'  # chromosome
        rep_2_read_1.ref_seq_len = 1000  # chromosome length
        rep_2_read_1.template_name = tests.util.template_name_str([2, 2])
        rep_2_read_2 = tests.bam.Read()
        error = tests.bam.set_read_pair_from_intervals(rep_2_read_1,
                                                       rep_2_read_2,
                                                       [[26, 100]],
                                                       [[1, 100], [401, 425]],
                                                       self._read_length)
        self.assertFalse(error)
        rep_2_bam.reads = [rep_2_read_1, rep_2_read_2]

        self._write_bams(sample_2_bams, sample_2_bams_path)
        return sample_2_bams

    def _cp_with_prefix(self, prefix, source_dir, dest_dir):
        source_paths = self._get_dot_rmats_paths(source_dir)
        command = [
            sys.executable, tests.test_config.CP_WITH_PREFIX, prefix, dest_dir
        ]
        command.extend(source_paths)
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)

    def _check_results(self):
        if self._sub_step == 'prep_1':
            self._check_results_prep_1()
        elif self._sub_step == 'inte_1_fail':
            self._check_results_inte_1_fail()
        elif self._sub_step == 'inte_1_pass':
            self._check_results_inte_1_pass()
        elif self._sub_step == 'prep_2':
            self._check_results_prep_2()
        elif self._sub_step == 'inte_2_fail':
            self._check_results_inte_2_fail()
        elif self._sub_step == 'inte_2_pass':
            self._check_results_inte_2_pass()
        elif self._sub_step == 'post':
            self._check_results_post()
        elif self._sub_step == 'duplicate_input_bam':
            self._check_results_dup_input_bam()
        elif self._sub_step == 'duplicate_prep_bam':
            self._check_results_dup_prep_bam()
        elif self._sub_step == 'missing_input_bam':
            self._check_results_miss_input_bam()
        elif self._sub_step == 'missing_prep_bam':
            self._check_results_miss_prep_bam()
        elif self._sub_step == 'extra_empty_dot_rmats':
            self._check_results_extra_empty_dot_rmats()
        else:
            self.fail('unexpected sub_step: {}'.format(self._sub_step))

    def _get_dot_rmats_paths(self, tmp_dir):
        dot_rmats_file_paths = glob.glob(os.path.join(tmp_dir, '*.rmats'))
        # filenames begin with a timestamp used for alphanumeric sort
        return sorted(dot_rmats_file_paths)

    def _check_results_prep_1(self):
        self._check_no_error_results()

        command_stdout_file_name = self._get_stdout_file_name()
        with open(command_stdout_file_name, 'rt') as out_f_h:
            out_lines = out_f_h.readlines()

        tests.util.assert_no_line_has(self, out_lines,
                                      'Processing count files')

        test_gene_id = tests.util.gene_id_str(1)
        quoted_test_gene_id = tests.util.double_quote(test_gene_id)
        dot_rmats_paths = self._get_dot_rmats_paths(self._prep_1_tmp_dir)
        self.assertEqual(len(dot_rmats_paths), 2)
        for dot_rmats_i in range(2):
            dot_rmats_contents, error = output_parser.parse_dot_rmats(
                dot_rmats_paths[dot_rmats_i])
            self.assertFalse(error)

            self.assertEqual(dot_rmats_contents['bams'],
                             [self._sample_1_bams[dot_rmats_i].path])
            self.assertEqual(dot_rmats_contents['read_length'],
                             self._read_length)

            novel_juncs = dot_rmats_contents['novel_juncs']
            self.assertEqual(novel_juncs, [dict()])

            exons = dot_rmats_contents['exons']
            if dot_rmats_i == 0:
                self.assertEqual(exons, [{
                    quoted_test_gene_id: [{
                        'start_box': [401, 499],
                        'end_box': [401, 499],
                        'count': 1
                    }]
                }])
            else:
                self.assertEqual(exons, [{
                    quoted_test_gene_id: [{
                        'start_box': [1, 99],
                        'end_box': [1, 99],
                        'count': 1
                    }]
                }])

            multis = dot_rmats_contents['multis']
            if dot_rmats_i == 0:
                self.assertEqual(multis, [{
                    quoted_test_gene_id: [{
                        'junction_pairs': [[1, 1], [100, 200], [299, 299]],
                        'count':
                        1
                    }]
                }])
            else:
                self.assertEqual(multis, [{
                    quoted_test_gene_id: [{
                        'junction_pairs': [[201, 201], [300, 400], [499, 499]],
                        'count':
                        1
                    }]
                }])

        self._cp_with_prefix('prep_1_', self._prep_1_tmp_dir,
                             self._post_tmp_dir)

    def _check_results_prep_2(self):
        self._check_no_error_results()

        command_stdout_file_name = self._get_stdout_file_name()
        with open(command_stdout_file_name, 'rt') as out_f_h:
            out_lines = out_f_h.readlines()

        tests.util.assert_no_line_has(self, out_lines,
                                      'Processing count files')

        test_gene_id = tests.util.gene_id_str(1)
        quoted_test_gene_id = tests.util.double_quote(test_gene_id)
        dot_rmats_paths = self._get_dot_rmats_paths(self._prep_2_tmp_dir)
        self.assertEqual(len(dot_rmats_paths), 2)
        for dot_rmats_i in range(2):
            dot_rmats_contents, error = output_parser.parse_dot_rmats(
                dot_rmats_paths[dot_rmats_i])
            self.assertFalse(error)

            self.assertEqual(dot_rmats_contents['bams'],
                             [self._sample_2_bams[dot_rmats_i].path])
            self.assertEqual(dot_rmats_contents['read_length'],
                             self._read_length)

            novel_juncs = dot_rmats_contents['novel_juncs']
            self.assertEqual(novel_juncs, [{quoted_test_gene_id: [[0, 0, 2]]}])

            exons = dot_rmats_contents['exons']
            if dot_rmats_i == 0:
                self.assertEqual(exons, [{
                    quoted_test_gene_id: [{
                        'start_box': [401, 499],
                        'end_box': [401, 499],
                        'count': 1
                    }]
                }])
            else:
                self.assertEqual(exons, [{
                    quoted_test_gene_id: [{
                        'start_box': [1, 99],
                        'end_box': [1, 99],
                        'count': 1
                    }]
                }])

        multis = dot_rmats_contents['multis']
        if dot_rmats_i == 0:
            self.assertEqual(multis, [{
                quoted_test_gene_id: [{
                    'junction_pairs': [[1, 1], [100, 400], [499, 499]],
                    'count':
                    1
                }]
            }])
        else:
            self.assertEqual(multis, [{
                quoted_test_gene_id: [{
                    'junction_pairs': [[1, 1], [100, 400], [499, 499]],
                    'count':
                    1
                }]
            }])

        self._cp_with_prefix('prep_2_', self._prep_2_tmp_dir,
                             self._post_tmp_dir)

    def _check_results_inte_1_fail(self):
        self.assertNotEqual(self._rmats_return_code, 0)

        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        tests.util.assert_some_line_has(
            self, err_lines, 'input bam files with no associated prep output')

    def _check_results_inte_1_pass(self):
        self._check_no_error_results()

    def _check_results_inte_2_fail(self):
        self.assertNotEqual(self._rmats_return_code, 0)

        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        tests.util.assert_some_line_has(
            self, err_lines,
            'bam files not in input but associated with prep output')

    def _check_results_inte_2_pass(self):
        self._check_no_error_results()

    def _check_results_post(self):
        self._check_no_error_results()
        self._check_results_post_shared(self._out_dir)

    def _check_results_post_shared(self, out_dir):
        command_stdout_file_name = self._get_stdout_file_name()
        with open(command_stdout_file_name, 'rt') as out_f_h:
            out_lines = out_f_h.readlines()

        tests.util.assert_some_line_has(self, out_lines,
                                        'Processing count files')

        from_gtf_se_path = os.path.join(out_dir, 'fromGTF.SE.txt')
        from_gtf_se_header, from_gtf_se_rows, error = output_parser.parse_from_gtf(
            from_gtf_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_se_rows), 1)
        from_gtf_se_row = from_gtf_se_rows[0]
        self.assertEqual(from_gtf_se_row['GeneID'],
                         tests.util.double_quote(tests.util.gene_id_str(1)))
        self.assertEqual(from_gtf_se_row['exonStart_0base'], '200')
        self.assertEqual(from_gtf_se_row['exonEnd'], '300')

        jc_raw_se_path = os.path.join(out_dir, 'JC.raw.input.SE.txt')
        jc_raw_se_header, jc_raw_se_rows, error = output_parser.parse_jc_raw(
            jc_raw_se_path)
        self.assertFalse(error)
        self.assertEqual(len(jc_raw_se_rows), 1)
        jc_raw_se_row = jc_raw_se_rows[0]
        self.assertEqual(jc_raw_se_row['ID'], from_gtf_se_row['ID'])
        self.assertEqual(jc_raw_se_row['IJC_SAMPLE_1'], '1,1')
        self.assertEqual(jc_raw_se_row['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(jc_raw_se_row['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(jc_raw_se_row['SJC_SAMPLE_2'], '1,1')

        se_mats_jc_path = os.path.join(out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 1)
        se_mats_jc_row = se_mats_jc_rows[0]
        pvalue = float(se_mats_jc_row['PValue'])
        tests.util.assert_within_bounds(self, pvalue, 0, 1)
        fdr = float(se_mats_jc_row['FDR'])
        tests.util.assert_within_bounds(self, fdr, 0, 1)
        inc_level_1_splits = se_mats_jc_row['IncLevel1'].split(',')
        self.assertEqual(len(inc_level_1_splits), 2)
        self.assertAlmostEqual(float(inc_level_1_splits[0]), 1)
        self.assertAlmostEqual(float(inc_level_1_splits[1]), 1)
        inc_level_2_splits = se_mats_jc_row['IncLevel2'].split(',')
        self.assertEqual(len(inc_level_2_splits), 2)
        self.assertAlmostEqual(float(inc_level_2_splits[0]), 0)
        self.assertAlmostEqual(float(inc_level_2_splits[1]), 0)
        self.assertAlmostEqual(float(se_mats_jc_row['IncLevelDifference']), 1)

    def _check_results_dup_input_bam(self):
        self.assertNotEqual(self._rmats_return_code, 0)

        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        dup_bam_path = self._sample_1_bams[0].path
        expected_error = '{} given 2 times'.format(dup_bam_path)
        tests.util.assert_some_line_has(self, err_lines, expected_error)

    def _check_results_dup_prep_bam(self):
        self.assertNotEqual(self._rmats_return_code, 0)

        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        for bam in self._sample_1_bams:
            dup_bam_path = bam.path
            expected_error = '{} found 2 times in .rmats'.format(dup_bam_path)
            tests.util.assert_some_line_has(self, err_lines, expected_error)

    def _check_results_miss_input_bam(self):
        self._check_no_error_results()

    def _check_results_miss_prep_bam(self):
        self.assertNotEqual(self._rmats_return_code, 0)

        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        for bam in self._sample_2_bams:
            miss_bam_path = bam.path
            expected_error = '{} not found in .rmats'.format(miss_bam_path)
            tests.util.assert_some_line_has(self, err_lines, expected_error)

    def _check_results_extra_empty_dot_rmats(self):
        self.assertEqual(self._rmats_return_code, 0)

        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        self.assertEqual(len(err_lines), 1)
        self.assertIn(
            'A .rmats file was found with no bams listed in it.'
            ' Ignoring that file', err_lines[0])

        # The extra empty file should not change the output
        self._check_results_post_shared(self._extra_empty_dot_rmats_out_dir)


if __name__ == '__main__':
    unittest.main(verbosity=2)
