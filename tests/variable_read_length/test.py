import os.path
import unittest

import tests.bam
import tests.base_test
import tests.gtf
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class VariableReadLengthBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir,
                                      'variable_read_length',
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
        self._task = 'both'
        self._read_length_1 = 50
        self._read_length_2 = 60

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
            '--tmp',
            self._tmp_dir,
            '--task',
            self._task,
        ]

    def _create_gtf(self, gtf_path):
        gtf = tests.gtf.GTF()
        gtf.path = gtf_path
        transcript_1 = tests.gtf.Transcript()
        transcript_1.chromosome = '1'
        transcript_1.strand = '+'
        transcript_1.gene_id = tests.util.gene_id_str(1)
        transcript_1.gene_name = tests.util.gene_name_str(1)
        transcript_1.transcript_id = tests.util.transcript_id_str(1)
        transcript_1.exons = [(1, 100), (201, 300), (401, 500), (601, 700)]
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
                                                       self._read_length_1)
        self.assertFalse(error)
        rep_1_bam.reads = [rep_1_read_1, rep_1_read_2]

        rep_2_read_1 = tests.bam.Read()
        rep_2_read_1.ref_seq_name = '1'  # chromosome
        rep_2_read_1.ref_seq_len = 1000  # chromosome length
        rep_2_read_1.template_name = tests.util.template_name_str([1, 2])
        rep_2_read_2 = tests.bam.Read()
        error = tests.bam.set_read_pair_from_intervals(
            rep_2_read_1, rep_2_read_2, [[226, 300]], [[401, 500], [601, 625]],
            self._read_length_2)
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
                                                       self._read_length_1)
        self.assertFalse(error)
        rep_1_bam.reads = [rep_1_read_1, rep_1_read_2]

        rep_2_read_1 = tests.bam.Read()
        rep_2_read_1.ref_seq_name = '1'  # chromosome
        rep_2_read_1.ref_seq_len = 1000  # chromosome length
        rep_2_read_1.template_name = tests.util.template_name_str([2, 2])
        rep_2_read_2 = tests.bam.Read()
        error = tests.bam.set_read_pair_from_intervals(
            rep_2_read_1, rep_2_read_2, [[226, 300]], [[201, 300], [601, 625]],
            self._read_length_2)
        self.assertFalse(error)
        rep_2_bam.reads = [rep_2_read_1, rep_2_read_2]

        self._write_bams(sample_2_bams, sample_2_bams_path)
        return sample_2_bams

    def _get_sorted_from_gtf_se_rows(self):
        from_gtf_se_path = os.path.join(self._out_dir, 'fromGTF.SE.txt')
        from_gtf_se_header, from_gtf_se_rows, error = output_parser.parse_from_gtf(
            from_gtf_se_path)
        self.assertFalse(error)
        return sorted(from_gtf_se_rows, key=lambda r: r['exonEnd'])

    def _get_matched_se_mats_jc_rows(self, gtf_rows):
        matched_rows = [None for _ in range(len(gtf_rows))]

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self.assertEqual(len(se_mats_jc_rows), len(gtf_rows))
        gtf_ids = [r['ID'] for r in gtf_rows]
        for mats_row in se_mats_jc_rows:
            idx = gtf_ids.index(mats_row['ID'])
            matched_rows[idx] = mats_row

        self.assertTrue(all(matched_rows))
        return matched_rows

    def _get_matched_se_mats_jcec_rows(self, gtf_rows):
        matched_rows = [None for _ in range(len(gtf_rows))]

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            se_mats_jcec_path)
        self.assertFalse(error)
        self.assertEqual(len(se_mats_jcec_rows), len(gtf_rows))
        gtf_ids = [r['ID'] for r in gtf_rows]
        for mats_row in se_mats_jcec_rows:
            idx = gtf_ids.index(mats_row['ID'])
            matched_rows[idx] = mats_row

        self.assertTrue(all(matched_rows))
        return matched_rows


class Length1Test(VariableReadLengthBaseTest):
    def _sub_test_name(self):
        return 'length_1'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend([
            '--readLength',
            str(self._read_length_1),
        ])
        return arguments

    def _check_results(self):
        self._check_no_error_results()

        sorted_from_gtf_se_rows = self._get_sorted_from_gtf_se_rows()
        self.assertEqual(len(sorted_from_gtf_se_rows), 1)

        from_gtf_se_row_1 = sorted_from_gtf_se_rows[0]
        self.assertEqual(from_gtf_se_row_1['exonStart_0base'], '200')
        self.assertEqual(from_gtf_se_row_1['exonEnd'], '300')

        matched_se_mats_jc_rows = self._get_matched_se_mats_jc_rows(
            sorted_from_gtf_se_rows)

        se_mats_jc_row_1 = matched_se_mats_jc_rows[0]
        self.assertEqual(se_mats_jc_row_1['ID'], from_gtf_se_row_1['ID'])
        self.assertEqual(se_mats_jc_row_1['IJC_SAMPLE_1'], '1,0')
        self.assertEqual(se_mats_jc_row_1['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jc_row_1['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(se_mats_jc_row_1['SJC_SAMPLE_2'], '1,0')
        self.assertEqual(se_mats_jc_row_1['IncFormLen'], '98')
        self.assertEqual(se_mats_jc_row_1['SkipFormLen'], '49')

        matched_se_mats_jcec_rows = self._get_matched_se_mats_jcec_rows(
            sorted_from_gtf_se_rows)

        se_mats_jcec_row_1 = matched_se_mats_jcec_rows[0]
        self.assertEqual(se_mats_jcec_row_1['ID'], from_gtf_se_row_1['ID'])
        self.assertEqual(se_mats_jcec_row_1['IJC_SAMPLE_1'], '1,0')
        self.assertEqual(se_mats_jcec_row_1['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jcec_row_1['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(se_mats_jcec_row_1['SJC_SAMPLE_2'], '1,0')
        self.assertEqual(se_mats_jcec_row_1['IncFormLen'], '149')
        self.assertEqual(se_mats_jcec_row_1['SkipFormLen'], '49')


class Length1VariableTest(VariableReadLengthBaseTest):
    def _sub_test_name(self):
        return 'length_1_variable'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend([
            '--readLength',
            str(self._read_length_1),
            '--variable-read-length',
        ])
        return arguments

    def _check_results(self):
        self._check_no_error_results()

        sorted_from_gtf_se_rows = self._get_sorted_from_gtf_se_rows()
        self.assertEqual(len(sorted_from_gtf_se_rows), 2)

        from_gtf_se_row_1 = sorted_from_gtf_se_rows[0]
        self.assertEqual(from_gtf_se_row_1['exonStart_0base'], '200')
        self.assertEqual(from_gtf_se_row_1['exonEnd'], '300')

        from_gtf_se_row_2 = sorted_from_gtf_se_rows[1]
        self.assertEqual(from_gtf_se_row_2['exonStart_0base'], '400')
        self.assertEqual(from_gtf_se_row_2['exonEnd'], '500')

        matched_se_mats_jc_rows = self._get_matched_se_mats_jc_rows(
            sorted_from_gtf_se_rows)

        se_mats_jc_row_1 = matched_se_mats_jc_rows[0]
        self.assertEqual(se_mats_jc_row_1['ID'], from_gtf_se_row_1['ID'])
        self.assertEqual(se_mats_jc_row_1['IJC_SAMPLE_1'], '1,0')
        self.assertEqual(se_mats_jc_row_1['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jc_row_1['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(se_mats_jc_row_1['SJC_SAMPLE_2'], '1,0')
        self.assertEqual(se_mats_jc_row_1['IncFormLen'], '98')
        self.assertEqual(se_mats_jc_row_1['SkipFormLen'], '49')

        se_mats_jc_row_2 = matched_se_mats_jc_rows[1]
        self.assertEqual(se_mats_jc_row_2['ID'], from_gtf_se_row_2['ID'])
        self.assertEqual(se_mats_jc_row_2['IJC_SAMPLE_1'], '0,1')
        self.assertEqual(se_mats_jc_row_2['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jc_row_2['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(se_mats_jc_row_2['SJC_SAMPLE_2'], '0,1')
        self.assertEqual(se_mats_jc_row_2['IncFormLen'], '98')
        self.assertEqual(se_mats_jc_row_2['SkipFormLen'], '49')

        matched_se_mats_jcec_rows = self._get_matched_se_mats_jcec_rows(
            sorted_from_gtf_se_rows)

        se_mats_jcec_row_1 = matched_se_mats_jcec_rows[0]
        self.assertEqual(se_mats_jcec_row_1['ID'], from_gtf_se_row_1['ID'])
        self.assertEqual(se_mats_jcec_row_1['IJC_SAMPLE_1'], '1,1')
        self.assertEqual(se_mats_jcec_row_1['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jcec_row_1['IJC_SAMPLE_2'], '0,1')
        self.assertEqual(se_mats_jcec_row_1['SJC_SAMPLE_2'], '1,0')
        self.assertEqual(se_mats_jcec_row_1['IncFormLen'], '149')
        self.assertEqual(se_mats_jcec_row_1['SkipFormLen'], '49')

        se_mats_jcec_row_2 = matched_se_mats_jcec_rows[1]
        self.assertEqual(se_mats_jcec_row_2['ID'], from_gtf_se_row_2['ID'])
        self.assertEqual(se_mats_jcec_row_2['IJC_SAMPLE_1'], '1,1')
        self.assertEqual(se_mats_jcec_row_2['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jcec_row_2['IJC_SAMPLE_2'], '1,0')
        self.assertEqual(se_mats_jcec_row_2['SJC_SAMPLE_2'], '0,1')
        self.assertEqual(se_mats_jcec_row_2['IncFormLen'], '149')
        self.assertEqual(se_mats_jcec_row_2['SkipFormLen'], '49')


class Length2Test(VariableReadLengthBaseTest):
    def _sub_test_name(self):
        return 'length_2'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend([
            '--readLength',
            str(self._read_length_2),
        ])
        return arguments

    def _check_results(self):
        self._check_no_error_results()

        sorted_from_gtf_se_rows = self._get_sorted_from_gtf_se_rows()
        self.assertEqual(len(sorted_from_gtf_se_rows), 1)

        from_gtf_se_row_1 = sorted_from_gtf_se_rows[0]
        self.assertEqual(from_gtf_se_row_1['exonStart_0base'], '400')
        self.assertEqual(from_gtf_se_row_1['exonEnd'], '500')

        matched_se_mats_jc_rows = self._get_matched_se_mats_jc_rows(
            sorted_from_gtf_se_rows)

        se_mats_jc_row_1 = matched_se_mats_jc_rows[0]
        self.assertEqual(se_mats_jc_row_1['ID'], from_gtf_se_row_1['ID'])
        self.assertEqual(se_mats_jc_row_1['IJC_SAMPLE_1'], '0,1')
        self.assertEqual(se_mats_jc_row_1['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jc_row_1['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(se_mats_jc_row_1['SJC_SAMPLE_2'], '0,1')
        self.assertEqual(se_mats_jc_row_1['IncFormLen'], '118')
        self.assertEqual(se_mats_jc_row_1['SkipFormLen'], '59')

        matched_se_mats_jcec_rows = self._get_matched_se_mats_jcec_rows(
            sorted_from_gtf_se_rows)

        se_mats_jcec_row_1 = matched_se_mats_jcec_rows[0]
        self.assertEqual(se_mats_jcec_row_1['ID'], from_gtf_se_row_1['ID'])
        self.assertEqual(se_mats_jcec_row_1['IJC_SAMPLE_1'], '0,1')
        self.assertEqual(se_mats_jcec_row_1['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jcec_row_1['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(se_mats_jcec_row_1['SJC_SAMPLE_2'], '0,1')
        self.assertEqual(se_mats_jcec_row_1['IncFormLen'], '159')
        self.assertEqual(se_mats_jcec_row_1['SkipFormLen'], '59')


class Length2VariableTest(VariableReadLengthBaseTest):
    def _sub_test_name(self):
        return 'length_2_variable'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.extend([
            '--readLength',
            str(self._read_length_2),
            '--variable-read-length',
        ])
        return arguments

    def _check_results(self):
        self._check_no_error_results()

        sorted_from_gtf_se_rows = self._get_sorted_from_gtf_se_rows()
        self.assertEqual(len(sorted_from_gtf_se_rows), 2)

        from_gtf_se_row_1 = sorted_from_gtf_se_rows[0]
        self.assertEqual(from_gtf_se_row_1['exonStart_0base'], '200')
        self.assertEqual(from_gtf_se_row_1['exonEnd'], '300')

        from_gtf_se_row_2 = sorted_from_gtf_se_rows[1]
        self.assertEqual(from_gtf_se_row_2['exonStart_0base'], '400')
        self.assertEqual(from_gtf_se_row_2['exonEnd'], '500')

        matched_se_mats_jc_rows = self._get_matched_se_mats_jc_rows(
            sorted_from_gtf_se_rows)

        se_mats_jc_row_1 = matched_se_mats_jc_rows[0]
        self.assertEqual(se_mats_jc_row_1['ID'], from_gtf_se_row_1['ID'])
        self.assertEqual(se_mats_jc_row_1['IJC_SAMPLE_1'], '1,0')
        self.assertEqual(se_mats_jc_row_1['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jc_row_1['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(se_mats_jc_row_1['SJC_SAMPLE_2'], '1,0')
        self.assertEqual(se_mats_jc_row_1['IncFormLen'], '118')
        self.assertEqual(se_mats_jc_row_1['SkipFormLen'], '59')

        se_mats_jc_row_2 = matched_se_mats_jc_rows[1]
        self.assertEqual(se_mats_jc_row_2['ID'], from_gtf_se_row_2['ID'])
        self.assertEqual(se_mats_jc_row_2['IJC_SAMPLE_1'], '0,1')
        self.assertEqual(se_mats_jc_row_2['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jc_row_2['IJC_SAMPLE_2'], '0,0')
        self.assertEqual(se_mats_jc_row_2['SJC_SAMPLE_2'], '0,1')
        self.assertEqual(se_mats_jc_row_1['IncFormLen'], '118')
        self.assertEqual(se_mats_jc_row_1['SkipFormLen'], '59')

        matched_se_mats_jcec_rows = self._get_matched_se_mats_jcec_rows(
            sorted_from_gtf_se_rows)

        se_mats_jcec_row_1 = matched_se_mats_jcec_rows[0]
        self.assertEqual(se_mats_jcec_row_1['ID'], from_gtf_se_row_1['ID'])
        self.assertEqual(se_mats_jcec_row_1['IJC_SAMPLE_1'], '1,1')
        self.assertEqual(se_mats_jcec_row_1['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jcec_row_1['IJC_SAMPLE_2'], '0,1')
        self.assertEqual(se_mats_jcec_row_1['SJC_SAMPLE_2'], '1,0')
        self.assertEqual(se_mats_jcec_row_1['IncFormLen'], '159')
        self.assertEqual(se_mats_jcec_row_1['SkipFormLen'], '59')

        se_mats_jcec_row_2 = matched_se_mats_jcec_rows[1]
        self.assertEqual(se_mats_jcec_row_2['ID'], from_gtf_se_row_2['ID'])
        self.assertEqual(se_mats_jcec_row_2['IJC_SAMPLE_1'], '1,1')
        self.assertEqual(se_mats_jcec_row_2['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(se_mats_jcec_row_2['IJC_SAMPLE_2'], '1,0')
        self.assertEqual(se_mats_jcec_row_2['SJC_SAMPLE_2'], '0,1')
        self.assertEqual(se_mats_jcec_row_2['IncFormLen'], '159')
        self.assertEqual(se_mats_jcec_row_2['SkipFormLen'], '59')


class NoLengthTest(VariableReadLengthBaseTest):
    def _sub_test_name(self):
        return 'no_length'

    def test(self):
        self._run_test()

    def _check_results(self):
        self.assertNotEqual(self._rmats_return_code, 0)

        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        tests.util.assert_some_line_has(self, err_lines,
                                        '--readLength is required')
        tests.util.assert_some_line_has(self, err_lines,
                                        '--variable-read-length')


if __name__ == '__main__':
    unittest.main(verbosity=2)
