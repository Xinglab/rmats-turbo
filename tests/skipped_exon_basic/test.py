import os.path
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
        self._test_dir = os.path.join(self._test_base_dir,
                                      'skipped_exon_basic')
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
        self._lib_type = 'fr-unstranded'
        self._task = 'both'

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

    def test(self):
        self._run_test()

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
            '--libType',
            self._lib_type,
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

    def _check_results(self):
        self._check_no_error_results()

        command_stdout_file_name = self._get_stdout_file_name()
        with open(command_stdout_file_name, 'rt') as out_f_h:
            out_lines = out_f_h.readlines()

        tests.util.assert_some_line_has(self, out_lines,
                                        'Found 1 exon skipping event')

        from_gtf_se_path = os.path.join(self._out_dir, 'fromGTF.SE.txt')
        from_gtf_se_header, from_gtf_se_rows, error = output_parser.parse_from_gtf(
            from_gtf_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_se_rows), 1)
        from_gtf_se_row = from_gtf_se_rows[0]
        self.assertEqual(from_gtf_se_row['GeneID'],
                         tests.util.double_quote(tests.util.gene_id_str(1)))
        self.assertEqual(from_gtf_se_row['geneSymbol'],
                         tests.util.double_quote(tests.util.gene_name_str(1)))
        self.assertEqual(from_gtf_se_row['chr'], 'chr1')
        self.assertEqual(from_gtf_se_row['strand'], '+')
        self.assertEqual(from_gtf_se_row['exonStart_0base'], '200')
        self.assertEqual(from_gtf_se_row['exonEnd'], '300')
        self.assertEqual(from_gtf_se_row['upstreamES'], '0')
        self.assertEqual(from_gtf_se_row['upstreamEE'], '100')
        self.assertEqual(from_gtf_se_row['downstreamES'], '400')
        self.assertEqual(from_gtf_se_row['downstreamEE'], '500')

        jc_raw_se_path = os.path.join(self._out_dir, 'JC.raw.input.SE.txt')
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


if __name__ == '__main__':
    unittest.main(verbosity=2)
