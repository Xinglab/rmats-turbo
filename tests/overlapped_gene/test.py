import os.path
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class Test(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir, 'overlapped_gene')
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

        self._sample_1_bams_path = os.path.join(self._generated_input_dir,
                                                'b1.txt')
        sample_1_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_1_rep_{}.bam')
        self._sample_1_bams = self._create_sample_1_bams(
            self._sample_1_bams_path, sample_1_bam_replicate_template)

        self._gtf_path = os.path.join(self._generated_input_dir, 'test.gtf')
        transcript_exons, transcript_genes = (
            self._exons_and_genes_by_transcript())
        self._gtf = self._create_gtf_from_transcripts(self._gtf_path,
                                                      transcript_exons,
                                                      genes=transcript_genes)

    def test(self):
        self._run_test()

    def _command_output_dir(self):
        return os.path.join(self._test_dir, 'command_output')

    def _rmats_arguments(self):
        return [
            '--gtf', self._gtf_path, '-t', self._read_type, '--readLength',
            str(self._read_length), '--od', self._out_dir, '--tmp',
            self._tmp_dir, '--b1', self._sample_1_bams_path, '--task', 'both',
            '--statoff'
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

    def _exons_and_genes_by_transcript(self):
        exons_and_genes = [
            ([(1, 100), (201, 300), (401, 500)], [1, 2]),  # SE 1 inc
            ([(1, 100), (401, 500)], [1, 2]),  # SE 1 skip
            ([(1001, 1100), (1201, 1300), (1401, 1500)], [1, 2]),  # SE 2 inc
            ([(1001, 1100), (1401, 1500)], [1]),  # SE 2 skip
        ]
        exons_by_transcript = list()
        genes_by_transcript = list()
        for exons, genes in exons_and_genes:
            for gene in genes:
                exons_by_transcript.append(exons)
                genes_by_transcript.append(gene)

        return exons_by_transcript, genes_by_transcript

    def _include_read_SE_1(self):
        return ([[81, 100], [201, 300]], [[201, 300]])

    def _skip_read_SE_1(self):
        return ([[81, 100], [401, 500]], [[401, 500]])

    def _include_read_SE_2(self):
        return ([[1081, 1100], [1201, 1300]], [[1201, 1300]])

    def _skip_read_SE_2(self):
        return ([[1081, 1100], [1401, 1500]], [[1401, 1500]])

    def _paired_read_coords_1_1(self):
        return [self._include_read_SE_1(), self._skip_read_SE_1()]

    def _paired_read_coords_1_2(self):
        return [self._include_read_SE_2(), self._skip_read_SE_2()]

    def _check_results(self):
        self._check_no_error_results()

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 2)
        expected_genes = [
            tests.util.double_quote(tests.util.gene_id_str(i)) for i in [1, 2]
        ]
        for row in se_mats_jc_rows:
            self.assertIn(row['exonStart_0base'], ['200', '1200'])
            if row['exonStart_0base'] == '200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,0')
                self.assertIn(row['GeneID'], expected_genes)
            elif row['exonStart_0base'] == '1200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,1')
                self.assertEqual(row['GeneID'], expected_genes[0])

        from_gtf_se_path = os.path.join(self._out_dir, 'fromGTF.SE.txt')
        from_gtf_se_header, from_gtf_se_rows, error = (
            output_parser.parse_from_gtf(from_gtf_se_path))
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_se_rows), 2)
        gtf_coords = list()
        for row in from_gtf_se_rows:
            gtf_coords.append(row['exonStart_0base'])

        self.assertEqual(sorted(gtf_coords), ['1200', '200'])

        from_gtf_nj_se_path = os.path.join(self._out_dir,
                                           'fromGTF.novelJunction.SE.txt')
        from_gtf_nj_se_header, from_gtf_nj_se_rows, error = (
            output_parser.parse_from_gtf_novel_junction(from_gtf_nj_se_path))
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_nj_se_rows), 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
