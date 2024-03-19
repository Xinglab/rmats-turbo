import os.path
import shutil
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class Test(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir, 'stat_large_file')
        self._generated_input_dir = os.path.join(self._test_dir,
                                                 'generated_input')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._tmp_dir = os.path.join(self._test_dir, 'tmp')

        tests.util.recreate_dirs([
            self._generated_input_dir, self._out_dir, self._tmp_dir,
            self._command_output_dir()
        ])

        self._gene_id = 'gene_id'
        self._gene_symbol = 'gene_symbol'
        self._chr = 'chr1'
        self._strand = '+'
        self._ids = [0, 1]
        self._base_exon_coords = [300, 400, 100, 200, 500, 600]
        self._num_sample_1_bams = int(1e4)
        self._num_sample_2_bams = int(1e4)
        self._inc_1_by_id = {0: 30, 1: 50}
        self._skip_1_by_id = {0: 70, 1: 50}
        self._inc_2_by_id = {0: 70, 1: 50}
        self._skip_2_by_id = {0: 30, 1: 50}
        self._inc_form_len = 200
        self._skip_form_len = 100

        self._from_gtf_se_path = os.path.join(self._generated_input_dir,
                                              'fromGTF.SE.txt')
        self._jc_raw_input_se_path = os.path.join(self._generated_input_dir,
                                                  'JC.raw.input.SE.txt')
        self._create_from_gtf(self._from_gtf_se_path)
        self._create_jc_input(self._jc_raw_input_se_path)
        shutil.copy(self._from_gtf_se_path, self._out_dir)
        shutil.copy(self._jc_raw_input_se_path, self._out_dir)

    def test(self):
        self._run_test()

    def _command_output_dir(self):
        return os.path.join(self._test_dir, 'command_output')

    def _rmats_arguments(self):
        return [
            '--od', self._out_dir, '--tmp', self._tmp_dir, '--task', 'stat'
        ]

    def _exon_coord_strs_for_i(self, i):
        exon_coords = [(coord + (i * 1000))
                       for coord in self._base_exon_coords]
        return [str(x) for x in exon_coords]

    def _write_tsv_line(self, handle, columns):
        handle.write('{}\n'.format('\t'.join(columns)))

    def _create_from_gtf(self, gtf_path):
        headers = [
            'ID', 'GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base',
            'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES',
            'downstreamEE'
        ]
        with open(gtf_path, 'wt') as handle:
            self._write_tsv_line(handle, headers)
            for i in self._ids:
                exon_coord_strs = self._exon_coord_strs_for_i(i)
                columns = [
                    str(i), self._gene_id, self._gene_symbol, self._chr,
                    self._strand
                ]
                columns.extend(exon_coord_strs)
                self._write_tsv_line(handle, columns)

    def _repeat_and_comma_separate(self, value, count):
        return ','.join([str(value)] * count)

    def _create_jc_input(self, jc_path):
        headers = [
            'ID', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2',
            'SJC_SAMPLE_2', 'IncFormLen', 'SkipFormLen'
        ]
        with open(jc_path, 'wt') as handle:
            self._write_tsv_line(handle, headers)
            for i in self._ids:
                inc_1 = self._inc_1_by_id[i]
                skip_1 = self._skip_1_by_id[i]
                inc_2 = self._inc_2_by_id[i]
                skip_2 = self._skip_2_by_id[i]
                inc_1_str = self._repeat_and_comma_separate(
                    inc_1, self._num_sample_1_bams)
                skip_1_str = self._repeat_and_comma_separate(
                    skip_1, self._num_sample_1_bams)
                inc_2_str = self._repeat_and_comma_separate(
                    inc_2, self._num_sample_2_bams)
                skip_2_str = self._repeat_and_comma_separate(
                    skip_2, self._num_sample_2_bams)
                columns = [
                    str(i), inc_1_str, skip_1_str, inc_2_str, skip_2_str,
                    str(self._inc_form_len),
                    str(self._skip_form_len)
                ]
                self._write_tsv_line(handle, columns)

    def _inc_level(self, inc, skip, inc_len, skip_len):
        normed_inc = (inc / inc_len)
        normed_skip = (skip / skip_len)
        inc_level = normed_inc / (normed_inc + normed_skip)
        return round(inc_level, 3)

    def _check_results(self):
        self.assertEqual(self._rmats_return_code, 0)
        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        # Only testing SE JC. Expect error line for SE JCEC and
        # all MXE, A3SS, A5SS, RI
        self.assertEqual(len(err_lines), 9)
        for err_line in err_lines:
            self.assertIn('Unable to produce final output files for', err_line)

        se_mats_jc_path = os.path.join(self._out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 2)
        for row_i, row in enumerate(se_mats_jc_rows):
            self.assertEqual(row['ID'], str(self._ids[row_i]))
            self.assertEqual(row['GeneID'], self._gene_id)
            self.assertEqual(row['geneSymbol'], self._gene_symbol)
            self.assertEqual(row['chr'], self._chr)
            self.assertEqual(row['strand'], self._strand)
            exon_coord_strs = self._exon_coord_strs_for_i(row_i)
            self.assertEqual(row['exonStart_0base'], exon_coord_strs[0])
            self.assertEqual(row['exonEnd'], exon_coord_strs[1])
            self.assertEqual(row['upstreamES'], exon_coord_strs[2])
            self.assertEqual(row['upstreamEE'], exon_coord_strs[3])
            self.assertEqual(row['downstreamES'], exon_coord_strs[4])
            self.assertEqual(row['downstreamEE'], exon_coord_strs[5])
            inc_1 = self._inc_1_by_id[row_i]
            skip_1 = self._skip_1_by_id[row_i]
            inc_2 = self._inc_2_by_id[row_i]
            skip_2 = self._skip_2_by_id[row_i]
            inc_1_str = self._repeat_and_comma_separate(
                inc_1, self._num_sample_1_bams)
            skip_1_str = self._repeat_and_comma_separate(
                skip_1, self._num_sample_1_bams)
            inc_2_str = self._repeat_and_comma_separate(
                inc_2, self._num_sample_2_bams)
            skip_2_str = self._repeat_and_comma_separate(
                skip_2, self._num_sample_2_bams)
            self.assertEqual(row['IJC_SAMPLE_1'], inc_1_str)
            self.assertEqual(row['SJC_SAMPLE_1'], skip_1_str)
            self.assertEqual(row['IJC_SAMPLE_2'], inc_2_str)
            self.assertEqual(row['SJC_SAMPLE_2'], skip_2_str)
            self.assertEqual(row['IncFormLen'], str(self._inc_form_len))
            self.assertEqual(row['SkipFormLen'], str(self._skip_form_len))
            inc_level_1 = self._inc_level(inc_1, skip_1, self._inc_form_len,
                                          self._skip_form_len)
            inc_level_2 = self._inc_level(inc_2, skip_2, self._inc_form_len,
                                          self._skip_form_len)
            inc_level_diff = inc_level_1 - inc_level_2
            inc_level_diff = round(inc_level_diff, 3)
            inc_level_1_str = self._repeat_and_comma_separate(
                inc_level_1, self._num_sample_1_bams)
            inc_level_2_str = self._repeat_and_comma_separate(
                inc_level_2, self._num_sample_2_bams)
            self.assertEqual(row['IncLevel1'], inc_level_1_str)
            self.assertEqual(row['IncLevel2'], inc_level_2_str)
            self.assertEqual(row['IncLevelDifference'], str(inc_level_diff))

        row_0 = se_mats_jc_rows[0]
        self.assertEqual(row_0['PValue'], '0')
        self.assertEqual(row_0['FDR'], '0.0')

        row_1 = se_mats_jc_rows[1]
        self.assertEqual(row_1['PValue'], '1')
        self.assertEqual(row_1['FDR'], '1.0')


if __name__ == '__main__':
    unittest.main(verbosity=2)
