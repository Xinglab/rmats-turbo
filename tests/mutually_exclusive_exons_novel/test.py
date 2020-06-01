import os.path
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class MutuallyExclusiveExonsNovelBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir,
                                      'mutually_exclusive_exons_novel',
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
        self._chromosome_length = 8000
        self._task = 'both'

        self._sample_1_bams_path = os.path.join(self._generated_input_dir,
                                                'b1.txt')
        sample_1_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_1_rep_{}.bam')
        self._sample_1_bams = self._create_sample_1_bams(
            self._sample_1_bams_path, sample_1_bam_replicate_template)
        self._gtf_path = os.path.join(self._generated_input_dir, 'test.gtf')
        self._gtf = self._create_gtf_from_transcripts(
            self._gtf_path, self._exons_by_transcript())

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

    def _create_sample_1_bams(self, sample_1_bams_path,
                              sample_1_replicate_template):
        rep_1_bam_path = sample_1_replicate_template.format(1)
        rep_1_bam = self._create_bam_from_paired_read_coords(
            rep_1_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords())
        sample_1_bams = [rep_1_bam]
        self._write_bams(sample_1_bams, sample_1_bams_path)
        return sample_1_bams

    def _sort_by_event_coords(self, rows):
        return sorted(rows,
                      key=lambda r:
                      (int(r['upstreamES']), int(r['upstreamEE']),
                       int(r['1stExonStart_0base']), int(r['1stExonEnd']),
                       int(r['2ndExonStart_0base']), int(r['2ndExonEnd']),
                       int(r['downstreamES']), int(r['downstreamEE'])))

    def _check_event_coords(self, row, u_start, u_end, f_start, f_end, s_start,
                            s_end, d_start, d_end):
        self.assertEqual(row['upstreamES'], u_start)
        self.assertEqual(row['upstreamEE'], u_end)
        self.assertEqual(row['1stExonStart_0base'], f_start)
        self.assertEqual(row['1stExonEnd'], f_end)
        self.assertEqual(row['2ndExonStart_0base'], s_start)
        self.assertEqual(row['2ndExonEnd'], s_end)
        self.assertEqual(row['downstreamES'], d_start)
        self.assertEqual(row['downstreamEE'], d_end)


class NovelJunction(MutuallyExclusiveExonsNovelBaseTest):
    def _sub_test_name(self):
        return 'novel_junction'

    def test(self):
        self._run_test()

    def _exons_by_transcript(self):
        return [
            [(1, 100), (201, 300), (601, 700)],
            [(1, 100), (401, 500), (601, 700)],
            [(1001, 1100), (1401, 1500)],
            [(801, 900), (1201, 1300), (1401, 1500)],
            [(1601, 1700), (1801, 1900), (2201, 2300)],
            [(1601, 1700), (1801, 1900), (2001, 2100), (2201, 2300)],
            [(2401, 2500), (2601, 2700)],
            [(2401, 2500), (2801, 2900), (3001, 3100)],
            [(3201, 3300), (3401, 3500), (4001, 4100)],
            [(3201, 3300), (3601, 3700), (3801, 3900), (4001, 4100)],
            [(4201, 4300), (4401, 4500), (4601, 4700), (5001, 5100)],
            [(4201, 4300), (4801, 4900), (5001, 5100)],
            [(5201, 5300), (5401, 5500), (5601, 5700), (5801, 5900)],
            [(5201, 5300), (5601, 5700), (5801, 5900)],
            [(6001, 6100), (6201, 6300), (6601, 6700)],
            [(6401, 6500), (6601, 6700)],
            [(6801, 6900), (7001, 7100), (7401, 7500)],
            [(6801, 6900), (7201, 7300)],
        ]

    def _paired_read_coords(self):
        return [
            ([[881, 900], [1001, 1100]], [[1001, 1100]]),
            ([[1681, 1700], [2001, 2100]], [[2001, 2100]]),
            ([[2681, 2700], [3001, 3100]], [[3001, 3100]]),
            ([[3681, 3700], [4001, 4100]], [[4001, 4100]]),
            ([[4281, 4300], [4601, 4700]], [[4601, 4700]]),
            ([[5481, 5500], [5801, 5900]], [[5801, 5900]]),
            ([[6081, 6100], [6401, 6500]], [[6401, 6500]]),
            ([[7281, 7300], [7401, 7500]], [[7401, 7500]]),
        ]

    def _check_results(self):
        self._check_no_error_results()

        from_gtf_mxe_path = os.path.join(self._out_dir, 'fromGTF.MXE.txt')
        from_gtf_mxe_header, from_gtf_mxe_rows, error = output_parser.parse_from_gtf(
            from_gtf_mxe_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_mxe_rows), 5)
        sorted_from_gtf_mxe_rows = self._sort_by_event_coords(
            from_gtf_mxe_rows)
        self._check_event_coords(sorted_from_gtf_mxe_rows[0], '0', '100',
                                 '200', '300', '400', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[1], '1600', '1700',
                                 '1800', '1900', '2000', '2100', '2200',
                                 '2300')
        self._check_event_coords(sorted_from_gtf_mxe_rows[2], '3200', '3300',
                                 '3400', '3500', '3600', '3700', '4000',
                                 '4100')
        self._check_event_coords(sorted_from_gtf_mxe_rows[3], '4200', '4300',
                                 '4600', '4700', '4800', '4900', '5000',
                                 '5100')
        self._check_event_coords(sorted_from_gtf_mxe_rows[4], '5200', '5300',
                                 '5400', '5500', '5600', '5700', '5800',
                                 '5900')

        from_gtf_novel_junction_mxe_path = os.path.join(
            self._out_dir, 'fromGTF.novelJunction.MXE.txt')
        from_gtf_novel_junction_mxe_header, from_gtf_novel_junction_mxe_rows, error = output_parser.parse_from_gtf_novel_junction(
            from_gtf_novel_junction_mxe_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_junction_mxe_rows), 4)
        sorted_from_gtf_novel_junction_mxe_rows = self._sort_by_event_coords(
            from_gtf_novel_junction_mxe_rows)
        self._check_event_coords(sorted_from_gtf_novel_junction_mxe_rows[0],
                                 '1600', '1700', '1800', '1900', '2000',
                                 '2100', '2200', '2300')
        self._check_event_coords(sorted_from_gtf_novel_junction_mxe_rows[1],
                                 '3200', '3300', '3400', '3500', '3600',
                                 '3700', '4000', '4100')
        self._check_event_coords(sorted_from_gtf_novel_junction_mxe_rows[2],
                                 '4200', '4300', '4600', '4700', '4800',
                                 '4900', '5000', '5100')
        self._check_event_coords(sorted_from_gtf_novel_junction_mxe_rows[3],
                                 '5200', '5300', '5400', '5500', '5600',
                                 '5700', '5800', '5900')

        from_gtf_novel_splice_site_mxe_path = os.path.join(
            self._out_dir, 'fromGTF.novelSpliceSite.MXE.txt')
        from_gtf_novel_splice_site_mxe_header, from_gtf_novel_splice_site_mxe_rows, error = output_parser.parse_from_gtf_novel_splice_site(
            from_gtf_novel_splice_site_mxe_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_splice_site_mxe_rows), 0)


class NovelSpliceSite(MutuallyExclusiveExonsNovelBaseTest):
    def _sub_test_name(self):
        return 'novel_splice_site'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.append('--novelSS')
        return arguments

    def _exons_by_transcript(self):
        return [
            [(1, 100), (201, 300), (601, 700)],
            [(1, 100), (401, 500), (601, 700)],
        ]

    def _paired_read_coords(self):
        return [
            ([[81, 98], [201, 300]], [[201, 300]]),
            ([[81, 96], [401, 500]], [[401, 500]]),
            ([[81, 94], [201, 300]], [[201, 300]]),
            ([[81, 94], [401, 500]], [[401, 500]]),
            ([[401, 500]], [[401, 500], [603, 620]]),
            ([[201, 300]], [[201, 300], [605, 620]]),
            ([[401, 500]], [[401, 500], [607, 620]]),
            ([[201, 300]], [[201, 300], [607, 620]]),
            ([[1, 100]], [[1, 100], [203, 220]]),
            ([[1, 100]], [[1, 100], [403, 420]]),
            ([[281, 298], [601, 700]], [[601, 700]]),
            ([[481, 498], [601, 700]], [[601, 700]]),
        ]

    def _check_results(self):
        self._check_no_error_results()

        from_gtf_mxe_path = os.path.join(self._out_dir, 'fromGTF.MXE.txt')
        from_gtf_mxe_header, from_gtf_mxe_rows, error = output_parser.parse_from_gtf(
            from_gtf_mxe_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_mxe_rows), 16)
        sorted_from_gtf_mxe_rows = self._sort_by_event_coords(
            from_gtf_mxe_rows)
        self._check_event_coords(sorted_from_gtf_mxe_rows[0], '0', '100',
                                 '200', '298', '400', '498', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[1], '0', '100',
                                 '200', '298', '400', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[2], '0', '100',
                                 '200', '298', '402', '498', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[3], '0', '100',
                                 '200', '298', '402', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[4], '0', '100',
                                 '200', '300', '400', '498', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[5], '0', '100',
                                 '200', '300', '400', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[6], '0', '100',
                                 '200', '300', '402', '498', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[7], '0', '100',
                                 '200', '300', '402', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[8], '0', '100',
                                 '202', '298', '400', '498', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[9], '0', '100',
                                 '202', '298', '400', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[10], '0', '100',
                                 '202', '298', '402', '498', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[11], '0', '100',
                                 '202', '298', '402', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[12], '0', '100',
                                 '202', '300', '400', '498', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[13], '0', '100',
                                 '202', '300', '400', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[14], '0', '100',
                                 '202', '300', '402', '498', '600', '700')
        self._check_event_coords(sorted_from_gtf_mxe_rows[15], '0', '100',
                                 '202', '300', '402', '500', '600', '700')

        from_gtf_novel_junction_mxe_path = os.path.join(
            self._out_dir, 'fromGTF.novelJunction.MXE.txt')
        from_gtf_novel_junction_mxe_header, from_gtf_novel_junction_mxe_rows, error = output_parser.parse_from_gtf_novel_junction(
            from_gtf_novel_junction_mxe_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_junction_mxe_rows), 0)

        from_gtf_novel_splice_site_mxe_path = os.path.join(
            self._out_dir, 'fromGTF.novelSpliceSite.MXE.txt')
        from_gtf_novel_splice_site_mxe_header, from_gtf_novel_splice_site_mxe_rows, error = output_parser.parse_from_gtf_novel_splice_site(
            from_gtf_novel_splice_site_mxe_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_splice_site_mxe_rows), 15)
        sorted_from_gtf_novel_splice_site_mxe_rows = self._sort_by_event_coords(
            from_gtf_novel_splice_site_mxe_rows)
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[0],
                                 '0', '100', '200', '298', '400', '498', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[1],
                                 '0', '100', '200', '298', '400', '500', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[2],
                                 '0', '100', '200', '298', '402', '498', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[3],
                                 '0', '100', '200', '298', '402', '500', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[4],
                                 '0', '100', '200', '300', '400', '498', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[5],
                                 '0', '100', '200', '300', '402', '498', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[6],
                                 '0', '100', '200', '300', '402', '500', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[7],
                                 '0', '100', '202', '298', '400', '498', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[8],
                                 '0', '100', '202', '298', '400', '500', '600',
                                 '700')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_mxe_rows[9],
                                 '0', '100', '202', '298', '402', '498', '600',
                                 '700')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_mxe_rows[10], '0', '100', '202',
            '298', '402', '500', '600', '700')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_mxe_rows[11], '0', '100', '202',
            '300', '400', '498', '600', '700')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_mxe_rows[12], '0', '100', '202',
            '300', '400', '500', '600', '700')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_mxe_rows[13], '0', '100', '202',
            '300', '402', '498', '600', '700')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_mxe_rows[14], '0', '100', '202',
            '300', '402', '500', '600', '700')


if __name__ == '__main__':
    unittest.main(verbosity=2)
