import os.path
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class Alternative3SpliceSiteNovelBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir,
                                      'alternative_3_splice_site_novel',
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
                      (int(r['flankingES']), int(r['flankingEE']),
                       int(r['longExonStart_0base']), int(r['shortES']),
                       int(r['shortEE']), int(r['longExonEnd'])))

    def _check_event_coords(self, row, f_start, f_end, l_start, s_start, end):
        self.assertEqual(row['flankingES'], f_start)
        self.assertEqual(row['flankingEE'], f_end)
        self.assertEqual(row['longExonStart_0base'], l_start)
        self.assertEqual(row['shortES'], s_start)
        self.assertEqual(row['shortEE'], end)
        self.assertEqual(row['longExonEnd'], end)


class NovelJunction(Alternative3SpliceSiteNovelBaseTest):
    def _sub_test_name(self):
        return 'novel_junction'

    def test(self):
        self._run_test()

    def _exons_by_transcript(self):
        return [
            [(1, 100), (201, 400)],
            [(1, 100), (301, 400)],
            [(801, 900)],
            [(501, 600), (701, 900)],
            [(1001, 1100), (1301, 1400)],
            [(1201, 1400)],
        ]

    def _paired_read_coords(self):
        return [
            ([[501, 600]], [[501, 600], [801, 820]]),
            ([[1001, 1100]], [[1001, 1100], [1201, 1220]]),
        ]

    def _check_results(self):
        self._check_no_error_results()

        from_gtf_a3ss_path = os.path.join(self._out_dir, 'fromGTF.A3SS.txt')
        from_gtf_a3ss_header, from_gtf_a3ss_rows, error = output_parser.parse_from_gtf(
            from_gtf_a3ss_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_a3ss_rows), 3)
        sorted_from_gtf_a3ss_rows = self._sort_by_event_coords(
            from_gtf_a3ss_rows)
        self._check_event_coords(sorted_from_gtf_a3ss_rows[0], '0', '100',
                                 '200', '300', '400')
        self._check_event_coords(sorted_from_gtf_a3ss_rows[1], '500', '600',
                                 '700', '800', '900')
        self._check_event_coords(sorted_from_gtf_a3ss_rows[2], '1000', '1100',
                                 '1200', '1300', '1400')

        from_gtf_novel_junction_a3ss_path = os.path.join(
            self._out_dir, 'fromGTF.novelJunction.A3SS.txt')
        from_gtf_novel_junction_a3ss_header, from_gtf_novel_junction_a3ss_rows, error = output_parser.parse_from_gtf_novel_junction(
            from_gtf_novel_junction_a3ss_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_junction_a3ss_rows), 2)
        sorted_from_gtf_novel_junction_a3ss_rows = self._sort_by_event_coords(
            from_gtf_novel_junction_a3ss_rows)
        self._check_event_coords(sorted_from_gtf_novel_junction_a3ss_rows[0],
                                 '500', '600', '700', '800', '900')
        self._check_event_coords(sorted_from_gtf_novel_junction_a3ss_rows[1],
                                 '1000', '1100', '1200', '1300', '1400')

        from_gtf_novel_splice_site_a3ss_path = os.path.join(
            self._out_dir, 'fromGTF.novelSpliceSite.A3SS.txt')
        from_gtf_novel_splice_site_a3ss_header, from_gtf_novel_splice_site_a3ss_rows, error = output_parser.parse_from_gtf_novel_splice_site(
            from_gtf_novel_splice_site_a3ss_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_splice_site_a3ss_rows), 0)


class NovelSpliceSite(Alternative3SpliceSiteNovelBaseTest):
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
            [(1, 100), (201, 400)],
            [(1, 100), (301, 400)],
        ]

    def _paired_read_coords(self):
        return [
            ([[81, 98], [201, 400]], [[201, 400]]),
            ([[81, 96], [301, 400]], [[301, 400]]),
            ([[81, 94], [201, 400]], [[201, 400]]),
            ([[81, 94], [301, 400]], [[301, 400]]),
            ([[1, 100]], [[1, 100], [303, 320]]),
            ([[1, 100]], [[1, 100], [203, 220]]),
        ]

    def _check_results(self):
        self._check_no_error_results()

        from_gtf_a3ss_path = os.path.join(self._out_dir, 'fromGTF.A3SS.txt')
        from_gtf_a3ss_header, from_gtf_a3ss_rows, error = output_parser.parse_from_gtf(
            from_gtf_a3ss_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_a3ss_rows), 7)
        sorted_from_gtf_a3ss_rows = self._sort_by_event_coords(
            from_gtf_a3ss_rows)
        self._check_event_coords(sorted_from_gtf_a3ss_rows[0], '0', '94',
                                 '200', '300', '400')
        self._check_event_coords(sorted_from_gtf_a3ss_rows[1], '0', '100',
                                 '200', '202', '400')
        self._check_event_coords(sorted_from_gtf_a3ss_rows[2], '0', '100',
                                 '200', '300', '400')
        self._check_event_coords(sorted_from_gtf_a3ss_rows[3], '0', '100',
                                 '200', '302', '400')
        self._check_event_coords(sorted_from_gtf_a3ss_rows[4], '0', '100',
                                 '202', '300', '400')
        self._check_event_coords(sorted_from_gtf_a3ss_rows[5], '0', '100',
                                 '202', '302', '400')
        self._check_event_coords(sorted_from_gtf_a3ss_rows[6], '0', '100',
                                 '300', '302', '400')

        from_gtf_novel_junction_a3ss_path = os.path.join(
            self._out_dir, 'fromGTF.novelJunction.A3SS.txt')
        from_gtf_novel_junction_a3ss_header, from_gtf_novel_junction_a3ss_rows, error = output_parser.parse_from_gtf_novel_junction(
            from_gtf_novel_junction_a3ss_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_junction_a3ss_rows), 0)

        from_gtf_novel_splice_site_a3ss_path = os.path.join(
            self._out_dir, 'fromGTF.novelSpliceSite.A3SS.txt')
        from_gtf_novel_splice_site_a3ss_header, from_gtf_novel_splice_site_a3ss_rows, error = output_parser.parse_from_gtf_novel_splice_site(
            from_gtf_novel_splice_site_a3ss_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_splice_site_a3ss_rows), 6)
        sorted_from_gtf_novel_splice_site_a3ss_rows = self._sort_by_event_coords(
            from_gtf_novel_splice_site_a3ss_rows)
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_a3ss_rows[0], '0', '94', '200',
            '300', '400')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_a3ss_rows[1], '0', '100', '200',
            '202', '400')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_a3ss_rows[2], '0', '100', '200',
            '302', '400')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_a3ss_rows[3], '0', '100', '202',
            '300', '400')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_a3ss_rows[4], '0', '100', '202',
            '302', '400')
        self._check_event_coords(
            sorted_from_gtf_novel_splice_site_a3ss_rows[5], '0', '100', '300',
            '302', '400')


if __name__ == '__main__':
    unittest.main(verbosity=2)
