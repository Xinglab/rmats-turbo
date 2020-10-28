import os.path
import unittest

import tests.base_test
import tests.output_parser as output_parser
import tests.test_config
import tests.util


class SkippedExonNovelBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()

        self._test_base_dir = tests.test_config.TEST_BASE_DIR
        self._test_dir = os.path.join(self._test_base_dir,
                                      'skipped_exon_novel',
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
                      (int(r['upstreamES']), int(r['upstreamEE']),
                       int(r['exonStart_0base']), int(r['exonEnd']),
                       int(r['downstreamES']), int(r['downstreamEE'])))

    def _check_event_coords(self, row, u_start, u_end, start, end, d_start,
                            d_end):
        self.assertEqual(row['upstreamES'], u_start)
        self.assertEqual(row['upstreamEE'], u_end)
        self.assertEqual(row['exonStart_0base'], start)
        self.assertEqual(row['exonEnd'], end)
        self.assertEqual(row['downstreamES'], d_start)
        self.assertEqual(row['downstreamEE'], d_end)


class NovelJunction(SkippedExonNovelBaseTest):
    def _sub_test_name(self):
        return 'novel_junction'

    def test(self):
        self._run_test()

    def _exons_by_transcript(self):
        return [
            [(1, 100)],
            [(201, 300), (401, 500), (601, 700), (801, 900)],
            [(401, 500), (801, 900)],
            [(601, 700), (1001, 1100)],
        ]

    def _paired_read_coords(self):
        return [
            ([[81, 100], [201, 300]], [[201, 300]]),
            ([[81, 100], [401, 500]], [[401, 500]]),
            ([[281, 300], [601, 700]], [[601, 700]]),
            ([[881, 900], [1001, 1100]], [[1001, 1100]]),
        ]

    def _check_results(self):
        self._check_no_error_results()

        from_gtf_se_path = os.path.join(self._out_dir, 'fromGTF.SE.txt')
        from_gtf_se_header, from_gtf_se_rows, error = output_parser.parse_from_gtf(
            from_gtf_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_se_rows), 4)
        sorted_from_gtf_se_rows = self._sort_by_event_coords(from_gtf_se_rows)
        self._check_event_coords(sorted_from_gtf_se_rows[0], '0', '100', '200',
                                 '300', '400', '500')
        self._check_event_coords(sorted_from_gtf_se_rows[1], '200', '300',
                                 '400', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_se_rows[2], '400', '500',
                                 '600', '700', '800', '900')
        self._check_event_coords(sorted_from_gtf_se_rows[3], '600', '700',
                                 '800', '900', '1000', '1100')

        from_gtf_novel_junction_se_path = os.path.join(
            self._out_dir, 'fromGTF.novelJunction.SE.txt')
        from_gtf_novel_junction_se_header, from_gtf_novel_junction_se_rows, error = output_parser.parse_from_gtf_novel_junction(
            from_gtf_novel_junction_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_junction_se_rows), 3)
        sorted_from_gtf_novel_junction_se_rows = self._sort_by_event_coords(
            from_gtf_novel_junction_se_rows)
        self._check_event_coords(sorted_from_gtf_novel_junction_se_rows[0],
                                 '0', '100', '200', '300', '400', '500')
        self._check_event_coords(sorted_from_gtf_novel_junction_se_rows[1],
                                 '200', '300', '400', '500', '600', '700')
        self._check_event_coords(sorted_from_gtf_novel_junction_se_rows[2],
                                 '600', '700', '800', '900', '1000', '1100')

        from_gtf_novel_splice_site_se_path = os.path.join(
            self._out_dir, 'fromGTF.novelSpliceSite.SE.txt')
        from_gtf_novel_splice_site_se_header, from_gtf_novel_splice_site_se_rows, error = output_parser.parse_from_gtf_novel_splice_site(
            from_gtf_novel_splice_site_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_splice_site_se_rows), 0)


class NovelSpliceSite(SkippedExonNovelBaseTest):
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
            [(1, 100), (201, 300), (401, 500)],
            [(1, 100), (401, 500)],
            [(601, 700), (801, 900), (1001, 1100)],
        ]

    def _paired_read_coords(self):
        return [
            ([[81, 98], [201, 300]], [[201, 300]]),
            ([[81, 96], [201, 300]], [[201, 300]]),
            ([[81, 96], [401, 500]], [[401, 500]]),
            ([[201, 300]], [[201, 300], [403, 420]]),
            ([[201, 300]], [[201, 300], [405, 420]]),
            ([[1, 100]], [[1, 100], [405, 420]]),
            ([[81, 96], [405, 500]], [[405, 500]]),
            ([[81, 100], [203, 300]], [[203, 300]]),
            ([[81, 96], [203, 300]], [[203, 300]]),
            ([[201, 298]], [[201, 298], [401, 420]]),
            ([[201, 298]], [[201, 298], [405, 420]]),
            ([[681, 700], [1001, 1100]], [[1001, 1100]]),
        ]

    def _check_results(self):
        self._check_no_error_results()

        from_gtf_se_path = os.path.join(self._out_dir, 'fromGTF.SE.txt')
        from_gtf_se_header, from_gtf_se_rows, error = output_parser.parse_from_gtf(
            from_gtf_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_se_rows), 7)
        sorted_from_gtf_se_rows = self._sort_by_event_coords(from_gtf_se_rows)
        self._check_event_coords(sorted_from_gtf_se_rows[0], '0', '96', '200',
                                 '300', '400', '500')
        self._check_event_coords(sorted_from_gtf_se_rows[1], '0', '100', '200',
                                 '298', '400', '500')
        self._check_event_coords(sorted_from_gtf_se_rows[2], '0', '100', '200',
                                 '300', '400', '500')
        self._check_event_coords(sorted_from_gtf_se_rows[3], '0', '100', '200',
                                 '300', '404', '500')
        self._check_event_coords(sorted_from_gtf_se_rows[4], '0', '100', '202',
                                 '298', '400', '500')
        self._check_event_coords(sorted_from_gtf_se_rows[5], '0', '100', '202',
                                 '300', '400', '500')
        self._check_event_coords(sorted_from_gtf_se_rows[6], '600', '700',
                                 '800', '900', '1000', '1100')

        from_gtf_novel_junction_se_path = os.path.join(
            self._out_dir, 'fromGTF.novelJunction.SE.txt')
        from_gtf_novel_junction_se_header, from_gtf_novel_junction_se_rows, error = output_parser.parse_from_gtf_novel_junction(
            from_gtf_novel_junction_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_junction_se_rows), 1)
        self._check_event_coords(from_gtf_novel_junction_se_rows[0], '600',
                                 '700', '800', '900', '1000', '1100')

        from_gtf_novel_splice_site_se_path = os.path.join(
            self._out_dir, 'fromGTF.novelSpliceSite.SE.txt')
        from_gtf_novel_splice_site_se_header, from_gtf_novel_splice_site_se_rows, error = output_parser.parse_from_gtf_novel_splice_site(
            from_gtf_novel_splice_site_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_splice_site_se_rows), 5)
        sorted_from_gtf_novel_splice_site_se_rows = self._sort_by_event_coords(
            from_gtf_novel_splice_site_se_rows)
        self._check_event_coords(sorted_from_gtf_novel_splice_site_se_rows[0],
                                 '0', '96', '200', '300', '400', '500')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_se_rows[1],
                                 '0', '100', '200', '298', '400', '500')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_se_rows[2],
                                 '0', '100', '200', '300', '404', '500')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_se_rows[3],
                                 '0', '100', '202', '298', '400', '500')
        self._check_event_coords(sorted_from_gtf_novel_splice_site_se_rows[4],
                                 '0', '100', '202', '300', '400', '500')


class NovelSpliceSiteNotNovelJunction(SkippedExonNovelBaseTest):
    '''
    There was a bug where --novelSS would label an event as
    novelJunction even though the event was fully annotated. This
    would happen if one of the exons in the event could have a novel
    splice site, but the novel splice site was not necessary to detect
    the event.
    '''
    def _sub_test_name(self):
        return 'nss_not_nj'

    def test(self):
        self._run_test()

    def _rmats_arguments(self):
        arguments = super()._rmats_arguments()
        arguments.append('--novelSS')
        return arguments

    def _exons_by_transcript(self):
        return [
            [(1, 100), (201, 300), (401, 500), (601, 700)],
            [(1, 100), (201, 300), (601, 700)],
        ]

    def _paired_read_coords(self):
        return [
            ([[81, 100], [205, 300]], [[205, 300]]),
            ([[281, 300], [401, 500]], [[401, 500]]),
            ([[281, 300], [601, 700]], [[601, 700]]),
            ([[401, 500]], [[401, 500], [601, 620]]),
            ([[201, 300]], [[201, 300], [601, 620]]),
        ]

    def _check_results(self):
        self._check_no_error_results()

        from_gtf_se_path = os.path.join(self._out_dir, 'fromGTF.SE.txt')
        from_gtf_se_header, from_gtf_se_rows, error = (
            output_parser.parse_from_gtf(from_gtf_se_path))
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_se_rows), 1)
        self._check_event_coords(from_gtf_se_rows[0], '200', '300', '400', '500',
                                 '600', '700')

        from_gtf_novel_junction_se_path = os.path.join(
            self._out_dir, 'fromGTF.novelJunction.SE.txt')
        (from_gtf_novel_junction_se_header, from_gtf_novel_junction_se_rows,
         error) = output_parser.parse_from_gtf_novel_junction(
             from_gtf_novel_junction_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_junction_se_rows), 0)

        from_gtf_novel_splice_site_se_path = os.path.join(
            self._out_dir, 'fromGTF.novelSpliceSite.SE.txt')
        (from_gtf_novel_splice_site_se_header,
         from_gtf_novel_splice_site_se_rows,
         error) = output_parser.parse_from_gtf_novel_splice_site(
             from_gtf_novel_splice_site_se_path)
        self.assertFalse(error)
        self.assertEqual(len(from_gtf_novel_splice_site_se_rows), 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
