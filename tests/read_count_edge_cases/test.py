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
        self._test_dir = os.path.join(self._test_base_dir,
                                      'read_count_edge_cases')
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
            '--individual-counts',
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

    # A read can have a novel splice site recorded as (-1,-1) as long as there is
    # another junction in the read with annotated splice sites.
    # Each event is in a different gene to test cases where a
    # read goes past the end of a gene.
    def _exons_and_genes_by_transcript(self):
        exons_and_genes = [
            # SE
            ([(1, 10), (21, 30), (41, 50)], [1]),
            ([(1, 10), (41, 50)], [1]),
            # A5SS (extra junction at beginning to allow novel junction)
            ([(181, 185), (190, 195), (201, 300), (501, 600)], [2]),
            ([(181, 185), (190, 195), (201, 400), (501, 600)], [2]),
            # MXE
            ([(701, 710), (721, 725), (741, 750)], [3]),
            ([(701, 710), (731, 735), (741, 750)], [3]),
            # RI (extra junction at beginning to allow novel junction)
            ([(791, 795), (801, 830), (901, 1000)], [4]),
            ([(791, 795), (801, 1000)], [4]),
        ]
        exons_by_transcript = list()
        genes_by_transcript = list()
        for exons, genes in exons_and_genes:
            for gene in genes:
                exons_by_transcript.append(exons)
                genes_by_transcript.append(gene)

        return exons_by_transcript, genes_by_transcript

    # Should not count as upstream to target because the read goes past the target end.
    def _se_upstream_to_target_past_end(self):
        return ([[8, 10], [21, 100]], [[21, 100]])

    # Should not count as across short because the novel junction starts in the long exon.
    def _a5ss_novel_junction_into_long(self):
        return ([[181, 185], [190, 195], [320, 400]], [[320, 400]])

    # Should not count as upstream to first because the read goes past the first exon.
    def _mxe_upstream_to_first_past_end(self):
        return ([[708, 710], [721, 800]], [[721, 800]])

    # Should not count as upstream to second because the read goes past the second exon.
    def _mxe_upstream_to_second_past_end(self):
        return ([[708, 710], [731, 800]], [[731, 800]])

    # Should not count as upstream to intron
    def _ri_novel_junction_to_end_after_upstream(self):
        return ([[791, 795], [801, 810], [915, 1000]], [[915, 1000]])

    def _paired_read_coords_1_1(self):
        return [
            self._se_upstream_to_target_past_end(),
            self._a5ss_novel_junction_into_long(),
            self._mxe_upstream_to_first_past_end(),
            self._mxe_upstream_to_second_past_end(),
            self._ri_novel_junction_to_end_after_upstream(),
        ]

    def _paired_read_coords_1_2(self):
        return [
            self._se_upstream_to_target_past_end(),
        ]

    def _check_results(self):
        self._check_no_error_results()

        se_mats_jcec_path = os.path.join(self._out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            se_mats_jcec_path)
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header,
                                        has_individual=True)
        self.assertEqual(len(se_mats_jcec_rows), 1)
        row = se_mats_jcec_rows[0]
        self.assertEqual(row['upstreamES'], '0')
        self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
        self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(row['IJC_SAMPLE_2'], '')
        self.assertEqual(row['SJC_SAMPLE_2'], '')
        self.assertEqual(row['upstream_to_target_count'], '0,0')
        self.assertEqual(row['target_to_downstream_count'], '0,0')
        self.assertEqual(row['target_count'], '0,0')
        self.assertEqual(row['upstream_to_downstream_count'], '0,0')

        a5ss_mats_jcec_path = os.path.join(self._out_dir, 'A5SS.MATS.JCEC.txt')
        a5ss_mats_jcec_header, a5ss_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            a5ss_mats_jcec_path)
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a5ss_mats_jcec_header,
                                           has_individual=True)
        self.assertEqual(len(a5ss_mats_jcec_rows), 1)
        row = a5ss_mats_jcec_rows[0]
        self.assertEqual(row['shortES'], '200')
        self.assertEqual(row['IJC_SAMPLE_1'], '1,0')
        self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(row['IJC_SAMPLE_2'], '')
        self.assertEqual(row['SJC_SAMPLE_2'], '')
        self.assertEqual(row['across_short_boundary_count'], '0,0')
        self.assertEqual(row['long_to_flanking_count'], '0,0')
        self.assertEqual(row['exclusive_to_long_count'], '1,0')
        self.assertEqual(row['short_to_flanking_count'], '0,0')

        mxe_mats_jcec_path = os.path.join(self._out_dir, 'MXE.MATS.JCEC.txt')
        mxe_mats_jcec_header, mxe_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            mxe_mats_jcec_path)
        self.assertFalse(error)
        self._check_mxe_mats_jcec_header(mxe_mats_jcec_header,
                                         has_individual=True)
        self.assertEqual(len(mxe_mats_jcec_rows), 1)
        row = mxe_mats_jcec_rows[0]
        self.assertEqual(row['upstreamES'], '700')
        self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
        self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(row['IJC_SAMPLE_2'], '')
        self.assertEqual(row['SJC_SAMPLE_2'], '')
        self.assertEqual(row['upstream_to_first_count'], '0,0')
        self.assertEqual(row['first_to_downstream_count'], '0,0')
        self.assertEqual(row['first_count'], '0,0')
        self.assertEqual(row['upstream_to_second_count'], '0,0')
        self.assertEqual(row['second_to_downstream_count'], '0,0')
        self.assertEqual(row['second_count'], '0,0')

        ri_mats_jcec_path = os.path.join(self._out_dir, 'RI.MATS.JCEC.txt')
        ri_mats_jcec_header, ri_mats_jcec_rows, error = output_parser.parse_mats_jcec(
            ri_mats_jcec_path)
        self.assertFalse(error)
        self._check_ri_mats_jcec_header(ri_mats_jcec_header,
                                        has_individual=True)
        self.assertEqual(len(ri_mats_jcec_rows), 1)
        row = ri_mats_jcec_rows[0]
        self.assertEqual(row['upstreamES'], '800')
        self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
        self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
        self.assertEqual(row['IJC_SAMPLE_2'], '')
        self.assertEqual(row['SJC_SAMPLE_2'], '')
        self.assertEqual(row['upstream_to_intron_count'], '0,0')
        self.assertEqual(row['intron_to_downstream_count'], '0,0')
        self.assertEqual(row['intron_count'], '0,0')
        self.assertEqual(row['upstream_to_downstream_count'], '0,0')


if __name__ == '__main__':
    unittest.main(verbosity=2)
