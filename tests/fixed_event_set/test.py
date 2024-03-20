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
        self._test_dir = os.path.join(self._test_base_dir, 'fixed_event_set')
        self._generated_input_dir = os.path.join(self._test_dir,
                                                 'generated_input')
        self._out_dir_orig = os.path.join(self._test_dir, 'out_orig')
        self._tmp_dir_orig = os.path.join(self._test_dir, 'tmp_orig')
        self._out_dir_fixed = os.path.join(self._test_dir, 'out_fixed')
        self._tmp_dir_fixed = os.path.join(self._test_dir, 'tmp_fixed')
        self._out_dir_id_strings = os.path.join(self._test_dir,
                                                'out_id_strings')
        self._tmp_dir_id_strings = os.path.join(self._test_dir,
                                                'tmp_id_strings')
        self._id_strings_dir = os.path.join(self._test_dir, 'id_strings')

        tests.util.recreate_dirs([
            self._generated_input_dir, self._out_dir_orig, self._tmp_dir_orig,
            self._out_dir_fixed, self._tmp_dir_fixed, self._out_dir_id_strings,
            self._tmp_dir_id_strings, self._id_strings_dir,
            self._command_output_dir()
        ])

        self._read_type = 'paired'
        self._read_length = 50
        self._chromosome_length = 16000

        self._sample_1_bams_path = os.path.join(self._generated_input_dir,
                                                'b1.txt')
        self._sample_2_bams_path = os.path.join(self._generated_input_dir,
                                                'b2.txt')
        self._sample_3_bams_path = os.path.join(self._generated_input_dir,
                                                'b3.txt')
        sample_1_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_1_rep_{}.bam')
        sample_2_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_2_rep_{}.bam')
        sample_3_bam_replicate_template = os.path.join(
            self._generated_input_dir, 'sample_3_rep_{}.bam')
        self._sample_1_bams = self._create_sample_1_bams(
            self._sample_1_bams_path, sample_1_bam_replicate_template)
        self._sample_2_bams = self._create_sample_2_bams(
            self._sample_2_bams_path, sample_2_bam_replicate_template)
        self._sample_3_bams = self._create_sample_3_bams(
            self._sample_3_bams_path, sample_3_bam_replicate_template)

        self._gtf_path = os.path.join(self._generated_input_dir, 'test.gtf')
        transcript_exons, transcript_genes, transcript_strands = (
            self._exons_genes_and_strands_by_transcript())
        self._gtf = self._create_gtf_from_transcripts(
            self._gtf_path,
            transcript_exons,
            genes=transcript_genes,
            strands=transcript_strands)
        self._sub_steps = [
            'original',
            'fixed',
            'id_strings',
        ]
        self._sub_step = None

    def test(self):
        for sub_step in self._sub_steps:
            self._sub_step = sub_step
            if sub_step == 'id_strings':
                self._setup_id_strings()

            self._run_test()

    def _command_output_dir(self):
        return os.path.join(self._test_dir, 'command_output')

    def _rmats_arguments(self):
        arguments = [
            '--gtf', self._gtf_path, '-t', self._read_type, '--readLength',
            str(self._read_length), '--task', 'both', '--statoff'
        ]
        if self._sub_step == 'original':
            arguments.extend([
                '--od', self._out_dir_orig, '--tmp', self._tmp_dir_orig,
                '--b1', self._sample_1_bams_path, '--b2',
                self._sample_2_bams_path
            ])
        if self._sub_step == 'fixed':
            arguments.extend([
                '--od', self._out_dir_fixed, '--tmp', self._tmp_dir_fixed,
                '--b1', self._sample_1_bams_path, '--b2',
                self._sample_3_bams_path, '--fixed-event-set',
                self._out_dir_orig
            ])
        if self._sub_step == 'id_strings':
            arguments.extend([
                '--od', self._out_dir_id_strings, '--tmp',
                self._tmp_dir_id_strings, '--b1', self._sample_1_bams_path,
                '--b2', self._sample_3_bams_path, '--fixed-event-set',
                self._id_strings_dir
            ])

        return arguments

    def _start_coord_to_id(self):
        return {
            '200': 'SE_1',
            '1200': 'SE_2',
            '3600': 'SE_MXE_1',
            '4600': 'SE_MXE_2',
            '3200': 'MXE_1',
            '4200': 'MXE_2',
            '6000': 'A5_1',
            '7000': 'A5_2',
            '9200': 'A3_1',
            '10200': 'A3_2',
            '12000': 'RI_1',
            '13000': 'RI_2',
            '15600': 'SE_MXE_minus',
            '15200': 'MXE_minus',
        }

    def _as_event_to_start_coord_header(self):
        return {
            'SE': 'exonStart_0base',
            'MXE': '1stExonStart_0base',
            'A3SS': 'longExonStart_0base',
            'A5SS': 'longExonStart_0base',
            'RI': 'riExonStart_0base'
        }

    def _setup_id_strings(self):
        coord_to_id = self._start_coord_to_id()
        event_to_coord_header = self._as_event_to_start_coord_header()
        for as_event, coord_header in event_to_coord_header.items():
            name = 'fromGTF.{}.txt'.format(as_event)
            orig_path = os.path.join(self._out_dir_orig, name)
            dest_path = os.path.join(self._id_strings_dir, name)
            with open(orig_path, 'rt') as in_f:
                with open(dest_path, 'wt') as out_f:
                    for i, line in enumerate(in_f):
                        values = line.strip().split('\t')
                        if i == 0:
                            headers = values
                            id_i = headers.index('ID')
                            coord_i = headers.index(coord_header)
                            out_f.write(line)
                            continue

                        coord_val = values[coord_i]
                        values[id_i] = coord_to_id[coord_val]
                        out_f.write('{}\n'.format('\t'.join(values)))

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

    def _create_sample_3_bams(self, sample_3_bams_path,
                              sample_3_replicate_template):
        rep_1_bam_path = sample_3_replicate_template.format(1)
        rep_1_bam = self._create_bam_from_paired_read_coords(
            rep_1_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords_3_1())

        rep_2_bam_path = sample_3_replicate_template.format(2)
        rep_2_bam = self._create_bam_from_paired_read_coords(
            rep_2_bam_path, self._chromosome_length, self._read_length,
            self._paired_read_coords_3_2())

        sample_3_bams = [rep_1_bam, rep_2_bam]
        self._write_bams(sample_3_bams, sample_3_bams_path)
        return sample_3_bams

    def _exons_genes_and_strands_by_transcript(self):
        # exons, genes, strands
        exons_genes_and_strands = [
            ([(1, 100), (201, 300), (401, 500)], [1], ['+']),  # SE 1
            ([(1001, 1100), (1201, 1300), (1401, 1500)], [1], ['+']),  # SE 2
            ([(2001, 2100), (2201, 2300), (2401, 2500)], [1], ['+']),  # SE 3
            # MXE inc 1
            ([(3001, 3100), (3201, 3300), (3801, 3900)], [1], ['+']),
            # MXE skip 1
            ([(3001, 3100), (3401, 3500), (3601, 3700),
              (3801, 3900)], [1], ['+']),
            # MXE inc 2
            ([(4001, 4100), (4201, 4300), (4801, 4900)], [1], ['+']),
            # MXE skip 2
            ([(4001, 4100), (4401, 4500), (4601, 4700),
              (4801, 4900)], [1], ['+']),
            # MXE inc 3
            ([(5001, 5100), (5201, 5300), (5801, 5900)], [1], ['+']),
            # MXE skip 3
            ([(5001, 5100), (5401, 5500), (5601, 5700),
              (5801, 5900)], [1], ['+']),
            ([(6001, 6200), (6301, 6400)], [1], ['+']),  # A5SS inc 1
            ([(6001, 6100)], [1], ['+']),  # A5SS skip 1
            ([(7001, 7200), (7301, 7400)], [1], ['+']),  # A5SS inc 2
            ([(7001, 7100)], [1], ['+']),  # A5SS skip 2
            ([(8001, 8200), (8301, 8400)], [1], ['+']),  # A5SS inc 3
            ([(8001, 8100)], [1], ['+']),  # A5SS skip 3
            ([(9001, 9100), (9201, 9400)], [1], ['+']),  # A3SS inc 1
            ([(9301, 9400)], [1], ['+']),  # A3SS skip 1
            ([(10001, 10100), (10201, 10400)], [1], ['+']),  # A3SS inc 2
            ([(10301, 10400)], [1], ['+']),  # A3SS skip 2
            ([(11001, 11100), (11201, 11400)], [1], ['+']),  # A3SS inc 3
            ([(11301, 11400)], [1], ['+']),  # A3SS skip 3
            ([(12001, 12400)], [1], ['+']),  # RI inc 1
            ([(12001, 12100)], [1], ['+']),  # RI skip 1_1
            ([(12301, 12400)], [1], ['+']),  # RI skip 1_2
            ([(13001, 13400)], [1], ['+']),  # RI inc 2
            ([(13001, 13100)], [1], ['+']),  # RI skip 2_1
            ([(13301, 13400)], [1], ['+']),  # RI skip 2_2
            ([(14001, 14400)], [1], ['+']),  # RI inc 3
            ([(14001, 14100)], [1], ['+']),  # RI skip 3_1
            ([(14301, 14400)], [1], ['+']),  # RI skip 3_2
            # MXE inc '-' strand
            ([(15001, 15100), (15201, 15300), (15801, 15900)], [2], ['-']),
            # MXE inc '-' strand
            ([(15001, 15100), (15401, 15430), (15601, 15700),
              (15801, 15900)], [2], ['-']),
        ]
        exons_by_transcript = list()
        genes_by_transcript = list()
        strands_by_transcript = list()
        for exons, genes, strands in exons_genes_and_strands:
            for gene in genes:
                for strand in strands:
                    exons_by_transcript.append(exons)
                    genes_by_transcript.append(gene)
                    strands_by_transcript.append(strand)

        return exons_by_transcript, genes_by_transcript, strands_by_transcript

    def _include_read_SE_1(self):
        return ([[81, 100], [201, 300]], [[201, 300]])

    def _skip_read_SE_1(self):
        return ([[81, 100], [401, 500]], [[401, 500]])

    def _include_read_SE_2(self):
        return ([[1081, 1100], [1201, 1300]], [[1201, 1300]])

    def _skip_read_SE_2(self):
        return ([[1081, 1100], [1401, 1500]], [[1401, 1500]])

    def _include_read_SE_3(self):
        return ([[2081, 2100], [2201, 2300]], [[2201, 2300]])

    def _skip_read_SE_3(self):
        return ([[2081, 2100], [2401, 2500]], [[2401, 2500]])

    def _include_read_MXE_1(self):
        return ([[3081, 3100], [3201, 3300]], [[3201, 3300]])

    def _skip_read_MXE_1(self):
        return ([[3481, 3500], [3801, 3900]], [[3801, 3900]])

    def _include_read_MXE_2(self):
        return ([[4081, 4100], [4201, 4300]], [[4201, 4300]])

    def _skip_read_MXE_2(self):
        return ([[4481, 4500], [4801, 4900]], [[4801, 4900]])

    def _include_read_MXE_3(self):
        return ([[5081, 5100], [5201, 5300]], [[5201, 5300]])

    def _skip_read_MXE_3(self):
        return ([[5481, 5500], [5801, 5900]], [[5801, 5900]])

    def _include_read_A5SS_1(self):
        return ([[6181, 6200], [6301, 6400]], [[6301, 6400]])

    def _skip_read_A5SS_1(self):
        return ([[6081, 6100], [6301, 6400]], [[6301, 6400]])

    def _include_read_A5SS_2(self):
        return ([[7181, 7200], [7301, 7400]], [[7301, 7400]])

    def _skip_read_A5SS_2(self):
        return ([[7081, 7100], [7301, 7400]], [[7301, 7400]])

    def _include_read_A5SS_3(self):
        return ([[8181, 8200], [8301, 8400]], [[8301, 8400]])

    def _skip_read_A5SS_3(self):
        return ([[8081, 8100], [8301, 8400]], [[8301, 8400]])

    def _include_read_A3SS_1(self):
        return ([[9001, 9100]], [[9001, 9100], [9201, 9220]])

    def _skip_read_A3SS_1(self):
        return ([[9001, 9100]], [[9001, 9100], [9301, 9320]])

    def _include_read_A3SS_2(self):
        return ([[10001, 10100]], [[10001, 10100], [10201, 10220]])

    def _skip_read_A3SS_2(self):
        return ([[10001, 10100]], [[10001, 10100], [10301, 10320]])

    def _include_read_A3SS_3(self):
        return ([[11001, 11100]], [[11001, 11100], [11201, 11220]])

    def _skip_read_A3SS_3(self):
        return ([[11001, 11100]], [[11001, 11100], [11301, 11320]])

    def _include_read_RI_1(self):
        return ([[12081, 12400]], [[12001, 12320]])

    def _skip_read_RI_1(self):
        return ([[12081, 12100], [12301, 12400]], [[12301, 12400]])

    def _include_read_RI_2(self):
        return ([[13081, 13400]], [[13001, 13320]])

    def _skip_read_RI_2(self):
        return ([[13081, 13100], [13301, 13400]], [[13301, 13400]])

    def _include_read_RI_3(self):
        return ([[14081, 14400]], [[14001, 14320]])

    def _skip_read_RI_3(self):
        return ([[14081, 14100], [14301, 14400]], [[14301, 14400]])

    # For '-' strand MXE the higher coordinate exon is the inclusion isoform
    def _include_read_MXE_minus_strand(self):
        return ([[15411, 15430], [15801, 15900]], [[15801, 15900]])

    def _skip_read_MXE_minus_strand(self):
        return ([[15081, 15100], [15201, 15300]], [[15201, 15300]])

    def _reads_by_count(self, i1, s1, i2, s2, i3, s3, minus_inc, minus_skip):
        reads = list()
        reads.extend(i1 * [
            self._include_read_SE_1(),
            self._include_read_MXE_1(),
            self._include_read_A5SS_1(),
            self._include_read_A3SS_1(),
            self._include_read_RI_1()
        ])
        reads.extend(s1 * [
            self._skip_read_SE_1(),
            self._skip_read_MXE_1(),
            self._skip_read_A5SS_1(),
            self._skip_read_A3SS_1(),
            self._skip_read_RI_1()
        ])
        reads.extend(i2 * [
            self._include_read_SE_2(),
            self._include_read_MXE_2(),
            self._include_read_A5SS_2(),
            self._include_read_A3SS_2(),
            self._include_read_RI_2()
        ])
        reads.extend(s2 * [
            self._skip_read_SE_2(),
            self._skip_read_MXE_2(),
            self._skip_read_A5SS_2(),
            self._skip_read_A3SS_2(),
            self._skip_read_RI_2()
        ])
        reads.extend(i3 * [
            self._include_read_SE_3(),
            self._include_read_MXE_3(),
            self._include_read_A5SS_3(),
            self._include_read_A3SS_3(),
            self._include_read_RI_3()
        ])
        reads.extend(s3 * [
            self._skip_read_SE_3(),
            self._skip_read_MXE_3(),
            self._skip_read_A5SS_3(),
            self._skip_read_A3SS_3(),
            self._skip_read_RI_3()
        ])
        reads.extend(minus_inc * [self._include_read_MXE_minus_strand()])
        reads.extend(minus_skip * [self._skip_read_MXE_minus_strand()])

        return reads

    def _paired_read_coords_1_1(self):
        return self._reads_by_count(1, 1, 0, 0, 0, 0, 1, 2)

    def _paired_read_coords_1_2(self):
        return self._reads_by_count(1, 1, 0, 0, 0, 0, 3, 4)

    def _paired_read_coords_2_1(self):
        return self._reads_by_count(0, 0, 1, 1, 0, 0, 5, 6)

    def _paired_read_coords_2_2(self):
        return self._reads_by_count(1, 1, 1, 1, 1, 0, 7, 8)

    def _paired_read_coords_3_1(self):
        return self._reads_by_count(0, 0, 0, 0, 1, 1, 9, 10)

    def _paired_read_coords_3_2(self):
        return self._reads_by_count(1, 1, 1, 0, 1, 1, 11, 12)

    def _check_results(self):
        if self._sub_step == 'original':
            self._check_results_original()
        elif self._sub_step == 'fixed':
            self._check_results_fixed()
        elif self._sub_step == 'id_strings':
            self._check_results_id_strings()
        else:
            self.fail('unexpected sub_step: {}'.format(self._sub_step))

    def _check_results_original(self):
        self._check_no_error_results()

        event_to_coord_header = self._as_event_to_start_coord_header()

        se_mats_jc_path = os.path.join(self._out_dir_orig, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 5)
        se_coord_header = event_to_coord_header['SE']
        for row in se_mats_jc_rows:
            # 3600, 4600 are found from the transcripts intended for MXE
            self.assertIn(row[se_coord_header],
                          ['200', '1200', '3600', '4600', '15600'])
            if row[se_coord_header] == '200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[se_coord_header] == '1200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')
            if row[se_coord_header] == '15600':
                self.assertEqual(row['strand'], '-')

        se_mats_jcec_path = os.path.join(self._out_dir_orig,
                                         'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(se_mats_jcec_path))
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 5)
        for row in se_mats_jcec_rows:
            self.assertIn(row[se_coord_header],
                          ['200', '1200', '3600', '4600', '15600'])
            if row[se_coord_header] == '200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '149')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[se_coord_header] == '1200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')
            if row[se_coord_header] == '15600':
                self.assertEqual(row['strand'], '-')

        mxe_mats_jc_path = os.path.join(self._out_dir_orig, 'MXE.MATS.JC.txt')
        mxe_mats_jc_header, mxe_mats_jc_rows, error = (
            output_parser.parse_mats_jc(mxe_mats_jc_path))
        self.assertFalse(error)
        self._check_mxe_mats_jc_header(mxe_mats_jc_header)
        self.assertEqual(len(mxe_mats_jc_rows), 3)
        mxe_coord_header = event_to_coord_header['MXE']
        for row in mxe_mats_jc_rows:
            self.assertIn(row[mxe_coord_header], ['3200', '4200', '15200'])
            if row[mxe_coord_header] == '3200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '98')
            if row[mxe_coord_header] == '4200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')
            if row[mxe_coord_header] == '15200':
                self.assertEqual(row['strand'], '-')
                self.assertEqual(row['IJC_SAMPLE_1'], '1,3')
                self.assertEqual(row['SJC_SAMPLE_1'], '2,4')
                self.assertEqual(row['IJC_SAMPLE_2'], '5,7')
                self.assertEqual(row['SJC_SAMPLE_2'], '6,8')
                self.assertEqual(row['IncFormLen'], '79')
                self.assertEqual(row['SkipFormLen'], '98')

        mxe_mats_jcec_path = os.path.join(self._out_dir_orig,
                                          'MXE.MATS.JCEC.txt')
        mxe_mats_jcec_header, mxe_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(mxe_mats_jcec_path))
        self.assertFalse(error)
        self._check_mxe_mats_jcec_header(mxe_mats_jcec_header)
        self.assertEqual(len(mxe_mats_jcec_rows), 3)
        for row in mxe_mats_jcec_rows:
            self.assertIn(row[mxe_coord_header], ['3200', '4200', '15200'])
            if row[mxe_coord_header] == '3200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '149')
                self.assertEqual(row['SkipFormLen'], '149')
            if row[mxe_coord_header] == '4200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')
            if row[mxe_coord_header] == '15200':
                self.assertEqual(row['strand'], '-')
                self.assertEqual(row['IJC_SAMPLE_1'], '1,3')
                self.assertEqual(row['SJC_SAMPLE_1'], '4,8')
                self.assertEqual(row['IJC_SAMPLE_2'], '5,7')
                self.assertEqual(row['SJC_SAMPLE_2'], '12,16')
                self.assertEqual(row['IncFormLen'], '79')
                self.assertEqual(row['SkipFormLen'], '149')

        a5ss_mats_jc_path = os.path.join(self._out_dir_orig,
                                         'A5SS.MATS.JC.txt')
        a5ss_mats_jc_header, a5ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a5ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a5ss_mats_jc_header)
        self.assertEqual(len(a5ss_mats_jc_rows), 2)
        a5ss_coord_header = event_to_coord_header['A5SS']
        for row in a5ss_mats_jc_rows:
            self.assertEqual(row['strand'], '+')
            self.assertIn(row[a5ss_coord_header], ['6000', '7000'])
            if row[a5ss_coord_header] == '6000':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[a5ss_coord_header] == '7000':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')

        a5ss_mats_jcec_path = os.path.join(self._out_dir_orig,
                                           'A5SS.MATS.JCEC.txt')
        a5ss_mats_jcec_header, a5ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a5ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a5ss_mats_jcec_header)
        self.assertEqual(len(a5ss_mats_jcec_rows), 2)
        for row in a5ss_mats_jcec_rows:
            self.assertEqual(row['strand'], '+')
            self.assertIn(row[a5ss_coord_header], ['6000', '7000'])
            if row[a5ss_coord_header] == '6000':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '149')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[a5ss_coord_header] == '7000':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')

        a3ss_mats_jc_path = os.path.join(self._out_dir_orig,
                                         'A3SS.MATS.JC.txt')
        a3ss_mats_jc_header, a3ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a3ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a3ss_mats_jc_header)
        self.assertEqual(len(a3ss_mats_jc_rows), 2)
        a3ss_coord_header = event_to_coord_header['A3SS']
        for row in a3ss_mats_jc_rows:
            self.assertEqual(row['strand'], '+')
            self.assertIn(row[a3ss_coord_header], ['9200', '10200'])
            if row[a3ss_coord_header] == '9200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[a3ss_coord_header] == '10200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')

        a3ss_mats_jcec_path = os.path.join(self._out_dir_orig,
                                           'A3SS.MATS.JCEC.txt')
        a3ss_mats_jcec_header, a3ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a3ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a3ss_mats_jcec_header)
        self.assertEqual(len(a3ss_mats_jcec_rows), 2)
        for row in a3ss_mats_jcec_rows:
            self.assertEqual(row['strand'], '+')
            self.assertIn(row[a3ss_coord_header], ['9200', '10200'])
            if row[a3ss_coord_header] == '9200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '149')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[a3ss_coord_header] == '10200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')

        ri_mats_jc_path = os.path.join(self._out_dir_orig, 'RI.MATS.JC.txt')
        ri_mats_jc_header, ri_mats_jc_rows, error = output_parser.parse_mats_jc(
            ri_mats_jc_path)
        self.assertFalse(error)
        self._check_ri_mats_jc_header(ri_mats_jc_header)
        self.assertEqual(len(ri_mats_jc_rows), 2)
        ri_coord_header = event_to_coord_header['RI']
        for row in ri_mats_jc_rows:
            self.assertEqual(row['strand'], '+')
            self.assertIn(row[ri_coord_header], ['12000', '13000'])
            if row[ri_coord_header] == '12000':
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[ri_coord_header] == '13000':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')

        ri_mats_jcec_path = os.path.join(self._out_dir_orig,
                                         'RI.MATS.JCEC.txt')
        ri_mats_jcec_header, ri_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(ri_mats_jcec_path))
        self.assertFalse(error)
        self._check_ri_mats_jcec_header(ri_mats_jcec_header)
        self.assertEqual(len(ri_mats_jcec_rows), 2)
        for row in ri_mats_jcec_rows:
            self.assertEqual(row['strand'], '+')
            self.assertIn(row[ri_coord_header], ['12000', '13000'])
            if row[ri_coord_header] == '12000':
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '249')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[ri_coord_header] == '13000':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '1,1')

    def _check_row_id(self, row_id, coord, coord_to_id):
        if coord_to_id is not None:
            self.assertEqual(row_id, coord_to_id.get(coord))
        else:
            self.assertTrue(row_id.isdigit())

    def _check_row_ids_match(self, as_event, out_dir, expected_ids):
        name_templates_and_parsers = [
            ('fromGTF.{}.txt', output_parser.parse_from_gtf),
            ('JC.raw.input.{}.txt', output_parser.parse_jc_raw),
            ('JCEC.raw.input.{}.txt', output_parser.parse_jcec_raw),
        ]
        for name_template, parser in name_templates_and_parsers:
            path = os.path.join(out_dir, name_template.format(as_event))
            header, rows, error = parser(path)
            self.assertFalse(error)
            for i, row in enumerate(rows):
                self.assertEqual(row['ID'], expected_ids[i])

    def _check_results_fixed_shared(self, out_dir, start_coord_to_id=None):
        self._check_no_error_results()

        event_to_coord_header = self._as_event_to_start_coord_header()

        se_mats_jc_path = os.path.join(out_dir, 'SE.MATS.JC.txt')
        se_mats_jc_header, se_mats_jc_rows, error = output_parser.parse_mats_jc(
            se_mats_jc_path)
        self.assertFalse(error)
        self._check_se_mats_jc_header(se_mats_jc_header)
        self.assertEqual(len(se_mats_jc_rows), 5)
        se_coord_header = event_to_coord_header['SE']
        se_mats_row_ids = list()
        for row in se_mats_jc_rows:
            se_mats_row_ids.append(row['ID'])
            self._check_row_id(row['ID'], row[se_coord_header],
                               start_coord_to_id)
            self.assertIn(row[se_coord_header],
                          ['200', '1200', '3600', '4600', '15600'])
            if row[se_coord_header] == '200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[se_coord_header] == '1200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
            if row[se_coord_header] == '15600':
                self.assertEqual(row['strand'], '-')

        se_mats_jcec_path = os.path.join(out_dir, 'SE.MATS.JCEC.txt')
        se_mats_jcec_header, se_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(se_mats_jcec_path))
        self.assertFalse(error)
        self._check_se_mats_jcec_header(se_mats_jcec_header)
        self.assertEqual(len(se_mats_jcec_rows), 5)
        for i, row in enumerate(se_mats_jcec_rows):
            self.assertEqual(row['ID'], se_mats_row_ids[i])
            self.assertIn(row[se_coord_header],
                          ['200', '1200', '3600', '4600', '15600'])
            if row[se_coord_header] == '200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '149')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[se_coord_header] == '1200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
            if row[se_coord_header] == '15600':
                self.assertEqual(row['strand'], '-')

        self._check_row_ids_match('SE', out_dir, se_mats_row_ids)

        mxe_mats_jc_path = os.path.join(out_dir, 'MXE.MATS.JC.txt')
        mxe_mats_jc_header, mxe_mats_jc_rows, error = (
            output_parser.parse_mats_jc(mxe_mats_jc_path))
        self.assertFalse(error)
        self._check_mxe_mats_jc_header(mxe_mats_jc_header)
        self.assertEqual(len(mxe_mats_jc_rows), 3)
        mxe_coord_header = event_to_coord_header['MXE']
        mxe_mats_row_ids = list()
        for row in mxe_mats_jc_rows:
            mxe_mats_row_ids.append(row['ID'])
            self._check_row_id(row['ID'], row[mxe_coord_header],
                               start_coord_to_id)
            self.assertIn(row[mxe_coord_header], ['3200', '4200', '15200'])
            if row[mxe_coord_header] == '3200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '98')
            if row[mxe_coord_header] == '4200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
            if row[mxe_coord_header] == '15200':
                self.assertEqual(row['strand'], '-')
                self.assertEqual(row['IJC_SAMPLE_1'], '1,3')
                self.assertEqual(row['SJC_SAMPLE_1'], '2,4')
                self.assertEqual(row['IJC_SAMPLE_2'], '9,11')
                self.assertEqual(row['SJC_SAMPLE_2'], '10,12')
                self.assertEqual(row['IncFormLen'], '79')
                self.assertEqual(row['SkipFormLen'], '98')

        mxe_mats_jcec_path = os.path.join(out_dir, 'MXE.MATS.JCEC.txt')
        mxe_mats_jcec_header, mxe_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(mxe_mats_jcec_path))
        self.assertFalse(error)
        self._check_mxe_mats_jcec_header(mxe_mats_jcec_header)
        self.assertEqual(len(mxe_mats_jcec_rows), 3)
        for i, row in enumerate(mxe_mats_jcec_rows):
            self.assertEqual(row['ID'], mxe_mats_row_ids[i])
            self.assertIn(row[mxe_coord_header], ['3200', '4200', '15200'])
            if row[mxe_coord_header] == '3200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '149')
                self.assertEqual(row['SkipFormLen'], '149')
            if row[mxe_coord_header] == '4200':
                self.assertEqual(row['strand'], '+')
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')
            if row[mxe_coord_header] == '15200':
                self.assertEqual(row['strand'], '-')
                self.assertEqual(row['IJC_SAMPLE_1'], '1,3')
                self.assertEqual(row['SJC_SAMPLE_1'], '4,8')
                self.assertEqual(row['IJC_SAMPLE_2'], '9,11')
                self.assertEqual(row['SJC_SAMPLE_2'], '20,24')
                self.assertEqual(row['IncFormLen'], '79')
                self.assertEqual(row['SkipFormLen'], '149')

        self._check_row_ids_match('MXE', out_dir, mxe_mats_row_ids)

        a5ss_mats_jc_path = os.path.join(out_dir, 'A5SS.MATS.JC.txt')
        a5ss_mats_jc_header, a5ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a5ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a5ss_mats_jc_header)
        self.assertEqual(len(a5ss_mats_jc_rows), 2)
        a5ss_coord_header = event_to_coord_header['A5SS']
        a5ss_mats_row_ids = list()
        for row in a5ss_mats_jc_rows:
            self.assertEqual(row['strand'], '+')
            a5ss_mats_row_ids.append(row['ID'])
            self._check_row_id(row['ID'], row[a5ss_coord_header],
                               start_coord_to_id)
            self.assertIn(row[a5ss_coord_header], ['6000', '7000'])
            if row[a5ss_coord_header] == '6000':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[a5ss_coord_header] == '7000':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')

        a5ss_mats_jcec_path = os.path.join(out_dir, 'A5SS.MATS.JCEC.txt')
        a5ss_mats_jcec_header, a5ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a5ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a5ss_mats_jcec_header)
        self.assertEqual(len(a5ss_mats_jcec_rows), 2)
        for i, row in enumerate(a5ss_mats_jcec_rows):
            self.assertEqual(row['strand'], '+')
            self.assertEqual(row['ID'], a5ss_mats_row_ids[i])
            self.assertIn(row[a5ss_coord_header], ['6000', '7000'])
            if row[a5ss_coord_header] == '6000':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '149')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[a5ss_coord_header] == '7000':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')

        self._check_row_ids_match('A5SS', out_dir, a5ss_mats_row_ids)

        a3ss_mats_jc_path = os.path.join(out_dir, 'A3SS.MATS.JC.txt')
        a3ss_mats_jc_header, a3ss_mats_jc_rows, error = (
            output_parser.parse_mats_jc(a3ss_mats_jc_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jc_header(a3ss_mats_jc_header)
        self.assertEqual(len(a3ss_mats_jc_rows), 2)
        a3ss_coord_header = event_to_coord_header['A3SS']
        a3ss_mats_row_ids = list()
        for row in a3ss_mats_jc_rows:
            self.assertEqual(row['strand'], '+')
            a3ss_mats_row_ids.append(row['ID'])
            self._check_row_id(row['ID'], row[a3ss_coord_header],
                               start_coord_to_id)
            self.assertIn(row[a3ss_coord_header], ['9200', '10200'])
            if row[a3ss_coord_header] == '9200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[a3ss_coord_header] == '10200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')

        a3ss_mats_jcec_path = os.path.join(out_dir, 'A3SS.MATS.JCEC.txt')
        a3ss_mats_jcec_header, a3ss_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(a3ss_mats_jcec_path))
        self.assertFalse(error)
        self._check_a35ss_mats_jcec_header(a3ss_mats_jcec_header)
        self.assertEqual(len(a3ss_mats_jcec_rows), 2)
        for i, row in enumerate(a3ss_mats_jcec_rows):
            self.assertEqual(row['strand'], '+')
            self.assertEqual(row['ID'], a3ss_mats_row_ids[i])
            self.assertIn(row[a3ss_coord_header], ['9200', '10200'])
            if row[a3ss_coord_header] == '9200':
                self.assertEqual(row['IJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '149')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[a3ss_coord_header] == '10200':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')

        self._check_row_ids_match('A3SS', out_dir, a3ss_mats_row_ids)

        ri_mats_jc_path = os.path.join(out_dir, 'RI.MATS.JC.txt')
        ri_mats_jc_header, ri_mats_jc_rows, error = output_parser.parse_mats_jc(
            ri_mats_jc_path)
        self.assertFalse(error)
        self._check_ri_mats_jc_header(ri_mats_jc_header)
        self.assertEqual(len(ri_mats_jc_rows), 2)
        ri_coord_header = event_to_coord_header['RI']
        ri_mats_row_ids = list()
        for row in ri_mats_jc_rows:
            self.assertEqual(row['strand'], '+')
            ri_mats_row_ids.append(row['ID'])
            self._check_row_id(row['ID'], row[ri_coord_header],
                               start_coord_to_id)
            self.assertIn(row[ri_coord_header], ['12000', '13000'])
            if row[ri_coord_header] == '12000':
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '98')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[ri_coord_header] == '13000':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')

        ri_mats_jcec_path = os.path.join(out_dir, 'RI.MATS.JCEC.txt')
        ri_mats_jcec_header, ri_mats_jcec_rows, error = (
            output_parser.parse_mats_jcec(ri_mats_jcec_path))
        self.assertFalse(error)
        self._check_ri_mats_jcec_header(ri_mats_jcec_header)
        self.assertEqual(len(ri_mats_jcec_rows), 2)
        for i, row in enumerate(ri_mats_jcec_rows):
            self.assertEqual(row['strand'], '+')
            self.assertEqual(row['ID'], ri_mats_row_ids[i])
            self.assertIn(row[ri_coord_header], ['12000', '13000'])
            if row[ri_coord_header] == '12000':
                self.assertEqual(row['IJC_SAMPLE_1'], '2,2')
                self.assertEqual(row['SJC_SAMPLE_1'], '1,1')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,1')
                self.assertEqual(row['IncFormLen'], '249')
                self.assertEqual(row['SkipFormLen'], '49')
            if row[ri_coord_header] == '13000':
                self.assertEqual(row['IJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['SJC_SAMPLE_1'], '0,0')
                self.assertEqual(row['IJC_SAMPLE_2'], '0,2')
                self.assertEqual(row['SJC_SAMPLE_2'], '0,0')

        self._check_row_ids_match('RI', out_dir, ri_mats_row_ids)

    def _check_results_fixed(self):
        self._check_results_fixed_shared(self._out_dir_fixed)

    def _check_results_id_strings(self):
        start_coord_to_id = self._start_coord_to_id()
        self._check_results_fixed_shared(self._out_dir_id_strings,
                                         start_coord_to_id=start_coord_to_id)


if __name__ == '__main__':
    unittest.main(verbosity=2)
