import os.path
import subprocess
import unittest

import tests.bam
import tests.gtf
import tests.test_config as test_config


class BaseTest(unittest.TestCase):
    def setUp(self):
        self._rmats_command = test_config.RMATS_COMMAND
        self._rmats_return_code = None

    def _run_test(self):
        self._run_rmats()
        self._check_results()

    def _command_output_dir(self):
        self.fail('derived classes need to override _command_output_dir')

    def _check_results(self):
        self.fail('derived classes need to override _check_results')

    def _rmats_arguments(self):
        self.fail('derived classes need to override _rmats_arguments')

    def _get_stdout_file_name(self):
        return os.path.join(self._command_output_dir(), 'rmats_stdout')

    def _get_stderr_file_name(self):
        return os.path.join(self._command_output_dir(), 'rmats_stderr')

    def _run_rmats(self):
        command = [self._rmats_command] + self._rmats_arguments()
        if self._command_output_dir():
            stdout_file = self._get_stdout_file_name()
            stderr_file = self._get_stderr_file_name()
            with open(stdout_file, 'wb') as out_f_handle:
                with open(stderr_file, 'wb') as err_f_handle:
                    result = subprocess.run(command,
                                            check=False,
                                            stdout=out_f_handle,
                                            stderr=err_f_handle)
        else:
            result = subprocess.run(command, check=False)

        self._rmats_return_code = result.returncode

    def _write_bams(self, bams, text_path):
        bam_paths = list()
        for bam in bams:
            error = bam.write()
            self.assertFalse(error)
            bam_paths.append(bam.path)

        with open(text_path, 'wt') as f_handle:
            f_handle.write(','.join(bam_paths))

    def _check_no_error_results(self):
        self.assertEqual(self._rmats_return_code, 0)

        command_stderr_file_name = self._get_stderr_file_name()
        with open(command_stderr_file_name, 'rt') as err_f_h:
            err_lines = err_f_h.readlines()

        self.assertEqual(err_lines, list())

    def _create_gtf_from_transcripts(self, gtf_path, exons_by_transcript):
        gtf = tests.gtf.GTF()
        gtf.path = gtf_path

        transcripts = list()
        for i, exons_for_transcript in enumerate(exons_by_transcript):
            transcript = tests.gtf.Transcript()
            transcript.chromosome = '1'
            transcript.strand = '+'
            transcript.gene_id = tests.util.gene_id_str(1)
            transcript.gene_name = tests.util.gene_name_str(1)
            transcript.transcript_id = tests.util.transcript_id_str(i)
            transcript.exons = exons_for_transcript
            transcripts.append(transcript)

        gtf.transcripts = transcripts
        error = gtf.write()
        self.assertFalse(error)
        return gtf

    def _create_bam_from_paired_read_coords(self,
                                            bam_path,
                                            chromosome_length,
                                            read_length,
                                            paired_read_coords,
                                            clip_length=None):
        bam = tests.bam.BAM()
        bam.path = bam_path

        bam_reads = list()
        for i, coord_pair in enumerate(paired_read_coords):
            read_1_coords, read_2_coords = coord_pair
            paired_read_1 = tests.bam.Read()
            paired_read_1.ref_seq_name = '1'  # chromosome
            paired_read_1.ref_seq_len = chromosome_length
            paired_read_1.template_name = tests.util.template_name_str([i])
            paired_read_2 = tests.bam.Read()
            error = tests.bam.set_read_pair_from_intervals(
                paired_read_1,
                paired_read_2,
                read_1_coords,
                read_2_coords,
                read_length,
                clip_length=clip_length)
            self.assertFalse(error)
            bam_reads.extend([paired_read_1, paired_read_2])

        bam.reads = bam_reads
        return bam

    def _check_se_mats_jc_header(self, header):
        self.assertEqual(header, [
            'ID', 'GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base',
            'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES',
            'downstreamEE', 'ID', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1',
            'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncFormLen', 'SkipFormLen',
            'PValue', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference'
        ])

    def _check_se_mats_jcec_header(self, header):
        self._check_se_mats_jc_header(header)
