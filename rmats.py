#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
# @file lite2.py
# @brief
# @author Zhijie Xie
# @date 2015-11-27

from __future__ import print_function

import sys
import os
import os.path
import argparse
import subprocess
import shutil
import time
from datetime import datetime
from rmatspipeline import run_pipe


VERSION = 'v4.3.0'
USAGE = '''%(prog)s [options]'''
pipe_tasks = set(['prep', 'post', 'both',])


def getstatusoutput(cmd):
    """ behave like commands.getstatusoutput which moved to
    subprocess.getstatusoutput in Python 3.
    Implmented with subprocess.check_output which is available
    in both Python 2 and 3.
    """
    status = 0
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        status = e.returncode
        output = e.output

    output = output.decode()
    if output[-1] == '\n':
        output = output[:-1]

    return (status, output)


def doSTARMapping(args): ## do STAR mapping
    fastqs = [args.s1, args.s2,]
    bams = [[], [],]

    for i in range(len(fastqs)):
        if fastqs[i] != '':
            sample = [pair.split(':') for pair in fastqs[i].split(',')]
            print("mapping the first sample")
            for rr, pair in enumerate(sample):
                map_folder = os.path.join(
                    args.tmp, '{}_bam{}_{}'.format(args.prep_prefix, i+1, rr+1))

                if os.path.exists(map_folder):
                    if os.path.isdir(map_folder):
                        os.rmdir(map_folder)
                    else:
                        os.unlink(map_folder)

                os.makedirs(map_folder)
                cmd = 'STAR --chimSegmentMin 2 --outFilterMismatchNmax 3'
                cmd += ' --twopassMode Basic'
                if not args.allow_clipping:
                    cmd += ' --alignEndsType EndToEnd'

                cmd += ' --runThreadN ' + str(max([4, args.nthread])) + ' --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate '
                cmd += '--alignSJDBoverhangMin ' + str(args.tophatAnchor) + ' --alignIntronMax 299999 --genomeDir ' + args.bIndex + ' --sjdbGTFfile ' + args.gtf
                cmd += ' --outFileNamePrefix ' + map_folder + '/ --readFilesIn '
                cmd += ' '.join(pair)
                if pair[0].endswith('.gz'):
                    cmd += ' --readFilesCommand zcat'
                status,output = getstatusoutput(cmd)
                print("mapping sample_%d, %s is done with status %s" % (i, ' '.join(pair), status))
                if (int(status)!=0): ## it did not go well
                    print("error in mapping sample_%d, %s: %s" % (i, ' '.join(pair),status))
                    print("error detail: %s" % output)
                    raise Exception()
                print(output)
                bams[i].append(os.path.join(map_folder, 'Aligned.sortedByCoord.out.bam'))

    return ','.join(bams[0]), ','.join(bams[1])
##### end of doSTARMapping ####


def get_args():
    """Supplies all the neccessary arguments to the argparse package, along with appropriate help, defaults, destinations, and choices.
    The function itself takes no arguments.
    Unless rMATS is called in stat mode, exits with appropriate errors if any of: sequence files, gtf, or readlength arguments are missing.
    If an output directory and/or a temporary directory aren't supplied, exits with appropriate errors.
    Creates output directory. Cleans user supplied bam or fastq filenames to remove trailing newlines, spaces, and commas.
    If FASTQs are supplied but not BAMs, aligns FASTQs using STAR and sets the resultant BAM file locations as BAM file arguments.
    If tstat is not supplied, sets the tstat argument equal to the nthread argument or the nthread default, if nthread is also not supplied..

    Returns an arg object.
    """
    parser = argparse.ArgumentParser(usage=USAGE)

    group1 = parser.add_mutually_exclusive_group()
    group2 = parser.add_mutually_exclusive_group()

    parser.add_argument('--version', action='version', version=VERSION)
    parser.add_argument('--gtf', action='store',
                        help='An annotation of genes and transcripts in GTF format', dest='gtf')

    group1.add_argument('--b1', action='store', default='',
                        help='A text file containing a comma separated list of the BAM files for sample_1. (Only if using BAM)', dest='b1')
    group2.add_argument('--b2', action='store', default='',
                        help='A text file containing a comma separated list of the BAM files for sample_2. (Only if using BAM)', dest='b2')
    group1.add_argument('--s1', action='store', default='',
                        help='A text file containing a comma separated list of the FASTQ files for sample_1. If using paired reads the format is ":" to separate pairs and "," to separate replicates. (Only if using fastq)', dest='s1')
    group2.add_argument('--s2', action='store', default='',
                        help='A text file containing a comma separated list of the FASTQ files for sample_2. If using paired reads the format is ":" to separate pairs and "," to separate replicates. (Only if using fastq)', dest='s2')

    parser.add_argument('--od', action='store',
                        help='The directory for final output from the post step', dest='od')
    parser.add_argument('--tmp', action='store',
                        help='The directory for intermediate output such as ".rmats" files from the prep step', dest='tmp')
    parser.add_argument('-t', action='store', default='paired',
                        choices=['paired', 'single'],
                        help='Type of read used in the analysis: either "paired" for paired-end data or "single" for single-end data. Default: %(default)s', dest='readtype')
    parser.add_argument('--libType', action='store', default='fr-unstranded',
                        choices=['fr-unstranded', 'fr-firststrand',
                                 'fr-secondstrand',],
                        help='Library type. Use fr-firststrand or fr-secondstrand for strand-specific data. Only relevant to the prep step, not the post step. Default: %(default)s', dest='dt')
    parser.add_argument('--readLength', action='store', type=int,
                        help='The length of each read. Required parameter, with the value set according to the RNA-seq read length', dest='readLength')
    parser.add_argument('--variable-read-length', action='store_true',
                        help='Allow reads with lengths that differ from --readLength to be processed. --readLength will still be used to determine IncFormLen and SkipFormLen',
                        dest='variable_read_length')
    parser.add_argument('--anchorLength', action='store', type=int, default=1,
                        help='The "anchor length" or "overhang length" used when counting the number of reads spanning splice junctions. A minimum number of "anchor length" nucleotides must be mapped to each end of a given splice junction. The minimum value is 1 and the default value is set to %(default)s to make use of all possible splice junction reads.', dest='anchorLength')
    parser.add_argument('--tophatAnchor', action='store', type=int, default=1,
                        help='The "anchor length" or "overhang length" used in the aligner. At least "anchor length" nucleotides must be mapped to each end of a given splice junction. The default is %(default)s. (Only if using fastq)', dest='tophatAnchor')
    parser.add_argument('--bi', action='store', default='',
                        help='The directory name of the STAR binary indices (name of the directory that contains the suffix array file). (Only if using fastq)', dest='bIndex')
    parser.add_argument('--nthread', action='store', type=int, default=1,
                        help='The number of threads. The optimal number of threads should be equal to the number of CPU cores. Default: %(default)s', dest='nthread')

    parser.add_argument('--tstat', action='store', type=int,
                        help='The number of threads for the statistical model. If not set then the value of --nthread is used', dest='tstat')
    parser.add_argument('--cstat', action='store', type=float, default=0.0001,
                        help='The cutoff splicing difference. The cutoff used in the null hypothesis test for differential alternative splicing. The default is 0.0001 for 0.01%% difference. Valid: 0 <= cutoff < 1. Does not apply to the paired stats model', dest='cstat')

    parser.add_argument('--task', action='store', default='both',
                        choices=['prep', 'post', 'both', 'inte', 'stat'],
                        help='Specify which step(s) of rMATS-turbo to run. Default: %(default)s. prep: preprocess BAM files and generate .rmats files. post: load .rmats files into memory, detect and count alternative splicing events, and calculate P value (if not --statoff). both: prep + post. inte (integrity): check that the BAM filenames recorded by the prep task(s) match the BAM filenames for the current command line. stat: run statistical test on existing output files', dest='task')
    parser.add_argument('--statoff', action='store_false',
                        help='Skip the statistical analysis', dest='stat')
    parser.add_argument('--paired-stats', action='store_true',
                        help='Use the paired stats model', dest='paired_stats')
    parser.add_argument('--darts-model', action='store_true',
                        help='Use the DARTS statistical model', dest='darts_model')
    parser.add_argument('--darts-cutoff', action='store', type=float, default=0.05,
                        help='The cutoff of delta-PSI in the DARTS model. The output posterior probability is P(abs(delta_psi) > cutoff). The default is %(default)s', dest='darts_cutoff')

    parser.add_argument('--novelSS', action='store_true',
                        help='Enable detection of novel splice sites (unannotated splice sites). Default is no detection of novel splice sites', dest='novelSS')
    parser.add_argument('--mil', action='store', type=int, default=50,
                        help='Minimum Intron Length. Only impacts --novelSS behavior. Default: %(default)s', dest='mil')
    parser.add_argument('--mel', action='store', type=int, default=500,
                        help='Maximum Exon Length. Only impacts --novelSS behavior. Default: %(default)s', dest='mel')
    parser.add_argument('--allow-clipping', action='store_true',
                        help='Allow alignments with soft or hard clipping to be used',
                        dest='allow_clipping')
    parser.add_argument('--fixed-event-set', action='store', help='A directory containing fromGTF.[AS].txt files to be used instead of detecting a new set of events')
    parser.add_argument('--individual-counts', action='store_true',
                        help='Output individualCounts.[AS_Event].txt files and add the individual count columns to [AS_Event].MATS.JC.txt',
                        dest='individual_counts')
    # The help text for --imbalance ratio is not added to the parser
    # since the parameter is only intended for internal use and it
    # defaults to no filtering.
    #
    # Filter events where the ratio of upstream junction reads to
    # downstream junction reads (or downstream to upstream) exceeds
    # --imbalance-ratio. The events are filtered before running the
    # stats model so that the FDR is based on the filtered
    # events. If not specified then no events are filtered.
    # Requires --individual-counts
    parser.add_argument('--imbalance-ratio', type=float,
                        help=argparse.SUPPRESS,
                        dest='imbalance_ratio')
    parser.add_argument('--drop-zero-read-replicates-for-stat', action='store_true', help='When calculating the pvalue, for each event ignore replicates that do not have any reads for that event')


    args = parser.parse_args()

    if args.task != 'stat':
        if args.b1 == '' and args.b2 == '' and args.s1 == '' and args.s2 == '':
            sys.exit('ERROR: BAM/FASTQ required. Please check b1, b2, s1 and s2.')
        if args.gtf == None:
            sys.exit('ERROR: GTF file required. Please check --gtf.')
        if args.readLength is None:
            sys.exit('ERROR: --readLength is required. An average or median'
                     ' --readLength can be used in combination with'
                     ' --variable-read-length when the reads do not have the'
                     ' same length.')
        args.junctionLength = 2 * (args.readLength - args.anchorLength)

    if args.od == None or args.tmp == None:
        sys.exit('ERROR: output folder and temporary folder required. Please check --od and --tmp.')
    if (args.s1 != '' or args.s2 != '') and args.bIndex == '':
        sys.exit('ERROR: STAR binary indexes required. Please check --bi.')
    if args.imbalance_ratio is not None and not args.individual_counts:
        sys.exit('ERROR: --imbalance-ratio requires --individual-counts')

    if len(args.b1) > 0:
        with open(args.b1, 'r') as fp:
            args.b1 = fp.read().strip().strip(',')
    if len(args.b2) > 0:
        with open(args.b2, 'r') as fp:
            args.b2 = fp.read().strip().strip(',')
    if len(args.s1) > 0:
        with open(args.s1, 'r') as fp:
            args.s1 = fp.read().strip().strip(',')
    if len(args.s2) > 0:
        with open(args.s2, 'r') as fp:
            args.s2 = fp.read().strip().strip(',')

    create_output_dirs(args)
    args.prep_prefix = claim_prep_prefix(args.task, args.tmp)

    if args.b1 == '' and args.b2 == '' and (args.s1 != '' or args.s2 != ''):
        args.b1, args.b2 = doSTARMapping(args)

    args.bams = ','.join([args.b1, args.b2]).strip(',')

    dt_map = {'fr-unstranded':0, 'fr-firststrand':1, 'fr-secondstrand':2}
    args.dt = dt_map[args.dt]

    if args.tstat is None:
        args.tstat = args.nthread

    return args


def check_integrity(input_bams_string, tmp_dir):
    """
    Purpose: Iterates over the supplied string of bam filenames and checks every file in tmp_dir
    to ensure every bam filename has exactly one prep file. Exits with appropriate errors if
    there are one or more of the following: duplicate bam files, bam files with no prep,
    bam files with multiple preps, or prep files with no corresponding bam. Otherwise, prints 'Ok.'

    Positional arguments:
    First: input_bams_string - A comma-delimited string of all the input bam filenames.
    Second: tmp_dir - The location of the current rMATS instance's temporary directory, ostensibly containing the prep files for the input bams.
    """

    input_bams = input_bams_string.split(',')
    duplicate_input_bams = list()
    prep_count_by_bam = dict()
    for input_bam in input_bams:
        if input_bam in prep_count_by_bam:
            duplicate_input_bams.append(input_bam)
        else:
            prep_count_by_bam[input_bam] = 0

    all_files = [os.path.join(root, name)
                 for root, dirs, files in os.walk(tmp_dir)
                 for name in files if name.endswith('.rmats')]

    bams_only_in_prep = list()
    for name in all_files:
        with open(name, 'r') as fp:
            prep_bams = fp.readline().strip().split(',')

        for prep_bam in prep_bams:
            if prep_bam not in prep_count_by_bam:
                bams_only_in_prep.append(prep_bam)
                prep_count_by_bam[prep_bam] = 1
            else:
                prep_count_by_bam[prep_bam] += 1

    missing_preps = list()
    duplicate_preps = list()
    for bam, prep_count in prep_count_by_bam.items():
        if prep_count == 0:
            missing_preps.append(bam)
        elif prep_count > 1:
            duplicate_preps.append(bam)

    def update_errors(errors, description, values):
        if not values:
            return

        errors.append('{}:\n\t{}'.format(description, '\n\t'.join(values)))

    errors = list()
    update_errors(errors, 'duplicate input bam files', duplicate_input_bams)
    update_errors(errors,
                  'bam files with multiple associations with prep output',
                  duplicate_preps)
    update_errors(errors,
                  'input bam files with no associated prep output',
                  missing_preps)
    update_errors(errors,
                  'bam files not in input but associated with prep output',
                  bams_only_in_prep)

    if errors:
        sys.exit('\n'.join(errors))

    print('Ok.')


def check_if_has_counts(counts_file_path):
    has_sample_1_counts = False
    has_sample_2_counts = False
    has_replicates = False
    with open(counts_file_path, 'rt') as f_handle:
        for i, line in enumerate(f_handle):
            if i == 0:
                continue  # skip header line

            values = line.strip().split('\t')
            inc_sample_1_vs = values[1]
            inc_sample_2_vs = values[3]
            sample_1_counts = inc_sample_1_vs.split(',')
            sample_2_counts = inc_sample_2_vs.split(',')
            has_sample_1_counts = inc_sample_1_vs != ''
            has_sample_2_counts = inc_sample_2_vs != ''
            has_replicates = ((len(sample_1_counts) > 1)
                              or (len(sample_2_counts) > 1))
            break  # only check the first row

    return {
        'has_sample_1_counts': has_sample_1_counts,
        'has_sample_2_counts': has_sample_2_counts,
        'has_replicates': has_replicates
    }


def filter_countfile(fn):
    """TODO: Docstring for filter_countfile.
    :returns: TODO

    """
    with open(fn, 'r+') as fp:
        data = fp.readlines()
        header = data[0]
        rest = [header,]

        for line in data[1:]:
            eles = line.split('\t')
            if len(eles) != 7:
                sys.exit('ERROR: corrupted read count file: %s' % (fn))
            sum_ic_1 = 0
            sum_sc_1 = 0
            sum_ic_2 = 0
            sum_sc_2 = 0
            incv1 = list(map(int, eles[1].split(',')))
            skpv1 = list(map(int, eles[2].split(',')))
            incv2 = list(map(int, eles[3].split(',')))
            skpv2 = list(map(int, eles[4].split(',')))
            inc_len = int(eles[5])
            skp_len = int(eles[6])

            for i in range(len(incv1)):
                sum_ic_1 += incv1[i]
                sum_sc_1 += skpv1[i]
            for i in range(len(incv2)):
                sum_ic_2 += incv2[i]
                sum_sc_2 += skpv2[i]

            if (sum_ic_1 + sum_sc_1 > 0 and\
                    sum_ic_2 + sum_sc_2 > 0 and\
                    (sum_ic_1 != 0 or sum_ic_2 != 0) and\
                    (sum_sc_1 != 0 or sum_sc_2 != 0) and\
                    inc_len != 0 and skp_len != 0):
                rest.append(line)

        fp.seek(0)
        fp.writelines(rest)
        fp.truncate()

    return


def process_counts(istat, tstat, counttype, ase, cstat, od, od_tmp, stat,
                   paired_stats, use_darts_model, darts_cutoff,
                   python_executable, root_dir, imbalance_ratio,
                   drop_zero_read_replicates_for_stat):
    """TODO: Docstring for process_counts.
    :returns: TODO

    """
    from_gtf_path = '%s/fromGTF.%s.txt' % (od, ase)
    if not os.path.exists(from_gtf_path):
        print('WARNING: Cannot find {}. Unable to produce final output files'
              ' for {} {}.'.format(from_gtf_path, ase, counttype),
              file=sys.stderr)
        return

    if not os.path.exists(istat):
        print('WARNING: Cannot find {}. Unable to produce final output files'
              ' for {} {}.'.format(istat, ase, counttype),
              file=sys.stderr)
        return

    indiv_counts_file_name = os.path.join(
        od, 'individualCounts.{}.txt'.format(ase))
    indiv_counts_temp_file_name = '{}.tmp'.format(indiv_counts_file_name)
    has_indiv_counts = os.path.exists(indiv_counts_file_name)

    if stat:
        has_counts_result = check_if_has_counts(istat)
        has_sample_1_counts = has_counts_result['has_sample_1_counts']
        has_sample_2_counts = has_counts_result['has_sample_2_counts']
        has_replicates = has_counts_result['has_replicates']
        if has_sample_1_counts and has_sample_2_counts:
            filter_countfile(istat)
        elif has_sample_1_counts or has_sample_2_counts:
            print('WARNING: Statistical step is skipped for {} {} because only'
                  ' one group is involved'.format(ase, counttype),
                  file=sys.stderr)
            stat = False

    if imbalance_ratio is not None:
        filter_countfile_by_imbalance_ratio(
            istat, indiv_counts_file_name, indiv_counts_temp_file_name,
            imbalance_ratio, ase)

    sec_tmp = os.path.join(od_tmp, '%s_%s' % (counttype, ase))
    if os.path.exists(sec_tmp):
        if os.path.isdir(sec_tmp):
            shutil.rmtree(sec_tmp)
        else:
            os.unlink(sec_tmp)
    os.mkdir(sec_tmp)
    ostat = os.path.join(sec_tmp, 'rMATS_result_%s.txt')

    FNULL = open(os.devnull, 'w')
    resfp = open(ostat % (''), 'w')
    finfn = os.path.join(od, '%s.MATS.%s.txt' % (ase, counttype))
    ostat_id = ostat % ('ID')
    ostat_inp = ostat % ('INP')
    ostat_il = ostat % ('I-L')
    ostat_pv = ostat % ('P-V')
    ostat_fdr = ostat % ('FDR')
    ostat_paired = ostat % ('paired')
    ostat_darts = ostat % ('darts')

    rmats_c = os.path.join(root_dir, 'rMATS_C/rMATSexe')
    paired_model = os.path.join(root_dir, 'rMATS_R/paired_model.R')
    darts_model = os.path.join(root_dir, 'rMATS_R/darts_model.R')
    pas_out = os.path.join(root_dir, 'rMATS_P/paste.py')
    inc_lvl = os.path.join(root_dir, 'rMATS_P/inclusion_level.py')
    fdr_cal = os.path.join(root_dir, 'rMATS_P/FDR.py')
    join_2f = os.path.join(root_dir, 'rMATS_P/joinFiles.py')

    # Calculate inclusion levels
    subprocess.call([python_executable, pas_out, '-i', istat, '--o1', ostat_id, '--o2', ostat_inp,], stdout=FNULL)
    subprocess.call([python_executable, inc_lvl, ostat_inp, ostat_il,], stdout=FNULL)

    # Calculate PValue and FDR
    if stat:
        if paired_stats:
            # PAIRADISE writes a status file to the working directory.
            # Use the cwd kwarg of subprocess.call to set which
            # directory the status file gets written to.
            # This requires passing absolute path names to paired_model.R
            abs_istat = os.path.abspath(istat)
            abs_ostat_fdr = os.path.abspath(ostat_fdr)
            paired_command = ['Rscript', paired_model, abs_istat, str(tstat),
                              abs_ostat_fdr]
            with open(ostat_paired, 'wb') as paired_fp:
                paired_model_return_code = subprocess.call(
                    paired_command, stdout=paired_fp, stderr=subprocess.STDOUT,
                    cwd=sec_tmp)

            if paired_model_return_code != 0:
                print('error in paired model', file=sys.stderr)
                print_file_to_stderr(ostat_paired)
        elif use_darts_model:
            has_replicates_str = 'true' if has_replicates else 'false'
            darts_command = ['Rscript', darts_model, istat, ostat_pv,
                             str(tstat), str(darts_cutoff), has_replicates_str]
            with open(ostat_darts, 'wb') as darts_fp:
                darts_return_code = subprocess.call(
                    darts_command, stdout=darts_fp, stderr=subprocess.STDOUT)

            if darts_return_code != 0:
                print('error running darts', file=sys.stderr)
                print_file_to_stderr(ostat_darts)
            else:
                write_fdr_file_from_darts_output(istat, ostat_pv, ostat_fdr)
        else:
            stat_in_file = istat
            stat_out_file = ostat_pv
            if drop_zero_read_replicates_for_stat:
                stat_in_file = os.path.join(sec_tmp, 'cleaned_stat_input.txt')
                stat_out_file = os.path.join(sec_tmp, 'cleaned_stat_output.txt')
                drop_zero_read_replicates(istat, stat_in_file)

            subprocess.call([rmats_c, '-i', stat_in_file, '-t', str(tstat), '-o', stat_out_file, '-c', str(cstat),], stdout=FNULL)
            if drop_zero_read_replicates_for_stat:
                copy_last_column_and_append(istat, stat_out_file, ostat_pv)

            subprocess.call([python_executable, fdr_cal, ostat_pv, ostat_fdr,], stdout=FNULL)
    else:
        append_columns_with_defaults(istat, ostat_fdr, ['PValue', 'FDR'], ['NA', 'NA'])

    # Combine into final output
    subprocess.call(['paste', ostat_fdr, ostat_il,], stdout=resfp)
    subprocess.call([python_executable, join_2f, from_gtf_path, resfp.name, '0', '0', finfn,], stdout=FNULL)
    if has_indiv_counts:
        append_individual_counts(indiv_counts_file_name, finfn, indiv_counts_temp_file_name)

    FNULL.close()
    resfp.close()


def drop_zero_read_replicates(orig_file_path, filtered_file_path):
    COUNT_HEADERS = ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2',
                     'SJC_SAMPLE_2']
    with open(orig_file_path, 'rt') as in_handle:
        with open(filtered_file_path, 'wt') as out_handle:
            for line_i, line in enumerate(in_handle):
                columns = line.rstrip('\n').split('\t')
                if line_i == 0:
                    headers = columns
                    out_handle.write(line)
                    continue

                row = dict(zip(headers, columns))
                for count_header in COUNT_HEADERS:
                    row[count_header] = parse_comma_ints(row[count_header])

                group_1_should_drop = find_zero_read_replicates(
                    row['IJC_SAMPLE_1'], row['SJC_SAMPLE_1'])
                group_2_should_drop = find_zero_read_replicates(
                    row['IJC_SAMPLE_2'], row['SJC_SAMPLE_2'])

                filtered_ijc_1, filtered_sjc_1 = drop_items_by_index(
                    group_1_should_drop, [row['IJC_SAMPLE_1'],
                                          row['SJC_SAMPLE_1']])
                filtered_ijc_2, filtered_sjc_2 = drop_items_by_index(
                    group_2_should_drop, [row['IJC_SAMPLE_2'],
                                          row['SJC_SAMPLE_2']])

                row['IJC_SAMPLE_1'] = int_list_to_string(filtered_ijc_1)
                row['SJC_SAMPLE_1'] = int_list_to_string(filtered_sjc_1)
                row['IJC_SAMPLE_2'] = int_list_to_string(filtered_ijc_2)
                row['SJC_SAMPLE_2'] = int_list_to_string(filtered_sjc_2)

                out_columns = list()
                for header in headers:
                    out_columns.append(row[header])

                out_handle.write('{}\n'.format('\t'.join(out_columns)))


def parse_comma_ints(string):
    ints = list()
    int_strings = string.split(',')
    for int_string in int_strings:
        parsed = int(int_string)
        ints.append(parsed)

    return ints


def int_list_to_string(int_list):
    return ','.join([str(x) for x in int_list])


def find_zero_read_replicates(ijc, sjc):
    is_zero = list()
    for i, inc_count in enumerate(ijc):
        skip_count = sjc[i]
        is_zero.append(inc_count == 0 and skip_count == 0)

    return is_zero


def drop_items_by_index(drop_list, item_lists):
    out_lists = [list() for _ in item_lists]
    for should_drop_i, should_drop in enumerate(drop_list):
        if should_drop:
            continue

        for list_i, item_list in enumerate(item_lists):
            out_lists[list_i].append(item_list[should_drop_i])

    return out_lists


def write_fdr_file_from_darts_output(istat, ostat_pv, ostat_fdr):
    darts_id_to_probability = dict()
    with open(ostat_pv, 'rt') as darts_out_handle:
        for line_i, line in enumerate(darts_out_handle):
            columns = line.rstrip('\n').split('\t')
            if line_i == 0:
                try:
                    darts_id_index = columns.index('ID')
                    darts_post_index = columns.index('post_pr')
                except ValueError:
                    print('error parsing column names in darts output: {}'
                          .format(line), file=sys.stderr)
                    return

                continue

            id_str = columns[darts_id_index]
            probability_str = columns[darts_post_index]
            darts_id_to_probability[id_str] = probability_str

    with open(istat, 'rt') as counts_handle:
        with open(ostat_fdr, 'wt') as fdr_out_handle:
            for line_i, line in enumerate(counts_handle):
                columns = line.rstrip('\n').split('\t')
                if line_i == 0:
                    columns.extend(['PValue', 'FDR'])
                else:
                    id_str = columns[0]
                    # DARTS will not produce an output value for an event if
                    # there are 0 counts for a replicate.
                    probability_str = darts_id_to_probability.get(id_str, 'NA')
                    columns.extend([probability_str, 'NA'])

                fdr_out_handle.write('{}\n'.format('\t'.join(columns)))


def print_file_to_stderr(filename):
    with open(filename, 'rt') as file_handle:
        for line in file_handle:
            print(line, file=sys.stderr)


def copy_last_column_and_append(base_file, column_file, out_file):
    with open(base_file, 'rt') as base_handle:
        with open(column_file, 'rt') as column_handle:
            with open(out_file, 'wt') as out_handle:
                for base_line in base_handle:
                    column_line = column_handle.readline()
                    base_line = base_line.rstrip('\n')
                    column_line = column_line.rstrip('\n')
                    all_columns = column_line.split('\t')
                    last_column = all_columns[-1]
                    out_line = '{}\t{}\n'.format(base_line, last_column)
                    out_handle.write(out_line)


def append_columns_with_defaults(in_file_name, out_file_name, column_names, default_values):
    column_names_addition = '\t{}\n'.format('\t'.join(column_names))
    default_values_addition = '\t{}\n'.format('\t'.join(default_values))

    with open(in_file_name, 'rt') as in_file:
        with open(out_file_name, 'wt') as out_file:
            for i, raw_line in enumerate(in_file):
                line = raw_line.rstrip('\n')
                out_file.write(line)
                if i == 0:
                    out_file.write(column_names_addition)
                else:
                    out_file.write(default_values_addition)


def filter_countfile_by_imbalance_ratio(
        counts_file_name, indiv_counts_file_name, temp_file_name,
        imbalance_ratio, splicing_event_type):
    with open(counts_file_name, 'rt') as counts_file_handle:
        with open(indiv_counts_file_name, 'rt') as indiv_counts_file_handle:
            with open(temp_file_name, 'wt') as temp_file_handle:
                error = filter_countfile_by_imbalance_ratio_with_handles(
                    counts_file_handle, indiv_counts_file_handle,
                    temp_file_handle, imbalance_ratio, splicing_event_type)

    if error:
        formatted_message = (
            'error in filter_countfile_by_imbalance_ratio({}, {}, {}, {}, {}): {}'
            .format(counts_file_name, indiv_counts_file_name, temp_file_name,
                    imbalance_ratio, splicing_event_type, error))
        print(formatted_message, file=sys.stderr)
    else:
        shutil.move(temp_file_name, counts_file_name)


def filter_countfile_by_imbalance_ratio_with_handles(
        counts_file_handle, indiv_counts_file_handle, temp_file_handle,
        imbalance_ratio, splicing_event_type):
    if imbalance_ratio == 0:
        return 'imbalance_ratio == 0'

    inverse_imbalance_ratio = 1 / imbalance_ratio
    calculate_ratios = get_calculate_ratios_for_event_type(splicing_event_type)
    if calculate_ratios is None:
        return 'unexpected splicing_event_type: {}'.format(splicing_event_type)

    for i, counts_line in enumerate(counts_file_handle):
        counts_line = counts_line.rstrip('\n')
        if i == 0:
            temp_file_handle.write('{}\n'.format(counts_line))
            indiv_counts_header_line = (
                indiv_counts_file_handle.readline().rstrip('\n'))
            indiv_counts_header_columns = indiv_counts_header_line.split('\t')
            continue

        counts_columns = counts_line.split('\t')
        counts_id = counts_columns[0]
        indiv_counts_columns = get_indiv_counts_columns_for_id(
            counts_id, indiv_counts_file_handle)
        if indiv_counts_columns is None:
            return 'ID: {} not found in individual count file'.format(counts_id)

        is_imbalanced, error = check_if_imbalanced(
            imbalance_ratio, inverse_imbalance_ratio,
            indiv_counts_header_columns, indiv_counts_columns, calculate_ratios)
        if error:
            return error

        if is_imbalanced:
            continue

        temp_file_handle.write('{}\n'.format(counts_line))

    return None


def append_individual_counts(counts_file_name, mats_file_name, temp_file_name):
    with open(counts_file_name, 'rt') as counts_file_handle:
        with open(mats_file_name, 'rt') as mats_file_handle:
            with open(temp_file_name, 'wt') as temp_file_handle:
                error = append_individual_counts_with_handles(
                    counts_file_handle, mats_file_handle, temp_file_handle)

    if error:
        formatted_message = (
            'error in append_individual_counts({}, {}, {}): {}'
            .format(counts_file_name, mats_file_name, temp_file_name, error))
        print(formatted_message, file=sys.stderr)
    else:
        shutil.move(temp_file_name, mats_file_name)


def append_individual_counts_with_handles(counts_file_handle, mats_file_handle,
                                          temp_file_handle):
    for i, mats_line in enumerate(mats_file_handle):
        mats_line = mats_line.rstrip('\n')
        mats_columns = mats_line.split('\t')
        if i == 0:
            counts_header_line = counts_file_handle.readline().rstrip('\n')
            counts_header_columns = counts_header_line.split('\t')
            counts_header_cols_without_id = counts_header_columns[1:]
            combined_line = '\t'.join(mats_columns
                                      + counts_header_cols_without_id)
            temp_file_handle.write('{}\n'.format(combined_line))
            continue

        mats_id = mats_columns[0]
        counts_columns = get_indiv_counts_columns_for_id(mats_id,
                                                         counts_file_handle)
        if counts_columns is None:
            return 'ID: {} not found in individual count file'.format(mats_id)

        counts_cols_without_id = counts_columns[1:]
        combined_line = '\t'.join(mats_columns + counts_cols_without_id)
        temp_file_handle.write('{}\n'.format(combined_line))

    return None


def try_parse_float(float_str):
    try:
        return float(float_str), None
    except ValueError as e:
        return None, 'error in try_parse_float({}): {}'.format(float_str, e)


def try_sum_floats(floats_string):
    float_strings = floats_string.split(',')
    float_sum = 0
    for float_string in float_strings:
        parsed_float, error = try_parse_float(float_string)
        if error:
            return None, 'try_sum_floats({}): {}'.format(floats_string, error)

        float_sum += parsed_float

    return float_sum, None


def calculate_ratio(a, b):
    if 0 in [a, b]:
        return 0

    return a / b


def calculate_ratios_se(headers, values):
    upstream_to_target_count = None
    target_to_downstream_count = None
    for i, header in enumerate(headers):
        if header == 'upstream_to_target_count':
            upstream_to_target_count, error = try_sum_floats(values[i])
            if error:
                return None, error

        if header == 'target_to_downstream_count':
            target_to_downstream_count, error = try_sum_floats(values[i])
            if error:
                return None, error

    if None in [upstream_to_target_count, target_to_downstream_count]:
        error = ('Missing expected headers in calculate_ratios_se({}, {})'
                 .format(headers, values))
        return None, error

    ratio = calculate_ratio(upstream_to_target_count, target_to_downstream_count)
    return [ratio], None


def calculate_ratios_mxe(headers, values):
    upstream_to_first_count = None
    first_to_downstream_count = None
    upstream_to_second_count = None
    second_to_downstream_count = None
    for i, header in enumerate(headers):
        if header == 'upstream_to_first_count':
            upstream_to_first_count, error = try_sum_floats(values[i])
            if error:
                return None, error

        if header == 'first_to_downstream_count':
            first_to_downstream_count, error = try_sum_floats(values[i])
            if error:
                return None, error

        if header == 'upstream_to_second_count':
            upstream_to_second_count, error = try_sum_floats(values[i])
            if error:
                return None, error

        if header == 'second_to_downstream_count':
            second_to_downstream_count, error = try_sum_floats(values[i])
            if error:
                return None, error

    if None in [upstream_to_first_count, first_to_downstream_count,
                upstream_to_second_count, second_to_downstream_count]:
        error = ('Missing expected headers in calculate_ratios_mxe({}, {})'
                 .format(headers, values))
        return None, error

    first_ratio = calculate_ratio(upstream_to_first_count,
                                  first_to_downstream_count)
    second_ratio = calculate_ratio(upstream_to_second_count,
                                   second_to_downstream_count)
    return [first_ratio, second_ratio], None


def calculate_ratios_alt_ss(headers, values):
    across_short_boundary_count = None
    long_to_flanking_count = None
    for i, header in enumerate(headers):
        if header == 'across_short_boundary_count':
            across_short_boundary_count, error = try_sum_floats(values[i])
            if error:
                return None, error

        if header == 'long_to_flanking_count':
            long_to_flanking_count, error = try_sum_floats(values[i])
            if error:
                return None, error

    if None in [across_short_boundary_count, long_to_flanking_count]:
        error = ('Missing expected headers in calculate_ratios_alt_ss({}, {})'
                 .format(headers, values))
        return None, error

    ratio = calculate_ratio(across_short_boundary_count, long_to_flanking_count)
    return [ratio], None


def calculate_ratios_ri(headers, values):
    upstream_to_intron_count = None
    intron_to_downstream_count = None
    for i, header in enumerate(headers):
        if header == 'upstream_to_intron_count':
            upstream_to_intron_count, error = try_sum_floats(values[i])
            if error:
                return None, error

        if header == 'intron_to_downstream_count':
            intron_to_downstream_count, error = try_sum_floats(values[i])
            if error:
                return None, error

    if None in [upstream_to_intron_count, intron_to_downstream_count]:
        error = ('Missing expected headers in calculate_ratios_ri({}, {})'
                 .format(headers, values))
        return None, error

    ratio = calculate_ratio(upstream_to_intron_count, intron_to_downstream_count)
    return [ratio], None


def get_calculate_ratios_for_event_type(event_type):
    if event_type == 'SE':
        return calculate_ratios_se
    if event_type == 'MXE':
        return calculate_ratios_mxe
    if event_type in ['A3SS', 'A5SS']:
        return calculate_ratios_alt_ss
    if event_type == 'RI':
        return calculate_ratios_ri

    return None


def check_if_imbalanced(imbalance_ratio, inverse_imbalance_ratio,
                        header_columns, counts_columns, calculate_ratios_func):
    ratios, error = calculate_ratios_func(header_columns,
                                          counts_columns)
    if error:
        return None, error

    for ratio in ratios:
        if ratio > imbalance_ratio or ratio < inverse_imbalance_ratio:
            return True, None

    return False, None


def get_indiv_counts_columns_for_id(mats_id, counts_file_handle):
    # The IDs have the same order in both files, but some IDs may only be
    # in the counts file.
    for counts_line in counts_file_handle:
        counts_line = counts_line.rstrip('\n')
        counts_columns = counts_line.split('\t')
        counts_id = counts_columns[0]
        if counts_id == mats_id:
            return counts_columns

    return None


def get_python_executable():
    # Try to get the absolute path of the executable for the running
    # Python interpreter.
    python_executable = sys.executable
    if not python_executable:
        # Fallback
        print('Absolute path for current Python interpreter not found.'
              ' Using "python" without a full path to run scripts',
              file=sys.stderr)
        python_executable = 'python'

    return python_executable


def generate_summary(python_executable, out_dir, root_dir):
    summary_script = os.path.join(root_dir, 'rMATS_P', 'summary.py')
    summary_out_file_path = os.path.join(out_dir, 'summary.txt')
    with open(summary_out_file_path, 'wb') as f_handle:
        subprocess.call([python_executable, summary_script, out_dir],
                        stdout=f_handle)


def claim_prep_prefix(task, tmp_dir):
    if task not in ['prep', 'both']:
        return None

    file_name_template = os.path.join(tmp_dir, '{}_0.rmats')
    prep_prefix = None
    while True:
        prep_prefix = datetime.fromtimestamp(time.time()).strftime(
            '%Y-%m-%d-%H_%M_%S_%f')
        file_path = file_name_template.format(prep_prefix)
        if not os.path.exists(file_path):
            with open(file_path, 'wt'):
                pass  # create the file to claim the timestamp

            break

    return prep_prefix


def create_output_dirs(args):
    args.out_tmp_sub_dir = os.path.join(args.od, 'tmp')
    for dir_path in [args.od, args.out_tmp_sub_dir, args.tmp]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)   # python2: makedirs() got an unexpected keyword argument 'exist_ok'


def apply_id_mapping(event_type, out_dir):
    id_to_orig = dict()
    mapping_path = os.path.join(out_dir, 'id_mapping.{}.txt'.format(event_type))
    with open(mapping_path, 'rt') as map_f:
        for i, line in enumerate(map_f):
            values = line.strip().split('\t')
            if i == 0:
                expected_mapping_headers = ['original_id', 'mapped_id']
                if values != expected_mapping_headers:
                    sys.exit('ERROR: expected headers in {} to be {}'
                             ' but found {}'.format(
                                 mapping_path, expected_mapping_headers, values))

                continue

            id_to_orig[values[1]] = values[0]

    file_templates = ['fromGTF.{}.txt', 'JC.raw.input.{}.txt',
                      'JCEC.raw.input.{}.txt', 'individualCounts.{}.txt']
    for file_template in file_templates:
        file_path = os.path.join(out_dir, file_template.format(event_type))
        if not os.path.exists(file_path):
            continue

        with open(file_path, 'rt') as in_f:
            # using mapping_path as a temporary file which is then removed
            with open(mapping_path, 'wt') as out_f:
                for i, line in enumerate(in_f):
                    values = line.strip().split('\t')
                    if i == 0:
                        id_i = values.index('ID')
                        out_f.write(line)
                        continue

                    mapped_id = values[id_i]
                    values[id_i] = id_to_orig[mapped_id]
                    out_f.write('{}\n'.format('\t'.join(values)))

        shutil.move(mapping_path, file_path)


def main():
    """Takes no arguments.
    Processes arguments supplied when rmats.py was called using get_args().
    If task argument is 'inte', checks BAM and prep file integrity.
    If task argument is 'prep', 'post', or 'both', runs pipeline using seperate module in rmatspipeline.
    If task argument is not valid, returns nothing.
    For each splicing event type, processes counts and outputs files.
    Generates output summary.
    Prints 'Done processing count files.' and returns nothing.
    """
    args = get_args()

    if args.task == 'inte':
        check_integrity(args.bams, args.tmp)
    if args.task in pipe_tasks:
        run_pipe(args)

    if args.task not in ['post', 'both', 'stat']:
        return

    jc_it = os.path.join(args.od, 'JC.raw.input.%s.txt')
    jcec_it = os.path.join(args.od, 'JCEC.raw.input.%s.txt')

    python_executable = get_python_executable()
    # realpath is used here to support the conda install which creates a
    # symlink for rmats.py
    root_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

    print('Processing count files.')
    for event_type in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']:
        if args.fixed_event_set:
            apply_id_mapping(event_type, args.od)

        process_counts(jc_it % (event_type), args.tstat, 'JC', event_type,
                       args.cstat, args.od, args.out_tmp_sub_dir, args.stat,
                       args.paired_stats, args.darts_model, args.darts_cutoff,
                       python_executable, root_dir, args.imbalance_ratio,
                       args.drop_zero_read_replicates_for_stat)
        process_counts(jcec_it % (event_type), args.tstat, 'JCEC', event_type,
                       args.cstat, args.od, args.out_tmp_sub_dir, args.stat,
                       args.paired_stats, args.darts_model, args.darts_cutoff,
                       python_executable, root_dir, args.imbalance_ratio,
                       args.drop_zero_read_replicates_for_stat)

    generate_summary(python_executable, args.od, root_dir)
    print('Done processing count files.')

    return


if __name__ == "__main__":
    main()
