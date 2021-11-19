from __future__ import print_function

import argparse
import math
import os.path
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description='summarize the rMATS results from the output files')
    parser.add_argument('output_dir',
                        help='path to the directory containing rMATS output')
    parser.add_argument(
        '--use-raw-p',
        action='store_true',
        help='check --p-cutoff against the raw PValue instead of against FDR')
    parser.add_argument(
        '--p-cutoff',
        type=float,
        help='cutoff value used for FDR (or PValue for --use-raw-p)'
        ' for determining significance. Default: %(default)s',
        default=0.05)
    parser.add_argument(
        '--inc-level-diff-cutoff',
        type=float,
        help='cutoff value used for the absolute value of IncLevelDifference'
        ' for determining significance. Default: %(default)s',
        default=0)

    output_file_group = parser.add_mutually_exclusive_group()
    output_file_group.add_argument(
        '--summary-path',
        help='path to write the generated summary to instead of stdout')
    output_file_group.add_argument(
        '--summary-prefix',
        help='prefix of path to write the generated summary to instead of'
        ' stdout. The cutoff values and ".txt" will be appended to the prefix'
        ' to obtain the actual path to write to')
    return parser.parse_args()


def get_output_file_path(args):
    if args.summary_path:
        return args.summary_path

    if args.summary_prefix:
        p_cutoff_col = 'PValue' if args.use_raw_p else 'FDR'
        return '{}_{}_{}_IncLevelDifference_{}.txt'.format(
            args.summary_prefix, p_cutoff_col, args.p_cutoff,
            args.inc_level_diff_cutoff)

    return None


def parse_float(s):
    try:
        return float(s)
    except ValueError:
        return float('nan')


def count_events(file_path, args):
    total = 0
    sig = 0
    sig_sample_1_higher = 0
    sig_sample_2_higher = 0
    if os.path.exists(file_path):
        with open(file_path, 'rt') as file_handle:
            for line_i, line in enumerate(file_handle):
                columns = line.strip().split('\t')
                if line_i == 0:
                    headers = columns
                    continue

                row = dict(zip(headers, columns))
                total += 1
                p_value = parse_float(row['PValue'])
                fdr = parse_float(row['FDR'])
                inc_level_diff = parse_float(row['IncLevelDifference'])
                check_p_value = p_value if args.use_raw_p else fdr
                if (math.isnan(check_p_value) or math.isnan(inc_level_diff)
                        or check_p_value > args.p_cutoff
                        or abs(inc_level_diff) < args.inc_level_diff_cutoff):
                    continue

                sig += 1
                if inc_level_diff > 0:
                    sig_sample_1_higher += 1
                else:
                    sig_sample_2_higher += 1

    return {
        'total': total,
        'sig': sig,
        'sig_sample_1_higher': sig_sample_1_higher,
        'sig_sample_2_higher': sig_sample_2_higher
    }


def summarize(args, output_file_handle):
    headers = [
        'EventType', 'TotalEventsJC', 'TotalEventsJCEC', 'SignificantEventsJC',
        'SigEventsJCSample1HigherInclusion',
        'SigEventsJCSample2HigherInclusion', 'SignificantEventsJCEC',
        'SigEventsJCECSample1HigherInclusion',
        'SigEventsJCECSample2HigherInclusion'
    ]

    print('\t'.join(headers), file=output_file_handle)
    for event in ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']:
        jc_path = os.path.join(args.output_dir, '{}.MATS.JC.txt'.format(event))
        jcec_path = os.path.join(args.output_dir,
                                 '{}.MATS.JCEC.txt'.format(event))
        jc_event_counts = count_events(jc_path, args)
        jcec_event_counts = count_events(jcec_path, args)
        jc_total = jc_event_counts['total']
        jcec_total = jcec_event_counts['total']

        values = [
            event, jc_total, jcec_total, jc_event_counts['sig'],
            jc_event_counts['sig_sample_1_higher'],
            jc_event_counts['sig_sample_2_higher'], jcec_event_counts['sig'],
            jcec_event_counts['sig_sample_1_higher'],
            jcec_event_counts['sig_sample_2_higher']
        ]
        print('\t'.join([str(v) for v in values]), file=output_file_handle)


def main():
    args = parse_args()
    output_file_path = get_output_file_path(args)
    if output_file_path:
        with open(output_file_path, 'wt') as output_file_handle:
            summarize(args, output_file_handle)
    else:
        summarize(args, sys.stdout)


if __name__ == '__main__':
    main()
