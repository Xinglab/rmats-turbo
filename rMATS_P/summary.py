from __future__ import print_function

import argparse
import csv
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
    return parser.parse_args()


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
    with open(file_path, 'rt') as file_handle:
        reader = csv.DictReader(file_handle, delimiter='\t')
        for row in reader:
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


def summarize(args):
    headers = [
        'EventType', 'TotalEvents', 'SignificantEventsJC',
        'SigEventsJCSample1HigherInclusion',
        'SigEventsJCSample2HigherInclusion', 'SignificantEventsJCEC',
        'SigEventsJCECSample1HigherInclusion',
        'SigEventsJCECSample2HigherInclusion'
    ]

    print('\t'.join(headers))
    for event in ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']:
        jc_path = os.path.join(args.output_dir, '{}.MATS.JC.txt'.format(event))
        jcec_path = os.path.join(args.output_dir,
                                 '{}.MATS.JCEC.txt'.format(event))
        jc_event_counts = count_events(jc_path, args)
        jcec_event_counts = count_events(jcec_path, args)
        jc_total = jc_event_counts['total']
        jcec_total = jcec_event_counts['total']
        if jc_total != jcec_total:
            print('Total {} event counts should match for JC and JCEC'
                  ' but saw: JC={}, JCEC={}'.format(event, jc_total,
                                                    jcec_total),
                  file=sys.stderr)

        values = [
            event, jc_total, jc_event_counts['sig'],
            jc_event_counts['sig_sample_1_higher'],
            jc_event_counts['sig_sample_2_higher'], jcec_event_counts['sig'],
            jcec_event_counts['sig_sample_1_higher'],
            jcec_event_counts['sig_sample_2_higher']
        ]
        print('\t'.join([str(v) for v in values]))


def main():
    args = parse_args()
    summarize(args)


if __name__ == '__main__':
    main()
