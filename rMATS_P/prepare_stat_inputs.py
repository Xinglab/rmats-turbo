from __future__ import print_function

import argparse
import os
import os.path
import shutil
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description='Modify the files from an rMATS output directory to prepare'
        ' a new --od directory for running --task stat')
    parser.add_argument(
        '--new-output-dir',
        required=True,
        help='Path to the directory to be used for running --task stat')
    parser.add_argument('--old-output-dir',
                        required=True,
                        help='Path to the directory of existing rMATS output')
    parser.add_argument(
        '--group-1-indices',
        required=True,
        help='Comma separated list of replicates to be used as group 1 in the'
        ' new statistical comparison. Each replicate is identified by its'
        ' 0-based index in the concatenation of replicates from groups 1 and 2'
        ' in the --old-output-dir. Example: If --old-output-dir had 2 reps in'
        ' group 1 and 3 in group 2 then --group-1-indices 0,1,4 would select'
        ' both reps from the old group 1 and the last rep from the old group'
        ' 2 to be the new group 1')
    parser.add_argument(
        '--group-2-indices',
        required=True,
        help='Comma separated list of replicates to be used as group 2 in the'
        ' new statistical comparison')

    return parser.parse_args()


def parse_indices(indices_str):
    indices_strs = indices_str.split(',')
    indices = list()
    for i, index_str in enumerate(indices_strs):
        try:
            index = int(index_str)
        except ValueError as e:
            sys.exit('could not parse index string {} ({}) in {}: {}'.format(
                i, index_str, indices_str, e))

        indices.append(index)

    return indices


def prepare_stat_inputs(args):
    group_1_indices = parse_indices(args.group_1_indices)
    group_2_indices = parse_indices(args.group_2_indices)
    if not os.path.exists(args.new_output_dir):
        print('creating {}'.format(args.new_output_dir))
        os.makedirs(args.new_output_dir)

    for event_type in ['SE', 'A3SS', 'A5SS', 'RI', 'MXE']:
        from_gtf_name = 'fromGTF.{}.txt'.format(event_type)
        old_from_gtf = os.path.join(args.old_output_dir, from_gtf_name)
        new_from_gtf = os.path.join(args.new_output_dir, from_gtf_name)
        print('copying {} to {}'.format(old_from_gtf, new_from_gtf))
        shutil.copy(old_from_gtf, new_from_gtf)

        for count_type in ['JC', 'JCEC']:
            count_name = '{}.raw.input.{}.txt'.format(count_type, event_type)
            old_counts = os.path.join(args.old_output_dir, count_name)
            new_counts = os.path.join(args.new_output_dir, count_name)
            print('creating {} based on {}'.format(new_counts, old_counts))
            with open(old_counts, 'rt') as old_counts_handle:
                with open(new_counts, 'wt') as new_counts_handle:
                    create_count_file(old_counts_handle, new_counts_handle,
                                      group_1_indices, group_2_indices)


def split_and_combine(vals_1, vals_2):
    split_1 = vals_1.split(',')
    split_2 = vals_2.split(',')
    combined = split_1 + split_2
    filtered = [x for x in combined if x != '']
    return filtered


def create_count_file(old_counts_handle, new_counts_handle, group_1_indices,
                      group_2_indices):
    header = old_counts_handle.readline()
    header_columns = header.strip().split('\t')
    expected_header_columns = [
        'ID', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2',
        'IncFormLen', 'SkipFormLen'
    ]
    if header_columns != expected_header_columns:
        sys.exit('Unexpected header_columns: {}. Expected: {}'.format(
            header_columns, expected_header_columns))

    new_counts_handle.write(header)
    ijc_1_index = expected_header_columns.index('IJC_SAMPLE_1')
    sjc_1_index = expected_header_columns.index('SJC_SAMPLE_1')
    ijc_2_index = expected_header_columns.index('IJC_SAMPLE_2')
    sjc_2_index = expected_header_columns.index('SJC_SAMPLE_2')

    for line in old_counts_handle:
        values = line.strip().split('\t')
        old_ijc_1 = values[ijc_1_index]
        old_sjc_1 = values[sjc_1_index]
        old_ijc_2 = values[ijc_2_index]
        old_sjc_2 = values[sjc_2_index]
        old_ijc_combined = split_and_combine(old_ijc_1, old_ijc_2)
        old_sjc_combined = split_and_combine(old_sjc_1, old_sjc_2)
        new_ijc_1 = [old_ijc_combined[i] for i in group_1_indices]
        new_sjc_1 = [old_sjc_combined[i] for i in group_1_indices]
        new_ijc_2 = [old_ijc_combined[i] for i in group_2_indices]
        new_sjc_2 = [old_sjc_combined[i] for i in group_2_indices]
        values[ijc_1_index] = ','.join(new_ijc_1)
        values[sjc_1_index] = ','.join(new_sjc_1)
        values[ijc_2_index] = ','.join(new_ijc_2)
        values[sjc_2_index] = ','.join(new_sjc_2)
        new_line = '\t'.join(values)
        new_counts_handle.write('{}\n'.format(new_line))


def main():
    args = parse_args()
    prepare_stat_inputs(args)


if __name__ == '__main__':
    main()
