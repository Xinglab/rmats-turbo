import argparse
import os
import os.path


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Extract PSI values from rMATS-turbo output'))
    parser.add_argument('--rmats-out-dir',
                        required=True,
                        help='the --od from rMATS')
    parser.add_argument('--out-tsv',
                        required=True,
                        help='where to write the PSI values')
    parser.add_argument('--average-read-count',
                        type=int,
                        default=10,
                        help=('Filter out events with lower average read count'
                              ' (default %(default)s)'))
    parser.add_argument(
        '--sample-names',
        required=False,
        help=('A comma separated list of sample names corresponding to'
              ' the order from --b1 and then --b2'))

    return parser.parse_args()


def write_tsv_line(columns, handle):
    handle.write('{}\n'.format('\t'.join([str(x) for x in columns])))


def read_tsv_line(line):
    columns = line.rstrip('\n').split('\t')
    return columns


def parse_comma_values(string, func):
    if string == '':
        return list()

    parts = string.split(',')
    result = list()
    for part in parts:
        try:
            value = func(part)
        except ValueError:
            value = None

        result.append(value)

    return result


def finalize_sample_names(sample_names, num_samples, file_name):
    if not sample_names:
        sample_names = list()
        highest_sample_str = str(num_samples - 1)
        num_digits = len(highest_sample_str)
        for i in range(num_samples):
            format_str = 'sample_{:0' + str(num_digits) + '}'
            sample_names.append(format_str.format(i))
    elif len(sample_names) != num_samples:
        raise Exception('Expected {} samples based on IncLevel columns'
                        ' in {}, but got {}'.format(num_samples, file_name,
                                                    sample_names))

    return sample_names


def write_psi_for_file(file_name, event_type, read_count_threshold,
                       is_first_file, out_headers, sample_names, in_handle,
                       out_handle):
    for line_i, line in enumerate(in_handle):
        in_columns = read_tsv_line(line)
        if line_i == 0:
            in_headers = in_columns
            continue

        row = dict(zip(in_headers, in_columns))
        event_id = row['ID']
        ijc_1_str = row['IJC_SAMPLE_1']
        sjc_1_str = row['SJC_SAMPLE_1']
        ijc_2_str = row['IJC_SAMPLE_2']
        sjc_2_str = row['SJC_SAMPLE_2']
        inc_1_str = row['IncLevel1']
        inc_2_str = row['IncLevel2']

        ijc_1 = parse_comma_values(ijc_1_str, int)
        sjc_1 = parse_comma_values(sjc_1_str, int)
        ijc_2 = parse_comma_values(ijc_2_str, int)
        sjc_2 = parse_comma_values(sjc_2_str, int)
        inc_1 = parse_comma_values(inc_1_str, float)
        inc_2 = parse_comma_values(inc_2_str, float)

        ijc = ijc_1 + ijc_2
        sjc = sjc_1 + sjc_2
        inc = inc_1 + inc_2
        num_samples = len(ijc)
        if (line_i == 1) and is_first_file:
            sample_names = finalize_sample_names(sample_names, num_samples,
                                                 file_name)
            out_headers.extend(sample_names)
            write_tsv_line(out_headers, out_handle)

        total_reads = sum(ijc) + sum(sjc)
        avg_reads = total_reads / num_samples
        if avg_reads < read_count_threshold:
            continue

        if any([x is None for x in inc]):
            continue

        write_tsv_line([event_type, event_id] + inc, out_handle)


def parse_sample_names(sample_names):
    if not sample_names:
        return None

    return sample_names.split(',')


def extract_psi_for_pca(rmats_dir, out_tsv, read_count_threshold,
                        sample_names):
    sample_names = parse_sample_names(sample_names)

    suffix = '.MATS.JC.txt'
    with open(out_tsv, 'wt') as out_handle:
        headers = ['event_type', 'id']
        is_first_file = True
        file_names = sorted(os.listdir(rmats_dir))
        for name in file_names:
            if not name.endswith(suffix):
                continue

            path = os.path.join(rmats_dir, name)
            event_type = name[:-len(suffix)]
            with open(path, 'rt') as in_handle:
                write_psi_for_file(name, event_type, read_count_threshold,
                                   is_first_file, headers, sample_names,
                                   in_handle, out_handle)
                is_first_file = False


def main():
    args = parse_args()
    extract_psi_for_pca(args.rmats_out_dir, args.out_tsv,
                        args.average_read_count, args.sample_names)


if __name__ == '__main__':
    main()
