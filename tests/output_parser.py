import os.path


def parse_tab_separated_value_file(file_path):
    header = list()
    rows = list()
    if not os.path.isfile(file_path):
        return header, rows, '{} is not a valid filepath'.format(file_path)

    with open(file_path, 'rt') as f_handle:
        for i, line in enumerate(f_handle):
            tokens = line.strip('\n').split('\t')
            if i == 0:
                header = tokens
                continue

            if len(tokens) != len(header):
                return header, rows, 'line {} has {} tokens but the header has {} in {}'.format(
                    i, len(tokens), len(header), file_path)

            rows.append(dict(zip(header, tokens)))

    return header, rows, None


def parse_from_gtf(from_gtf_path):
    return parse_tab_separated_value_file(from_gtf_path)


def parse_from_gtf_novel_junction(from_gtf_novel_junction_path):
    return parse_tab_separated_value_file(from_gtf_novel_junction_path)


def parse_from_gtf_novel_splice_site(from_gtf_novel_splice_site_path):
    return parse_tab_separated_value_file(from_gtf_novel_splice_site_path)


def parse_jc_raw(jc_raw_path):
    return parse_tab_separated_value_file(jc_raw_path)


def parse_jcec_raw(jcec_raw_path):
    return parse_tab_separated_value_file(jcec_raw_path)


def parse_mats_jc(mats_jc_path):
    return parse_tab_separated_value_file(mats_jc_path)


def parse_mats_jcec(mats_jcec_path):
    return parse_tab_separated_value_file(mats_jcec_path)


def _parse_int(s):
    try:
        return int(s), None
    except ValueError as e:
        return None, 'could not parse {} as int: {}'.format(s, e)


def _parse_int_line(f_handle):
    int_line = f_handle.readline()
    return _parse_int(int_line.strip())


def _parse_by_gene_by_bam_section(bams, f_handle, value_parser):
    values_by_gene_by_bam = list()
    for _ in range(len(bams)):
        values_by_gene = dict()
        values_by_gene_by_bam.append(values_by_gene)
        num_genes_for_bam, error = _parse_int_line(f_handle)
        if error:
            return None, error

        for _ in range(num_genes_for_bam):
            values_line = f_handle.readline()
            semi_splits = values_line.strip().split(';')
            num_values_for_gene = len(semi_splits) - 1
            if num_values_for_gene < 1:
                return None, 'no values for gene in line: {}'.format(
                    values_line)

            gene_id = semi_splits[0]
            values_for_gene = list()
            values_by_gene[gene_id] = values_for_gene
            for value_str in semi_splits[1:]:
                value, error = value_parser(value_str)
                if error:
                    return None, error

                values_for_gene.append(value)

    return values_by_gene_by_bam, None


def _parse_n_separated_ints(n, sep, string):
    int_strs = string.split(sep)
    if len(int_strs) != n:
        return None, 'expected {} ints: {}'.format(n, string)

    ints = list()
    for int_str in int_strs:
        int_v, error = _parse_int(int_str)
        if error:
            return None, error

        ints.append(int_v)

    return ints, None


def _parse_novel_junc(novel_junc_str):
    return _parse_n_separated_ints(3, ',', novel_junc_str)


def _parse_exons(exon_str):
    exon_ints, error = _parse_n_separated_ints(5, ',', exon_str)
    if error:
        return None, error

    return {
        'start_box': exon_ints[0:2],
        'end_box': exon_ints[2:4],
        'count': exon_ints[4]
    }, None


def _parse_multis(multi_str):
    comma_splits = multi_str.split(',')
    if len(comma_splits) != 2:
        return None, 'expected a read string and a count in {}'.format(
            multi_str)

    read_str, count_str = comma_splits
    junction_pairs = list()
    junction_pair_strs = read_str.split('=')
    for junction_pair_str in junction_pair_strs:
        int_pair, error = _parse_n_separated_ints(2, ':', junction_pair_str)
        if error:
            return None, error

        junction_pairs.append(int_pair)

    count, error = _parse_int(count_str)
    if error:
        return None, error

    return {'junction_pairs': junction_pairs, 'count': count}, None


def parse_dot_rmats(dot_rmats_path):
    if not os.path.isfile(dot_rmats_path):
        return None, '{} is not a valid filepath'.format(dot_rmats_path)

    with open(dot_rmats_path, 'rt') as f_handle:
        bam_line = f_handle.readline()
        bams = bam_line.strip().split(',')

        read_length, error = _parse_int_line(f_handle)
        if error:
            return None, error

        novel_juncs, error = _parse_by_gene_by_bam_section(
            bams, f_handle, _parse_novel_junc)
        if error:
            return None, error

        exons, error = _parse_by_gene_by_bam_section(bams, f_handle,
                                                     _parse_exons)
        if error:
            return None, error

        multis, error = _parse_by_gene_by_bam_section(bams, f_handle,
                                                      _parse_multis)
        if error:
            return None, error

    return {
        'bams': bams,
        'read_length': read_length,
        'novel_juncs': novel_juncs,
        'exons': exons,
        'multis': multis
    }, None
