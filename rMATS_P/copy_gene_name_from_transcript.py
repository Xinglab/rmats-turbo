import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Copy the gene_name attribute in a gtf from'
                     ' transcript lines to exon lines'))
    parser.add_argument('--in-gtf', help='input gtf', required=True)
    parser.add_argument('--out-gtf', help='output gtf', required=True)

    return parser.parse_args()


def parse_attribute_str(attribute_str):
    attribute_str = attribute_str.strip()
    attribute_str = attribute_str.strip(';')
    key_value_strings = attribute_str.split('; ')
    attributes = dict()
    for key_value_string in key_value_strings:
        first_space = key_value_string.find(' ')
        if first_space < 0:
            continue  # no space found

        key = key_value_string[:first_space]
        value = key_value_string[first_space + 1:]
        attributes[key] = value

    return attributes


def format_attribute_str(attributes):
    parts = list()
    for key, value in attributes.items():
        part = '{} {};'.format(key, value)
        parts.append(part)

    return ' '.join(parts)


def parse_gtf_line(line):
    line = line.strip()
    if line.startswith('#'):
        # ignore comment lines
        return None

    columns = line.split('\t')
    seq_name = columns[0]
    source = columns[1]
    feature = columns[2]
    start = columns[3]
    end = columns[4]
    score = columns[5]
    strand = columns[6]
    frame = columns[7]
    attribute_str = columns[8]
    attributes = parse_attribute_str(attribute_str)

    parsed = {
        'seq': seq_name,
        'source': source,
        'feature': feature,
        'start': start,
        'end': end,
        'score': score,
        'strand': strand,
        'frame': frame,
        'attributes': attributes
    }

    return parsed


def write_gtf_line(parsed, handle):
    attribute_str = format_attribute_str(parsed['attributes'])
    formatted = '\t'.join([
        parsed['seq'], parsed['source'], parsed['feature'],
        parsed['start'], parsed['end'], parsed['score'], parsed['strand'],
        parsed['frame'], attribute_str
    ])
    handle.write('{}\n'.format(formatted))


def get_gene_names(in_gtf):
    transcript_id_to_gene_name = dict()
    with open(in_gtf, 'rt') as in_handle:
        for line in in_handle:
            parsed_line = parse_gtf_line(line)
            if parsed_line is None:
                continue

            transcript_id = parsed_line['attributes'].get('transcript_id')
            gene_name = parsed_line['attributes'].get('gene_name')
            if None in [transcript_id, gene_name]:
                continue

            transcript_id_to_gene_name[transcript_id] = gene_name

    return transcript_id_to_gene_name


def add_gene_names(in_gtf, out_gtf, transcript_id_to_gene_name):
    with open(in_gtf, 'rt') as in_handle:
        with open(out_gtf, 'wt') as out_handle:
            for line in in_handle:
                parsed_line = parse_gtf_line(line)
                if parsed_line is None:
                    out_handle.write(line)
                    continue

                transcript_id = parsed_line['attributes'].get('transcript_id')
                gene_name = parsed_line['attributes'].get('gene_name')
                gene_name_lookup = transcript_id_to_gene_name.get(transcript_id)
                if (gene_name is not None) or gene_name_lookup is None:
                    out_handle.write(line)
                    continue

                parsed_line['attributes']['gene_name'] = gene_name_lookup
                write_gtf_line(parsed_line, out_handle)


def main():
    args = parse_args()
    transcript_id_to_gene_name = get_gene_names(args.in_gtf)
    add_gene_names(args.in_gtf, args.out_gtf, transcript_id_to_gene_name)


if __name__ == '__main__':
    main()
