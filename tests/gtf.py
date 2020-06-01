import tests.util


def _format_attributes(attrs):
    formatted_pairs = list()
    for k, v in attrs:
        formatted_pairs.append(' '.join([k, tests.util.double_quote(v)]))

    return '; '.join(formatted_pairs)


def _group_transcripts_by_gene_id(transcripts):
    transcripts_by_gene_id = dict()
    for transcript in transcripts:
        gene_id = transcript.gene_id
        transcripts_for_gene_id = transcripts_by_gene_id.get(gene_id)
        if not transcripts_for_gene_id:
            transcripts_for_gene_id = list()
            transcripts_by_gene_id[gene_id] = transcripts_for_gene_id

        transcripts_for_gene_id.append(transcript)

    return transcripts_by_gene_id


def _get_lines_for_transcript(transcript, gene_id, gene_name, chromosome,
                              strand, score_placeholder, frame_placeholder):
    lines_for_transcript = list()
    trans_min_start, trans_max_end = transcript.exons[0]
    for exon_start, exon_end in transcript.exons:
        if exon_start < trans_min_start:
            trans_min_start = exon_start

        if exon_end > trans_max_end:
            trans_max_end = exon_end

        attributes = [('gene_id', gene_id),
                      ('transcript_id', transcript.transcript_id),
                      ('gene_name', gene_name)]
        formatted_attributes = _format_attributes(attributes)
        lines_for_transcript.append('\t'.join([
            chromosome,
            'processed_transcript',
            'exon',
            str(exon_start),
            str(exon_end),
            score_placeholder,
            strand,
            frame_placeholder,
            formatted_attributes,
        ]))

    attributes = [('gene_id', gene_id),
                  ('transcript_id', transcript.transcript_id),
                  ('gene_name', gene_name)]
    formatted_attributes = _format_attributes(attributes)
    lines_for_transcript.insert(
        0, '\t'.join([
            chromosome,
            'processed_transcript',
            'transcript',
            str(trans_min_start),
            str(trans_max_end),
            score_placeholder,
            strand,
            frame_placeholder,
            formatted_attributes,
        ]))

    return {
        'lines_for_transcript': lines_for_transcript,
        'trans_min_start': trans_min_start,
        'trans_max_end': trans_max_end
    }


def _write_lines_for_gene_id(f_handle, gene_id, transcripts):
    score_placeholder = '.'
    frame_placeholder = '.'

    lines_for_gene = list()
    gene_min_start = None
    gene_max_end = None
    gene_name = None
    chromosome = None
    strand = None
    for transcript in transcripts:
        gene_name = transcript.gene_name
        chromosome = transcript.chromosome
        strand = transcript.strand

        if not transcript.exons:
            return 'no exons in transcript'

        lines_for_transcript_result = _get_lines_for_transcript(
            transcript, gene_id, gene_name, chromosome, strand,
            score_placeholder, frame_placeholder)
        lines_for_transcript = lines_for_transcript_result[
            'lines_for_transcript']
        trans_min_start = lines_for_transcript_result['trans_min_start']
        trans_max_end = lines_for_transcript_result['trans_max_end']

        lines_for_gene.extend(lines_for_transcript)
        if gene_min_start is None or trans_min_start < gene_min_start:
            gene_min_start = trans_min_start

        if gene_max_end is None or trans_max_end > gene_max_end:
            gene_max_end = trans_max_end

    attributes = [('gene_id', gene_id), ('gene_name', gene_name)]
    formatted_attributes = _format_attributes(attributes)
    lines_for_gene.insert(
        0, '\t'.join([
            chromosome,
            'pseudogene',
            'gene',
            str(gene_min_start),
            str(gene_max_end),
            score_placeholder,
            strand,
            frame_placeholder,
            formatted_attributes,
        ]))
    for line in lines_for_gene:
        f_handle.write('{}\n'.format(line))

    return None


class GTF(object):
    def __init__(self):
        self.path = None
        self.transcripts = None

    def write(self):
        if not self.path:
            return 'path not valid: {}'.format(self.path)

        if not self.transcripts:
            return 'no transcripts to write'

        transcripts_by_gene_id = _group_transcripts_by_gene_id(
            self.transcripts)
        with open(self.path, 'wt') as f_handle:
            for gene_id, transcripts in transcripts_by_gene_id.items():
                error = _write_lines_for_gene_id(f_handle, gene_id,
                                                 transcripts)
                if error:
                    return 'error in gene_id: {}: {}'.format(gene_id, error)

        return None


class Transcript(object):
    def __init__(self):
        self.chromosome = None
        self.strand = None
        self.gene_id = None
        self.gene_name = None
        self.transcript_id = None
        self.exons = None
