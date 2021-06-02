import os
import os.path
import subprocess
import tempfile


def _write_read_line(read, f_handle):
    line_values = [
        read.template_name,
        str(read.flags_int()),
        read.ref_seq_name,
        str(read.start_coord),
        str(read.map_quality),
        read.cigar,
        read.next_ref_seq_name,
        str(read.next_start_coord),
        str(read.template_len),
        read.read_sequence,
        read.read_quality,
        'NH:i:{}'.format(read.number_of_alignments),
    ]
    f_handle.write('{}\n'.format('\t'.join(line_values)))


class BAM(object):
    def __init__(self):
        self.path = None
        self.reads = None
        self._sorted_reads_by_ref_seq = None
        self._sorted_ref_seqs = None

    def write(self):
        if not self.path:
            return 'path not valid: {}'.format(self.path)

        if not self.reads:
            return 'no reads to write'

        self._sort_reads()

        temp_sam_file_path = None
        try:
            with tempfile.NamedTemporaryFile(mode='wt',
                                             prefix='tmp',
                                             suffix='.sam',
                                             delete=False) as temp_f_handle:
                temp_sam_file_path = temp_f_handle.name

                self._write_header(temp_f_handle)
                for ref_seq in self._sorted_ref_seqs:
                    for read in self._sorted_reads_by_ref_seq[ref_seq]:
                        _write_read_line(read, temp_f_handle)

            sam_to_bam_command = [
                'samtools', 'view', '-b', temp_sam_file_path, '-o', self.path
            ]
            result = subprocess.run(sam_to_bam_command,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=False)
            if result.returncode != 0:
                return 'could not convert sam to bam: {}\nstdout:\n{}\nstderr:\n{}'.format(
                    result.returncode, result.stdout.decode(),
                    result.stderr.decode())

        finally:
            if temp_sam_file_path and os.path.exists(temp_sam_file_path):
                os.remove(temp_sam_file_path)

        return None

    def _sort_reads(self):
        self._sorted_reads_by_ref_seq = dict()
        for read in self.reads:
            reads = self._sorted_reads_by_ref_seq.get(read.ref_seq_name)
            if not reads:
                reads = list()
                self._sorted_reads_by_ref_seq[read.ref_seq_name] = reads

            reads.append(read)

        self._sorted_ref_seqs = sorted(self._sorted_reads_by_ref_seq.keys())
        for reads in self._sorted_reads_by_ref_seq.values():
            reads.sort(key=lambda read: read.start_coord)

    def _write_header(self, f_handle):
        hd_values = ['@HD', 'VN:1.4', 'SO:coordinate']
        f_handle.write('{}\n'.format('\t'.join(hd_values)))
        for ref_seq in self._sorted_ref_seqs:
            reads = self._sorted_reads_by_ref_seq[ref_seq]
            first_read = reads[0]
            sq_values = [
                '@SQ', 'SN:{}'.format(ref_seq),
                'LN:{}'.format(first_read.ref_seq_len)
            ]
            f_handle.write('{}\n'.format('\t'.join(sq_values)))


def _clip_length_from_cigar_ends(clip_length, cigar):
    first_op = cigar[0]
    if first_op[1] != 'M':
        return 'cannot clip without initial M cigar operation: {}'.format(
            cigar)
    if first_op[0] <= clip_length:
        return 'not enough initial M to clip: {} {}'.format(clip_length, cigar)
    cigar[0][0] -= clip_length
    cigar.insert(0, [clip_length, 'H'])

    last_op = cigar[-1]
    if last_op[1] != 'M':
        return 'cannot clip without trailing M cigar operation: {}'.format(
            cigar)
    if last_op[0] <= clip_length:
        return 'not enough trailing M to clip: {} {}'.format(
            clip_length, cigar)

    cigar[-1][0] -= clip_length
    cigar.append([clip_length, 'H'])

    return None


def _string_from_cigar_ops(cigar_ops):
    cigar_s = ''
    for length, op in cigar_ops:
        cigar_s += '{}{}'.format(length, op)

    return cigar_s


def set_read_pair_from_intervals(read_1,
                                 read_2,
                                 intervals_1,
                                 intervals_2,
                                 read_length,
                                 clip_length=None,
                                 is_reversed_1=False,
                                 is_reversed_2=True):
    cigar_1 = list()
    remaining_length = read_length
    prev_end = None
    for start, end in intervals_1:
        if remaining_length == 0:
            return 'not enough read_length for all read_1 intervals'

        if prev_end is not None:
            skip_len = (start - prev_end) - 1
            cigar_1.append([skip_len, 'N'])

        length = (end - start) + 1
        if length >= remaining_length:
            cigar_1.append([remaining_length, 'M'])
            remaining_length = 0
        else:
            cigar_1.append([length, 'M'])
            remaining_length -= length

        prev_end = end

    read_1_start = intervals_1[0][0]
    if clip_length:
        read_1_start += clip_length
        clip_error = _clip_length_from_cigar_ends(clip_length, cigar_1)
        if clip_error:
            return clip_error

    read_1.start_coord = read_1_start
    read_1.cigar = _string_from_cigar_ops(cigar_1)

    read_2_start = None
    cigar_2 = list()
    remaining_length = read_length
    prev_start = None
    for start, end in reversed(intervals_2):
        if remaining_length == 0:
            return 'not enough read_length for all read_2 intervals'

        if prev_start is not None:
            skip_len = (prev_start - end) - 1
            cigar_2.append([skip_len, 'N'])

        length = (end - start) + 1
        if length >= remaining_length:
            read_2_start = end - (remaining_length - 1)
            cigar_2.append([remaining_length, 'M'])
            remaining_length = 0
        else:
            read_2_start = start
            cigar_2.append([length, 'M'])
            remaining_length -= length

        prev_start = start

    cigar_2.reverse()
    read_2_end = intervals_2[-1][-1]
    if clip_length:
        read_2_start += clip_length
        read_2_end -= clip_length
        clip_error = _clip_length_from_cigar_ends(clip_length, cigar_2)
        if clip_error:
            return clip_error

    read_2.start_coord = read_2_start
    read_2.cigar = _string_from_cigar_ops(cigar_2)

    read_1.template_len = (read_2_end - read_1_start) + 1
    read_1.is_reversed = is_reversed_1
    read_2.is_reversed = is_reversed_2
    make_read_pair(read_1, read_2)
    return None


def make_read_pair(read_1, read_2):
    read_2.ref_seq_name = read_1.ref_seq_name
    read_2.ref_seq_len = read_1.ref_seq_len
    read_2.template_name = read_1.template_name
    read_1.next_ref_seq_name = '='
    read_2.next_ref_seq_name = '='
    read_1.next_start_coord = read_2.start_coord
    read_2.next_start_coord = read_1.start_coord
    read_2.template_len = -read_1.template_len
    read_1.has_multiple_segments = True
    read_2.has_multiple_segments = True
    read_1.all_segments_aligned = ((not read_1.segment_unmapped)
                                   and (not read_2.segment_unmapped))
    read_2.all_segments_aligned = read_1.all_segments_aligned
    read_1.next_segment_unmapped = read_2.segment_unmapped
    read_2.next_segment_unmapped = read_1.segment_unmapped
    read_1.next_segment_is_reversed = read_2.is_reversed
    read_2.next_segment_is_reversed = read_1.is_reversed
    read_1.is_first_segment = True
    read_2.is_first_segment = False
    read_1.is_last_segment = False
    read_2.is_last_segment = True


def set_single_end_read_from_intervals(read,
                                       intervals,
                                       read_length,
                                       clip_length=None,
                                       is_reversed=False):
    cigar = list()
    remaining_length = read_length
    prev_end = None
    for start, end in intervals:
        if remaining_length == 0:
            return 'not enough read_length for all read intervals'

        if prev_end is not None:
            skip_len = (start - prev_end) - 1
            cigar.append([skip_len, 'N'])

        length = (end - start) + 1
        if length >= remaining_length:
            cigar.append([remaining_length, 'M'])
            remaining_length = 0
        else:
            cigar.append([length, 'M'])
            remaining_length -= length

        prev_end = end

    read_start = intervals[0][0]
    read_end = prev_end
    if clip_length:
        read_start += clip_length
        read_end -= clip_length
        clip_error = _clip_length_from_cigar_ends(clip_length, cigar)
        if clip_error:
            return clip_error

    read.start_coord = read_start
    read.cigar = _string_from_cigar_ops(cigar)
    read.template_len = 0  # only set for paired reads
    read.is_reversed = is_reversed
    read.next_ref_seq_name = '*'
    read.next_start_coord = 0
    read.has_multiple_segments = False
    read.all_segments_aligned = not read.segment_unmapped
    read.next_segment_unmapped = False
    read.next_segment_is_reversed = False
    read.is_first_segment = True
    read.is_last_segment = True
    return None


class Read(object):
    def __init__(self):
        self.ref_seq_name = None  # chromosome
        self.ref_seq_len = None  # chromosome length
        self.template_name = None
        self.start_coord = None
        self.cigar = None
        self.map_quality = 0
        self.next_ref_seq_name = '*'
        self.next_start_coord = 0
        self.template_len = 0
        self.read_sequence = '*'
        self.read_quality = '*'
        self.has_multiple_segments = False
        self.all_segments_aligned = True
        self.segment_unmapped = False
        self.next_segment_unmapped = False
        self.is_reversed = False
        self.next_segment_is_reversed = False
        self.is_first_segment = True
        self.is_last_segment = False
        self.number_of_alignments = 1

    def flags_int(self):
        flags_int = 0
        if self.has_multiple_segments:
            flags_int += 1
        if self.all_segments_aligned:
            flags_int += 2
        if self.segment_unmapped:
            flags_int += 4
        if self.next_segment_unmapped:
            flags_int += 8
        if self.is_reversed:
            flags_int += 16
        if self.next_segment_is_reversed:
            flags_int += 32
        if self.is_first_segment:
            flags_int += 64
        if self.is_last_segment:
            flags_int += 128

        return flags_int
