#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
# @file lite2.py
# @brief 
# @author Zhijie Xie
# @date 2015-11-27
# cython: language_level=2
# cython: c_string_type=str, c_string_encoding=ascii


import shutil
import sys
import time
from nogilbam cimport *
from os import mkdir, walk
from os.path import basename, splitext, join, exists
from datetime import datetime
from cython.parallel import prange
from cython import boundscheck, wraparound
from cython.operator cimport dereference as deref, preincrement as inc

from rmatspipeline_declarations cimport *


USAGE = '''usage: %(prog)s [options] arg1 arg2'''


cdef:
    size_t refer_len = 1000
    int sg_mode = 0, read_mode = 1, both_mode = 2
    char ccolon = ':', cdot = '.'
    char plus_mark = '+'
    char minus_mark = '-'
    char equal_mark = '='
    char empty_mark = ' '
    string plus_mark_str = "+"
    string minus_mark_str = "-"
    string equal_mark_str = "="
    string empty_mark_str = " "
    string NH = 'NH'
    string CHR = "chr"
    string dot_str = '.'
    string ntxID = '_novel_'
    cbool GTF_TX = False
    cbool BAM_TX = True
    int FRUNSTRANDED = 0
    int FRFIRSTSTRAND = 1
    int FRSECONDSTRAND = 2

    # enum for counting the outcome of each read from the BAMs
    int READ_USED = 0
    int READ_NOT_PAIRED = 1
    int READ_NOT_NH_1 = 2
    int READ_NOT_EXPECTED_CIGAR = 3
    int READ_NOT_EXPECTED_READ_LENGTH = 4
    int READ_NOT_EXPECTED_STRAND = 5
    int READ_EXON_NOT_MATCHED_TO_ANNOTATION = 6
    int READ_JUNCTION_NOT_MATCHED_TO_ANNOTATION = 7
    int READ_CLIPPED = 8
    int READ_ENUM_VALUE_COUNT = 9
    vector[string] READ_ENUM_NAMES = [
        "USED",
        "NOT_PAIRED",
        "NOT_NH_1",
        "NOT_EXPECTED_CIGAR",
        "NOT_EXPECTED_READ_LENGTH",
        "NOT_EXPECTED_STRAND",
        "EXON_NOT_MATCHED_TO_ANNOTATION",
        "JUNCTION_NOT_MATCHED_TO_ANNOTATION",
        "CLIPPED"
    ]

    char* se_template = '%ld\t%s\t%s\t%s\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n'
    char* mxe_template = '%ld\t%s\t%s\t%s\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n'
    char* alt35_template = '%ld\t%s\t%s\t%s\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n'
    char* ri_template = '%ld\t%s\t%s\t%s\t%c\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n'

    char* ceHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n"
    char* mxeHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n"
    char* altSSHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\n"
    char* riHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\triExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\n"

    char* count_header = 'ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n'
    char* count_tmp = '%d\t%s\t%s\t%s\t%s\t%d\t%d\n'

    char* se_count_header = 'ID\tupstream_to_target_count\ttarget_to_downstream_count\ttarget_count\tupstream_to_downstream_count\n'
    char* se_count_template = '%d\t%s\t%s\t%s\t%s\n'
    char* mxe_count_header = 'ID\tupstream_to_first_count\tfirst_to_downstream_count\tfirst_count\tupstream_to_second_count\tsecond_to_downstream_count\tsecond_count\n'
    char* mxe_count_template = '%d\t%s\t%s\t%s\t%s\t%s\t%s\n'
    char* alt35_count_header = 'ID\tacross_short_boundary_count\tlong_to_flanking_count\texclusive_to_long_count\tshort_to_flanking_count\n'
    char* alt35_count_template = '%d\t%s\t%s\t%s\t%s\n'
    char* ri_count_header = 'ID\tupstream_to_intron_count\tintron_to_downstream_count\tintron_count\tupstream_to_downstream_count\n'
    char* ri_count_template = '%d\t%s\t%s\t%s\t%s\n'

    inline int c_max(long a, long b) nogil: return a if a >= b else b
    inline int c_min(long a, long b) nogil: return a if a <= b else b

    string tmptmp = "chr1"
    string tmpgene = '\"ENSG00000124508\"'


@boundscheck(False)
@wraparound(False)
cdef void parse_gtf(str gtff, unordered_map[int,cset[string]]& geneGroup,
                    unordered_map[string,Gene]& genes,
                    unordered_map[string,SupInfo]& supple) except *:
    """TODO: Docstring for parse_gtf.
    :returns: TODO

    """
    cdef:
        size_t i, idx
        str striped, chrom, rtype, strand, eled
        str gID, tID, txID, gName, dName, dVal
        list ele, desc, group
        tuple exon, brange
        Transcript tx
        unordered_map[string,Transcript] vts
        cmap[pair[long,long],size_t].iterator imap

    gtf = open(gtff, 'r')
    for striped in (line.strip() for line in gtf): ## for each line
        if striped[0] == '#': ## comments, skip this line
            continue

        ele = striped.split('\t')
        desc = ele[8].split(';')
        gID, txID, gName = '', '', 'NA' ## init..

        for eled in desc: ## for each element of description
            eleds = eled.strip().split(' ')
            if len(eled.strip()) < 2 or len(eleds) < 2:
                continue ## probably the last description

            dName = eleds[0].upper()
            dVal = eleds[1]
            if dName == 'GENE_ID': ## it is a description for gene_id
                gID = dVal
            elif dName == 'TRANSCRIPT_ID': ## it is a description for transcript_id
                txID = dVal
            elif dName == 'GENE_NAME': ## it is a description for gene_name
                gName = dVal

        if gID == '' or txID == '': ## wrong one..
            continue ## process next line

        chrom = ele[0]

        # TODO what if chrom[0:3] != 'chr' ?
        if chrom[0:3] != 'chr':
            chrom = 'chr' + chrom

        rtype = ele[2] ## exon, intron, CDS, start_codon, stop_codon..
        brange = (int(ele[3]), int(ele[4])) # start coord, end coord, 1-base
        # group = range(brange[0]/refer_len-1, brange[1]/refer_len+1) ## groups this line could belong to
        strand = ele[6]

        # for i in group: ## for each possible group
        #     geneGroup[i].insert(gID)

        if rtype == 'exon':  ## process exon
            if genes.find(gID) != genes.end():
                ## this transcript is added already
                if genes[gID].trans.find(txID) != genes[gID].trans.end():
                    ## add exon to the existing Tx
                    genes[gID].trans[txID].exons.push_back(brange)
                else: ## first time processing this Tx
                    tx.exons = [brange,]
                    genes[gID].trans[txID] = tx ## add first exon
            else:
                tx.exons = [brange,]
                genes[gID].trans[txID] = tx
                supple[gID].set_info(gName, chrom, strand)

    cdef:
        unordered_map[string,Gene].iterator igs
        unordered_map[string,Transcript].iterator itx

    igs = genes.begin()
    while igs != genes.end():
        # TODO should we skip this?
        # TODO No. However, I don't know why it is necessary.
        # TODO If we skip it, SE will lose some count. It can't be...
        # if deref(igs).second.trans.size() == 1:
        #     inc(igs)
        #     continue
        itx = deref(igs).second.trans.begin()
        while itx != deref(igs).second.trans.end():
            # TODO should we skip this tx?
            # TODO No. RI need it.
            # if deref(itx).second.exons.size() == 1:
            #     inc(itx)
            #     continue
            sort(deref(itx).second.exons.begin(), deref(itx).second.exons.end())
            for idx in range(deref(itx).second.exons.size()):
                deref(igs).second.exon_idx[deref(itx).second.exons[idx]] = 0
                # len(first) < len(exons) and len(second) < len(exons)
                # first: long one preserved. second: short one preserved.
                deref(itx).second.first[deref(itx).second.exons[idx].first] = idx
                deref(itx).second.second[deref(itx).second.exons[idx].second] = idx

            inc(itx)

        i = 0
        deref(igs).second.idx_exon.resize(deref(igs).second.exon_idx.size())
        imap = deref(igs).second.exon_idx.begin()
        while imap != deref(igs).second.exon_idx.end():
            deref(imap).second = i
            deref(igs).second.idx_exon[i] = deref(imap).first
            i += 1
            inc(imap)

        ## for each possible group this line could belong to
        for i in range(deref(igs).second.idx_exon.front().first/refer_len,
                       deref(igs).second.idx_exon.back().second/refer_len+1):
            geneGroup[i].insert(deref(igs).first)

        inc(igs)

    return


@boundscheck(False)
@wraparound(False)
cdef cbool read_has_nh_tag_1(const BamAlignment& bread) nogil:
    cdef:
        char ttype = '0'
        int8_t int8 = 0
        int16_t int16 = 0
        int32_t int32 = 0
        uint8_t uint8 = 0
        uint16_t uint16 = 0
        uint32_t uint32 = 0
        float ff = 0

    if bread.GetTagType(NH, ttype):
        if ttype == BAM_TAG_TYPE_INT8:
            if bread.GetTag(NH, int8) == False:
                return False
            elif int8 != 1:
                return False
        elif ttype == BAM_TAG_TYPE_UINT8:
            if bread.GetTag(NH, uint8) == False:
                return False
            elif uint8 != 1:
                return False
        elif ttype == BAM_TAG_TYPE_INT16:
            if bread.GetTag(NH, int16) == False:
                return False
            elif int16 != 1:
                return False
        elif ttype == BAM_TAG_TYPE_UINT16:
            if bread.GetTag(NH, uint16) == False:
                return False
            elif uint16 != 1:
                return False
        elif ttype == BAM_TAG_TYPE_INT32:
            if bread.GetTag(NH, int32) == False:
                return False
            elif int32 != 1:
                return False
        elif ttype == BAM_TAG_TYPE_UINT32:
            if bread.GetTag(NH, uint32) == False:
                return False
            elif uint32 != 1:
                return False
        elif ttype == BAM_TAG_TYPE_FLOAT:
            if bread.GetTag(NH, ff) == False:
                return False
            elif ff != 1:
                return False
    else:
        if bread.GetTag(NH, uint8) == False:
            return False
        elif uint8 != 1:
            return False

    return True


@boundscheck(False)
@wraparound(False)
cdef int filter_read(const BamAlignment& bread, const cbool& ispaired) nogil:
    if ispaired and not bread.IsProperPair():
        return READ_NOT_PAIRED

    if not read_has_nh_tag_1(bread):
        return READ_NOT_NH_1

    return READ_USED


@boundscheck(False)
@wraparound(False)
cdef int is_bam_exonread(const vector[CigarOp]& cigar_data,
                         const int amount_clipped,
                         const int& readLength,
                         const cbool variable_read_length) nogil:
    cdef:
        int length

    if cigar_data.size() != 1 or cigar_data[0].Type != 'M':
        return READ_NOT_EXPECTED_CIGAR

    length = amount_clipped + cigar_data[0].Length
    if (not variable_read_length) and length != readLength:
        return READ_NOT_EXPECTED_READ_LENGTH

    return READ_USED


@boundscheck(False)
@wraparound(False)
cdef int is_bam_multijunc(const vector[CigarOp]& cigar_data,
                          const int amount_clipped,
                          const int readLength,
                          const cbool variable_read_length) nogil:
    """TODO: Docstring for is_bam_multijunc.
    :returns: TODO

    """
    cdef:
        size_t i
        size_t num = cigar_data.size()
        int length = amount_clipped

    if num < 3 or num % 2 == 0:
        return READ_NOT_EXPECTED_CIGAR

    for i in range(num):
        if i % 2 == 1:
            if cigar_data[i].Type != 'N':
                return READ_NOT_EXPECTED_CIGAR
        else:
            if cigar_data[i].Type != 'M':
                return READ_NOT_EXPECTED_CIGAR

            length += cigar_data[i].Length

    if (not variable_read_length) and length != readLength:
        return READ_NOT_EXPECTED_READ_LENGTH

    return READ_USED


# jS 1-based
# jE 0-based
@boundscheck(False)
@wraparound(False)
cdef cbool locate_novel(long& jS, long& jE, cset[Triad]& ntx,
                        Triad& triad, Gene& gene, cbool& novelSS, long& mil) nogil:
    cdef:
        cbool connected = False
        int leftanchored = -1, rightanchored = -1
        size_t nexon = 0, i = 0
        size_t ui, di # if ui/di == nexon, ui/di == -1.
        cset[Triad].iterator intx
        unordered_map[long,size_t].iterator uInd, dInd
        unordered_map[string,Transcript].iterator ctx

    ntx.clear()

    ctx = gene.trans.begin()
    while ctx != gene.trans.end():
        ## for each transcript in the candidate gene
        uInd = deref(ctx).second.second.find(jS)
        if uInd == deref(ctx).second.second.end():
            inc(ctx)
            continue
        ui = deref(uInd).second
        dInd = deref(ctx).second.first.find(jE+1)
        if dInd == deref(ctx).second.first.end():
            inc(ctx)
            continue
        di = deref(dInd).second
        if di - ui == 1:
            ntx.clear()
            connected = True
            break
        elif di - ui > 1:
            # TODO 0-1-2-3-4-5 1-2-4-5 1,4 0-1-4-5 1-4-5
            # TODO rbegin()
            connected = True
            triad.left = gene.exon_idx[deref(ctx).second.exons[ui]]
            triad.mid = triad.left
            triad.right = gene.exon_idx[deref(ctx).second.exons[di]]
            ntx.insert(triad)
            if ui > 0:
                triad.left = gene.exon_idx[deref(ctx).second.exons[ui-1]]
                ntx.insert(triad)
            if di < deref(ctx).second.exons.size() - 1:
                ## not the last exon
                triad.left = triad.mid
                triad.mid = triad.right
                triad.right = gene.exon_idx[deref(ctx).second.exons[di+1]]
                ntx.insert(triad)
        else: ### should not be here..
            pass

        inc(ctx)

    if not connected:
        for i in range(gene.idx_exon.size()):
            if gene.idx_exon[i].second == jS:
                leftanchored = i
            elif gene.idx_exon[i].first == jE+1:
                rightanchored = i

        if leftanchored != -1 and rightanchored != -1:
            connected = True
            triad.left = leftanchored
            triad.mid = triad.left
            triad.right = rightanchored
            ntx.insert(triad)
        elif novelSS and leftanchored != -1 and jE - jS >= mil:
            connected = True
            triad.left = leftanchored
            triad.mid = triad.left
            triad.right = -(jE+1)
            ntx.insert(triad)
        elif novelSS and rightanchored != -1 and jE - jS >= mil:
            connected = True
            triad.left = -(jS)
            triad.mid = triad.left
            triad.right = rightanchored
            ntx.insert(triad)

    return connected


@boundscheck(False)
@wraparound(False)
cdef cbool locate_exon(long mc, long mec, int& rl_jl,
                       Tetrad& tetrad, Gene& gene) nogil:
    cdef:
        int ilen
        long idxi, idxj
        cbool enclosed

    enclosed = False
    tetrad.set(-1,-1,-1,-1)
    ilen = gene.idx_exon.size()

    # TODO Truly epic optimization task.
    for i in range(ilen):
        # 0-based
        idxi = gene.idx_exon[i].first - 1

        if mec > idxi:
            # 1-based
            idxj = gene.idx_exon[i].second

            if mc > idxi:
                if mec <= idxj:
                    enclosed = True
                if mc > idxj:
                    if tetrad.first < idxj + 1:
                        tetrad.first = idxj + 1
                elif mc < idxj:
                    if tetrad.first < idxi + 1:
                        tetrad.first = idxi + 1
                    if tetrad.second > idxj - 1 or tetrad.second == -1:
                        tetrad.second = idxj - 1
                    if tetrad.second > idxj - rl_jl:
                        if mc == idxj - rl_jl + 1:
                            tetrad.first = idxj - rl_jl + 1
                            tetrad.second = tetrad.first
                        elif mc < idxj - rl_jl + 1:
                            tetrad.second = idxj - rl_jl
                else:
                    tetrad.first = idxj
                    tetrad.second = idxj
            elif mc < idxi:
                if tetrad.second > idxi - 1 or tetrad.second == -1:
                    tetrad.second = idxi - 1
                if tetrad.second > idxi - rl_jl:
                    if mc == idxi - rl_jl + 1:
                        tetrad.first = idxi - rl_jl + 1
                        tetrad.second = tetrad.first
                    elif mc < idxi - rl_jl + 1:
                        tetrad.second = idxi - rl_jl
            else:
                if mec <= idxj:
                    enclosed = True
                tetrad.first = idxi
                tetrad.second = idxi

            if mec < idxj:
                if tetrad.fourth > idxj - 1 or tetrad.fourth == -1:
                    tetrad.fourth = idxj - 1
                if tetrad.third < idxi + 1:
                    tetrad.third = idxi + 1
                if tetrad.third < idxi + rl_jl:
                    if mec == idxi + rl_jl:
                        tetrad.third = idxi + rl_jl
                        tetrad.fourth = tetrad.third
                    elif mec > idxi + rl_jl:
                        tetrad.third = idxi + rl_jl + 1
            elif mec > idxj:
                if tetrad.third < idxj + 1:
                    tetrad.third = idxj + 1
                if tetrad.third < idxj + rl_jl:
                    if mec == idxj + rl_jl:
                        tetrad.third = idxj + rl_jl
                        tetrad.fourth = tetrad.third
                    elif mec > idxj + rl_jl:
                        tetrad.third = idxj + rl_jl + 1
            else:
                tetrad.third = idxj
                tetrad.fourth = idxj
        elif mec < idxi:
            if tetrad.fourth > idxi - 1 or tetrad.fourth == -1:
                tetrad.fourth = idxi - 1
            if tetrad.second > idxi - 1 or tetrad.second == -1:
                tetrad.second = idxi - 1
            if tetrad.second > idxi - rl_jl:
                if mc == idxi - rl_jl + 1:
                    tetrad.first = idxi - rl_jl + 1
                    tetrad.second = tetrad.first
                elif mc < idxi - rl_jl + 1:
                    tetrad.second = idxi - rl_jl # + 1
            break
        else:
            tetrad.third = idxi
            tetrad.fourth = idxi

    return enclosed


# mc 0-based.
@boundscheck(False)
@wraparound(False)
cdef void locate_multi(long& mc, vector[CigarOp]& cigars, int& rl_jl,
                       cset[Triad]& ntx, string& multiread,
                       Gene& gene, char* numstr, cbool& valid, cbool& novelSS,
                       long& mil, long& mel) nogil:
    cdef:
        size_t i = 0
        cbool enclosed = False
        Triad triad
        Triad exon_junc
        Tetrad tetrad
        cbool flag = True, connected = False
        long jstart = -1, jend = -1
        cset[Triad] local_ntx

    ntx.clear()
    multiread.clear()
    exon_junc.set(mc, -1, -1) # 0-based
    for i in range(0, cigars.size()):

        if i%2 == 0:
            exon_junc.mid = exon_junc.left + cigars[i].Length # 1-based
        elif i%2 == 1:
            exon_junc.left = exon_junc.right # 0-based
            # exon_junc.left = exon_junc.mid + cigars[i].Length # 0-based
            continue

        if i < cigars.size()-1:
            exon_junc.right = exon_junc.mid + cigars[i+1].Length # 0-based
            connected = locate_novel(exon_junc.mid, exon_junc.right, local_ntx, triad, gene, novelSS, mil)
            insert_set[Triad](ntx, local_ntx.begin(), local_ntx.end())

            if valid:
                if connected:
                    jstart = exon_junc.mid # 1-based
                    jend = exon_junc.right # 0-based
                    flag = False
                else:
                    jstart = -1
                    jend = -1
                if multiread.length() == 0:
                    # TODO +1 or not.
                    enclosed = locate_exon(exon_junc.left, exon_junc.mid, rl_jl, tetrad, gene)
                    if not novelSS: # or exon_junc.left - tetrad.first < mel:
                        multiread.append(join_pair(numstr, tetrad.first, tetrad.first, ccolon))
                    else:
                        multiread.append(join_pair(numstr, exon_junc.left, exon_junc.left, ccolon))
                multiread.append(equal_mark_str+join_pair(numstr, jstart, jend, ccolon))

        elif i == cigars.size() - 1:
            # TODO +1 or not.
            enclosed = locate_exon(exon_junc.left, exon_junc.mid, rl_jl, tetrad, gene)
            if not novelSS: # or tetrad.fourth - exon_junc.mid < mel:
                multiread.append(equal_mark_str+join_pair(numstr, tetrad.fourth, tetrad.fourth, ccolon))
            else:
                multiread.append(equal_mark_str+join_pair(numstr, exon_junc.mid, exon_junc.mid, ccolon))

    if flag:
        multiread.clear()


# In an fr-firststrand paired library, the second read in the pair is
# aligned to the original strand and the first read is aligned to the
# opposite strand.
# In an fr-secondstrand paired library, the first read in the pair is
# aligned to the original strand and the second read is aligned to the
# opposite strand.
@boundscheck(False)
@wraparound(False)
cdef char check_strand_paired_fr_first_strand(const BamAlignment& bread) nogil:
    cdef:
        cbool is_rev = bread.IsReverseStrand()
        cbool is_mate_rev = bread.IsMateReverseStrand()
        cbool is_first = bread.IsFirstMate()
        cbool is_second = bread.IsSecondMate()
        cbool one_mate_rev = is_rev ^ is_mate_rev
        cbool either_first_or_second = is_first ^ is_second

    if not (one_mate_rev and either_first_or_second):
        # strand check failed
        return cdot

    if is_first:
        if is_rev:
            return plus_mark
        else:
            return minus_mark
    else:
        if is_rev:
            return minus_mark
        else:
            return plus_mark


@boundscheck(False)
@wraparound(False)
cdef char check_strand_paired_fr_second_strand(const BamAlignment& bread) nogil:
    cdef:
        cbool is_rev = bread.IsReverseStrand()
        cbool is_mate_rev = bread.IsMateReverseStrand()
        cbool is_first = bread.IsFirstMate()
        cbool is_second = bread.IsSecondMate()
        cbool one_mate_rev = is_rev ^ is_mate_rev
        cbool either_first_or_second = is_first ^ is_second

    if not (one_mate_rev and either_first_or_second):
        # strand check failed
        return cdot

    if is_first:
        if is_rev:
            return minus_mark
        else:
            return plus_mark
    else:
        if is_rev:
            return plus_mark
        else:
            return minus_mark


# The single end checks are simplified since there is just one read.
@boundscheck(False)
@wraparound(False)
cdef char check_strand_single_end_fr_first_strand(const BamAlignment& bread) nogil:
    cdef:
        cbool is_rev = bread.IsReverseStrand()

    if is_rev:
        return plus_mark
    else:
        return minus_mark


@boundscheck(False)
@wraparound(False)
cdef char check_strand_single_end_fr_second_strand(const BamAlignment& bread) nogil:
    cdef:
        cbool is_rev = bread.IsReverseStrand()

    if is_rev:
        return minus_mark
    else:
        return plus_mark


@boundscheck(False)
@wraparound(False)
cdef char check_strand(const BamAlignment& bread, const cbool& ispaired, const int& dt) nogil:
    if dt == FRFIRSTSTRAND:
        if ispaired:
            return check_strand_paired_fr_first_strand(bread)

        return check_strand_single_end_fr_first_strand(bread)

    if dt == FRSECONDSTRAND:
        if ispaired:
            return check_strand_paired_fr_second_strand(bread)

        return check_strand_single_end_fr_second_strand(bread)

    # FRUNSTRANDED is the last expected case.
    # If the library type is anything else (should never happen) the read will
    # be marked as READ_NOT_EXPECTED_STRAND.
    return cdot


@boundscheck(False)
@wraparound(False)
cdef void drop_clipping_info(const BamAlignment& bread,
                             vector[CigarOp]* no_clip_cigar_data,
                             int* amount_clipped,
                             cbool* was_clipped,
                             cbool* invalid_cigar) nogil:
    cdef:
        int i = 0
        int first_non_clip_i = 0
        int first_trailing_clip_i = 0

    amount_clipped[0] = 0
    invalid_cigar[0] = False
    no_clip_cigar_data[0].clear()

    for i in range(bread.CigarData.size()):
        if bread.CigarData[i].Type == 'S' or bread.CigarData[i].Type == 'H':
            amount_clipped[0] += bread.CigarData[i].Length
            first_non_clip_i += 1
            continue
        else:
            break

    first_trailing_clip_i = first_non_clip_i
    for i in range(first_non_clip_i, bread.CigarData.size()):
        if bread.CigarData[i].Type == 'S' or bread.CigarData[i].Type == 'H':
            break
        else:
            no_clip_cigar_data[0].push_back(bread.CigarData[i])
            first_trailing_clip_i += 1

    for i in range(first_trailing_clip_i, bread.CigarData.size()):
        if bread.CigarData[i].Type == 'S' or bread.CigarData[i].Type == 'H':
            amount_clipped[0] += bread.CigarData[i].Length
            continue
        else:
            invalid_cigar[0] = True
            break

    was_clipped[0] = no_clip_cigar_data[0].size() != bread.CigarData.size()


@boundscheck(False)
@wraparound(False)
cdef void parse_bam(long fidx, string bam,
                    unordered_map[int,cset[string]]& geneGroup,
                    unordered_map[string,Gene]& genes,
                    unordered_map[string,SupInfo]& supple,
                    unordered_map[string,vector[Triad]]& novel_juncs,
                    unordered_map[string,cmap[Tetrad,int]]& exons,
                    unordered_map[string,cmap[string,int]]& multis,
                    cbool issingle, int jld2, int readLength,
                    cbool variable_read_length, int dt, cbool& novelSS,
                    long& mil, long& mel, cbool allow_clipping,
                    vector[int64_t]& read_outcome_counts) nogil:
    """TODO: Docstring for parse_bam.
    :returns: TODO

    """
    cdef:
        cset[string].iterator cg
        cset[string] visited
        cbool ispaired = not issingle
        long mc, mec
        size_t i, j
        int rl_jl = readLength-jld2
        int estart = 0, eend = 0, jend = 0
        cbool enclosed
        Triad triad
        cset[Triad] ntx
        Transcript tx, trans
        Tetrad tetrad
        string bref_name
        unordered_map[int,string] refid2str
        int minAnchor, ilen = 0
        string multiread
        char[1000] numstr
        char strand
        cbool valid

    cdef:
        BamReader br
        BamAlignment bread
        RefVector refv
        vector[CigarOp] cigar_data_after_clipping
        cbool bread_was_clipped
        cbool drop_clipping_info_invalid_cigar
        int amount_clipped

    if not br.Open(bam):
        with gil:
            print('Fail to open {}: {}'.format(bam, br.GetErrorString()))
        return

    refv = br.GetReferenceData()
    for i in range(refv.size()):
        # TODO what if chr[0:3] != 'chr' ?
        if refv[i].RefName.substr(0,3).compare(CHR) != 0:
            refid2str[br.GetReferenceID(refv[i].RefName)] = CHR + refv[i].RefName
        else:
            refid2str[br.GetReferenceID(refv[i].RefName)] = refv[i].RefName

    while br.GetNextAlignment(bread):
        drop_clipping_info(bread, &cigar_data_after_clipping, &amount_clipped,
                           &bread_was_clipped, &drop_clipping_info_invalid_cigar)
        if drop_clipping_info_invalid_cigar:
            read_outcome_counts[READ_NOT_EXPECTED_CIGAR] += 1
            continue

        if bread_was_clipped and not allow_clipping:
            read_outcome_counts[READ_CLIPPED] += 1
            continue

        filter_outcome = filter_read(bread, ispaired)
        if filter_outcome != READ_USED:
            read_outcome_counts[filter_outcome] += 1
            continue

        strand = check_strand(bread, ispaired, dt)
        if dt != FRUNSTRANDED and strand == cdot:
            read_outcome_counts[READ_NOT_EXPECTED_STRAND] += 1
            continue

        any_exon_match = False
        any_multijunc_match = False
        exon_outcome = is_bam_exonread(
            cigar_data_after_clipping, amount_clipped, readLength,
            variable_read_length)
        multijunc_outcome = is_bam_multijunc(
            cigar_data_after_clipping, amount_clipped, readLength,
            variable_read_length)

        if exon_outcome == READ_USED:
            mc = bread.Position + 1 # position (1-based) where alignment starts
            mec = mc + cigar_data_after_clipping[0].Length - 1
            bref_name = refid2str[bread.RefID]

            visited.clear()
            for i in range(mc/refer_len, mec/refer_len+1):
                if geneGroup.find(i) == geneGroup.end():
                    continue

                cg = geneGroup[i].begin()
                while cg != geneGroup[i].end():
                    ## for each candidate gene
                    if ((supple[deref(cg)].chrom != bref_name
                         or (supple[deref(cg)].strand != strand
                             and dt != FRUNSTRANDED)
                         or visited.find(deref(cg)) != visited.end())):
                        inc(cg)
                        continue

                    visited.insert(deref(cg))
                    enclosed = locate_exon(mc, mec, rl_jl, tetrad, genes[deref(cg)])

                    if (tetrad.first != -1 and tetrad.second != -1) or\
                            (tetrad.third != -1 and tetrad.fourth != -1):
                        any_exon_match = True
                        exons[deref(cg)][tetrad] += 1

                    inc(cg)

        elif multijunc_outcome == READ_USED:
            mc = bread.Position # position (0-based) where alignment starts

            minAnchor = c_min(cigar_data_after_clipping[0].Length,
                              cigar_data_after_clipping.back().Length)
            valid = (minAnchor >= rl_jl)
            bref_name = refid2str[bread.RefID]

            visited.clear()
            estart = mc
            for i in range(0, cigar_data_after_clipping.size()):
                if i%2 == 0:
                    eend = estart + cigar_data_after_clipping[i].Length # 1-based, end of exonic part
                elif i%2 == 1:
                    estart = eend + cigar_data_after_clipping[i].Length# 0-based, start of exonic part
                    continue

                for j in range(estart/refer_len, eend/refer_len+1):
                    if geneGroup.find(j) == geneGroup.end():
                        continue

                    cg = geneGroup[j].begin()
                    while cg != geneGroup[j].end():
                        ## for each candidate gene
                        if ((supple[deref(cg)].chrom != bref_name
                             or (supple[deref(cg)].strand != strand
                                 and dt != FRUNSTRANDED)
                             or visited.find(deref(cg)) != visited.end())):
                            inc(cg)
                            continue

                        visited.insert(deref(cg))
                        locate_multi(mc, cigar_data_after_clipping, rl_jl, ntx, multiread,
                                     genes[deref(cg)], numstr, valid, novelSS,
                                     mil, mel)

                        if multiread.length() > 0:
                            any_multijunc_match = True
                            multis[deref(cg)][multiread] += 1
                        if ntx.size() != 0:
                            any_multijunc_match = True
                            novel_juncs[deref(cg)].insert(novel_juncs[deref(cg)].begin(),
                                                          ntx.begin(), ntx.end())

                        inc(cg)

        if any_exon_match or any_multijunc_match:
            read_outcome_counts[READ_USED] += 1
        elif exon_outcome == READ_USED:
            read_outcome_counts[READ_EXON_NOT_MATCHED_TO_ANNOTATION] += 1
        elif multijunc_outcome == READ_USED:
            read_outcome_counts[READ_JUNCTION_NOT_MATCHED_TO_ANNOTATION] += 1
        elif (exon_outcome == READ_NOT_EXPECTED_READ_LENGTH
              or multijunc_outcome == READ_NOT_EXPECTED_READ_LENGTH):
            read_outcome_counts[READ_NOT_EXPECTED_READ_LENGTH] += 1
        else:
            read_outcome_counts[READ_NOT_EXPECTED_CIGAR] += 1

    br.Close()

    return

@boundscheck(False)
@wraparound(False)
cdef void output_read_outcomes(const vector[vector[int64_t]]& read_outcome_counts,
                               const vector[string]& vbams, str tmp_dir,
                               str prep_prefix):
    cdef:
        vector[int64_t] aggregated_read_outcome_counts
        int64_t total_for_bam
        int64_t total

    # initialize counts to zero
    aggregated_read_outcome_counts.resize(READ_ENUM_VALUE_COUNT)

    f_name = join(tmp_dir, '{}_read_outcomes_by_bam.txt'.format(prep_prefix))
    with open(f_name, 'wt') as f_handle:
        for i in range(vbams.size()):
            f_handle.write('{}\n'.format(vbams[i]))
            total_for_bam = 0
            for j in range(READ_ENUM_VALUE_COUNT):
                f_handle.write('{}: {}\n'.format(READ_ENUM_NAMES[j],
                                                 read_outcome_counts[i][j]))
                aggregated_read_outcome_counts[j] += read_outcome_counts[i][j]
                total_for_bam += read_outcome_counts[i][j]

            f_handle.write('TOTAL_FOR_BAM: {}\n'.format(total_for_bam))

    print('')
    print('read outcome totals across all BAMs')
    total = 0
    for i in range(READ_ENUM_VALUE_COUNT):
        print('{}: {}'.format(READ_ENUM_NAMES[i],
                              aggregated_read_outcome_counts[i]))
        total += aggregated_read_outcome_counts[i]

    print('total: {}'.format(total))
    print('outcomes by BAM written to: {}'.format(f_name))
    print('')


@boundscheck(False)
@wraparound(False)
cdef void detect_novel(str bams, unordered_map[int,cset[string]]& geneGroup,
                       unordered_map[string,Gene]& genes,
                       unordered_map[string,SupInfo]& supple,
                       vector[unordered_map[string,vector[Triad]]]& novel_juncs,
                       vector[unordered_map[string,cmap[Tetrad,int]]]& exons,
                       vector[unordered_map[string,cmap[string,int]]]& multis, args):
    """TODO: Docstring for detect_novel.
    :returns: TODO

    """
    cdef:
        int fidx, dt
        cbool issingle = False, novelSS = args.novelSS
        cbool variable_read_length = args.variable_read_length
        long vlen, count = 0, nthread = args.nthread
        int readLength = args.readLength, jld2 = args.junctionLength/2
        vector[string] vbams = bams.split(',')
        long mil = args.mil
        long mel = args.mel
        cbool allow_clipping = args.allow_clipping
        vector[vector[int64_t]] read_outcome_counts

    dt = args.dt
    vlen = vbams.size()
    novel_juncs.resize(vlen)
    exons.resize(vlen)
    multis.resize(vlen)
    read_outcome_counts.resize(vlen)

    # initialize outcome counts to zero
    for i in range(vlen):
        read_outcome_counts[i].resize(READ_ENUM_VALUE_COUNT)

    if args.readtype == 'single':
        issingle = True

    for fidx in prange(vlen, schedule='static', num_threads=nthread, nogil=True):
        parse_bam(fidx, vbams[fidx], geneGroup, genes, supple, novel_juncs[fidx],
                  exons[fidx], multis[fidx], issingle, jld2, readLength,
                  variable_read_length, dt, novelSS, mil, mel, allow_clipping,
                  read_outcome_counts[fidx])

    output_read_outcomes(read_outcome_counts, vbams, args.tmp, args.prep_prefix)


@boundscheck(False)
@wraparound(False)
cdef void buildup_graph(const string& gID, Gene& gene,
                        vector[vector[Triad]]& novel_j, cbool& novelSS):
    """TODO: Docstring for buildup_graphs.
    :returns: TODO

    """
    cdef:
        size_t i, j, num_exon = 0, left, mid, right
        unordered_map[string,Transcript].iterator itx 
        vector[pair[long,long]].iterator iexon
        Transcript tx

    num_exon = gene.exon_idx.size()
    if num_exon == 0:
        return

    gene.sg.resize(num_exon)
    for i in range(num_exon):
        gene.sg[i].resize(num_exon)

    itx = gene.trans.begin()
    while itx != gene.trans.end():
        num_exon = deref(itx).second.exons.size()

        if num_exon >= 2:
            left = gene.exon_idx[deref(itx).second.exons[0]]
            right = gene.exon_idx[deref(itx).second.exons[1]]
            gene.sg[left][right].insert(pair[size_t,cbool](left, GTF_TX))

            if num_exon > 2:
                for i in range(deref(itx).second.exons.size()-2):
                    left = gene.exon_idx[deref(itx).second.exons[i]]
                    mid = gene.exon_idx[deref(itx).second.exons[i+1]]
                    right = gene.exon_idx[deref(itx).second.exons[i+2]]
                    gene.sg[mid][right].insert(pair[size_t,cbool](left, GTF_TX))

        inc(itx)

    for i in range(novel_j.size()):
        for j in range(novel_j[i].size()):
            if novelSS and novel_j[i][j].left == novel_j[i][j].mid:
                tx.exons.clear()
                tx.exons.push_back(gene.idx_exon[novel_j[i][j].left])
                tx.exons.push_back(gene.idx_exon[novel_j[i][j].right])
                gene.trans[gID+ntxID+num2str(novel_j[i][j].left)+\
                        num2str(novel_j[i][j].right)] = tx
            elif novelSS and novel_j[i][j].left != novel_j[i][j].mid:
                tx.exons.clear()
                tx.exons.push_back(gene.idx_exon[novel_j[i][j].left])
                tx.exons.push_back(gene.idx_exon[novel_j[i][j].mid])
                tx.exons.push_back(gene.idx_exon[novel_j[i][j].right])
                gene.trans[gID+ntxID+num2str(novel_j[i][j].left)+\
                        num2str(novel_j[i][j].mid)+num2str(novel_j[i][j].right)] = tx
            gene.sg[novel_j[i][j].mid]\
                    [novel_j[i][j].right].insert(
                            pair[size_t,cbool](novel_j[i][j].left, BAM_TX))

    return


@boundscheck(False)
@wraparound(False)
cdef size_t convert_idx(vector[pair[size_t, size_t]]& offset, const size_t idx):
    cdef:
        size_t k
        vector[pair[size_t, size_t]].reverse_iterator i

    i = offset.rbegin()
    while i != offset.rend():
        if idx >= deref(i).first:
            return idx + deref(i).second
        inc(i)

    return idx


@boundscheck(False)
@wraparound(False)
cdef cbool expand_graph(const string& gID, Gene& gene,
                        vector[vector[Triad]]& novel_ss, cbool& novelSS,
                        long& mel):
    """TODO: Update idx_exon, exon_idx, sg, but not transcript.
    :returns: TODO

    """
    cdef:
        size_t i, j, k, curr = 0, pos = 0
        size_t left, mid, right, num_novel = 0, num_known = 0
        vector[pair[size_t, size_t]] offset
        vector[Triad] leftanchor, rightanchor
        cmap[Triad, vector[pair[size_t, size_t]]].iterator c
        cmap[Triad, vector[pair[size_t, size_t]]] candidates
        vector[Tetrad] both
        Tetrad tetrad
        cset[pair[size_t,cbool]] new_set
        cset[pair[size_t,cbool]].iterator iset
        cbool updated = False
        pair[size_t, size_t] tran
        pair[long, long] exon
        pair[size_t, cbool] sgele
        vector[pair[long,long]] local_idx_exon
        cmap[size_t, cset[pair[long, long]]] bucket # bucket for novel exons.
        cmap[size_t, cset[pair[long, long]]].iterator ibucket
        cmap[size_t, cset[pair[long, long]]].reverse_iterator ribucket
        vector[pair[long,long]].iterator iexon
        unordered_map[string,Transcript].iterator ctx
        long e_length = 0
        Transcript tx
        size_t numtx = 0

    if not novelSS:
        return False

    numtx = gene.trans.size()
    num_known = gene.idx_exon.size()
    testdic = gene.exon_idx

    # TODO note: leftanchor, rightanhor: left,mid,right could < 0.
    for i in range(novel_ss.size()):
        for j in range(novel_ss[i].size()):
            if novel_ss[i][j].left < -1: # TODO Should it be zero?
                rightanchor.push_back(novel_ss[i][j])
            elif novel_ss[i][j].right < -1: # TODO Should it be zero?
                leftanchor.push_back(novel_ss[i][j])

    # TODO note: leftanchor, rightanhor: left,mid,right could < 0.
    ctx = gene.trans.begin()
    while ctx != gene.trans.end():
        for i in range(deref(ctx).second.exons.size()):
            # leftanchor.
            for j in range(leftanchor.size()):
                e_length = deref(ctx).second.exons[i].second + leftanchor[j].right + 1

                if e_length > 0 and e_length < mel:
                    tran.first = gene.exon_idx[deref(ctx).second.exons[i]]
                    if i < deref(ctx).second.exons.size()-1:
                        tran.second = gene.exon_idx[deref(ctx).second.exons[i+1]]
                    else:
                        tran.second = -1
                    candidates[leftanchor[j]].push_back(tran)
                    exon.first = -leftanchor[j].right
                    exon.second = deref(ctx).second.exons[i].second
                    iexon = lower_bound(gene.idx_exon.begin(), gene.idx_exon.end(), exon)
                    j = iexon - gene.idx_exon.begin()
                    if gene.exon_idx.find(exon) == gene.exon_idx.end():
                        bucket[j].insert(exon)

            # rightanchor.
            for j in range(rightanchor.size()):
                e_length = -rightanchor[j].left - deref(ctx).second.exons[i].first + 1

                if e_length > 0 and e_length < mel:
                    tran.second = gene.exon_idx[deref(ctx).second.exons[i]]
                    if i > 0:
                        tran.first = gene.exon_idx[deref(ctx).second.exons[i-1]]
                    else:
                        tran.first = -1
                    candidates[rightanchor[j]].push_back(tran)
                    exon.first = deref(ctx).second.exons[i].first
                    exon.second = -rightanchor[j].left
                    iexon = lower_bound(gene.idx_exon.begin(), gene.idx_exon.end(), exon)
                    j = iexon - gene.idx_exon.begin()
                    if gene.exon_idx.find(exon) == gene.exon_idx.end():
                        bucket[j].insert(exon)

        inc(ctx)

    # TODO note: candidates: left, mid, right could < 0.
    for i in range(leftanchor.size()):
        for j in range(rightanchor.size()):
            e_length = -rightanchor[j].left + leftanchor[i].right
            if e_length > 0 and e_length < mel:
                tetrad.set(leftanchor[i].left, -leftanchor[i].right,
                           -rightanchor[j].left, rightanchor[j].right)
                both.push_back(tetrad)

    for i in range(both.size()):
        exon.first = both[i].second
        exon.second = both[i].third
        iexon = lower_bound(gene.idx_exon.begin(), gene.idx_exon.end(), exon)
        j = iexon - gene.idx_exon.begin()
        if gene.exon_idx.find(exon) == gene.exon_idx.end():
            bucket[j].insert(exon)

    num_novel = 0
    ibucket = bucket.begin()
    while ibucket != bucket.end():
        num_novel += deref(ibucket).second.size()
        inc(ibucket)

    # Backward. Updating idx_exon.
    ribucket = bucket.rbegin()
    while ribucket != bucket.rend():
        pos = deref(ribucket).first

        # Expanding idx_exon.
        gene.idx_exon.insert(gene.idx_exon.begin() + pos,
                             deref(ribucket).second.begin(),
                             deref(ribucket).second.end())

        # Expanding sg.
        for i in range(gene.sg.size()):
            gene.sg[i].insert(gene.sg[i].begin() + pos, deref(ribucket).second.size(),
                              cset[pair[size_t,cbool]]())

        inc(ribucket)

    # Expanding sg.
    ribucket = bucket.rbegin()
    while ribucket != bucket.rend():
        pos = deref(ribucket).first
        gene.sg.insert(gene.sg.begin() + pos, deref(ribucket).second.size(),
                       vector[cset[pair[size_t,cbool]]](num_novel + num_known))
        inc(ribucket)

    # Updating exon_idx.
    for i in range(gene.idx_exon.size()):
        gene.exon_idx[gene.idx_exon[i]] = i

    # Updating sg.
    curr = 0
    ibucket = bucket.begin()
    while ibucket != bucket.end():
        curr += deref(ibucket).second.size()
        offset.push_back(pair[size_t,size_t](deref(ibucket).first, curr))
        inc(ibucket)

    for i in range(gene.sg.size()-1):
        for j in range(i+1, gene.sg.size()):
            iset = gene.sg[i][j].begin()
            updated = False
            new_set.clear()

            # workaround to erase iset.
            while iset != gene.sg[i][j].end():
                pos = convert_idx(offset, deref(iset).first)
                if pos != deref(iset).first:
                    updated = True
                    sgele.first = pos
                    sgele.second = deref(iset).second
                    new_set.insert(sgele)
                else:
                    new_set.insert(deref(iset))
                inc(iset)
            if updated:
                gene.sg[i][j] = new_set

    # TODO update using exon_idx.
    c = candidates.begin()
    while c != candidates.end():
        # leftanchor
        if deref(c).first.right < 0:
            left = convert_idx(offset, deref(c).first.left)
            for i in range(deref(c).second.size()):
                mid = convert_idx(offset, deref(c).second[i].first)
                exon.first = -deref(c).first.right
                exon.second = gene.idx_exon[mid].second
                mid = gene.exon_idx[exon]
                if deref(c).second[i].second == -1:
                    right = mid
                    mid = left
                    tx.exons.clear()
                    tx.exons.push_back(gene.idx_exon[left])
                    tx.exons.push_back(gene.idx_exon[right])
                    gene.trans[gID+ntxID+num2str(left)+num2str(right)] = tx
                else:
                    right = convert_idx(offset, deref(c).second[i].second)
                    gene.sg[left][mid].insert(pair[size_t,cbool](left, BAM_TX))
                    tx.exons.clear()
                    tx.exons.push_back(gene.idx_exon[left])
                    tx.exons.push_back(gene.idx_exon[mid])
                    tx.exons.push_back(gene.idx_exon[right])
                    gene.trans[gID+ntxID+num2str(left)+num2str(mid)+num2str(right)] = tx
                gene.sg[mid][right].insert(pair[size_t,cbool](left, BAM_TX))

        # rightanchor
        elif deref(c).first.left < 0:
            right = convert_idx(offset, deref(c).first.right)
            for i in range(deref(c).second.size()):
                mid = convert_idx(offset, deref(c).second[i].second)
                exon.first = gene.idx_exon[mid].first
                exon.second = -deref(c).first.left
                mid = gene.exon_idx[exon]
                if deref(c).second[i].first == -1:
                    left = mid
                    tx.exons.clear()
                    tx.exons.push_back(gene.idx_exon[left])
                    tx.exons.push_back(gene.idx_exon[right])
                    gene.trans[gID+ntxID+num2str(left)+num2str(right)] = tx
                else:
                    left = convert_idx(offset, deref(c).second[i].first)
                    gene.sg[left][mid].insert(pair[size_t,cbool](left, BAM_TX))
                    tx.exons.clear()
                    tx.exons.push_back(gene.idx_exon[left])
                    tx.exons.push_back(gene.idx_exon[mid])
                    tx.exons.push_back(gene.idx_exon[right])
                    gene.trans[gID+ntxID+num2str(left)+num2str(mid)+num2str(right)] = tx
                gene.sg[mid][right].insert(pair[size_t,cbool](left, BAM_TX))

        inc(c)

    for i in range(both.size()):
        left = convert_idx(offset, both[i].first)
        exon.first = both[i].second
        exon.second = both[i].third
        mid = gene.exon_idx[exon]
        right = convert_idx(offset, both[i].fourth)
        gene.sg[left][mid].insert(pair[size_t,cbool](left, BAM_TX))
        gene.sg[mid][right].insert(pair[size_t,cbool](left, BAM_TX))

        tx.exons.clear()
        tx.exons.push_back(gene.idx_exon[left])
        tx.exons.push_back(gene.idx_exon[mid])
        tx.exons.push_back(gene.idx_exon[right])
        gene.trans[gID+ntxID+num2str(left)+num2str(mid)+num2str(right)] = tx

    if numtx < gene.trans.size():
        return True
    else:
        return False


@boundscheck(False)
@wraparound(False)
cdef void check_edges(const cset[pair[size_t, cbool]]& edges, cbool* found_any, cbool* all_novel):
    cdef cset[pair[size_t, cbool]].iterator iset

    # [0] is used because cython does not support *found_any = False
    found_any[0] = False
    all_novel[0] = True
    iset = edges.begin()
    while iset != edges.end():
        found_any[0] = True
        if deref(iset).second != BAM_TX:
            all_novel[0] = False
            return

        inc(iset)


@boundscheck(False)
@wraparound(False)
cdef void sm_inclen(const long& ts, const long& te, const long& us,
                    const long& ue, const long& ds, const long& de,
                    Tetrad *res, const int& jld2, const int& rl, const int& rl_jl):
    cdef:
        int tlen = te-ts
        int ulen = c_min(ue-us, jld2)
        int dlen = c_min(de-ds, jld2)

    deref(res).first = rl-2*rl_jl+1+c_min(tlen, rl-2*rl_jl+1)
    deref(res).second = rl-2*rl_jl+1
    deref(res).third = deref(res).first+c_max(0, tlen-rl+1)
    deref(res).fourth = deref(res).second


@boundscheck(False)
@wraparound(False)
cdef void ms_inclen(const long& ts, const long& te,
                    const long& ss, const long& se,
                    const long& us, const long& ue,
                    const long& ds, const long& de,
                    Tetrad *res, const int& jld2, const int& rl, const int& rl_jl):
    cdef:
        int tlen = te-ts
        int slen = se-ss
        int ulen = c_min(ue-us, jld2)
        int dlen = c_min(de-ds, jld2)

    deref(res).first = rl-2*rl_jl+1+c_min(tlen, rl-2*rl_jl+1)
    deref(res).second = rl-2*rl_jl+1+c_min(slen, rl-2*rl_jl+1)
    deref(res).third = deref(res).first+c_max(0, tlen-rl+1)
    deref(res).fourth = deref(res).second+c_max(0, slen-rl+1)


# TODO multi-thread.
@boundscheck(False)
@wraparound(False)
cdef void detect_se_mxe(const string& gID, Gene& gene, SupInfo& supInfo,
                        pair[long,long]& exon, size_t idx,
                        cmap[SE_key,SE_info]& junction_se,
                        cmap[MXE_key,MXE_info]& junction_mxe,
                        int& jld2, int& rl, int& rl_jl,
                        const cbool includes_novel_ss):
    """TODO: Docstring for detect_se_mxe.

    :arg1: TODO
    :returns: TODO

    """
    cdef:
        int iid
        SE_key se_key
        MXE_key mxe_key
        pair[size_t,cbool] gtf_pair, bam_pair
        size_t left, right, mid, skip_left, skip_right
        Tetrad len_pair
        cmap[SE_key,SE_info].iterator ise
        cmap[MXE_key,MXE_info].iterator imxe
        cbool found_any, all_novel_left, all_novel_right, all_novel_skip
        cbool event_tx_type, increased_inc_skp_len, converts_to_novel_ss
        cbool converts_to_annotated, maintains_tx_type
        cbool found_mid_gtf_transcript, found_mid_bam_transcript
        cbool found_mid_transcript, found_idx_gtf_transcript
        cbool found_idx_bam_transcript, found_idx_transcript

    gtf_pair.second = GTF_TX
    bam_pair.second = BAM_TX
    se_key.third = exon.first
    se_key.fourth = exon.second
    se_key.chrom = supInfo.chrom
    mxe_key.third = exon.first
    mxe_key.fourth = exon.second
    mxe_key.chrom = supInfo.chrom

    # TODO longest ASE.

    for left in range(idx):
        check_edges(gene.sg[left][idx], &found_any, &all_novel_left)
        if not found_any:
            continue

        gtf_pair.first = left
        bam_pair.first = left
        se_key.first = gene.idx_exon[left].second
        mxe_key.first = gene.idx_exon[left].second
        for right in range(idx+1, gene.exon_idx.size()):
            check_edges(gene.sg[idx][right], &found_any, &all_novel_right)
            if not found_any:
                continue

            # TODO Should we make sure that left-mid-right is in the same transcript.
            # According to group meeting, no.
            # if gene.sg[idx][right].find(gtf_pair) != gene.sg[idx][right].end() or\
            #         gene.sg[idx][right].find(bam_pair) != gene.sg[idx][right].end():
            se_key.second = gene.idx_exon[right].first
            mxe_key.second = gene.idx_exon[right].first

            sm_inclen(exon.first-1, exon.second, gene.idx_exon[left].first-1,
                      se_key.first, se_key.second-1, gene.idx_exon[right].second,
                      &len_pair, jld2, rl, rl_jl)

            for skip_left in range(idx):
                if gene.idx_exon[skip_left].second != gene.idx_exon[left].second:
                    continue

                for skip_right in range(idx+1, gene.idx_exon.size()):
                    if gene.idx_exon[skip_right].first != gene.idx_exon[right].first:
                        continue

                    check_edges(gene.sg[skip_left][skip_right], &found_any, &all_novel_skip)
                    if not found_any:
                        continue

                    event_tx_type = GTF_TX
                    if ((not includes_novel_ss)
                        and (all_novel_left
                             or all_novel_right
                             or all_novel_skip)):
                        event_tx_type = BAM_TX

                    ise = junction_se.find(se_key)
                    if ise == junction_se.end():
                        iid = junction_se.size()
                        junction_se[se_key].set(iid, gID, supInfo,
                                                exon.first-1, exon.second,
                                                gene.idx_exon[left].first-1,
                                                se_key.first, se_key.second-1,
                                                gene.idx_exon[right].second,
                                                len_pair.first, len_pair.second,
                                                len_pair.third, len_pair.fourth,
                                                event_tx_type,
                                                includes_novel_ss)
                    else:
                        increased_inc_skp_len = (
                            len_pair.first > deref(ise).second.inc_len
                            or (len_pair.first == deref(ise).second.inc_len
                                and len_pair.second > deref(ise).second.skp_len)
                        )
                        converts_to_novel_ss = (
                            includes_novel_ss
                            and (not deref(ise).second.includes_novel_ss)
                        )
                        converts_to_annotated = (
                            event_tx_type == GTF_TX
                            and deref(ise).second.txtype == BAM_TX
                        )
                        maintains_tx_type = (
                            event_tx_type == deref(ise).second.txtype
                        )
                        if ((not converts_to_novel_ss)
                            and (converts_to_annotated
                                 or (increased_inc_skp_len
                                     and maintains_tx_type))):
                            deref(ise).second.set(-1, gID, supInfo,
                                                  exon.first-1, exon.second,
                                                  gene.idx_exon[left].first-1,
                                                  se_key.first, se_key.second-1,
                                                  gene.idx_exon[right].second,
                                                  len_pair.first, len_pair.second,
                                                  len_pair.third, len_pair.fourth,
                                                  event_tx_type,
                                                  includes_novel_ss)

                    break
                else:
                    continue

                break

            # TODO junction set.
            for mid in range(idx, gene.idx_exon.size()):
                # TODO mutually exclusive?
                if gene.idx_exon[mid].first <= exon.second:
                    continue

                mxe_key.mxe_first = gene.idx_exon[mid].first
                mxe_key.mxe_second = gene.idx_exon[mid].second
                ms_inclen(exon.first-1, exon.second, gene.idx_exon[mid].first-1, gene.idx_exon[mid].second,
                          gene.idx_exon[left].first-1, mxe_key.first, mxe_key.second-1, gene.idx_exon[right].second,
                          &len_pair, jld2, rl, rl_jl)
                found_mid_gtf_transcript = (gene.sg[mid][right].find(gtf_pair)
                                            != gene.sg[mid][right].end())
                found_mid_bam_transcript = (gene.sg[mid][right].find(bam_pair)
                                            != gene.sg[mid][right].end())
                found_mid_transcript = (found_mid_gtf_transcript
                                        or found_mid_bam_transcript)
                found_idx_gtf_transcript = (gene.sg[idx][right].find(gtf_pair)
                                            != gene.sg[idx][right].end())
                found_idx_bam_transcript = (gene.sg[idx][right].find(bam_pair)
                                            != gene.sg[idx][right].end())
                found_idx_transcript = (found_idx_gtf_transcript
                                        or found_idx_bam_transcript)
                if not (found_mid_transcript and found_idx_transcript):
                    continue

                # all_novel_left and all_novel_right do not need to be
                # checked here since found_idx_bam_transcript would be
                # True if either all_novel_left or all_novel_right was.
                event_tx_type = GTF_TX
                if ((not includes_novel_ss)
                    and (found_mid_bam_transcript
                         or found_idx_bam_transcript)):
                    event_tx_type = BAM_TX

                imxe = junction_mxe.find(mxe_key)
                # TODO plus_mark: se_key.first,second == mxe_key.first.second
                if supInfo.strand == minus_mark:
                    ms_inclen(gene.idx_exon[mid].first-1, gene.idx_exon[mid].second, exon.first-1, exon.second,
                              gene.idx_exon[left].first-1, mxe_key.first, mxe_key.second-1, gene.idx_exon[right].second,
                              &len_pair, jld2, rl, rl_jl)
                # TODO mark.
                if imxe == junction_mxe.end():
                    iid = junction_mxe.size()
                    junction_mxe[mxe_key].set(iid, gID, supInfo,
                                              exon.first-1, exon.second,
                                              gene.idx_exon[mid].first-1, gene.idx_exon[mid].second,
                                              gene.idx_exon[left].first-1, mxe_key.first,
                                              mxe_key.second-1, gene.idx_exon[right].second,
                                              len_pair.first, len_pair.second,
                                              len_pair.third, len_pair.fourth,
                                              event_tx_type, includes_novel_ss)
                else:
                    increased_inc_skp_len = (
                        len_pair.first > deref(imxe).second.inc_len
                        or (len_pair.first == deref(imxe).second.inc_len
                            and len_pair.second > deref(imxe).second.skp_len)
                    )
                    converts_to_novel_ss = (
                        includes_novel_ss
                        and (not deref(imxe).second.includes_novel_ss)
                    )
                    converts_to_annotated = (
                        event_tx_type == GTF_TX
                        and deref(imxe).second.txtype == BAM_TX
                    )
                    maintains_tx_type = (
                        event_tx_type == deref(imxe).second.txtype
                    )
                    if ((not converts_to_novel_ss)
                        and (converts_to_annotated
                             or (increased_inc_skp_len
                                 and maintains_tx_type))):
                        deref(imxe).second.set(-1, gID, supInfo,
                                               exon.first-1, exon.second,
                                               gene.idx_exon[mid].first-1, gene.idx_exon[mid].second,
                                               gene.idx_exon[left].first-1, mxe_key.first,
                                               mxe_key.second-1, gene.idx_exon[right].second,
                                               len_pair.first, len_pair.second,
                                               len_pair.third, len_pair.fourth,
                                               event_tx_type,
                                               includes_novel_ss)


@boundscheck(False)
@wraparound(False)
cdef void alt_inclen(const long& ls, const long& le, const long& ss,
                     const long& se, const long& fs, const long& fe,
                     Tetrad *res, const int& jld2, const int& rl, const int& rl_jl):
    cdef:
        int llen = c_min(le-ls, jld2)
        int slen = c_min(se-ss, jld2)
        int flen = c_min(fe-fs, jld2)
        int alen = le-ls-(se-ss)

    deref(res).first = rl-2*rl_jl+1+c_min(alen, rl-2*rl_jl+1)
    deref(res).second = rl-2*rl_jl+1
    deref(res).third = deref(res).first+c_max(0, alen-rl+1)
    deref(res).fourth = deref(res).second


@boundscheck(False)
@wraparound(False)
cdef void update_alt35_right_flank_event(const string& gID, const Gene& gene,
                                         const SupInfo& supInfo,
                                         const pair[long, long]& exon,
                                         const size_t idx, const size_t i,
                                         const size_t j, const Tetrad& len_pair,
                                         const cbool all_novel_i,
                                         const cbool all_novel_j,
                                         const cbool includes_novel_ss,
                                         const ALT35_key& alt35_key,
                                         cmap[ALT35_key, ALT35_info]& junction_35):
    cdef:
        cmap[ALT35_key, ALT35_info].iterator ialt35
        int iid
        cbool event_tx_type, increased_inc_skp_len, converts_to_novel_ss
        cbool converts_to_annotated, maintains_tx_type

    event_tx_type = GTF_TX
    if ((not includes_novel_ss)
        and (all_novel_i or all_novel_j)):
        event_tx_type = BAM_TX

    ialt35 = junction_35.find(alt35_key)
    if ialt35 == junction_35.end():
        iid = junction_35.size()
        junction_35[alt35_key].set(iid, gID, supInfo,
                                   gene.idx_exon[i].first-1, alt35_key.second,
                                   gene.idx_exon[i].first-1, alt35_key.first,
                                   exon.first-1, exon.second,
                                   len_pair.first, len_pair.second,
                                   len_pair.third, len_pair.fourth,
                                   event_tx_type,
                                   includes_novel_ss)
    else:
        increased_inc_skp_len = (
            len_pair.first > deref(ialt35).second.inc_len
            or (len_pair.first == deref(ialt35).second.inc_len
                and len_pair.second > deref(ialt35).second.skp_len)
        )
        converts_to_novel_ss = (
            includes_novel_ss
            and (not deref(ialt35).second.includes_novel_ss)
        )
        converts_to_annotated = (
            event_tx_type == GTF_TX
            and deref(ialt35).second.txtype == BAM_TX
        )
        maintains_tx_type = (
            event_tx_type == deref(ialt35).second.txtype
        )
        if ((not converts_to_novel_ss)
            and (converts_to_annotated
                 or (increased_inc_skp_len and maintains_tx_type))):
            deref(ialt35).second.set(-1, gID, supInfo,
                                     gene.idx_exon[i].first-1, alt35_key.second,
                                     gene.idx_exon[i].first-1, alt35_key.first,
                                     exon.first-1, exon.second,
                                     len_pair.first, len_pair.second,
                                     len_pair.third, len_pair.fourth,
                                     event_tx_type,
                                     includes_novel_ss)


@boundscheck(False)
@wraparound(False)
cdef void update_alt35_left_flank_event(const string& gID, const Gene& gene,
                                        const SupInfo& supInfo,
                                        const pair[long, long]& exon,
                                        const size_t idx, const size_t i,
                                        const size_t j, const Tetrad& len_pair,
                                        const cbool all_novel_i,
                                        const cbool all_novel_j,
                                        const cbool includes_novel_ss,
                                        const ALT35_key& alt35_key,
                                        cmap[ALT35_key, ALT35_info]& junction_35):
    cdef:
        cmap[ALT35_key, ALT35_info].iterator ialt35
        int iid
        cbool event_tx_type, increased_inc_skp_len, converts_to_novel_ss
        cbool converts_to_annotated, maintains_tx_type

    event_tx_type = GTF_TX
    if ((not includes_novel_ss)
        and (all_novel_i or all_novel_j)):
        event_tx_type = BAM_TX

    ialt35 = junction_35.find(alt35_key)
    if ialt35 == junction_35.end():
        iid = junction_35.size()
        junction_35[alt35_key].set(iid, gID, supInfo,
                                   alt35_key.second, gene.idx_exon[i].second,
                                   alt35_key.third, gene.idx_exon[i].second,
                                   exon.first-1, exon.second,
                                   len_pair.first, len_pair.second,
                                   len_pair.third, len_pair.fourth,
                                   event_tx_type,
                                   includes_novel_ss)
    else:
        increased_inc_skp_len = (
            len_pair.first > deref(ialt35).second.inc_len
            or (len_pair.first == deref(ialt35).second.inc_len
                and len_pair.second > deref(ialt35).second.skp_len)
        )
        converts_to_novel_ss = (
            includes_novel_ss
            and (not deref(ialt35).second.includes_novel_ss)
        )
        converts_to_annotated = (
            event_tx_type == GTF_TX
            and deref(ialt35).second.txtype == BAM_TX
        )
        maintains_tx_type = (
            event_tx_type == deref(ialt35).second.txtype
        )
        if ((not converts_to_novel_ss)
            and (converts_to_annotated
                 or (increased_inc_skp_len and maintains_tx_type))):
            deref(ialt35).second.set(-1, gID, supInfo,
                                     alt35_key.second, gene.idx_exon[i].second,
                                     alt35_key.third, gene.idx_exon[i].second,
                                     exon.first-1, exon.second,
                                     len_pair.first, len_pair.second,
                                     len_pair.third, len_pair.fourth,
                                     event_tx_type,
                                     includes_novel_ss)


@boundscheck(False)
@wraparound(False)
cdef void detect_alt35(const string& gID, Gene& gene, SupInfo& supInfo,
                       pair[long,long]& exon, size_t idx,
                       cmap[ALT35_key,ALT35_info]& junction_3,
                       cmap[ALT35_key,ALT35_info]& junction_5,
                       int& jld2, int& rl, int& rl_jl,
                       const cbool includes_novel_ss):
    """TODO: Docstring for detect_alt35.

    :arg1: TODO
    :returns: TODO

    """
    cdef:
        size_t i, j
        ALT35_key alt35_key
        Tetrad len_pair
        cmap[ALT35_key,ALT35_info].iterator ialt35
        cbool found_any, all_novel_i, all_novel_j

    alt35_key.chrom = supInfo.chrom

    # TODO global lookup.
    # TODO longest ASE.

    for i in range(idx):
        check_edges(gene.sg[i][idx], &found_any, &all_novel_i)
        if not found_any:
            continue

        for j in range(i+1, idx):
            check_edges(gene.sg[j][idx], &found_any, &all_novel_j)
            same_start = gene.idx_exon[i].first == gene.idx_exon[j].first
            if (not found_any) or (not same_start):
                continue

            alt35_key.first = gene.idx_exon[i].second
            alt35_key.second = gene.idx_exon[j].second
            alt35_key.third = exon.first-1

            alt_inclen(gene.idx_exon[i].first-1, alt35_key.second,
                       gene.idx_exon[i].first-1, alt35_key.first,
                       exon.first-1, exon.second, &len_pair, jld2, rl, rl_jl)

            if supInfo.strand == plus_mark:
                update_alt35_right_flank_event(gID, gene, supInfo, exon, idx, i,
                                               j, len_pair, all_novel_i,
                                               all_novel_j, includes_novel_ss,
                                               alt35_key, junction_5)
            elif supInfo.strand == minus_mark:
                update_alt35_right_flank_event(gID, gene, supInfo, exon, idx, i,
                                               j, len_pair, all_novel_i,
                                               all_novel_j, includes_novel_ss,
                                               alt35_key, junction_3)

    for i in range(idx+1, gene.exon_idx.size()):
        check_edges(gene.sg[idx][i], &found_any, &all_novel_i)
        if not found_any:
            continue

        for j in range(i+1, gene.exon_idx.size()):
            check_edges(gene.sg[idx][j], &found_any, &all_novel_j)
            same_end = gene.idx_exon[i].second == gene.idx_exon[j].second
            if (not found_any) or (not same_end):
                continue

            alt35_key.first = exon.second
            alt35_key.second = gene.idx_exon[i].first-1
            alt35_key.third = gene.idx_exon[j].first-1

            alt_inclen(alt35_key.second, gene.idx_exon[i].second,
                       alt35_key.third, gene.idx_exon[i].second,
                       exon.first-1, exon.second, &len_pair, jld2, rl, rl_jl)

            if supInfo.strand == plus_mark:
                update_alt35_left_flank_event(gID, gene, supInfo, exon, idx, i,
                                              j, len_pair, all_novel_i,
                                              all_novel_j, includes_novel_ss,
                                              alt35_key, junction_3)
            elif supInfo.strand == minus_mark:
                update_alt35_left_flank_event(gID, gene, supInfo, exon, idx, i,
                                              j, len_pair, all_novel_i,
                                              all_novel_j, includes_novel_ss,
                                              alt35_key, junction_5)


@boundscheck(False)
@wraparound(False)
cdef void ri_inclen(const long& rs, const long& re, const long& us,
                    const long& ue, const long& ds, const long& de, Tetrad *res,
                    const int& jld2, const int& rl, const int& rl_jl):
    cdef:
        int ulen = c_min(ue-us, jld2)
        int dlen = c_min(de-ds, jld2)
        int rilen = ds-ue

    deref(res).first = rl-2*rl_jl+1+c_min(rilen, rl-2*rl_jl+1)
    deref(res).second = rl-2*rl_jl+1
    deref(res).third = deref(res).first+c_max(0, rilen-rl+1)
    deref(res).fourth = deref(res).second


@boundscheck(False)
@wraparound(False)
cdef void detect_ri(const string& gID, Gene& gene, SupInfo& supInfo,
                    pair[long,long]& exon, size_t idx,
                    cmap[RI_key,RI_info]& junction_ri,
                    int& jld2, int& rl, int& rl_jl,
                    const cbool includes_novel_ss):
    """TODO: Docstring for detect_ri.

    :arg1: TODO
    :returns: TODO

    """
    cdef:
        int iid
        size_t i
        pair[long,long] tmp_pair
        RI_key ri_key
        Tetrad len_pair
        cmap[RI_key,RI_info].iterator iri
        cbool found_any, all_novel, event_tx_type, increased_inc_skp_len
        cbool converts_to_novel_ss, converts_to_annotated, maintains_tx_type

    ri_key.chrom = supInfo.chrom

    # TODO longest ASE.

    for i in range(idx):
        check_edges(gene.sg[i][idx], &found_any, &all_novel)
        if not found_any:
            continue

        tmp_pair.first = gene.idx_exon[i].first
        tmp_pair.second = exon.second
        ri_key.first = gene.idx_exon[i].second
        ri_key.second = exon.first-1

        if gene.exon_idx.find(tmp_pair) == gene.exon_idx.end():
            continue

        event_tx_type = GTF_TX
        if (not includes_novel_ss) and all_novel:
            event_tx_type = BAM_TX

        iri = junction_ri.find(ri_key)
        ri_inclen(gene.idx_exon[i].first-1, exon.second,
                  gene.idx_exon[i].first-1, ri_key.first,
                  ri_key.second, exon.second, &len_pair, jld2, rl, rl_jl)
        if iri == junction_ri.end():
            iid = junction_ri.size()
            junction_ri[ri_key].set(iid, gID, supInfo,
                                    gene.idx_exon[i].first-1, exon.second,
                                    gene.idx_exon[i].first-1, ri_key.first,
                                    ri_key.second, exon.second,
                                    len_pair.first, len_pair.second,
                                    len_pair.third, len_pair.fourth,
                                    event_tx_type,
                                    includes_novel_ss)
        else:
            increased_inc_skp_len = (
                len_pair.first > deref(iri).second.inc_len
                or (len_pair.first == deref(iri).second.inc_len
                    and len_pair.second > deref(iri).second.skp_len)
            )
            converts_to_novel_ss = (
                includes_novel_ss
                and (not deref(iri).second.includes_novel_ss)
            )
            converts_to_annotated = (
                event_tx_type == GTF_TX
                and deref(iri).second.txtype == BAM_TX
            )
            maintains_tx_type = (
                event_tx_type == deref(iri).second.txtype
            )
            if ((not converts_to_novel_ss)
                and (converts_to_annotated
                     or (increased_inc_skp_len and maintains_tx_type))):
                deref(iri).second.set(-1, gID, supInfo,
                                      gene.idx_exon[i].first-1, exon.second,
                                      gene.idx_exon[i].first-1, ri_key.first,
                                      ri_key.second, exon.second,
                                      len_pair.first, len_pair.second,
                                      len_pair.third, len_pair.fourth,
                                      event_tx_type,
                                      includes_novel_ss)


@boundscheck(False)
@wraparound(False)
cdef void count_se_junction(const vector[pair[long,long]]& junction_read,
                            int read_count, const SE_info& se_event,
                            SE_counts_for_event& se_counts) nogil:
    cdef:
        size_t i

    for i in range(1, junction_read.size()-1):
        if ((i == 1
             and junction_read[i].first == se_event.te
             and junction_read[i].second == se_event.ds
             and junction_read[0].second >= se_event.ts)):
            se_counts.target_to_downstream_count += read_count
            se_counts.jc_counts.inc_count += read_count
            se_counts.jcec_counts.inc_count += read_count
            return
        if ((i == junction_read.size() - 2
             and junction_read[i].first == se_event.ue
             and junction_read[i].second == se_event.ts
             and junction_read[junction_read.size()-1].first != -1
             and junction_read[junction_read.size()-1].first <= se_event.te)):
            se_counts.upstream_to_target_count += read_count
            se_counts.jc_counts.inc_count += read_count
            se_counts.jcec_counts.inc_count += read_count
            return
        if ((i < junction_read.size()-2
             and junction_read[i].first == se_event.ue
             and junction_read[i].second == se_event.ts
             and junction_read[i+1].first == se_event.te
             and junction_read[i+1].second == se_event.ds)):
            se_counts.upstream_to_target_count += read_count
            se_counts.target_to_downstream_count += read_count
            se_counts.jc_counts.inc_count += read_count
            se_counts.jcec_counts.inc_count += read_count
            return
        if ((junction_read[i].first == se_event.ue
             and junction_read[i].second == se_event.ds)):
            se_counts.upstream_to_downstream_count += read_count
            se_counts.jc_counts.skp_count += read_count
            se_counts.jcec_counts.skp_count += read_count
            return


@boundscheck(False)
@wraparound(False)
cdef void count_se_exon(const Tetrad& exon_read,
                        const int count,
                        const SE_info& se_event,
                        SE_counts_for_event& se_counts) nogil:
    if ((exon_read.first > se_event.ts
         and exon_read.fourth != -1
         and exon_read.fourth <= se_event.te)):
        se_counts.target_count += count
        se_counts.jcec_counts.inc_count += count


@boundscheck(False)
@wraparound(False)
cdef void count_se(cset[SE_info]& junction_se,
                   unordered_map[string,cmap[Tetrad,int]]& exons,
                   unordered_map[string,cmap[vector[pair[long,long]],int]]& juncs,
                   vector[SE_counts_for_event_by_bam]& se_counts,
                   const int bam_i) nogil:
    cdef:
        int idx
        cmap[Tetrad,int].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,int]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[SE_info].iterator ise = junction_se.begin()

    while ise != junction_se.end():
        idx = deref(ise).iid

        iunmap = juncs.find(deref(ise).gID)
        if iunmap != juncs.end():
            imap = deref(iunmap).second.begin()
            while imap != deref(iunmap).second.end():
                count_se_junction(deref(imap).first, deref(imap).second,
                                  deref(ise), se_counts[idx].counts[bam_i])
                inc(imap)

        iunmap2 = exons.find(deref(ise).gID)
        if iunmap2 != exons.end():
            imap2 = deref(iunmap2).second.begin()
            while imap2 != deref(iunmap2).second.end():
                count_se_exon(deref(imap2).first, deref(imap2).second,
                              deref(ise), se_counts[idx].counts[bam_i])
                inc(imap2)

        inc(ise)


@boundscheck(False)
@wraparound(False)
cdef void count_mxe_junction(const vector[pair[long,long]]& junction_read,
                             int read_count, const MXE_info& mxe_event,
                             MXE_counts_for_event& mxe_counts) nogil:
    cdef:
        size_t i

    for i in range(1, junction_read.size()-1):
        if ((i == 1
             and junction_read[i].first == mxe_event.te
             and junction_read[i].second == mxe_event.ds
             and junction_read[0].second >= mxe_event.ts)):
            mxe_counts.first_to_downstream_count += read_count
            mxe_counts.jc_counts.inc_count += read_count
            mxe_counts.jcec_counts.inc_count += read_count
            return
        if ((i == junction_read.size()-2
             and junction_read[i].first == mxe_event.ue
             and junction_read[i].second == mxe_event.ts
             and junction_read[junction_read.size()-1].first != -1
             and junction_read[junction_read.size()-1].first <= mxe_event.te)):
            mxe_counts.upstream_to_first_count += read_count
            mxe_counts.jc_counts.inc_count += read_count
            mxe_counts.jcec_counts.inc_count += read_count
            return
        if ((i < junction_read.size()-2
             and junction_read[i].first == mxe_event.ue
             and junction_read[i].second == mxe_event.ts
             and junction_read[i+1].first == mxe_event.te
             and junction_read[i+1].second == mxe_event.ds)):
            mxe_counts.upstream_to_first_count += read_count
            mxe_counts.first_to_downstream_count += read_count
            mxe_counts.jc_counts.inc_count += read_count
            mxe_counts.jcec_counts.inc_count += read_count
            return
        if ((i == 1
             and junction_read[i].first == mxe_event.se
             and junction_read[i].second == mxe_event.ds
             and junction_read[0].second >= mxe_event.ss)):
            mxe_counts.second_to_downstream_count += read_count
            mxe_counts.jc_counts.skp_count += read_count
            mxe_counts.jcec_counts.skp_count += read_count
            return
        if ((i == junction_read.size()-2
             and junction_read[i].first == mxe_event.ue
             and junction_read[i].second == mxe_event.ss
             and junction_read[junction_read.size()-1].first != -1
             and junction_read[junction_read.size()-1].first <= mxe_event.se)):
            mxe_counts.upstream_to_second_count += read_count
            mxe_counts.jc_counts.skp_count += read_count
            mxe_counts.jcec_counts.skp_count += read_count
            return
        if ((i < junction_read.size()-2
             and junction_read[i].first == mxe_event.ue
             and junction_read[i].second == mxe_event.ss
             and junction_read[i+1].first == mxe_event.se
             and junction_read[i+1].second == mxe_event.ds)):
            mxe_counts.upstream_to_second_count += read_count
            mxe_counts.second_to_downstream_count += read_count
            mxe_counts.jc_counts.skp_count += read_count
            mxe_counts.jcec_counts.skp_count += read_count
            return


@boundscheck(False)
@wraparound(False)
cdef void count_mxe_exon(const Tetrad& exon_read,
                         const int count,
                         const MXE_info& mxe_event,
                         MXE_counts_for_event& mxe_counts) nogil:
    if ((exon_read.first > mxe_event.ts
         and exon_read.fourth != -1
         and exon_read.fourth <= mxe_event.te)):
        mxe_counts.first_count += count
        mxe_counts.jcec_counts.inc_count += count
        return
    if ((exon_read.first > mxe_event.ss
           and exon_read.fourth != -1
           and exon_read.fourth <= mxe_event.se)):
        mxe_counts.second_count += count
        mxe_counts.jcec_counts.skp_count += count


@boundscheck(False)
@wraparound(False)
cdef void count_mxe(cset[MXE_info]& junction_mxe,
                    unordered_map[string,cmap[Tetrad,int]]& exons,
                    unordered_map[string,cmap[vector[pair[long,long]],int]]& juncs,
                    vector[MXE_counts_for_event_by_bam]& mxe_counts,
                    const int bam_i) nogil:
    cdef:
        int idx
        cmap[Tetrad,int].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,int]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[MXE_info].iterator imxe = junction_mxe.begin()

    while imxe != junction_mxe.end():
        idx = deref(imxe).iid

        iunmap = juncs.find(deref(imxe).gID)
        if iunmap != juncs.end():
            imap = deref(iunmap).second.begin()
            while imap != deref(iunmap).second.end():
                count_mxe_junction(deref(imap).first, deref(imap).second,
                                   deref(imxe), mxe_counts[idx].counts[bam_i])
                inc(imap)

        iunmap2 = exons.find(deref(imxe).gID)
        if iunmap2 != exons.end():
            imap2 = deref(iunmap2).second.begin()
            while imap2 != deref(iunmap2).second.end():
                count_mxe_exon(deref(imap2).first, deref(imap2).second,
                               deref(imxe), mxe_counts[idx].counts[bam_i])
                inc(imap2)

        inc(imxe)


@boundscheck(False)
@wraparound(False)
cdef void count_alt35_right_flank_junction(
    const vector[pair[long,long]]& junction_read,
    int read_count, const ALT35_info& alt35_event, int rl_jl,
    ALT35_counts_for_event& alt35_counts) nogil:

    cdef:
        size_t i

    for i in range(1, junction_read.size()-1):
        if ((junction_read[i].first == alt35_event.le
             and junction_read[i].second == alt35_event.fs)):
            alt35_counts.long_to_flanking_count += read_count
            alt35_counts.jc_counts.inc_count += read_count
            alt35_counts.jcec_counts.inc_count += read_count
            return
        if ((junction_read[i].first == alt35_event.se
               and junction_read[i].second == alt35_event.fs)):
            alt35_counts.short_to_flanking_count += read_count
            alt35_counts.jc_counts.skp_count += read_count
            alt35_counts.jcec_counts.skp_count += read_count
            return
        if ((i == junction_read.size()-2
             and junction_read[i].second != -1
             and junction_read[i].second <= alt35_event.se-rl_jl
             and junction_read[i+1].first >= alt35_event.se+rl_jl
             and junction_read[i+1].first <= alt35_event.le)):
            alt35_counts.across_short_boundary_count += read_count
            alt35_counts.jc_counts.inc_count += read_count
            alt35_counts.jcec_counts.inc_count += read_count
            return


@boundscheck(False)
@wraparound(False)
cdef void count_alt35_right_flank_exon(
    const Tetrad& exon_read,
    const int count,
    const ALT35_info& alt35_event, int rl_jl,
    ALT35_counts_for_event& alt35_counts) nogil:

    if ((exon_read.second <= alt35_event.se-rl_jl+1
         and exon_read.second != -1
         and exon_read.fourth != -1
         and exon_read.fourth <= alt35_event.le
         and exon_read.third >= alt35_event.se+rl_jl)):
        alt35_counts.across_short_boundary_count += count
        alt35_counts.jc_counts.inc_count += count
        alt35_counts.jcec_counts.inc_count += count
    if ((exon_read.first > alt35_event.se
         and exon_read.fourth != -1
         and exon_read.fourth <= alt35_event.le)):
        alt35_counts.exclusive_to_long_count += count
        alt35_counts.jcec_counts.inc_count += count


@boundscheck(False)
@wraparound(False)
cdef void count_alt35_left_flank_junction(
    const vector[pair[long,long]]& junction_read,
    int read_count, const ALT35_info& alt35_event, int rl_jl,
    ALT35_counts_for_event& alt35_counts) nogil:

    cdef:
        size_t i

    for i in range(1, junction_read.size()-1):
        if ((junction_read[i].first == alt35_event.fe
             and junction_read[i].second == alt35_event.ls)):
            alt35_counts.long_to_flanking_count += read_count
            alt35_counts.jc_counts.inc_count += read_count
            alt35_counts.jcec_counts.inc_count += read_count
            return
        if ((junction_read[i].first == alt35_event.fe
             and junction_read[i].second == alt35_event.ss)):
            alt35_counts.short_to_flanking_count += read_count
            alt35_counts.jc_counts.skp_count += read_count
            alt35_counts.jcec_counts.skp_count += read_count
            return
        if ((i == 1
             and junction_read[0].second >= alt35_event.ls
             and junction_read[0].second <= alt35_event.ss-rl_jl
             and junction_read[1].first >= alt35_event.ss+rl_jl)):
            alt35_counts.across_short_boundary_count += read_count
            alt35_counts.jc_counts.inc_count += read_count
            alt35_counts.jcec_counts.inc_count += read_count
            return


@boundscheck(False)
@wraparound(False)
cdef void count_alt35_left_flank_exon(
    const Tetrad& exon_read,
    const int count,
    const ALT35_info& alt35_event, int rl_jl,
    ALT35_counts_for_event& alt35_counts) nogil:

    if ((exon_read.first > alt35_event.ls
         and exon_read.second != -1
         and exon_read.second <= alt35_event.ss-rl_jl+1
         and exon_read.third >= alt35_event.ss+rl_jl)):
        alt35_counts.across_short_boundary_count += count
        alt35_counts.jc_counts.inc_count += count
        alt35_counts.jcec_counts.inc_count += count
    if ((exon_read.first > alt35_event.ls
         and exon_read.fourth != -1
         and exon_read.fourth <= alt35_event.ss)):
        alt35_counts.exclusive_to_long_count += count
        alt35_counts.jcec_counts.inc_count += count


@boundscheck(False)
@wraparound(False)
cdef void count_alt35(cset[ALT35_info]& junction_35,
                      unordered_map[string,cmap[Tetrad,int]]& exons,
                      unordered_map[string,cmap[vector[pair[long,long]],int]]& juncs,
                      vector[ALT35_counts_for_event_by_bam]& alt35_counts,
                      int& jld2, int& rl, const int bam_i) nogil:
    cdef:
        int idx
        int rl_jl = rl - jld2
        cmap[Tetrad,int].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,int]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[ALT35_info].iterator ialt35 = junction_35.begin()

    while ialt35 != junction_35.end():
        idx = deref(ialt35).iid

        if deref(ialt35).fs > deref(ialt35).le:
            iunmap = juncs.find(deref(ialt35).gID)
            if iunmap != juncs.end():
                imap = deref(iunmap).second.begin()
                while imap != deref(iunmap).second.end():
                    count_alt35_right_flank_junction(
                        deref(imap).first, deref(imap).second,
                        deref(ialt35), rl_jl,
                        alt35_counts[idx].counts[bam_i])
                    inc(imap)

            iunmap2 = exons.find(deref(ialt35).gID)
            if iunmap2 != exons.end():
                imap2 = deref(iunmap2).second.begin()
                while imap2 != deref(iunmap2).second.end():
                    count_alt35_right_flank_exon(
                        deref(imap2).first, deref(imap2).second,
                        deref(ialt35), rl_jl,
                        alt35_counts[idx].counts[bam_i])
                    inc(imap2)
        else:
            iunmap = juncs.find(deref(ialt35).gID)
            if iunmap != juncs.end():
                imap = deref(iunmap).second.begin()
                while imap != deref(iunmap).second.end():
                    count_alt35_left_flank_junction(
                        deref(imap).first, deref(imap).second,
                        deref(ialt35), rl_jl,
                        alt35_counts[idx].counts[bam_i])
                    inc(imap)

            iunmap2 = exons.find(deref(ialt35).gID)
            if iunmap2 != exons.end():
                imap2 = deref(iunmap2).second.begin()
                while imap2 != deref(iunmap2).second.end():
                    count_alt35_left_flank_exon(
                        deref(imap2).first, deref(imap2).second,
                        deref(ialt35), rl_jl,
                        alt35_counts[idx].counts[bam_i])
                    inc(imap2)

        inc(ialt35)


@boundscheck(False)
@wraparound(False)
cdef void count_ri_junction(const vector[pair[long,long]]& junction_read,
                            int read_count, const RI_info& ri_event,
                            int rl_jl,
                            RI_counts_for_event& ri_counts) nogil:
    cdef:
        size_t i

    for i in range(1, junction_read.size()-1):
        if ((junction_read[i].first == ri_event.ue
             and junction_read[i].second == ri_event.ds)):
            ri_counts.upstream_to_downstream_count += read_count
            ri_counts.jc_counts.skp_count += read_count
            ri_counts.jcec_counts.skp_count += read_count
            return
        if ((i == 1
             and junction_read[0].second != -1
             and junction_read[0].second <= ri_event.ds-rl_jl
             and junction_read[1].first >= ri_event.ds+rl_jl)):
            ri_counts.intron_to_downstream_count += read_count
            ri_counts.jc_counts.inc_count += read_count
            ri_counts.jcec_counts.inc_count += read_count
            return
        if ((i == junction_read.size()-2
             and junction_read[i].second != -1
             and junction_read[i].second <= ri_event.ue-rl_jl
             and junction_read[i+1].first >= ri_event.ue+rl_jl)):
            ri_counts.upstream_to_intron_count += read_count
            ri_counts.jc_counts.inc_count += read_count
            ri_counts.jcec_counts.inc_count += read_count
            return


@boundscheck(False)
@wraparound(False)
cdef void count_ri_exon(const Tetrad& exon_read,
                        const int count,
                        const RI_info& ri_event, int rl_jl,
                        RI_counts_for_event& ri_counts) nogil:
    cdef:
        cbool counted = False

    if ((exon_read.second <= ri_event.ue-rl_jl+1
         and exon_read.second != -1
         and exon_read.third >= ri_event.ue+rl_jl)):
        ri_counts.upstream_to_intron_count += count
        ri_counts.jc_counts.inc_count += count
        ri_counts.jcec_counts.inc_count += count
        counted = True
    if ((exon_read.second != -1
         and exon_read.second <= ri_event.ds-rl_jl+1
         and exon_read.third >= ri_event.ds+rl_jl)):
        ri_counts.intron_to_downstream_count += count
        if not counted:
            ri_counts.jc_counts.inc_count += count
            ri_counts.jcec_counts.inc_count += count
    if ((exon_read.first > ri_event.ue
         and exon_read.fourth != -1
         and exon_read.fourth <= ri_event.ds)):
        ri_counts.intron_count += count
        ri_counts.jcec_counts.inc_count += count


@boundscheck(False)
@wraparound(False)
cdef void count_ri(cset[RI_info]& junction_ri,
                   unordered_map[string,cmap[Tetrad,int]]& exons,
                   unordered_map[string,cmap[vector[pair[long,long]],int]]& juncs,
                   vector[RI_counts_for_event_by_bam]& ri_counts,
                   int& jld2, int& rl, const int bam_i) nogil:
    cdef:
        int idx
        long rl_jl = rl - jld2
        cmap[Tetrad,int].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,int]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[RI_info].iterator iri = junction_ri.begin()

    while iri != junction_ri.end():
        idx = deref(iri).iid

        iunmap = juncs.find(deref(iri).gID)
        if iunmap != juncs.end():
            imap = deref(iunmap).second.begin()
            while imap != deref(iunmap).second.end():
                count_ri_junction(deref(imap).first, deref(imap).second,
                                  deref(iri), rl_jl,
                                  ri_counts[idx].counts[bam_i])
                inc(imap)

        iunmap2 = exons.find(deref(iri).gID)
        if iunmap2 != exons.end():
            imap2 = deref(iunmap2).second.begin()

            while imap2 != deref(iunmap2).second.end():
                count_ri_exon(deref(imap2).first, deref(imap2).second,
                              deref(iri), rl_jl,
                              ri_counts[idx].counts[bam_i])
                inc(imap2)

        inc(iri)


@boundscheck(False)
@wraparound(False)
cdef count_occurrence(str bams, list dot_rmats_paths, str od,
                      cset[SE_info]& se, cset[MXE_info]& mxe,
                      cset[ALT35_info]& alt3, cset[ALT35_info]& alt5,
                      cset[RI_info]& ri, int sam1len, int& jld2, int& rl,
                      int& nthread, cbool stat, cbool individual_counts):
    cdef:
        size_t idx = 0
        list vbams = bams.split(',')
        int num = len(vbams), vlen, fidx, bam_i
        cset[SE_info].iterator ise = se.begin()
        cset[MXE_info].iterator imxe = mxe.begin()
        cset[ALT35_info].iterator ialt3 = alt3.begin()
        cset[ALT35_info].iterator ialt5 = alt5.begin()
        cset[RI_info].iterator iri = ri.begin()
        vector[SE_counts_for_event_by_bam] se_counts = (
            vector[SE_counts_for_event_by_bam](se.size()))
        vector[MXE_counts_for_event_by_bam] mxe_counts = (
            vector[MXE_counts_for_event_by_bam](mxe.size()))
        vector[ALT35_counts_for_event_by_bam] alt3_counts = (
            vector[ALT35_counts_for_event_by_bam](alt3.size()))
        vector[ALT35_counts_for_event_by_bam] alt5_counts = (
            vector[ALT35_counts_for_event_by_bam](alt5.size()))
        vector[RI_counts_for_event_by_bam] ri_counts = (
            vector[RI_counts_for_event_by_bam](ri.size()))
        vector[unordered_map[string,cmap[Tetrad,int]]] exons
        vector[unordered_map[string,cmap[vector[pair[long,long]],int]]] juncs

    while ise != se.end():
        idx = deref(ise).iid
        se_counts[idx].counts = vector[SE_counts_for_event](num)
        se_counts[idx].jc_lengths.inc_len = deref(ise).inc_len
        se_counts[idx].jc_lengths.skp_len = deref(ise).skp_len
        se_counts[idx].jcec_lengths.inc_len = deref(ise).inc_len_jcec
        se_counts[idx].jcec_lengths.skp_len = deref(ise).skp_len_jcec
        se_counts[idx].strand = deref(ise).supInfo.strand
        inc(ise)

    while imxe != mxe.end():
        idx = deref(imxe).iid
        mxe_counts[idx].counts = vector[MXE_counts_for_event](num)
        mxe_counts[idx].jc_lengths.inc_len = deref(imxe).inc_len
        mxe_counts[idx].jc_lengths.skp_len = deref(imxe).skp_len
        mxe_counts[idx].jcec_lengths.inc_len = deref(imxe).inc_len_jcec
        mxe_counts[idx].jcec_lengths.skp_len = deref(imxe).skp_len_jcec
        mxe_counts[idx].strand = deref(imxe).supInfo.strand
        inc(imxe)

    while ialt3 != alt3.end():
        idx = deref(ialt3).iid
        alt3_counts[idx].counts = vector[ALT35_counts_for_event](num)
        alt3_counts[idx].jc_lengths.inc_len = deref(ialt3).inc_len
        alt3_counts[idx].jc_lengths.skp_len = deref(ialt3).skp_len
        alt3_counts[idx].jcec_lengths.inc_len = deref(ialt3).inc_len_jcec
        alt3_counts[idx].jcec_lengths.skp_len = deref(ialt3).skp_len_jcec
        alt3_counts[idx].strand = deref(ialt3).supInfo.strand
        inc(ialt3)

    while ialt5 != alt5.end():
        idx = deref(ialt5).iid
        alt5_counts[idx].counts = vector[ALT35_counts_for_event](num)
        alt5_counts[idx].jc_lengths.inc_len = deref(ialt5).inc_len
        alt5_counts[idx].jc_lengths.skp_len = deref(ialt5).skp_len
        alt5_counts[idx].jcec_lengths.inc_len = deref(ialt5).inc_len_jcec
        alt5_counts[idx].jcec_lengths.skp_len = deref(ialt5).skp_len_jcec
        alt5_counts[idx].strand = deref(ialt5).supInfo.strand
        inc(ialt5)

    while iri != ri.end():
        idx = deref(iri).iid
        ri_counts[idx].counts = vector[RI_counts_for_event](num)
        ri_counts[idx].jc_lengths.inc_len = deref(iri).inc_len
        ri_counts[idx].jc_lengths.skp_len = deref(iri).skp_len
        ri_counts[idx].jcec_lengths.inc_len = deref(iri).inc_len_jcec
        ri_counts[idx].jcec_lengths.skp_len = deref(iri).skp_len_jcec
        ri_counts[idx].strand = deref(iri).supInfo.strand
        inc(iri)

    vlen = len(dot_rmats_paths)
    for fidx in prange(vlen, schedule='static', num_threads=nthread, nogil=True):
        with gil:
            bam_i = load_read(bams, dot_rmats_paths[fidx], exons, juncs)

        count_se(se, exons[bam_i], juncs[bam_i], se_counts, bam_i)
        count_mxe(mxe, exons[bam_i], juncs[bam_i], mxe_counts, bam_i)
        count_alt35(alt3, exons[bam_i], juncs[bam_i], alt3_counts, jld2, rl, bam_i)
        count_alt35(alt5, exons[bam_i], juncs[bam_i], alt5_counts, jld2, rl, bam_i)
        count_ri(ri, exons[bam_i], juncs[bam_i], ri_counts, jld2, rl, bam_i)

        exons[bam_i].clear()
        juncs[bam_i].clear()

    save_ct(od, se_counts, mxe_counts, alt3_counts, alt5_counts, ri_counts,
            sam1len, stat, individual_counts)


@boundscheck(False)
@wraparound(False)
cdef void write_se_count_rows(const vector[SE_counts_for_event_by_bam]& se_counts,
                              int sam1len, FILE* se_fp,
                              FILE* se_fp_n, FILE* se_fp_individual,
                              cbool individual_counts):
    cdef:
        SE_joined_count_strings joined_strings

    for i in range(se_counts.size()):
        se_counts[i].join_counts_across_bams(sam1len, &joined_strings)
        fprintf(se_fp, count_tmp, i,
                joined_strings.jc_inc_1.c_str(),
                joined_strings.jc_skp_1.c_str(),
                joined_strings.jc_inc_2.c_str(),
                joined_strings.jc_skp_2.c_str(),
                se_counts[i].jc_lengths.inc_len,
                se_counts[i].jc_lengths.skp_len)
        fprintf(se_fp_n, count_tmp, i,
                joined_strings.jcec_inc_1.c_str(),
                joined_strings.jcec_skp_1.c_str(),
                joined_strings.jcec_inc_2.c_str(),
                joined_strings.jcec_skp_2.c_str(),
                se_counts[i].jcec_lengths.inc_len,
                se_counts[i].jcec_lengths.skp_len)
        if individual_counts:
            fprintf(se_fp_individual, se_count_template, i,
                    joined_strings.upstream_to_target.c_str(),
                    joined_strings.target_to_downstream.c_str(),
                    joined_strings.target.c_str(),
                    joined_strings.upstream_to_downstream.c_str())


@boundscheck(False)
@wraparound(False)
cdef void write_mxe_count_rows(const vector[MXE_counts_for_event_by_bam]& mxe_counts,
                               int sam1len, FILE* mxe_fp,
                               FILE* mxe_fp_n, FILE* mxe_fp_individual,
                               cbool individual_counts):
    cdef:
        MXE_joined_count_strings joined_strings

    for i in range(mxe_counts.size()):
        mxe_counts[i].join_counts_across_bams(sam1len, &joined_strings)
        if mxe_counts[i].strand == plus_mark:
            fprintf(mxe_fp, count_tmp, i,
                    joined_strings.jc_inc_1.c_str(),
                    joined_strings.jc_skp_1.c_str(),
                    joined_strings.jc_inc_2.c_str(),
                    joined_strings.jc_skp_2.c_str(),
                    mxe_counts[i].jc_lengths.inc_len,
                    mxe_counts[i].jc_lengths.skp_len)
            fprintf(mxe_fp_n, count_tmp, i,
                    joined_strings.jcec_inc_1.c_str(),
                    joined_strings.jcec_skp_1.c_str(),
                    joined_strings.jcec_inc_2.c_str(),
                    joined_strings.jcec_skp_2.c_str(),
                    mxe_counts[i].jcec_lengths.inc_len,
                    mxe_counts[i].jcec_lengths.skp_len)
            if individual_counts:
                fprintf(mxe_fp_individual, mxe_count_template, i,
                        joined_strings.upstream_to_first.c_str(),
                        joined_strings.first_to_downstream.c_str(),
                        joined_strings.first.c_str(),
                        joined_strings.upstream_to_second.c_str(),
                        joined_strings.second_to_downstream.c_str(),
                        joined_strings.second.c_str())
        else:
            fprintf(mxe_fp, count_tmp, i,
                    joined_strings.jc_skp_1.c_str(),
                    joined_strings.jc_inc_1.c_str(),
                    joined_strings.jc_skp_2.c_str(),
                    joined_strings.jc_inc_2.c_str(),
                    mxe_counts[i].jc_lengths.inc_len,
                    mxe_counts[i].jc_lengths.skp_len)
            fprintf(mxe_fp_n, count_tmp, i,
                    joined_strings.jcec_skp_1.c_str(),
                    joined_strings.jcec_inc_1.c_str(),
                    joined_strings.jcec_skp_2.c_str(),
                    joined_strings.jcec_inc_2.c_str(),
                    mxe_counts[i].jcec_lengths.inc_len,
                    mxe_counts[i].jcec_lengths.skp_len)
            if individual_counts:
                fprintf(mxe_fp_individual, mxe_count_template, i,
                        joined_strings.upstream_to_second.c_str(),
                        joined_strings.second_to_downstream.c_str(),
                        joined_strings.second.c_str(),
                        joined_strings.upstream_to_first.c_str(),
                        joined_strings.first_to_downstream.c_str(),
                        joined_strings.first.c_str())


@boundscheck(False)
@wraparound(False)
cdef void write_alt35_count_rows(const vector[ALT35_counts_for_event_by_bam]& alt35_counts,
                                 int sam1len, FILE* alt35_fp,
                                 FILE* alt35_fp_n, FILE* alt35_fp_individual,
                                 cbool individual_counts):
    cdef:
        ALT35_joined_count_strings joined_strings

    for i in range(alt35_counts.size()):
        alt35_counts[i].join_counts_across_bams(sam1len, &joined_strings)
        fprintf(alt35_fp, count_tmp, i,
                joined_strings.jc_inc_1.c_str(),
                joined_strings.jc_skp_1.c_str(),
                joined_strings.jc_inc_2.c_str(),
                joined_strings.jc_skp_2.c_str(),
                alt35_counts[i].jc_lengths.inc_len,
                alt35_counts[i].jc_lengths.skp_len)
        fprintf(alt35_fp_n, count_tmp, i,
                joined_strings.jcec_inc_1.c_str(),
                joined_strings.jcec_skp_1.c_str(),
                joined_strings.jcec_inc_2.c_str(),
                joined_strings.jcec_skp_2.c_str(),
                alt35_counts[i].jcec_lengths.inc_len,
                alt35_counts[i].jcec_lengths.skp_len)
        if individual_counts:
            fprintf(alt35_fp_individual, alt35_count_template, i,
                    joined_strings.across_short_boundary.c_str(),
                    joined_strings.long_to_flanking.c_str(),
                    joined_strings.exclusive_to_long.c_str(),
                    joined_strings.short_to_flanking.c_str())


@boundscheck(False)
@wraparound(False)
cdef void write_ri_count_rows(const vector[RI_counts_for_event_by_bam]& ri_counts,
                              int sam1len, FILE* ri_fp,
                              FILE* ri_fp_n, FILE* ri_fp_individual,
                              cbool individual_counts):
    cdef:
        RI_joined_count_strings joined_strings

    for i in range(ri_counts.size()):
        ri_counts[i].join_counts_across_bams(sam1len, &joined_strings)
        fprintf(ri_fp, count_tmp, i,
                joined_strings.jc_inc_1.c_str(),
                joined_strings.jc_skp_1.c_str(),
                joined_strings.jc_inc_2.c_str(),
                joined_strings.jc_skp_2.c_str(),
                ri_counts[i].jc_lengths.inc_len,
                ri_counts[i].jc_lengths.skp_len)
        fprintf(ri_fp_n, count_tmp, i,
                joined_strings.jcec_inc_1.c_str(),
                joined_strings.jcec_skp_1.c_str(),
                joined_strings.jcec_inc_2.c_str(),
                joined_strings.jcec_skp_2.c_str(),
                ri_counts[i].jcec_lengths.inc_len,
                ri_counts[i].jcec_lengths.skp_len)
        if individual_counts:
            fprintf(ri_fp_individual, ri_count_template, i,
                    joined_strings.upstream_to_intron.c_str(),
                    joined_strings.intron_to_downstream.c_str(),
                    joined_strings.intron.c_str(),
                    joined_strings.upstream_to_downstream.c_str())


@boundscheck(False)
@wraparound(False)
cdef save_ct(str od, vector[SE_counts_for_event_by_bam]& se_counts,
             vector[MXE_counts_for_event_by_bam]& mxe_counts,
             vector[ALT35_counts_for_event_by_bam]& alt3_counts,
             vector[ALT35_counts_for_event_by_bam]& alt5_counts,
             vector[RI_counts_for_event_by_bam]& ri_counts,
             int sam1len, cbool stat, cbool individual_counts):
    cdef:
        size_t i
        FILE *se_fp
        FILE *mxe_fp
        FILE *alt3_fp
        FILE *alt5_fp
        FILE *ri_fp
        FILE *se_fp_n
        FILE *mxe_fp_n
        FILE *alt3_fp_n
        FILE *alt5_fp_n
        FILE *ri_fp_n
        FILE *se_fp_individual
        FILE *mxe_fp_individual
        FILE *alt3_fp_individual
        FILE *alt5_fp_individual
        FILE *ri_fp_individual
        int total

    se_fp = fopen('%s/JC.raw.input.SE.txt' % (od), 'w')
    se_fp_n = fopen('%s/JCEC.raw.input.SE.txt' % (od), 'w')
    mxe_fp = fopen('%s/JC.raw.input.MXE.txt' % (od), 'w')
    mxe_fp_n = fopen('%s/JCEC.raw.input.MXE.txt' % (od), 'w')
    alt3_fp = fopen('%s/JC.raw.input.A3SS.txt' % (od), 'w')
    alt3_fp_n = fopen('%s/JCEC.raw.input.A3SS.txt' % (od), 'w')
    alt5_fp = fopen('%s/JC.raw.input.A5SS.txt' % (od), 'w')
    alt5_fp_n = fopen('%s/JCEC.raw.input.A5SS.txt' % (od), 'w')
    ri_fp = fopen('%s/JC.raw.input.RI.txt' % (od), 'w')
    ri_fp_n = fopen('%s/JCEC.raw.input.RI.txt' % (od), 'w')
    if individual_counts:
        se_fp_individual = fopen('%s/individualCounts.SE.txt' % (od), 'w')
        mxe_fp_individual = fopen('%s/individualCounts.MXE.txt' % (od), 'w')
        alt3_fp_individual = fopen('%s/individualCounts.A3SS.txt' % (od), 'w')
        alt5_fp_individual = fopen('%s/individualCounts.A5SS.txt' % (od), 'w')
        ri_fp_individual = fopen('%s/individualCounts.RI.txt' % (od), 'w')

    fprintf(se_fp, count_header)
    fprintf(se_fp_n, count_header)
    fprintf(mxe_fp, count_header)
    fprintf(mxe_fp_n, count_header)
    fprintf(alt3_fp, count_header)
    fprintf(alt3_fp_n, count_header)
    fprintf(alt5_fp, count_header)
    fprintf(alt5_fp_n, count_header)
    fprintf(ri_fp, count_header)
    fprintf(ri_fp_n, count_header)

    if individual_counts:
        fprintf(se_fp_individual, se_count_header)
        fprintf(mxe_fp_individual, mxe_count_header)
        fprintf(alt3_fp_individual, alt35_count_header)
        fprintf(alt5_fp_individual, alt35_count_header)
        fprintf(ri_fp_individual, ri_count_header)

    write_se_count_rows(se_counts, sam1len, se_fp, se_fp_n,
                        se_fp_individual, individual_counts)
    write_mxe_count_rows(mxe_counts, sam1len, mxe_fp, mxe_fp_n,
                         mxe_fp_individual, individual_counts)
    write_alt35_count_rows(alt3_counts, sam1len, alt3_fp, alt3_fp_n,
                           alt3_fp_individual, individual_counts)
    write_alt35_count_rows(alt5_counts, sam1len, alt5_fp, alt5_fp_n,
                           alt5_fp_individual, individual_counts)
    write_ri_count_rows(ri_counts, sam1len, ri_fp, ri_fp_n,
                        ri_fp_individual, individual_counts)

    fclose(se_fp)
    fclose(se_fp_n)
    fclose(mxe_fp)
    fclose(mxe_fp_n)
    fclose(alt3_fp)
    fclose(alt3_fp_n)
    fclose(alt5_fp)
    fclose(alt5_fp_n)
    fclose(ri_fp)
    fclose(ri_fp_n)

    if individual_counts:
        fclose(se_fp_individual)
        fclose(mxe_fp_individual)
        fclose(alt3_fp_individual)
        fclose(alt5_fp_individual)
        fclose(ri_fp_individual)


@boundscheck(False)
@wraparound(False)
cdef void save_se_ase(FILE* fp, const SE_info& se):
    fprintf(fp, se_template, se.iid,
            se.gID.c_str(),
            se.supInfo.g_name.c_str(),
            se.supInfo.chrom.c_str(),
            se.supInfo.strand,
            se.ts, se.te,
            se.us, se.ue,
            se.ds, se.de)


@boundscheck(False)
@wraparound(False)
cdef void save_mxe_ase(FILE* fp, const MXE_info& mxe):
    fprintf(fp, mxe_template, mxe.iid,
            mxe.gID.c_str(),
            mxe.supInfo.g_name.c_str(),
            mxe.supInfo.chrom.c_str(),
            mxe.supInfo.strand,
            mxe.ts, mxe.te,
            mxe.ss, mxe.se,
            mxe.us, mxe.ue,
            mxe.ds, mxe.de)


@boundscheck(False)
@wraparound(False)
cdef void save_alt35_ase(FILE* fp, const ALT35_info& alt35):
    fprintf(fp, alt35_template, alt35.iid,
            alt35.gID.c_str(),
            alt35.supInfo.g_name.c_str(),
            alt35.supInfo.chrom.c_str(),
            alt35.supInfo.strand,
            alt35.ls, alt35.le,
            alt35.ss, alt35.se,
            alt35.fs, alt35.fe)


@boundscheck(False)
@wraparound(False)
cdef void save_ri_ase(FILE* fp, const RI_info& ri):
    fprintf(fp, ri_template, ri.iid,
            ri.gID.c_str(),
            ri.supInfo.g_name.c_str(),
            ri.supInfo.chrom.c_str(),
            ri.supInfo.strand,
            ri.rs, ri.re,
            ri.us, ri.ue,
            ri.ds, ri.de)


@boundscheck(False)
@wraparound(False)
cdef save_ase(string& od, cset[SE_info]& juncs_se,
              cset[MXE_info]& juncs_mxe,
              cset[ALT35_info]& juncs_alt3,
              cset[ALT35_info]& juncs_alt5,
              cset[RI_info]& juncs_ri):
    cdef:
        FILE *se_fp
        FILE *mxe_fp
        FILE *alt3_fp
        FILE *alt5_fp
        FILE *ri_fp
        FILE *se_fp_nj
        FILE *mxe_fp_nj
        FILE *alt3_fp_nj
        FILE *alt5_fp_nj
        FILE *ri_fp_nj
        FILE *se_fp_nss
        FILE *mxe_fp_nss
        FILE *alt3_fp_nss
        FILE *alt5_fp_nss
        FILE *ri_fp_nss
        cset[SE_info].iterator ise = juncs_se.begin()
        cset[MXE_info].iterator imxe = juncs_mxe.begin()
        cset[ALT35_info].iterator ialt3 = juncs_alt3.begin()
        cset[ALT35_info].iterator ialt5 = juncs_alt5.begin()
        cset[RI_info].iterator iri = juncs_ri.begin()

    se_fp = fopen('%s/fromGTF.SE.txt' % (od), 'w')
    se_fp_nj = fopen('%s/fromGTF.novelJunction.SE.txt' % (od), 'w')
    se_fp_nss = fopen('%s/fromGTF.novelSpliceSite.SE.txt' % (od), 'w')
    mxe_fp = fopen('%s/fromGTF.MXE.txt' % (od), 'w')
    mxe_fp_nj = fopen('%s/fromGTF.novelJunction.MXE.txt' % (od), 'w')
    mxe_fp_nss = fopen('%s/fromGTF.novelSpliceSite.MXE.txt' % (od), 'w')
    alt3_fp = fopen('%s/fromGTF.A3SS.txt' % (od), 'w')
    alt3_fp_nj = fopen('%s/fromGTF.novelJunction.A3SS.txt' % (od), 'w')
    alt3_fp_nss = fopen('%s/fromGTF.novelSpliceSite.A3SS.txt' % (od), 'w')
    alt5_fp = fopen('%s/fromGTF.A5SS.txt' % (od), 'w')
    alt5_fp_nj = fopen('%s/fromGTF.novelJunction.A5SS.txt' % (od), 'w')
    alt5_fp_nss = fopen('%s/fromGTF.novelSpliceSite.A5SS.txt' % (od), 'w')
    ri_fp = fopen('%s/fromGTF.RI.txt' % (od), 'w')
    ri_fp_nj = fopen('%s/fromGTF.novelJunction.RI.txt' % (od), 'w')
    ri_fp_nss = fopen('%s/fromGTF.novelSpliceSite.RI.txt' % (od), 'w')

    fprintf(se_fp, ceHeader)
    fprintf(se_fp_nj, ceHeader)
    fprintf(se_fp_nss, ceHeader)
    fprintf(mxe_fp, mxeHeader)
    fprintf(mxe_fp_nj, mxeHeader)
    fprintf(mxe_fp_nss, mxeHeader)
    fprintf(alt3_fp, altSSHeader)
    fprintf(alt3_fp_nj, altSSHeader)
    fprintf(alt3_fp_nss, altSSHeader)
    fprintf(alt5_fp, altSSHeader)
    fprintf(alt5_fp_nj, altSSHeader)
    fprintf(alt5_fp_nss, altSSHeader)
    fprintf(ri_fp, riHeader)
    fprintf(ri_fp_nj, riHeader)
    fprintf(ri_fp_nss, riHeader)

    while ise != juncs_se.end():
        save_se_ase(se_fp, deref(ise))

        if deref(ise).includes_novel_ss:
            save_se_ase(se_fp_nss, deref(ise))
        elif deref(ise).txtype == BAM_TX:
            save_se_ase(se_fp_nj, deref(ise))

        inc(ise)

    while imxe != juncs_mxe.end():
        save_mxe_ase(mxe_fp, deref(imxe))

        if deref(imxe).includes_novel_ss:
            save_mxe_ase(mxe_fp_nss, deref(imxe))
        elif deref(imxe).txtype == BAM_TX:
            save_mxe_ase(mxe_fp_nj, deref(imxe))

        inc(imxe)

    while ialt3 != juncs_alt3.end():
        save_alt35_ase(alt3_fp, deref(ialt3))

        if deref(ialt3).includes_novel_ss:
            save_alt35_ase(alt3_fp_nss, deref(ialt3))
        elif deref(ialt3).txtype == BAM_TX:
            save_alt35_ase(alt3_fp_nj, deref(ialt3))

        inc(ialt3)

    while ialt5 != juncs_alt5.end():
        save_alt35_ase(alt5_fp, deref(ialt5))

        if deref(ialt5).includes_novel_ss:
            save_alt35_ase(alt5_fp_nss, deref(ialt5))
        elif deref(ialt5).txtype == BAM_TX:
            save_alt35_ase(alt5_fp_nj, deref(ialt5))

        inc(ialt5)

    while iri != juncs_ri.end():
        save_ri_ase(ri_fp, deref(iri))

        if deref(iri).includes_novel_ss:
            save_ri_ase(ri_fp_nss, deref(iri))
        elif deref(iri).txtype == BAM_TX:
            save_ri_ase(ri_fp_nj, deref(iri))

        inc(iri)

    fclose(se_fp)
    fclose(se_fp_nj)
    fclose(se_fp_nss)
    fclose(mxe_fp)
    fclose(mxe_fp_nj)
    fclose(mxe_fp_nss)
    fclose(alt3_fp)
    fclose(alt3_fp_nj)
    fclose(alt3_fp_nss)
    fclose(alt5_fp)
    fclose(alt5_fp_nj)
    fclose(alt5_fp_nss)
    fclose(ri_fp)
    fclose(ri_fp_nj)
    fclose(ri_fp_nss)


@boundscheck(False)
@wraparound(False)
cdef split_nj(const string& gID,
              vector[vector[Triad]]& novel_j,
              vector[vector[Triad]]& novel_ss,
              vector[unordered_map[string,vector[Triad]]]& novel_juncs,
              cbool& novelSS):
    cdef:
        size_t i, j
        int a, b
        unordered_map[string,vector[Triad]].iterator ibam

    for i in range(novel_juncs.size()):
        ibam = novel_juncs[i].find(gID)
        novel_j[i].clear()
        novel_ss[i].clear()
        if ibam != novel_juncs[i].end():
            for j in range(deref(ibam).second.size()):
                if deref(ibam).second[j].left < 0 or\
                        deref(ibam).second[j].mid < 0 or\
                        deref(ibam).second[j].right < 0:
                    if novelSS:
                        novel_ss[i].push_back(deref(ibam).second[j])
                else:
                    novel_j[i].push_back(deref(ibam).second[j])

    return


@boundscheck(False)
@wraparound(False)
cdef detect_ase(unordered_map[string,Gene]& genes,
                unordered_map[string,SupInfo]& supple, string& od,
                vector[unordered_map[string,vector[Triad]]]& novel_juncs,
                cset[SE_info]& se, cset[MXE_info]& mxe,
                cset[ALT35_info]& alt3, cset[ALT35_info]& alt5,
                cset[RI_info]& ri, int& jld2, int& rl, cbool& novelSS,
                long& mel):
    cdef:
        size_t idx
        cmap[pair[long,long],int].iterator iexon
        unordered_map[string,Gene].iterator igs = genes.begin()
        long numse = 0, nummxe = 0, num3 = 0, num5 = 0, numri = 0
        int rl_jl = rl-jld2
        cmap[SE_key,SE_info] junction_se
        cmap[MXE_key,MXE_info] junction_mxe
        cmap[ALT35_key,ALT35_info] junction_3
        cmap[ALT35_key,ALT35_info] junction_5
        cmap[RI_key,RI_info] junction_ri
        cmap[SE_key,SE_info].iterator ise
        cmap[MXE_key,MXE_info].iterator imxe
        cmap[ALT35_key,ALT35_info].iterator ialt3
        cmap[ALT35_key,ALT35_info].iterator ialt5
        cmap[RI_key,RI_info].iterator iri
        vector[vector[Triad]] novel_j
        vector[vector[Triad]] novel_ss
        int maxiter = 10

    novel_j.resize(novel_juncs.size())
    novel_ss.resize(novel_juncs.size())

    includes_novel_ss = False
    while igs != genes.end():
        split_nj(deref(igs).first, novel_j, novel_ss, novel_juncs, novelSS)
        buildup_graph(deref(igs).first, deref(igs).second, novel_j, novelSS)

        for idx in range(deref(igs).second.idx_exon.size()):
            detect_se_mxe(deref(igs).first, deref(igs).second, supple[deref(igs).first],
                          deref(igs).second.idx_exon[idx], idx,
                          junction_se, junction_mxe, jld2, rl, rl_jl,
                          includes_novel_ss)
            detect_alt35(deref(igs).first, deref(igs).second, supple[deref(igs).first],
                         deref(igs).second.idx_exon[idx], idx,
                         junction_3, junction_5, jld2, rl, rl_jl,
                         includes_novel_ss)
            detect_ri(deref(igs).first, deref(igs).second, supple[deref(igs).first],
                      deref(igs).second.idx_exon[idx], idx,
                      junction_ri, jld2, rl, rl_jl, includes_novel_ss)

        # Release the memory because we don't need it any more.
        if novelSS:
            deref(igs).second.sg.clear()
        else:
            deref(igs).second = Gene()

        inc(igs)

    if novelSS:
        includes_novel_ss = True
        igs = genes.begin()
        while igs != genes.end():
            split_nj(deref(igs).first, novel_j, novel_ss, novel_juncs, novelSS)
            buildup_graph(deref(igs).first, deref(igs).second, novel_j, novelSS)
            expand_graph(deref(igs).first, deref(igs).second, novel_ss, novelSS, mel)
            # maxiter = 10
            # while maxiter > 0 and expand_graph(deref(igs).first, deref(igs).second,
            #                                    novel_ss, novelSS):
            #     maxiter -= 1

            for idx in range(deref(igs).second.idx_exon.size()):
                detect_se_mxe(deref(igs).first, deref(igs).second, supple[deref(igs).first],
                              deref(igs).second.idx_exon[idx], idx,
                              junction_se, junction_mxe, jld2, rl, rl_jl,
                              includes_novel_ss)
                detect_alt35(deref(igs).first, deref(igs).second, supple[deref(igs).first],
                             deref(igs).second.idx_exon[idx], idx,
                             junction_3, junction_5, jld2, rl, rl_jl,
                             includes_novel_ss)
                detect_ri(deref(igs).first, deref(igs).second, supple[deref(igs).first],
                          deref(igs).second.idx_exon[idx], idx,
                          junction_ri, jld2, rl, rl_jl, includes_novel_ss)

            deref(igs).second = Gene()

            inc(igs)

    # TODO release memory.
    novel_juncs.clear()

    # TODO
    ise = junction_se.begin()
    while ise != junction_se.end():
        se.insert(deref(ise).second)
        inc(ise)
    junction_se.clear()

    imxe = junction_mxe.begin()
    while imxe != junction_mxe.end():
        mxe.insert(deref(imxe).second)
        inc(imxe)
    junction_mxe.clear()

    ialt3 = junction_3.begin()
    while ialt3 != junction_3.end():
        alt3.insert(deref(ialt3).second)
        inc(ialt3)
    junction_3.clear()

    ialt5 = junction_5.begin()
    while ialt5 != junction_5.end():
        alt5.insert(deref(ialt5).second)
        inc(ialt5)
    junction_5.clear()

    iri = junction_ri.begin()
    while iri != junction_ri.end():
        ri.insert(deref(iri).second)
        inc(iri)
    junction_ri.clear()

    # TODO if the key involve two end of current exon, we can use local lookup table,
    # or, we should use global lookup table.
    # update: No, there is a complicated situation. We should always use a global one.
    numse = se.size()
    nummxe = mxe.size()
    num3 = alt3.size()
    num5 = alt5.size()
    numri = ri.size()

    print '\n=========='
    print 'Done processing each gene from dictionary to compile AS events'
    print 'Found %d exon skipping events' % (numse)
    print 'Found %d exon MX events' % (nummxe)
    print 'Found %d alt SS events' % (num3+num5)
    print 'There are %d alt 3 SS events and %d alt 5 SS events.' % (num3, num5)
    print 'Found %d RI events' % (numri)
    print '==========\n'

    save_ase(od, se, mxe, alt3, alt5, ri)


@boundscheck(False)
@wraparound(False)
cdef str copy_from_gtf(str src_dir, str dest_dir, str event):
    cdef:
        str from_gtf_base, src_gtf_path, path_of_copy, mapping_path
        str mapping_path_template = 'id_mapping.{}.txt'
        str from_gtf_template = 'fromGTF.{}.txt'
        int i, id_i
        str line, mapped_line, orig_id, mapped_id
        list values

    from_gtf_base = from_gtf_template.format(event)
    src_gtf_path = join(src_dir, from_gtf_base)
    path_of_copy = join(dest_dir, from_gtf_base)
    mapping_path = join(dest_dir, mapping_path_template.format(event))
    with open(src_gtf_path, 'rt') as src_f:
        with open(path_of_copy, 'wt') as dest_f:
            with open(mapping_path, 'wt') as map_f:
                for i, line in enumerate(src_f):
                    values = line.strip().split('\t')
                    if i == 0:
                        id_i = values.index('ID')
                        dest_f.write(line)
                        mapped_line = '\t'.join(['original_id', 'mapped_id'])
                        map_f.write('{}\n'.format(mapped_line))
                        continue

                    orig_id = values[id_i]
                    # the mapped_id is (i-1) so that the first id is 0
                    mapped_id = str(i - 1)
                    values[id_i] = mapped_id
                    dest_f.write('{}\n'.format('\t'.join(values)))
                    map_f.write('{}\n'.format('\t'.join([orig_id, mapped_id])))

    return path_of_copy


@boundscheck(False)
@wraparound(False)
cdef read_event_sets(str fixed_event_set_dir, str out_dir, cset[SE_info]& se,
                     cset[MXE_info]& mxe, cset[ALT35_info]& alt3,
                     cset[ALT35_info]& alt5, cset[RI_info]& ri,
                     const int jld2, const int rl):
    cdef:
        str copied_from_gtf_path
        int rl_jl = rl-jld2

    copied_from_gtf_path = copy_from_gtf(fixed_event_set_dir, out_dir, 'SE')
    read_se_event_set(copied_from_gtf_path, jld2, rl, rl_jl, se)

    copied_from_gtf_path = copy_from_gtf(fixed_event_set_dir, out_dir, 'MXE')
    read_mxe_event_set(copied_from_gtf_path, jld2, rl, rl_jl, mxe)

    copied_from_gtf_path = copy_from_gtf(fixed_event_set_dir, out_dir, 'A3SS')
    read_alt35_event_set(copied_from_gtf_path, jld2, rl, rl_jl, alt3)

    copied_from_gtf_path = copy_from_gtf(fixed_event_set_dir, out_dir, 'A5SS')
    read_alt35_event_set(copied_from_gtf_path, jld2, rl, rl_jl, alt5)

    copied_from_gtf_path = copy_from_gtf(fixed_event_set_dir, out_dir, 'RI')
    read_ri_event_set(copied_from_gtf_path, jld2, rl, rl_jl, ri)


cdef struct FromGtfSharedColIndices:
    int event_id_index
    int g_id_index
    int g_sym_index
    int chrom_index
    int strand_index


cdef struct FromGtfSharedColValues:
    int event_id
    string g_id
    string g_sym
    string chrom
    string strand


@boundscheck(False)
@wraparound(False)
cdef find_shared_col_indices(list expected_headers,
                             FromGtfSharedColIndices* shared_col_indices):
    shared_col_indices[0].event_id_index = expected_headers.index('ID')
    shared_col_indices[0].g_id_index = expected_headers.index('GeneID')
    shared_col_indices[0].g_sym_index = expected_headers.index('geneSymbol')
    shared_col_indices[0].chrom_index = expected_headers.index('chr')
    shared_col_indices[0].strand_index = expected_headers.index('strand')


@boundscheck(False)
@wraparound(False)
cdef parse_shared_col_values(list col_vals,
                             const FromGtfSharedColIndices& shared_col_indices,
                             FromGtfSharedColValues* shared_col_values):
    shared_col_values[0].event_id = int(
        col_vals[shared_col_indices.event_id_index])
    shared_col_values[0].g_id = col_vals[shared_col_indices.g_id_index]
    shared_col_values[0].g_sym = col_vals[shared_col_indices.g_sym_index]
    shared_col_values[0].chrom = col_vals[shared_col_indices.chrom_index]
    shared_col_values[0].strand = col_vals[shared_col_indices.strand_index]


@boundscheck(False)
@wraparound(False)
cdef read_se_event_set(str from_gtf_path, const int jld2, const int rl,
                       const int rl_jl, cset[SE_info]& se):
    cdef:
        list expected_headers, col_vals
        FromGtfSharedColIndices shared_col_indices
        FromGtfSharedColValues shared_col_values
        int ex_start_index, ex_end_index, up_start_index, up_end_index
        int down_start_index, down_end_index
        int ex_start, ex_end, up_start, up_end, down_start, down_end
        int exon_i, up_i, down_i
        int line_i
        str line
        cbool is_novel_junc, is_novel_ss
        Tetrad inc_skip_lens
        SupInfo sup_info
        SE_info se_info

    expected_headers = ['ID', 'GeneID', 'geneSymbol', 'chr', 'strand',
                        'exonStart_0base', 'exonEnd', 'upstreamES',
                        'upstreamEE', 'downstreamES', 'downstreamEE']
    find_shared_col_indices(expected_headers, &shared_col_indices)
    ex_start_index = expected_headers.index('exonStart_0base')
    ex_end_index = expected_headers.index('exonEnd')
    up_start_index = expected_headers.index('upstreamES')
    up_end_index = expected_headers.index('upstreamEE')
    down_start_index = expected_headers.index('downstreamES')
    down_end_index = expected_headers.index('downstreamEE')

    with open(from_gtf_path, 'rt') as f_handle:
        for line_i, line in enumerate(f_handle):
            col_vals = line.strip().split('\t')
            if line_i == 0:
                if col_vals != expected_headers:
                    sys.exit('ERROR: unable to read event set from {}.'
                             ' Expected headers to be {}, but saw {}'.format(
                                 from_gtf_path, expected_headers, col_vals))

                continue

            parse_shared_col_values(col_vals, shared_col_indices,
                                    &shared_col_values)
            ex_start = int(col_vals[ex_start_index])
            ex_end = int(col_vals[ex_end_index])
            up_start = int(col_vals[up_start_index])
            up_end = int(col_vals[up_end_index])
            down_start = int(col_vals[down_start_index])
            down_end = int(col_vals[down_end_index])

            sup_info.set_info(shared_col_values.g_sym, shared_col_values.chrom,
                              shared_col_values.strand)
            sm_inclen(ex_start, ex_end, up_start, up_end, down_start, down_end,
                      &inc_skip_lens, jld2, rl, rl_jl)
            # The events are provided as input so are not considered novel here
            is_novel_junc = False
            is_novel_ss = False
            se_info.set(shared_col_values.event_id, shared_col_values.g_id,
                        sup_info, ex_start, ex_end, up_start, up_end,
                        down_start, down_end,
                        inc_skip_lens.first, inc_skip_lens.second,
                        inc_skip_lens.third, inc_skip_lens.fourth,
                        is_novel_junc, is_novel_ss)
            se.insert(se_info)


@boundscheck(False)
@wraparound(False)
cdef read_mxe_event_set(str from_gtf_path, const int jld2, const int rl,
                        const int rl_jl, cset[MXE_info]& mxe):
    cdef:
        list expected_headers, col_vals
        FromGtfSharedColIndices shared_col_indices
        FromGtfSharedColValues shared_col_values
        int first_ex_start_index, first_ex_end_index, second_ex_start_index
        int second_ex_end_index, up_start_index, up_end_index
        int down_start_index, down_end_index
        int first_ex_start, first_ex_end, second_ex_start, second_ex_end,
        int up_start, up_end, down_start, down_end
        int first_exon_i, second_exon_i, up_i, down_i
        int line_i
        str line
        cbool is_novel_junc, is_novel_ss
        Tetrad inc_skip_lens
        SupInfo sup_info
        MXE_info mxe_info

    expected_headers = ['ID', 'GeneID', 'geneSymbol', 'chr', 'strand',
                        '1stExonStart_0base', '1stExonEnd',
                        '2ndExonStart_0base', '2ndExonEnd', 'upstreamES',
                        'upstreamEE', 'downstreamES', 'downstreamEE']
    find_shared_col_indices(expected_headers, &shared_col_indices)
    first_ex_start_index = expected_headers.index('1stExonStart_0base')
    first_ex_end_index = expected_headers.index('1stExonEnd')
    second_ex_start_index = expected_headers.index('2ndExonStart_0base')
    second_ex_end_index = expected_headers.index('2ndExonEnd')
    up_start_index = expected_headers.index('upstreamES')
    up_end_index = expected_headers.index('upstreamEE')
    down_start_index = expected_headers.index('downstreamES')
    down_end_index = expected_headers.index('downstreamEE')

    with open(from_gtf_path, 'rt') as f_handle:
        for line_i, line in enumerate(f_handle):
            col_vals = line.strip().split('\t')
            if line_i == 0:
                if col_vals != expected_headers:
                    sys.exit('ERROR: unable to read event set from {}.'
                             ' Expected headers to be {}, but saw {}'.format(
                                 from_gtf_path, expected_headers, col_vals))

                continue

            parse_shared_col_values(col_vals, shared_col_indices,
                                    &shared_col_values)
            first_ex_start = int(col_vals[first_ex_start_index])
            first_ex_end = int(col_vals[first_ex_end_index])
            second_ex_start = int(col_vals[second_ex_start_index])
            second_ex_end = int(col_vals[second_ex_end_index])
            up_start = int(col_vals[up_start_index])
            up_end = int(col_vals[up_end_index])
            down_start = int(col_vals[down_start_index])
            down_end = int(col_vals[down_end_index])

            sup_info.set_info(shared_col_values.g_sym, shared_col_values.chrom,
                              shared_col_values.strand)
            if shared_col_values.strand == '-':
                ms_inclen(second_ex_start, second_ex_end, first_ex_start,
                          first_ex_end, up_start, up_end, down_start, down_end,
                          &inc_skip_lens, jld2, rl, rl_jl)
            else:
                ms_inclen(first_ex_start, first_ex_end, second_ex_start,
                          second_ex_end, up_start, up_end, down_start, down_end,
                          &inc_skip_lens, jld2, rl, rl_jl)

            is_novel_junc = False
            is_novel_ss = False
            mxe_info.set(shared_col_values.event_id, shared_col_values.g_id,
                         sup_info, first_ex_start, first_ex_end,
                         second_ex_start, second_ex_end, up_start, up_end,
                         down_start, down_end, inc_skip_lens.first,
                         inc_skip_lens.second, inc_skip_lens.third,
                         inc_skip_lens.fourth, is_novel_junc, is_novel_ss)
            mxe.insert(mxe_info)


@boundscheck(False)
@wraparound(False)
cdef read_alt35_event_set(str from_gtf_path, const int jld2, const int rl,
                          const int rl_jl, cset[ALT35_info]& alt35):
    cdef:
        list expected_headers, col_vals
        FromGtfSharedColIndices shared_col_indices
        FromGtfSharedColValues shared_col_values
        int long_start_index, long_end_index, short_start_index
        int short_end_index, flank_start_index, flank_end_index
        int long_start, long_end, short_start, short_end, flank_start, flank_end
        int long_i, short_i, flank_i
        int line_i
        str line
        cbool is_novel_junc, is_novel_ss
        Tetrad inc_skip_lens
        SupInfo sup_info
        ALT35_info alt35_info

    expected_headers = ['ID', 'GeneID', 'geneSymbol', 'chr', 'strand',
                        'longExonStart_0base', 'longExonEnd', 'shortES',
                        'shortEE', 'flankingES', 'flankingEE']
    find_shared_col_indices(expected_headers, &shared_col_indices)
    long_start_index = expected_headers.index('longExonStart_0base')
    long_end_index = expected_headers.index('longExonEnd')
    short_start_index = expected_headers.index('shortES')
    short_end_index = expected_headers.index('shortEE')
    flank_start_index = expected_headers.index('flankingES')
    flank_end_index = expected_headers.index('flankingEE')

    with open(from_gtf_path, 'rt') as f_handle:
        for line_i, line in enumerate(f_handle):
            col_vals = line.strip().split('\t')
            if line_i == 0:
                if col_vals != expected_headers:
                    sys.exit('ERROR: unable to read event set from {}.'
                             ' Expected headers to be {}, but saw {}'.format(
                                 from_gtf_path, expected_headers, col_vals))

                continue

            parse_shared_col_values(col_vals, shared_col_indices,
                                    &shared_col_values)
            long_start = int(col_vals[long_start_index])
            long_end = int(col_vals[long_end_index])
            short_start = int(col_vals[short_start_index])
            short_end = int(col_vals[short_end_index])
            flank_start = int(col_vals[flank_start_index])
            flank_end = int(col_vals[flank_end_index])

            sup_info.set_info(shared_col_values.g_sym, shared_col_values.chrom,
                              shared_col_values.strand)
            alt_inclen(long_start, long_end, short_start, short_end,
                       flank_start, flank_end, &inc_skip_lens, jld2, rl, rl_jl)
            is_novel_junc = False
            is_novel_ss = False
            alt35_info.set(shared_col_values.event_id, shared_col_values.g_id,
                           sup_info, long_start, long_end, short_start,
                           short_end, flank_start, flank_end,
                           inc_skip_lens.first, inc_skip_lens.second,
                           inc_skip_lens.third, inc_skip_lens.fourth,
                           is_novel_junc, is_novel_ss)
            alt35.insert(alt35_info)


@boundscheck(False)
@wraparound(False)
cdef read_ri_event_set(str from_gtf_path, const int jld2, const int rl,
                       const int rl_jl, cset[RI_info]& ri):
    cdef:
        list expected_headers, col_vals
        FromGtfSharedColIndices shared_col_indices
        FromGtfSharedColValues shared_col_values
        int ri_start_index, ri_end_index, up_start_index, up_end_index
        int down_start_index, down_end_index
        int ri_start, ri_end, up_start, up_end, down_start, down_end
        int ri_i, up_i, down_i
        int line_i
        str line
        cbool is_novel_junc, is_novel_ss
        Tetrad inc_skip_lens
        SupInfo sup_info
        RI_info ri_info

    expected_headers = ['ID', 'GeneID', 'geneSymbol', 'chr', 'strand',
                        'riExonStart_0base', 'riExonEnd', 'upstreamES',
                        'upstreamEE', 'downstreamES', 'downstreamEE']
    find_shared_col_indices(expected_headers, &shared_col_indices)
    ri_start_index = expected_headers.index('riExonStart_0base')
    ri_end_index = expected_headers.index('riExonEnd')
    up_start_index = expected_headers.index('upstreamES')
    up_end_index = expected_headers.index('upstreamEE')
    down_start_index = expected_headers.index('downstreamES')
    down_end_index = expected_headers.index('downstreamEE')

    with open(from_gtf_path, 'rt') as f_handle:
        for line_i, line in enumerate(f_handle):
            col_vals = line.strip().split('\t')
            if line_i == 0:
                if col_vals != expected_headers:
                    sys.exit('ERROR: unable to read event set from {}.'
                             ' Expected headers to be {}, but saw {}'.format(
                                 from_gtf_path, expected_headers, col_vals))

                continue

            parse_shared_col_values(col_vals, shared_col_indices,
                                    &shared_col_values)
            ri_start = int(col_vals[ri_start_index])
            ri_end = int(col_vals[ri_end_index])
            up_start = int(col_vals[up_start_index])
            up_end = int(col_vals[up_end_index])
            down_start = int(col_vals[down_start_index])
            down_end = int(col_vals[down_end_index])

            sup_info.set_info(shared_col_values.g_sym, shared_col_values.chrom,
                              shared_col_values.strand)
            ri_inclen(ri_start, ri_end, up_start, up_end, down_start, down_end,
                      &inc_skip_lens, jld2, rl, rl_jl)
            is_novel_junc = False
            is_novel_ss = False
            ri_info.set(shared_col_values.event_id, shared_col_values.g_id,
                        sup_info, ri_start, ri_end, up_start, up_end,
                        down_start, down_end,
                        inc_skip_lens.first, inc_skip_lens.second,
                        inc_skip_lens.third, inc_skip_lens.fourth,
                        is_novel_junc, is_novel_ss)
            ri.insert(ri_info)


@boundscheck(False)
@wraparound(False)
cdef void statistic(unordered_map[string,Gene]& genes,
                    unordered_map[int,cset[string]]& geneGroup):
    """TODO: Docstring for statistic.
    :returns: TODO

    """
    cdef:
        long nGene = genes.size() ## number of genes in genes dict
        long nTx = 0 ## number of transcripts
        long oneTx = 0 ## number of one-tx genes
        long nExon = 0 ## number of exons
        long oneExon = 0 ## number of one-exon transcripts
        long count_gene_per_group = 0
        long oneTxOneExon = 0 ## number of one-tx genes with only one exon
        unordered_map[string,Gene].iterator igs
        unordered_map[string,Transcript].iterator itx
        unordered_map[int,cset[string]].iterator igg

    igs = genes.begin()
    while igs != genes.end():
        nTx += deref(igs).second.trans.size() 
        if deref(igs).second.trans.size() == 1:
            oneTx += 1
        itx = deref(igs).second.trans.begin()
        while itx != deref(igs).second.trans.end():
            nExon += deref(itx).second.exons.size()
            if deref(itx).second.exons.size() == 1:
                oneExon += 1
                if deref(igs).second.trans.size() == 1:
                    oneTxOneExon += 1

            inc(itx)

        inc(igs)

    igg = geneGroup.begin()
    while igg != geneGroup.end():
        count_gene_per_group += deref(igg).second.size()
        inc(igg)

    print "There are %d distinct gene ID in the gtf file" % nGene
    print "There are %d distinct transcript ID in the gtf file" % nTx
    print "There are %d one-transcript genes in the gtf file" % oneTx
    print "There are %d exons in the gtf file" % nExon
    print "There are %d one-exon transcripts in the gtf file" % oneExon
    print "There are %d one-transcript genes with only one exon in the transcript" % oneTxOneExon
    if (nGene>0): ## to avoid divided by zero exception
        print "Average number of transcripts per gene is %f" % (float(nTx)/nGene)
    if (nTx>0): ## to avoid divided by zero exception
        print "Average number of exons per transcript is %f" % (float(nExon)/nTx)
    if (nTx-oneExon)>0: ## to avoid divided by zero exception
        print "Average number of exons per transcript excluding one-exon tx is %f" % (float(nExon-oneExon)/(nTx-oneExon))
    print "Average number of gene per geneGroup is %f" % (float(count_gene_per_group)/geneGroup.size())


@boundscheck(False)
@wraparound(False)
cdef save_nj(fp, unordered_map[string,vector[Triad]]& novel_juncs):
    cdef:
        string line
        list formated
        size_t i = 0
        unordered_map[string,vector[Triad]].iterator inj
        cset[Triad] tmp_set

    line = '%d\n' % (novel_juncs.size())
    fp.write(line)

    inj = novel_juncs.begin()
    while inj != novel_juncs.end():
        tmp_set.clear()
        for i in range(deref(inj).second.size()):
            tmp_set.insert(deref(inj).second[i])

        deref(inj).second.assign(tmp_set.begin(), tmp_set.end())
        formated = ['%d,%d,%d' % (deref(inj).second[i].left,
                    deref(inj).second[i].mid, deref(inj).second[i].right,)
                    for i in range(deref(inj).second.size())]
        formated.insert(0,deref(inj).first)
        line = ';'.join(formated)
        fp.write(line)
        fp.write('\n')
        inc(inj)


@boundscheck(False)
@wraparound(False)
cdef save_exons(fp, unordered_map[string,cmap[Tetrad,int]]& exons):
    cdef:
        string line
        unordered_map[string,cmap[Tetrad,int]].iterator iexons
        cmap[Tetrad,int].iterator imap

    line = '%d\n' % (exons.size())
    fp.write(line)

    iexons = exons.begin()
    while iexons != exons.end():
        line = deref(iexons).first
        imap = deref(iexons).second.begin()
        while imap != deref(iexons).second.end():
            line = '%s;%d,%d,%d,%d,%d' % (line, deref(imap).first.first,
                    deref(imap).first.second, deref(imap).first.third,
                    deref(imap).first.fourth, deref(imap).second)
            inc(imap)

        fp.write(line)
        fp.write('\n')
        inc(iexons)


@boundscheck(False)
@wraparound(False)
cdef save_multis(fp, unordered_map[string,cmap[string,int]]& multis):
    cdef:
        string line
        unordered_map[string,cmap[string,int]].iterator imultis
        cmap[string,int].iterator imap

    line = '%d\n' % (multis.size())
    fp.write(line)

    imultis = multis.begin()
    while imultis != multis.end():
        line = deref(imultis).first
        imap = deref(imultis).second.begin()
        while imap != deref(imultis).second.end():
            line = '%s;%s,%d' % (line, deref(imap).first, deref(imap).second)
            inc(imap)

        fp.write(line)
        fp.write('\n')
        inc(imultis)


@boundscheck(False)
@wraparound(False)
cdef save_job(str bams, const string& tmp_dir, str prep_prefix,
              const int& readLength,
              vector[unordered_map[string,vector[Triad]]]& novel_juncs,
              vector[unordered_map[string,cmap[Tetrad,int]]]& exons,
              vector[unordered_map[string,cmap[string,int]]]& multis):
    cdef:
        int bam_i
        str bam
        str file_name_template
        str file_path
        vector[string] vbams = bams.split(',')

    file_name_template = join(tmp_dir, '{}_{}.rmats')
    for bam_i, bam in enumerate(vbams):
        file_path = file_name_template.format(prep_prefix, bam_i)
        with open(file_path, 'wt') as fp:
            fp.write('{}\n'.format(bam))
            fp.write('{}\n'.format(readLength))

            save_nj(fp, novel_juncs[bam_i])
            save_exons(fp, exons[bam_i])
            save_multis(fp, multis[bam_i])

    print('The splicing graph and candidate read have been saved into {}'
          .format(file_name_template.format(prep_prefix, '*')))


@boundscheck(False)
@wraparound(False)
cdef int try_get_index(list values, object value, cbool* found):
    cdef int idx = 0
    try:
        idx = values.index(value)
    except ValueError:
        found[0] = False
        return 0

    found[0] = True
    return idx


@boundscheck(False)
@wraparound(False)
cdef int _load_job(str rmatsf, list vbams,
                   vector[unordered_map[string,vector[Triad]]]& novel_juncs,
                   vector[unordered_map[string,cmap[Tetrad,int]]]& exons,
                   vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs,
                   int mode):
    cdef:
        int i = 0, j = 0, num = 0, idx = 0
        size_t vlen
        Triad triad
        Tetrad tetrad
        vector[pair[long,long]] vp
        str line, gene_id, bam
        list bams, ele, eles, aligns, coords
        cbool index_found

    vlen = len(vbams)
    if novel_juncs.size() != vlen:
        novel_juncs.resize(vlen)
    if exons.size() != vlen:
        exons.resize(vlen)
    if juncs.size() != vlen:
        juncs.resize(vlen)

    with open(rmatsf, 'r') as fp:
        bams = fp.readline().strip().split(',')
        if len(bams) != 1:
            sys.exit('Expected 1 bam per .rmats but got: {}\n'.format(bams))

        bam = bams[0]
        idx = try_get_index(vbams, bam, &index_found)
        if not index_found:
            sys.exit('Found data for unexpected bam in .rmats: {}\n'
                     .format(bam))

        novel_juncs[idx].clear()
        exons[idx].clear()
        juncs[idx].clear()

        # Skip over read length line. Already handled in split_sg_files_by_bam.
        fp.readline()

        # processing novel junctions
        num = int(fp.readline())
        for i in range(num):
            line = fp.readline().strip()
            if mode == read_mode:
                continue

            eles = line.split(';')
            gene_id = eles[0]
            for line in eles[1:]:
                ele = [int(s) for s in line.split(',')]
                triad.set(ele[0], ele[1], ele[2])
                novel_juncs[idx][gene_id].push_back(triad)

        if mode == sg_mode:
            return idx

        # processing exonic reads
        num = int(fp.readline())
        for i in range(num):
            line = fp.readline().strip()
            eles = line.split(';')
            gene_id = eles[0]
            for line in eles[1:]:
                ele = [int(s) for s in line.split(',')]
                tetrad.set(ele[0], ele[1], ele[2], ele[3])
                exons[idx][gene_id][tetrad] = ele[4]

        # processing junction reads
        num = int(fp.readline())
        for i in range(num):
            line = fp.readline().strip()
            eles = line.split(';')
            gene_id = eles[0]
            for line in eles[1:]:
                ele = [s for s in line.split(',')]
                aligns = ele[0].split('=')
                vp = vector[pair[long,long]](len(aligns))
                for j in range(len(aligns)):
                    coords = [int(s) for s in aligns[j].split(':')]
                    vp[j].first = coords[0]
                    vp[j].second = coords[1]
                juncs[idx][gene_id][vp] = int(ele[1])

        return idx


@boundscheck(False)
@wraparound(False)
cdef load_sg(str bams, list dot_rmats_paths,
             vector[unordered_map[string,vector[Triad]]]& novel_juncs):
    cdef:
        int num = 0, bam_i
        list vbams = bams.split(',')
        list prep_counts_by_bam
        vector[unordered_map[string,cmap[Tetrad,int]]] exons
        vector[unordered_map[string,cmap[vector[pair[long,long]],int]]] juncs

    num = len(vbams)
    prep_counts_by_bam = [0 for i in range(num)]

    for name in dot_rmats_paths:
        bam_i = _load_job(name, vbams, novel_juncs, exons, juncs, sg_mode)
        prep_counts_by_bam[bam_i] += 1

    prep_counts_by_bam_name = {bam_name: 0 for bam_name in vbams}
    input_counts_by_bam_name = {bam_name: 0 for bam_name in vbams}
    for i, prep_count in enumerate(prep_counts_by_bam):
        bam_name = vbams[i]
        input_counts_by_bam_name[bam_name] += 1
        prep_counts_by_bam_name[bam_name] += prep_count

    any_error = False
    for bam_name, input_count in input_counts_by_bam_name.items():
        if input_count != 1:
            sys.stderr.write('{} given {} times in input\n'.format(
                bam_name, input_count))
            any_error = True

    for bam_name, prep_count in prep_counts_by_bam_name.items():
        if prep_count == 0:
            sys.stderr.write('{} not found in .rmats files\n'.format(
                bam_name))
            any_error = True
        elif prep_count > 1:
            sys.stderr.write('{} found {} times in .rmats files\n'.format(
                bam_name, prep_count))
            any_error = True

    if any_error:
        sys.exit(1)


@boundscheck(False)
@wraparound(False)
cdef dict split_sg_files_by_bam(str bams, str tmp_dir, str out_dir,
                                const int read_length):
    cdef:
        dict result
        int gene_i
        int num_genes
        int orig_i
        int read_length_from_file
        int section_i
        int split_i
        list all_orig_dot_rmats
        list bams_from_file
        list dot_rmats_file_paths
        list in_input_by_splits
        list splits_for_this_orig
        set input_bam_paths
        str bam_from_file
        str gene_line
        str num_genes_line
        str orig_basename
        str orig_dot_rmats
        str split_dot_rmats_dir_path
        str split_path

    split_dot_rmats_dir_path = join(out_dir, 'split_dot_rmats')
    if exists(split_dot_rmats_dir_path):
        shutil.rmtree(split_dot_rmats_dir_path)

    mkdir(split_dot_rmats_dir_path)

    input_bam_paths = set(bams.split(','))
    all_orig_dot_rmats = [join(root, name)
                          for root, dirs, files in walk(tmp_dir)
                          for name in files if name.endswith('.rmats')]
    dot_rmats_file_paths = list()
    for orig_i, orig_dot_rmats in enumerate(all_orig_dot_rmats):
        with open(orig_dot_rmats, 'r') as orig_handle:
            bams_from_file = orig_handle.readline().strip().split(',')
            if bams_from_file == ['']:
                sys.stderr.write(
                    'WARNING: A .rmats file was found with no bams listed in'
                    ' it. Ignoring that file: {}\n'.format(orig_dot_rmats))
                continue

            read_length_from_file = int(orig_handle.readline().strip())
            if read_length_from_file != read_length:
                print('WARNING: The post step should use the same read length'
                      ' as the prep step.'
                      '\n         The prep step\'s read length: {}'
                      '\n         The post step\'s read length: {}'
                      '\n         Please check {}'
                      .format(read_length_from_file, read_length,
                              orig_dot_rmats))

            # Already a single bam
            if len(bams_from_file) == 1:
                # Only track files with data from input bams
                if bams_from_file[0] in input_bam_paths:
                    dot_rmats_file_paths.append(orig_dot_rmats)

                continue

            splits_for_this_orig = list()
            in_input_by_splits = list()
            for split_i, bam_from_file in enumerate(bams_from_file):
                in_input_by_splits.append(bam_from_file in input_bam_paths)
                orig_basename = basename(orig_dot_rmats)
                split_path = join(
                    split_dot_rmats_dir_path,
                    '{}_{}_{}'.format(orig_i, split_i, orig_basename))
                splits_for_this_orig.append(split_path)

                # Only track files with data from input bams
                if not in_input_by_splits[split_i]:
                    continue

                dot_rmats_file_paths.append(split_path)
                with open(split_path, 'wt') as split_handle:
                    split_handle.write('{}\n'.format(bam_from_file))
                    split_handle.write('{}\n'.format(read_length_from_file))

            # copy the novel junctions, exon reads, and junction reads sections
            for section_i in range(3):
                for split_i, split_path in enumerate(splits_for_this_orig):
                    num_genes_line = orig_handle.readline()
                    num_genes = int(num_genes_line)
                    # Only write files with data from input bams
                    if in_input_by_splits[split_i]:
                        with open(split_path, 'at') as split_handle:
                            split_handle.write(num_genes_line)

                            for gene_i in range(num_genes):
                                gene_line = orig_handle.readline()
                                split_handle.write(gene_line)
                    else:
                        for gene_i in range(num_genes):
                            # skip over without writing
                            orig_handle.readline()

    result = dict()
    result['file_paths'] = dot_rmats_file_paths
    result['dir_path'] = split_dot_rmats_dir_path
    return result


@boundscheck(False)
@wraparound(False)
cdef int load_read(str bams, str fn,
                   vector[unordered_map[string,cmap[Tetrad,int]]]& exons,
                   vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs):
    cdef:
        int bam_i
        list vbams = bams.split(',')
        vector[unordered_map[string,vector[Triad]]] novel_juncs

    bam_i = _load_job(fn, vbams, novel_juncs, exons, juncs, read_mode)
    return bam_i


def run_pipe(args):
    """TODO: Docstring for run_pipe.
    :returns: TODO

    """
    cdef:
        unordered_map[int,cset[string]] geneGroup
        unordered_map[string,Gene] genes
        unordered_map[string,SupInfo] supple
        vector[unordered_map[string,vector[Triad]]] novel_juncs
        vector[unordered_map[string,cmap[Tetrad,int]]] exons
        vector[unordered_map[string,cmap[string,int]]] multis
        cset[SE_info] se,
        cset[MXE_info] mxe,
        cset[ALT35_info] alt3,
        cset[ALT35_info] alt5,
        cset[RI_info] ri,
        int sam1len = len(args.b1.split(','))
        int jld2

    start = time.time()
    try:
        parse_gtf(args.gtf, geneGroup, genes, supple)
    except Exception:
        print('unable to parse the gtf: {}'.format(args.gtf))
        print('please check that the --gtf argument is a valid'
              ' .gtf file that is not compressed')
        raise

    print 'gtf:', time.time() - start

    start = time.time()
    statistic(genes, geneGroup)
    print 'statistic:', time.time() - start

    prep_tasks = set(['prep', 'both',])
    post_tasks = set(['post', 'both',])

    if args.task in prep_tasks:

        start = time.time()
        detect_novel(args.bams, geneGroup, genes, supple,
                     novel_juncs, exons, multis, args)
        print 'novel:', time.time() - start

        start = time.time()
        save_job(args.bams, args.tmp, args.prep_prefix, args.readLength,
                 novel_juncs, exons, multis)
        print 'save:', time.time() - start

    if args.task in post_tasks:

        multis.clear()
        start = time.time()
        split_sg_files_result = split_sg_files_by_bam(
            args.bams, args.tmp, args.od, args.readLength)
        dot_rmats_file_paths = split_sg_files_result['file_paths']
        split_dot_rmats_dir_path = split_sg_files_result['dir_path']

        load_sg(args.bams, dot_rmats_file_paths, novel_juncs)
        print 'loadsg:', time.time() - start

        start = time.time()
        jld2 = args.junctionLength/2
        if args.fixed_event_set:
            read_event_sets(args.fixed_event_set, args.od, se, mxe, alt3, alt5,
                            ri, jld2, args.readLength)
        else:
            detect_ase(genes, supple, args.od, novel_juncs,
                       se, mxe, alt3, alt5, ri,
                       jld2, args.readLength, args.novelSS,
                       args.mel)

        # release memory
        genes.clear()
        supple.clear()
        novel_juncs.clear()
        print 'ase:', time.time() - start

        start = time.time()
        count_occurrence(args.bams, dot_rmats_file_paths, args.od, se, mxe,
                         alt3, alt5, ri, sam1len, jld2,
                         args.readLength, args.nthread, args.stat,
                         args.individual_counts)
        print 'count:', time.time() - start

        shutil.rmtree(split_dot_rmats_dir_path)
