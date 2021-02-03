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

    inline int c_max(long a, long b) nogil: return a if a >= b else b
    inline int c_min(long a, long b) nogil: return a if a <= b else b

    string tmptmp = "chr1"
    string tmpgene = '\"ENSG00000124508\"'


@boundscheck(False)
@wraparound(False)
cdef void parse_gtf(str gtff, unordered_map[int,cset[string]]& geneGroup,
                    unordered_map[string,Gene]& genes,
                    unordered_map[string,SupInfo]& supple):
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


@boundscheck(False)
@wraparound(False)
cdef char check_strand(const BamAlignment& bread, const cbool& ispaired, const int& dt) nogil:
    if dt == FRFIRSTSTRAND:
        if ispaired:
            if (bread.IsFirstMate() and bread.IsReverseStrand()) or\
                    (bread.IsSecondMate() and bread.IsMateReverseStrand()):
                return plus_mark
            elif (bread.IsFirstMate() and bread.IsMateReverseStrand()) or\
                    (bread.IsSecondMate() and bread.IsReverseStrand()):
                return minus_mark
        else:
            if bread.IsReverseStrand():
                return plus_mark
            elif bread.IsMateReverseStrand():
                return minus_mark
    elif dt == FRSECONDSTRAND:
        if ispaired:
            if (bread.IsFirstMate() and bread.IsMateReverseStrand()) or\
                    (bread.IsSecondMate() and bread.IsReverseStrand()):
                return plus_mark
            elif (bread.IsFirstMate() and bread.IsReverseStrand()) or\
                    (bread.IsSecondMate() and bread.IsMateReverseStrand()):
                return minus_mark
        else:
            if bread.IsMateReverseStrand():
                return plus_mark
            elif bread.IsReverseStrand():
                return minus_mark
    elif dt == FRUNSTRANDED:
        pass

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
                    unordered_map[string,cmap[Tetrad,pair[int,int]]]& exons,
                    unordered_map[string,cmap[string,int]]& multis,
                    cbool issingle, int jld2, int readLength,
                    cbool variable_read_length, int dt, cbool& novelSS,
                    long& mil, long& mel, cbool allow_clipping,
                    vector[int]& read_outcome_counts) nogil:
    """TODO: Docstring for parse_bam.
    :returns: TODO

    """
    cdef:
        cset[pair[int,int]] tmp_set
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
            print 'Fail to open %s' % (bam)
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

            strand = check_strand(bread, ispaired, dt)
            if dt != FRUNSTRANDED and strand == cdot:
                read_outcome_counts[READ_NOT_EXPECTED_STRAND] += 1
                continue

            mec = mc + cigar_data_after_clipping[0].Length - 1
            bref_name = refid2str[bread.RefID]

            visited.clear()
            for i in range(mc/refer_len, mec/refer_len+1):
                cg = geneGroup[i].begin()
                while cg != geneGroup[i].end():
                    ## for each candidate gene

                    if supple[deref(cg)].chrom != bref_name or\
                        visited.find(deref(cg)) != visited.end():
                        inc(cg)
                        continue

                    visited.insert(deref(cg))
                    enclosed = locate_exon(mc, mec, rl_jl, tetrad, genes[deref(cg)])

                    if (tetrad.first != -1 and tetrad.second != -1) or\
                            (tetrad.third != -1 and tetrad.fourth != -1):
                        any_exon_match = True
                        if dt == FRUNSTRANDED:
                            exons[deref(cg)][tetrad].first = exons[deref(cg)][tetrad].first + 1
                        elif strand == plus_mark:
                            exons[deref(cg)][tetrad].first = exons[deref(cg)][tetrad].first + 1
                        elif strand == minus_mark:
                            exons[deref(cg)][tetrad].second = exons[deref(cg)][tetrad].second + 1

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
                    cg = geneGroup[j].begin()
                    while cg != geneGroup[j].end():
                        ## for each candidate gene
                        if supple[deref(cg)].chrom != bref_name or\
                                visited.find(deref(cg)) != visited.end():
                            inc(cg)
                            continue

                        visited.insert(deref(cg))
                        locate_multi(mc, cigar_data_after_clipping, rl_jl, ntx, multiread,
                                     genes[deref(cg)], numstr, valid, novelSS,
                                     mil, mel)

                        if multiread.length() > 0:
                            any_multijunc_match = True
                            multis[deref(cg)][multiread] = multis[deref(cg)][multiread] + 1
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
cdef void output_read_outcomes(const vector[vector[int]]& read_outcome_counts,
                               const vector[string]& vbams, str tmp_dir,
                               str prep_prefix):
    cdef:
        vector[int] aggregated_read_outcome_counts
        int total_for_bam
        int total

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
                       vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
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
        vector[vector[int]] read_outcome_counts

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

                    ise = junction_se.find(se_key)
                    if ise == junction_se.end():
                        iid = junction_se.size()
                        junction_se[se_key].set(iid, gID, supInfo,
                                                exon.first-1, exon.second,
                                                gene.idx_exon[left].first-1,
                                                se_key.first, se_key.second-1,
                                                gene.idx_exon[right].second,
                                                idx, left, right,
                                                len_pair.first, len_pair.second,
                                                len_pair.third, len_pair.fourth,
                                                GTF_TX, includes_novel_ss)
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
                        if increased_inc_skp_len and not converts_to_novel_ss:
                            deref(ise).second.set(-1, gID, supInfo,
                                                  exon.first-1, exon.second,
                                                  gene.idx_exon[left].first-1,
                                                  se_key.first, se_key.second-1,
                                                  gene.idx_exon[right].second,
                                                  idx, left, right,
                                                  len_pair.first, len_pair.second,
                                                  len_pair.third, len_pair.fourth,
                                                  deref(ise).second.txtype,
                                                  deref(ise).second.includes_novel_ss)

                    # Do not update the novelJunction status when
                    # handling novelSS
                    if not includes_novel_ss:
                        ise = junction_se.find(se_key)
                        if (deref(ise).second.txtype == GTF_TX
                            and (all_novel_left
                                 or all_novel_right
                                 or all_novel_skip)):
                            deref(ise).second.txtype = BAM_TX

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
                                              idx, mid, left, right,
                                              len_pair.first, len_pair.second,
                                              len_pair.third, len_pair.fourth,
                                              GTF_TX, includes_novel_ss)
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
                    if increased_inc_skp_len and not converts_to_novel_ss:
                        deref(imxe).second.set(-1, gID, supInfo,
                                               exon.first-1, exon.second,
                                               gene.idx_exon[mid].first-1, gene.idx_exon[mid].second,
                                               gene.idx_exon[left].first-1, mxe_key.first,
                                               mxe_key.second-1, gene.idx_exon[right].second,
                                               idx, mid, left, right,
                                               len_pair.first, len_pair.second,
                                               len_pair.third, len_pair.fourth,
                                               deref(imxe).second.txtype,
                                               deref(imxe).second.includes_novel_ss)
                imxe = junction_mxe.find(mxe_key)
                # all_novel_left and all_novel_right do not need to be
                # checked here since found_idx_bam_transcript would be
                # True if either all_novel_left or all_novel_right was.
                #
                # Do not update the novelJunction status when handling novelSS
                if not includes_novel_ss:
                    if (deref(imxe).second.txtype == GTF_TX
                        and (found_mid_bam_transcript
                             or found_idx_bam_transcript)):
                        deref(imxe).second.txtype = BAM_TX


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

    ialt35 = junction_35.find(alt35_key)
    if ialt35 == junction_35.end():
        iid = junction_35.size()
        junction_35[alt35_key].set(iid, gID, supInfo,
                                   gene.idx_exon[i].first-1, alt35_key.second,
                                   gene.idx_exon[i].first-1, alt35_key.first,
                                   exon.first-1, exon.second,
                                   j, i, idx,
                                   len_pair.first, len_pair.second,
                                   len_pair.third, len_pair.fourth, GTF_TX,
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
        if increased_inc_skp_len and not converts_to_novel_ss:
            deref(ialt35).second.set(-1, gID, supInfo,
                                     gene.idx_exon[i].first-1, alt35_key.second,
                                     gene.idx_exon[i].first-1, alt35_key.first,
                                     exon.first-1, exon.second,
                                     j, i, idx,
                                     len_pair.first, len_pair.second,
                                     len_pair.third, len_pair.fourth,
                                     deref(ialt35).second.txtype,
                                     deref(ialt35).second.includes_novel_ss)

    # Do not update the novelJunction status when handling novelSS
    if not includes_novel_ss:
        ialt35 = junction_35.find(alt35_key)
        if (deref(ialt35).second.txtype == GTF_TX
            and (all_novel_i or all_novel_j)):
            deref(ialt35).second.txtype = BAM_TX


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

    ialt35 = junction_35.find(alt35_key)
    if ialt35 == junction_35.end():
        iid = junction_35.size()
        junction_35[alt35_key].set(iid, gID, supInfo,
                                   alt35_key.second, gene.idx_exon[i].second,
                                   alt35_key.third, gene.idx_exon[i].second,
                                   exon.first-1, exon.second,
                                   i, j, idx,
                                   len_pair.first, len_pair.second,
                                   len_pair.third, len_pair.fourth, GTF_TX,
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
        if increased_inc_skp_len and not converts_to_novel_ss:
            deref(ialt35).second.set(-1, gID, supInfo,
                                     alt35_key.second, gene.idx_exon[i].second,
                                     alt35_key.third, gene.idx_exon[i].second,
                                     exon.first-1, exon.second,
                                     i, j, idx,
                                     len_pair.first, len_pair.second,
                                     len_pair.third, len_pair.fourth,
                                     deref(ialt35).second.txtype,
                                     deref(ialt35).second.includes_novel_ss)

    # Do not update the novelJunction status when handling novelSS
    if not includes_novel_ss:
        ialt35 = junction_35.find(alt35_key)
        if (deref(ialt35).second.txtype == GTF_TX
            and (all_novel_i or all_novel_j)):
            deref(ialt35).second.txtype = BAM_TX


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
        cbool found_any, all_novel

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
                                    gene.exon_idx[tmp_pair], i, idx,
                                    len_pair.first, len_pair.second,
                                    len_pair.third, len_pair.fourth, GTF_TX,
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
            if increased_inc_skp_len and not converts_to_novel_ss:
                deref(iri).second.set(-1, gID, supInfo,
                                      gene.idx_exon[i].first-1, exon.second,
                                      gene.idx_exon[i].first-1, ri_key.first,
                                      ri_key.second, exon.second,
                                      gene.exon_idx[tmp_pair], i, idx,
                                      len_pair.first, len_pair.second,
                                      len_pair.third, len_pair.fourth,
                                      deref(iri).second.txtype,
                                      deref(iri).second.includes_novel_ss)

        # Do not update the novelJunction status when handling novelSS
        if not includes_novel_ss:
            iri = junction_ri.find(ri_key)
            if deref(iri).second.txtype == GTF_TX and all_novel:
                deref(iri).second.txtype = BAM_TX


@boundscheck(False)
@wraparound(False)
cdef void count_se(cset[SE_info]& junction_se,
                   vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
                   vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs,
                   vector[Read_count_table]& jc_se, vector[Read_count_table]& jcec_se,
                   vector[size_t]& residx, int& dt) nogil:
    cdef:
        int idx, count = 0
        size_t i, j, 
        Tetrad tetrad
        cmap[Tetrad,pair[int,int]].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,pair[int,int]]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[SE_info].iterator ise = junction_se.begin()

    while ise != junction_se.end():
        idx = deref(ise).iid

        for i in residx:
            iunmap = juncs[i].find(deref(ise).gID)
            if iunmap != juncs[i].end():
                imap = deref(iunmap).second.begin()
                while imap != deref(iunmap).second.end():
                    for j in range(1, deref(imap).first.size()-1):
                        if (j == 1 and\
                                deref(imap).first[j].first == deref(ise).te and\
                                deref(imap).first[j].second == deref(ise).ds and\
                                deref(imap).first[0].second >= deref(ise).ts) or\
                           (j == deref(imap).first.size() - 2 and\
                                deref(imap).first[j].first == deref(ise).ue and\
                                deref(imap).first[j].second == deref(ise).ts and\
                                deref(imap).first[deref(imap).first.size()-1].first <= deref(ise).te) or\
                           (j < deref(imap).first.size()-2 and\
                                deref(imap).first[j].first == deref(ise).ue and\
                                deref(imap).first[j].second == deref(ise).ts and\
                                deref(imap).first[j+1].first == deref(ise).te and\
                                deref(imap).first[j+1].second == deref(ise).ds):
                            jc_se[idx].incv[i] += deref(imap).second
                            jcec_se[idx].incv[i] += deref(imap).second
                            break
                        elif deref(imap).first[j].first == deref(ise).ue and\
                                deref(imap).first[j].second == deref(ise).ds:
                            jc_se[idx].skpv[i] += deref(imap).second
                            jcec_se[idx].skpv[i] += deref(imap).second
                            break

                    inc(imap)

            iunmap2 = exons[i].find(deref(ise).gID)
            if iunmap2 != exons[i].end():
                imap2 = deref(iunmap2).second.begin()
                while imap2 != deref(iunmap2).second.end():
                    if dt == FRUNSTRANDED:
                        count = deref(imap2).second.first
                    elif deref(ise).supInfo.strand == plus_mark:
                        count = deref(imap2).second.first
                    elif deref(ise).supInfo.strand == minus_mark:
                        count = deref(imap2).second.second
                    if deref(imap2).first.first > deref(ise).ts and\
                            deref(imap2).first.fourth != -1 and\
                            deref(imap2).first.fourth <= deref(ise).te:
                        jcec_se[idx].incv[i] += count
                    inc(imap2)

        inc(ise)


@boundscheck(False)
@wraparound(False)
cdef void count_mxe(cset[MXE_info]& junction_mxe,
                    vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
                    vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs,
                    vector[Read_count_table]& jc_mxe, vector[Read_count_table]& jcec_mxe,
                    vector[size_t]& residx, int& dt) nogil:
    cdef:
        int idx, count = 0
        size_t i, j
        Tetrad tetrad
        cmap[Tetrad,pair[int,int]].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,pair[int,int]]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[MXE_info].iterator imxe = junction_mxe.begin()

    while imxe != junction_mxe.end():
        idx = deref(imxe).iid

        for i in residx:
            iunmap = juncs[i].find(deref(imxe).gID)
            if iunmap != juncs[i].end():
                imap = deref(iunmap).second.begin()
                while imap != deref(iunmap).second.end():
                    for j in range(1, deref(imap).first.size()-1):
                        if (j == 1 and\
                                deref(imap).first[j].first == deref(imxe).te and\
                                deref(imap).first[j].second == deref(imxe).ds and\
                                deref(imap).first[0].second >= deref(imxe).ts) or\
                           (j == deref(imap).first.size()-2 and\
                                deref(imap).first[j].first == deref(imxe).ue and\
                                deref(imap).first[j].second == deref(imxe).ts and\
                                deref(imap).first[deref(imap).first.size()-1].first <= deref(imxe).te) or\
                           (j < deref(imap).first.size()-2 and\
                                deref(imap).first[j].first == deref(imxe).ue and\
                                deref(imap).first[j].second == deref(imxe).ts and\
                                deref(imap).first[j+1].first == deref(imxe).te and\
                                deref(imap).first[j+1].second == deref(imxe).ds):
                            jc_mxe[idx].incv[i] += deref(imap).second
                            jcec_mxe[idx].incv[i] += deref(imap).second
                            break
                        elif (j == 1 and\
                                deref(imap).first[j].first == deref(imxe).se and\
                                deref(imap).first[j].second == deref(imxe).ds and\
                                deref(imap).first[0].second >= deref(imxe).ss) or\
                           (j == deref(imap).first.size()-2 and\
                                deref(imap).first[j].first == deref(imxe).ue and\
                                deref(imap).first[j].second == deref(imxe).ss and\
                                deref(imap).first[deref(imap).first.size()-1].first <= deref(imxe).se) or\
                           (j < deref(imap).first.size()-2 and\
                                deref(imap).first[j].first == deref(imxe).ue and\
                                deref(imap).first[j].second == deref(imxe).ss and\
                                deref(imap).first[j+1].first == deref(imxe).se and\
                                deref(imap).first[j+1].second == deref(imxe).ds):
                            jc_mxe[idx].skpv[i] += deref(imap).second
                            jcec_mxe[idx].skpv[i] += deref(imap).second
                            break

                    inc(imap)

            iunmap2 = exons[i].find(deref(imxe).gID)
            if iunmap2 != exons[i].end():
                imap2 = deref(iunmap2).second.begin()
                while imap2 != deref(iunmap2).second.end():
                    if dt == FRUNSTRANDED:
                        count = deref(imap2).second.first
                    elif deref(imxe).supInfo.strand == plus_mark:
                        count = deref(imap2).second.first
                    elif deref(imxe).supInfo.strand == minus_mark:
                        count = deref(imap2).second.second
                    if deref(imap2).first.first > deref(imxe).ts and\
                            deref(imap2).first.fourth != -1 and\
                            deref(imap2).first.fourth <= deref(imxe).te:
                        jcec_mxe[idx].incv[i] += count
                    elif deref(imap2).first.first > deref(imxe).ss and\
                            deref(imap2).first.fourth != -1 and\
                            deref(imap2).first.fourth <= deref(imxe).se:
                        jcec_mxe[idx].skpv[i] += count

                    inc(imap2)

        inc(imxe)


@boundscheck(False)
@wraparound(False)
cdef void count_alt3(cset[ALT35_info]& junction_3,
                     vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
                     vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs,
                     vector[Read_count_table]& jc_alt3, vector[Read_count_table]& jcec_alt3,
                     int& jld2, int& rl, vector[size_t]& residx, int& dt) nogil:
    cdef:
        int idx, count = 0
        size_t i, j
        int rl_jl = rl - jld2
        Tetrad tetrad
        cmap[Tetrad,pair[int,int]].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,pair[int,int]]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[ALT35_info].iterator ialt3 = junction_3.begin()

    while ialt3 != junction_3.end():
        idx = deref(ialt3).iid

        if deref(ialt3).fs > deref(ialt3).le:
            for i in residx:
                iunmap = juncs[i].find(deref(ialt3).gID)
                if iunmap != juncs[i].end():
                    imap = deref(iunmap).second.begin()
                    while imap != deref(iunmap).second.end():
                        for j in range(1, deref(imap).first.size()-1):
                            if deref(imap).first[j].first == deref(ialt3).le and\
                                    deref(imap).first[j].second == deref(ialt3).fs:
                                jc_alt3[idx].incv[i] += deref(imap).second
                                jcec_alt3[idx].incv[i] += deref(imap).second
                                break
                            elif deref(imap).first[j].first == deref(ialt3).se and\
                                    deref(imap).first[j].second == deref(ialt3).fs:
                                jc_alt3[idx].skpv[i] += deref(imap).second
                                jcec_alt3[idx].skpv[i] += deref(imap).second
                                break
                            elif j == deref(imap).first.size()-2 and\
                                    deref(imap).first[j].second <= deref(ialt3).se-rl_jl and\
                                    deref(imap).first[j+1].first >= deref(ialt3).se+rl_jl and\
                                    deref(imap).first[j+1].first <= deref(ialt3).le:
                                jc_alt3[idx].incv[i] += deref(imap).second
                                jcec_alt3[idx].incv[i] += deref(imap).second
                                break

                        inc(imap)

                iunmap2 = exons[i].find(deref(ialt3).gID)
                if iunmap2 != exons[i].end():
                    imap2 = deref(iunmap2).second.begin()
                    while imap2 != deref(iunmap2).second.end():
                        if dt == FRUNSTRANDED:
                            count = deref(imap2).second.first
                        elif deref(ialt3).supInfo.strand == plus_mark:
                            count = deref(imap2).second.first
                        elif deref(ialt3).supInfo.strand == minus_mark:
                            count = deref(imap2).second.second
                        if deref(imap2).first.second <= deref(ialt3).se-rl_jl+1 and\
                                deref(imap2).first.second != -1 and\
                                deref(imap2).first.fourth != -1 and\
                                deref(imap2).first.fourth <= deref(ialt3).le and\
                                deref(imap2).first.third >= deref(ialt3).se+rl_jl:
                            jc_alt3[idx].incv[i] += count
                            jcec_alt3[idx].incv[i] += count
                        if deref(imap2).first.first > deref(ialt3).se and\
                                deref(imap2).first.fourth != -1 and\
                                deref(imap2).first.fourth <= deref(ialt3).le:
                            jcec_alt3[idx].incv[i] += count

                        inc(imap2)

        else:
            for i in residx:
                iunmap = juncs[i].find(deref(ialt3).gID)
                if iunmap != juncs[i].end():
                    imap = deref(iunmap).second.begin()
                    while imap != deref(iunmap).second.end():
                        for j in range(1, deref(imap).first.size()-1):
                            if deref(imap).first[j].first == deref(ialt3).fe and\
                                    deref(imap).first[j].second == deref(ialt3).ls:
                                jc_alt3[idx].incv[i] += deref(imap).second
                                jcec_alt3[idx].incv[i] += deref(imap).second
                                break
                            elif deref(imap).first[j].first == deref(ialt3).fe and\
                                    deref(imap).first[j].second == deref(ialt3).ss:
                                jc_alt3[idx].skpv[i] += deref(imap).second
                                jcec_alt3[idx].skpv[i] += deref(imap).second
                                break
                            elif j == 1 and\
                                    deref(imap).first[0].second >= deref(ialt3).ls and\
                                    deref(imap).first[0].second <= deref(ialt3).ss-rl_jl and\
                                    deref(imap).first[1].first >= deref(ialt3).ss+rl_jl:
                                jc_alt3[idx].incv[i] += deref(imap).second
                                jcec_alt3[idx].incv[i] += deref(imap).second
                                break

                        inc(imap)

                iunmap2 = exons[i].find(deref(ialt3).gID)
                if iunmap2 != exons[i].end():
                    imap2 = deref(iunmap2).second.begin()
                    while imap2 != deref(iunmap2).second.end():
                        if dt == FRUNSTRANDED:
                            count = deref(imap2).second.first
                        elif deref(ialt3).supInfo.strand == plus_mark:
                            count = deref(imap2).second.first
                        elif deref(ialt3).supInfo.strand == minus_mark:
                            count = deref(imap2).second.second
                        if deref(imap2).first.first > deref(ialt3).ls and\
                                deref(imap2).first.second != -1 and\
                                deref(imap2).first.second <= deref(ialt3).ss-rl_jl+1 and\
                                deref(imap2).first.third >= deref(ialt3).ss+rl_jl:
                            jc_alt3[idx].incv[i] += count
                            jcec_alt3[idx].incv[i] += count
                        if deref(imap2).first.first > deref(ialt3).ls and\
                                deref(imap2).first.fourth != -1 and\
                                deref(imap2).first.fourth <= deref(ialt3).ss:
                            jcec_alt3[idx].incv[i] += count

                        inc(imap2)

        inc(ialt3)


@boundscheck(False)
@wraparound(False)
cdef void count_alt5(cset[ALT35_info]& junction_5,
                     vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
                     vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs,
                     vector[Read_count_table]& jc_alt5, vector[Read_count_table]& jcec_alt5,
                     int& jld2, int& rl, vector[size_t]& residx, int& dt) nogil:
    cdef:
        int idx, count = 0
        size_t i, j
        long rl_jl = rl - jld2
        Tetrad tetrad
        cmap[Tetrad,pair[int,int]].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,pair[int,int]]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[ALT35_info].iterator ialt5 = junction_5.begin()

    while ialt5 != junction_5.end():
        idx = deref(ialt5).iid

        if deref(ialt5).fs > deref(ialt5).le:
            for i in residx:
                iunmap = juncs[i].find(deref(ialt5).gID)
                if iunmap != juncs[i].end():
                    imap = deref(iunmap).second.begin()
                    while imap != deref(iunmap).second.end():
                        for j in range(1, deref(imap).first.size()-1):
                            if deref(imap).first[j].first == deref(ialt5).le and\
                                    deref(imap).first[j].second == deref(ialt5).fs:
                                jc_alt5[idx].incv[i] += deref(imap).second
                                jcec_alt5[idx].incv[i] += deref(imap).second
                                break
                            elif deref(imap).first[j].first == deref(ialt5).se and\
                                    deref(imap).first[j].second == deref(ialt5).fs:
                                jc_alt5[idx].skpv[i] += deref(imap).second
                                jcec_alt5[idx].skpv[i] += deref(imap).second
                                break
                            elif j == deref(imap).first.size()-2 and\
                                    deref(imap).first[j].second <= deref(ialt5).se-rl_jl and\
                                    deref(imap).first[j+1].first >= deref(ialt5).se+rl_jl and\
                                    deref(imap).first[j+1].first <= deref(ialt5).le:
                                jc_alt5[idx].incv[i] += deref(imap).second
                                jcec_alt5[idx].incv[i] += deref(imap).second
                                break

                        inc(imap)

                iunmap2 = exons[i].find(deref(ialt5).gID)
                if iunmap2 != exons[i].end():
                    imap2 = deref(iunmap2).second.begin()
                    while imap2 != deref(iunmap2).second.end():
                        if dt == FRUNSTRANDED:
                            count = deref(imap2).second.first
                        elif deref(ialt5).supInfo.strand == plus_mark:
                            count = deref(imap2).second.first
                        elif deref(ialt5).supInfo.strand == minus_mark:
                            count = deref(imap2).second.second
                        if deref(imap2).first.second <= deref(ialt5).se-rl_jl+1 and\
                                deref(imap2).first.second != -1 and\
                                deref(imap2).first.fourth != -1 and\
                                deref(imap2).first.fourth <= deref(ialt5).le and\
                                deref(imap2).first.third >= deref(ialt5).se+rl_jl:
                            jc_alt5[idx].incv[i] += count
                            jcec_alt5[idx].incv[i] += count
                        if deref(imap2).first.first > deref(ialt5).se and\
                                deref(imap2).first.fourth != -1 and\
                                deref(imap2).first.fourth <= deref(ialt5).le:
                            jcec_alt5[idx].incv[i] += count

                        inc(imap2)

        else:
            for i in residx:
                iunmap = juncs[i].find(deref(ialt5).gID)
                if iunmap != juncs[i].end():
                    imap = deref(iunmap).second.begin()
                    while imap != deref(iunmap).second.end():
                        for j in range(1, deref(imap).first.size()-1):
                            if deref(imap).first[j].first == deref(ialt5).fe and\
                                    deref(imap).first[j].second == deref(ialt5).ls:
                                jc_alt5[idx].incv[i] += deref(imap).second
                                jcec_alt5[idx].incv[i] += deref(imap).second
                                break
                            elif deref(imap).first[j].first == deref(ialt5).fe and\
                                    deref(imap).first[j].second == deref(ialt5).ss:
                                jc_alt5[idx].skpv[i] += deref(imap).second
                                jcec_alt5[idx].skpv[i] += deref(imap).second
                                break
                            elif j == 1 and\
                                    deref(imap).first[0].second >= deref(ialt5).ls and\
                                    deref(imap).first[0].second <= deref(ialt5).ss-rl_jl and\
                                    deref(imap).first[1].first >= deref(ialt5).ss+rl_jl:
                                jc_alt5[idx].incv[i] += deref(imap).second
                                jcec_alt5[idx].incv[i] += deref(imap).second
                                break

                        inc(imap)

                iunmap2 = exons[i].find(deref(ialt5).gID)
                if iunmap2 != exons[i].end():
                    imap2 = deref(iunmap2).second.begin()
                    while imap2 != deref(iunmap2).second.end():
                        if dt == FRUNSTRANDED:
                            count = deref(imap2).second.first
                        elif deref(ialt5).supInfo.strand == plus_mark:
                            count = deref(imap2).second.first
                        elif deref(ialt5).supInfo.strand == minus_mark:
                            count = deref(imap2).second.second
                        if deref(imap2).first.first > deref(ialt5).ls and\
                                deref(imap2).first.second != -1 and\
                                deref(imap2).first.second <= deref(ialt5).ss-rl_jl+1 and\
                                deref(imap2).first.third >= deref(ialt5).ss+rl_jl:
                            jc_alt5[idx].incv[i] += count
                            jcec_alt5[idx].incv[i] += count
                        if deref(imap2).first.first > deref(ialt5).ls and\
                                deref(imap2).first.fourth != -1 and\
                                deref(imap2).first.fourth <= deref(ialt5).ss:
                            jcec_alt5[idx].incv[i] += count

                        inc(imap2)

        inc(ialt5)


@boundscheck(False)
@wraparound(False)
cdef void count_ri(cset[RI_info]& junction_ri,
                   vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
                   vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs,
                   vector[Read_count_table]& jc_ri, vector[Read_count_table]& jcec_ri,
                   int& jld2, int& rl, vector[size_t]& residx, int& dt) nogil:
    cdef:
        int idx, count = 0
        size_t i, j
        long rl_jl = rl - jld2
        Tetrad tetrad
        cmap[Tetrad,pair[int,int]].iterator imap2
        cmap[vector[pair[long,long]],int].iterator imap
        unordered_map[string,cmap[Tetrad,pair[int,int]]].iterator iunmap2
        unordered_map[string,cmap[vector[pair[long,long]],int]].iterator iunmap
        cset[RI_info].iterator iri = junction_ri.begin()

    while iri != junction_ri.end():
        idx = deref(iri).iid

        for i in residx:
            iunmap = juncs[i].find(deref(iri).gID)
            if iunmap != juncs[i].end():
                imap = deref(iunmap).second.begin()
                while imap != deref(iunmap).second.end():
                    for j in range(1, deref(imap).first.size()-1):
                        if deref(imap).first[j].first == deref(iri).ue and\
                                deref(imap).first[j].second == deref(iri).ds:
                            jc_ri[idx].skpv[i] += deref(imap).second
                            jcec_ri[idx].skpv[i] += deref(imap).second
                            break
                        elif j == 1 and\
                                deref(imap).first[0].second <= deref(iri).ds-rl_jl and\
                                deref(imap).first[1].first >= deref(iri).ds+rl_jl:
                            jc_ri[idx].incv[i] += deref(imap).second
                            jcec_ri[idx].incv[i] += deref(imap).second
                            break
                        elif j == deref(imap).first.size()-2 and\
                                deref(imap).first[j].second <= deref(iri).ue-rl_jl and\
                                deref(imap).first[j+1].first >= deref(iri).ue+rl_jl:
                            jc_ri[idx].incv[i] += deref(imap).second
                            jcec_ri[idx].incv[i] += deref(imap).second
                            break

                    inc(imap)

            iunmap2 = exons[i].find(deref(iri).gID)
            if iunmap2 != exons[i].end():
                imap2 = deref(iunmap2).second.begin()

                while imap2 != deref(iunmap2).second.end():
                    if dt == FRUNSTRANDED:
                        count = deref(imap2).second.first
                    elif deref(iri).supInfo.strand == plus_mark:
                        count = deref(imap2).second.first
                    elif deref(iri).supInfo.strand == minus_mark:
                        count = deref(imap2).second.second
                    if (deref(imap2).first.second <= deref(iri).ue-rl_jl+1 and\
                            deref(imap2).first.second != -1 and\
                            deref(imap2).first.third >= deref(iri).ue+rl_jl)\
                        or\
                       (deref(imap2).first.second != -1 and\
                            deref(imap2).first.second <= deref(iri).ds-rl_jl+1 and\
                            deref(imap2).first.third >= deref(iri).ds+rl_jl):
                        jc_ri[idx].incv[i] += count
                        jcec_ri[idx].incv[i] += count
                    if deref(imap2).first.first > deref(iri).ue and\
                            deref(imap2).first.fourth != -1 and\
                            deref(imap2).first.fourth <= deref(iri).ds:
                        jcec_ri[idx].incv[i] += count

                    inc(imap2)

        inc(iri)


@boundscheck(False)
@wraparound(False)
cdef count_occurrence(str bams, list dot_rmats_paths, str od,
                      cset[SE_info]& se, cset[MXE_info]& mxe,
                      cset[ALT35_info]& alt3, cset[ALT35_info]& alt5,
                      cset[RI_info]& ri, int sam1len, int& jld2, int& rl,
                      int& nthread, int& dt, cbool stat):
    cdef:
        size_t idx = 0
        list vbams = bams.split(',')
        int num = len(vbams), vlen, fidx
        cset[SE_info].iterator ise = se.begin()
        cset[MXE_info].iterator imxe = mxe.begin()
        cset[ALT35_info].iterator ialt3 = alt3.begin()
        cset[ALT35_info].iterator ialt5 = alt5.begin()
        cset[RI_info].iterator iri = ri.begin()
        vector[Read_count_table] jc_se = vector[Read_count_table](se.size())
        vector[Read_count_table] jc_mxe = vector[Read_count_table](mxe.size())
        vector[Read_count_table] jc_alt3 = vector[Read_count_table](alt3.size())
        vector[Read_count_table] jc_alt5 = vector[Read_count_table](alt5.size())
        vector[Read_count_table] jc_ri = vector[Read_count_table](ri.size())
        vector[Read_count_table] jcec_se = vector[Read_count_table](se.size())
        vector[Read_count_table] jcec_mxe = vector[Read_count_table](mxe.size())
        vector[Read_count_table] jcec_alt3 = vector[Read_count_table](alt3.size())
        vector[Read_count_table] jcec_alt5 = vector[Read_count_table](alt5.size())
        vector[Read_count_table] jcec_ri = vector[Read_count_table](ri.size())
        vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]] exons
        vector[unordered_map[string,cmap[vector[pair[long,long]],int]]] juncs
        vector[vector[size_t]] resindice

    while ise != se.end():
        idx = deref(ise).iid
        jc_se[idx].incv = vector[int](num, 0)
        jc_se[idx].skpv = vector[int](num, 0)
        jc_se[idx].inc_len = deref(ise).inc_len
        jc_se[idx].skp_len = deref(ise).skp_len
        jc_se[idx].strand = deref(ise).supInfo.strand
        jcec_se[idx].incv = vector[int](num, 0)
        jcec_se[idx].skpv = vector[int](num, 0)
        jcec_se[idx].inc_len = deref(ise).inc_len_jcec
        jcec_se[idx].skp_len = deref(ise).skp_len_jcec
        jcec_se[idx].strand = deref(ise).supInfo.strand
        inc(ise)

    while imxe != mxe.end():
        idx = deref(imxe).iid
        jc_mxe[idx].incv = vector[int](num, 0)
        jc_mxe[idx].skpv = vector[int](num, 0)
        jc_mxe[idx].inc_len = deref(imxe).inc_len
        jc_mxe[idx].skp_len = deref(imxe).skp_len
        jc_mxe[idx].strand = deref(imxe).supInfo.strand
        jcec_mxe[idx].incv = vector[int](num, 0)
        jcec_mxe[idx].skpv = vector[int](num, 0)
        jcec_mxe[idx].inc_len = deref(imxe).inc_len_jcec
        jcec_mxe[idx].skp_len = deref(imxe).skp_len_jcec
        jcec_mxe[idx].strand = deref(imxe).supInfo.strand
        inc(imxe)

    while ialt3 != alt3.end():
        idx = deref(ialt3).iid
        jc_alt3[idx].incv = vector[int](num, 0)
        jc_alt3[idx].skpv = vector[int](num, 0)
        jc_alt3[idx].inc_len = deref(ialt3).inc_len
        jc_alt3[idx].skp_len = deref(ialt3).skp_len
        jc_alt3[idx].strand = deref(ialt3).supInfo.strand
        jcec_alt3[idx].incv = vector[int](num, 0)
        jcec_alt3[idx].skpv = vector[int](num, 0)
        jcec_alt3[idx].inc_len = deref(ialt3).inc_len_jcec
        jcec_alt3[idx].skp_len = deref(ialt3).skp_len_jcec
        jcec_alt3[idx].strand = deref(ialt3).supInfo.strand
        inc(ialt3)

    while ialt5 != alt5.end():
        idx = deref(ialt5).iid
        jc_alt5[idx].incv = vector[int](num, 0)
        jc_alt5[idx].skpv = vector[int](num, 0)
        jc_alt5[idx].inc_len = deref(ialt5).inc_len
        jc_alt5[idx].skp_len = deref(ialt5).skp_len
        jc_alt5[idx].strand = deref(ialt5).supInfo.strand
        jcec_alt5[idx].incv = vector[int](num, 0)
        jcec_alt5[idx].skpv = vector[int](num, 0)
        jcec_alt5[idx].inc_len = deref(ialt5).inc_len_jcec
        jcec_alt5[idx].skp_len = deref(ialt5).skp_len_jcec
        jcec_alt5[idx].strand = deref(ialt5).supInfo.strand
        inc(ialt5)

    while iri != ri.end():
        idx = deref(iri).iid
        jc_ri[idx].incv = vector[int](num, 0)
        jc_ri[idx].skpv = vector[int](num, 0)
        jc_ri[idx].inc_len = deref(iri).inc_len
        jc_ri[idx].skp_len = deref(iri).skp_len
        jc_ri[idx].strand = deref(iri).supInfo.strand
        jcec_ri[idx].incv = vector[int](num, 0)
        jcec_ri[idx].skpv = vector[int](num, 0)
        jcec_ri[idx].inc_len = deref(iri).inc_len_jcec
        jcec_ri[idx].skp_len = deref(iri).skp_len_jcec
        jcec_ri[idx].strand = deref(iri).supInfo.strand
        inc(iri)

    vlen = len(dot_rmats_paths)
    resindice.resize(vlen)
    for fidx in prange(vlen, schedule='static', num_threads=nthread, nogil=True):
        with gil:
            resindice[fidx] = load_read(bams, dot_rmats_paths[fidx], exons, juncs)

        count_se(se, exons, juncs, jc_se, jcec_se, resindice[fidx], dt)
        count_mxe(mxe, exons, juncs, jc_mxe, jcec_mxe, resindice[fidx], dt)
        count_alt3(alt3, exons, juncs, jc_alt3, jcec_alt3, jld2, rl, resindice[fidx], dt)
        count_alt5(alt5, exons, juncs, jc_alt5, jcec_alt5, jld2, rl, resindice[fidx], dt)
        count_ri(ri, exons, juncs, jc_ri, jcec_ri, jld2, rl, resindice[fidx], dt)

        for idx in resindice[fidx]:
            exons[idx].clear()
            juncs[idx].clear()

    save_ct(od, jc_se, jc_mxe, jc_alt3, jc_alt5, jc_ri,
            jcec_se, jcec_mxe, jcec_alt3, jcec_alt5, jcec_ri, sam1len, stat)


@boundscheck(False)
@wraparound(False)
cdef save_ct(str od, vector[Read_count_table]& jc_se,
             vector[Read_count_table]& jc_mxe, vector[Read_count_table]& jc_alt3,
             vector[Read_count_table]& jc_alt5, vector[Read_count_table]& jc_ri,
             vector[Read_count_table]& jcec_se,
             vector[Read_count_table]& jcec_mxe, vector[Read_count_table]& jcec_alt3,
             vector[Read_count_table]& jcec_alt5, vector[Read_count_table]& jcec_ri,
             int sam1len, cbool stat):
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

    for i in range(jc_se.size()):
        fprintf(se_fp, count_tmp, i,
                cjoin[vector[int].iterator](jc_se[i].incv.begin(), jc_se[i].incv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jc_se[i].skpv.begin(), jc_se[i].skpv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jc_se[i].incv.begin()+sam1len, jc_se[i].incv.end(), ',').c_str(),
                cjoin[vector[int].iterator](jc_se[i].skpv.begin()+sam1len, jc_se[i].skpv.end(), ',').c_str(),
                jc_se[i].inc_len, jc_se[i].skp_len)
        fprintf(se_fp_n, count_tmp, i,
                cjoin[vector[int].iterator](jcec_se[i].incv.begin(), jcec_se[i].incv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jcec_se[i].skpv.begin(), jcec_se[i].skpv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jcec_se[i].incv.begin()+sam1len, jcec_se[i].incv.end(), ',').c_str(),
                cjoin[vector[int].iterator](jcec_se[i].skpv.begin()+sam1len, jcec_se[i].skpv.end(), ',').c_str(),
                jcec_se[i].inc_len, jcec_se[i].skp_len)

    for i in range(jc_mxe.size()):
        if jc_mxe[i].strand == plus_mark:
            fprintf(mxe_fp, count_tmp, i,
                    cjoin[vector[int].iterator](jc_mxe[i].incv.begin(), jc_mxe[i].incv.begin()+sam1len, ',').c_str(),
                    cjoin[vector[int].iterator](jc_mxe[i].skpv.begin(), jc_mxe[i].skpv.begin()+sam1len, ',').c_str(),
                    cjoin[vector[int].iterator](jc_mxe[i].incv.begin()+sam1len, jc_mxe[i].incv.end(), ',').c_str(),
                    cjoin[vector[int].iterator](jc_mxe[i].skpv.begin()+sam1len, jc_mxe[i].skpv.end(), ',').c_str(),
                    jc_mxe[i].inc_len, jc_mxe[i].skp_len)
            fprintf(mxe_fp_n, count_tmp, i,
                    cjoin[vector[int].iterator](jcec_mxe[i].incv.begin(), jcec_mxe[i].incv.begin()+sam1len, ',').c_str(),
                    cjoin[vector[int].iterator](jcec_mxe[i].skpv.begin(), jcec_mxe[i].skpv.begin()+sam1len, ',').c_str(),
                    cjoin[vector[int].iterator](jcec_mxe[i].incv.begin()+sam1len, jcec_mxe[i].incv.end(), ',').c_str(),
                    cjoin[vector[int].iterator](jcec_mxe[i].skpv.begin()+sam1len, jcec_mxe[i].skpv.end(), ',').c_str(),
                    jcec_mxe[i].inc_len, jcec_mxe[i].skp_len)
        else:
            fprintf(mxe_fp, count_tmp, i,
                    cjoin[vector[int].iterator](jc_mxe[i].skpv.begin(), jc_mxe[i].skpv.begin()+sam1len, ',').c_str(),
                    cjoin[vector[int].iterator](jc_mxe[i].incv.begin(), jc_mxe[i].incv.begin()+sam1len, ',').c_str(),
                    cjoin[vector[int].iterator](jc_mxe[i].skpv.begin()+sam1len, jc_mxe[i].skpv.end(), ',').c_str(),
                    cjoin[vector[int].iterator](jc_mxe[i].incv.begin()+sam1len, jc_mxe[i].incv.end(), ',').c_str(),
                    jc_mxe[i].inc_len, jc_mxe[i].skp_len)
            fprintf(mxe_fp_n, count_tmp, i,
                    cjoin[vector[int].iterator](jcec_mxe[i].skpv.begin(), jcec_mxe[i].skpv.begin()+sam1len, ',').c_str(),
                    cjoin[vector[int].iterator](jcec_mxe[i].incv.begin(), jcec_mxe[i].incv.begin()+sam1len, ',').c_str(),
                    cjoin[vector[int].iterator](jcec_mxe[i].skpv.begin()+sam1len, jcec_mxe[i].skpv.end(), ',').c_str(),
                    cjoin[vector[int].iterator](jcec_mxe[i].incv.begin()+sam1len, jcec_mxe[i].incv.end(), ',').c_str(),
                    jcec_mxe[i].inc_len, jcec_mxe[i].skp_len)

    for i in range(jc_alt3.size()):
        fprintf(alt3_fp, count_tmp, i,
                cjoin[vector[int].iterator](jc_alt3[i].incv.begin(), jc_alt3[i].incv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jc_alt3[i].skpv.begin(), jc_alt3[i].skpv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jc_alt3[i].incv.begin()+sam1len, jc_alt3[i].incv.end(), ',').c_str(),
                cjoin[vector[int].iterator](jc_alt3[i].skpv.begin()+sam1len, jc_alt3[i].skpv.end(), ',').c_str(),
                jc_alt3[i].inc_len, jc_alt3[i].skp_len)
        fprintf(alt3_fp_n, count_tmp, i,
                cjoin[vector[int].iterator](jcec_alt3[i].incv.begin(), jcec_alt3[i].incv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jcec_alt3[i].skpv.begin(), jcec_alt3[i].skpv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jcec_alt3[i].incv.begin()+sam1len, jcec_alt3[i].incv.end(), ',').c_str(),
                cjoin[vector[int].iterator](jcec_alt3[i].skpv.begin()+sam1len, jcec_alt3[i].skpv.end(), ',').c_str(),
                jcec_alt3[i].inc_len, jcec_alt3[i].skp_len)

    for i in range(jc_alt5.size()):
        fprintf(alt5_fp, count_tmp, i,
                cjoin[vector[int].iterator](jc_alt5[i].incv.begin(), jc_alt5[i].incv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jc_alt5[i].skpv.begin(), jc_alt5[i].skpv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jc_alt5[i].incv.begin()+sam1len, jc_alt5[i].incv.end(), ',').c_str(),
                cjoin[vector[int].iterator](jc_alt5[i].skpv.begin()+sam1len, jc_alt5[i].skpv.end(), ',').c_str(),
                jc_alt5[i].inc_len, jc_alt5[i].skp_len)
        fprintf(alt5_fp_n, count_tmp, i,
                cjoin[vector[int].iterator](jcec_alt5[i].incv.begin(), jcec_alt5[i].incv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jcec_alt5[i].skpv.begin(), jcec_alt5[i].skpv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jcec_alt5[i].incv.begin()+sam1len, jcec_alt5[i].incv.end(), ',').c_str(),
                cjoin[vector[int].iterator](jcec_alt5[i].skpv.begin()+sam1len, jcec_alt5[i].skpv.end(), ',').c_str(),
                jcec_alt5[i].inc_len, jcec_alt5[i].skp_len)

    for i in range(jc_ri.size()):
        fprintf(ri_fp, count_tmp, i,
                cjoin[vector[int].iterator](jc_ri[i].incv.begin(), jc_ri[i].incv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jc_ri[i].skpv.begin(), jc_ri[i].skpv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jc_ri[i].incv.begin()+sam1len, jc_ri[i].incv.end(), ',').c_str(),
                cjoin[vector[int].iterator](jc_ri[i].skpv.begin()+sam1len, jc_ri[i].skpv.end(), ',').c_str(),
                jc_ri[i].inc_len, jc_ri[i].skp_len)
        fprintf(ri_fp_n, count_tmp, i,
                cjoin[vector[int].iterator](jcec_ri[i].incv.begin(), jcec_ri[i].incv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jcec_ri[i].skpv.begin(), jcec_ri[i].skpv.begin()+sam1len, ',').c_str(),
                cjoin[vector[int].iterator](jcec_ri[i].incv.begin()+sam1len, jcec_ri[i].incv.end(), ',').c_str(),
                cjoin[vector[int].iterator](jcec_ri[i].skpv.begin()+sam1len, jcec_ri[i].skpv.end(), ',').c_str(),
                jcec_ri[i].inc_len, jcec_ri[i].skp_len)

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
cdef save_exons(fp, unordered_map[string,cmap[Tetrad,pair[int,int]]]& exons):
    cdef:
        string line
        unordered_map[string,cmap[Tetrad,pair[int,int]]].iterator iexons
        cmap[Tetrad,pair[int,int]].iterator imap

    line = '%d\n' % (exons.size())
    fp.write(line)

    iexons = exons.begin()
    while iexons != exons.end():
        line = deref(iexons).first
        imap = deref(iexons).second.begin()
        while imap != deref(iexons).second.end():
            line = '%s;%d,%d,%d,%d,%d,%d' % (line, deref(imap).first.first,
                    deref(imap).first.second, deref(imap).first.third,
                    deref(imap).first.fourth, deref(imap).second.first,
                    deref(imap).second.second)
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
              vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
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
cdef size_t _load_job(str rmatsf, list vbams, list prep_counts_by_bam,
                      vector[unordered_map[string,vector[Triad]]]& novel_juncs,
                      vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
                      vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs,
                      int mode):
    cdef:
        int i = 0, j = 0, k = 0, num = 0, idx = 0
        size_t vlen
        Triad triad
        Tetrad tetrad
        vector[pair[long,long]] vp
        str line, gene_id
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
        # Skip over read length line. Already handled in split_sg_files_by_bam.
        fp.readline()
        for i in range(len(bams)):
            idx = try_get_index(vbams, bams[i], &index_found)
            if not index_found:
            # It's ok to have data for bams besides vbams
                continue

            novel_juncs[idx].clear()
            exons[idx].clear()
            juncs[idx].clear()
            prep_counts_by_bam[idx] += 1

        # processing novel junctions
        for i in range(len(bams)):
            num = int(fp.readline())
            idx = try_get_index(vbams, bams[i], &index_found)

            for j in range(num):
                line = fp.readline().strip()

                # Still need to read past the lines for this bam
                # even if not index_found
                if mode == read_mode or not index_found:
                    continue

                eles = line.split(';')
                gene_id = eles[0]

                for line in eles[1:]:
                    ele = [int(s) for s in line.split(',')]
                    triad.set(ele[0], ele[1], ele[2])
                    novel_juncs[idx][gene_id].push_back(triad)

        if mode == sg_mode:
            return len(bams)

        # processing exonic reads
        for i in range(len(bams)):
            num = int(fp.readline())
            idx = try_get_index(vbams, bams[i], &index_found)
            for j in range(num):
                line = fp.readline().strip()

                if not index_found:
                    continue

                eles = line.split(';')
                gene_id = eles[0]

                for line in eles[1:]:
                    ele = [int(s) for s in line.split(',')]
                    tetrad.set(ele[0], ele[1], ele[2], ele[3])
                    exons[idx][gene_id][tetrad].first = ele[4]
                    exons[idx][gene_id][tetrad].second = ele[5]

        # processing junction reads
        for i in range(len(bams)):
            num = int(fp.readline())
            idx = try_get_index(vbams, bams[i], &index_found)

            for j in range(num):
                line = fp.readline().strip()

                if not index_found:
                    continue

                eles = line.split(';')
                gene_id = eles[0]

                for line in eles[1:]:
                    ele = [s for s in line.split(',')]
                    aligns = ele[0].split('=')
                    vp = vector[pair[long,long]](len(aligns))
                    for k in range(len(aligns)):
                        coords = [int(s) for s in aligns[k].split(':')]
                        vp[k].first = coords[0]
                        vp[k].second = coords[1]
                    juncs[idx][gene_id][vp] = int(ele[1])

        return len(bams)


@boundscheck(False)
@wraparound(False)
cdef load_sg(str bams, list dot_rmats_paths,
             vector[unordered_map[string,vector[Triad]]]& novel_juncs):
    cdef:
        int num = 0, num_file = 0
        list vbams = bams.split(',')
        list prep_counts_by_bam
        vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]] exons
        vector[unordered_map[string,cmap[vector[pair[long,long]],int]]] juncs

    num = len(vbams)
    prep_counts_by_bam = [0 for i in range(num)]

    for name in dot_rmats_paths:
        _load_job(name, vbams, prep_counts_by_bam, novel_juncs, exons, juncs, sg_mode)

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
cdef vector[size_t] load_read(str bams, str fn,
                              vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]]& exons,
                              vector[unordered_map[string,cmap[vector[pair[long,long]],int]]]& juncs):
    cdef:
        int num = 0, num_file = 0
        list vbams = bams.split(',')
        list prep_counts_by_bam
        vector[unordered_map[string,vector[Triad]]] novel_juncs
        vector[size_t] residx

    num = len(vbams)
    prep_counts_by_bam = [0 for i in range(num)]

    _load_job(fn, vbams, prep_counts_by_bam, novel_juncs, exons, juncs, read_mode)
    residx = [i for i in range(num) if prep_counts_by_bam[i] == 1]

    return residx


def run_pipe(args):
    """TODO: Docstring for run_pipe.
    :returns: TODO

    """
    cdef:
        unordered_map[int,cset[string]] geneGroup
        unordered_map[string,Gene] genes
        unordered_map[string,SupInfo] supple
        vector[unordered_map[string,vector[Triad]]] novel_juncs
        vector[unordered_map[string,cmap[Tetrad,pair[int,int]]]] exons
        vector[unordered_map[string,cmap[string,int]]] multis
        cset[SE_info] se,
        cset[MXE_info] mxe,
        cset[ALT35_info] alt3,
        cset[ALT35_info] alt5,
        cset[RI_info] ri,
        int sam1len = len(args.b1.split(','))

    start = time.time()
    parse_gtf(args.gtf, geneGroup, genes, supple)
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
        detect_ase(genes, supple, args.od, novel_juncs,
                   se, mxe, alt3, alt5, ri,
                   args.junctionLength/2, args.readLength, args.novelSS, args.mel)
        print 'ase:', time.time() - start

        start = time.time()
        count_occurrence(args.bams, dot_rmats_file_paths, args.od, se, mxe,
                         alt3, alt5, ri, sam1len, args.junctionLength/2,
                         args.readLength, args.nthread, args.dt, args.stat)
        print 'count:', time.time() - start

        shutil.rmtree(split_dot_rmats_dir_path)
