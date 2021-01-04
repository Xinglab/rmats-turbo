from libcpp cimport bool as cbool
from libcpp.set cimport set as cset
from libcpp.map cimport map as cmap
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_set cimport unordered_set
from libcpp.unordered_map cimport unordered_map
from libcpp.algorithm cimport sort, lower_bound
from libc.string cimport strcmp, strlen
from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE, fprintf, fopen, fclose


# ctypedef unordered_map.iterator map_iter


cdef extern from '<algorithm>' namespace 'std' nogil:
    void fill[ForwardIt,T](ForwardIt first, ForwardIt last, const T& value)
    OutputIt set_union[InputIt1, InputIt2, OutputIt](InputIt1 first1,
            InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt d_first)


cdef extern from '<string>' namespace 'std' nogil:
    long stol(const string&, size_t*, int)
    long stol(const string&, size_t*)
    long stol(const string&)


cdef extern from 'tcx.h' namespace 'rmats' nogil:
    cdef cppclass SE_key:
        long first
        long second
        long third
        long fourth
        string chrom

    cdef cppclass SE_info:
        int iid
        string gID
        SupInfo supInfo
        long ts
        long te
        long us
        long ue
        long ds
        long de
        size_t tidx
        size_t uidx
        size_t didx
        int inc_len
        int skp_len
        int inc_len_jcec
        int skp_len_jcec
        cbool txtype
        cbool includes_novel_ss
        void set(const int& iiid, const string& igID, const SupInfo& isupInfo,
                 const long& its, const long& ite, const long& ius, const long& iue,
                 const long& ids, const long& ide, const size_t& itidx,
                 const size_t& iuidx, const size_t& ididx, const int& il,
                 const int& sl, const int& iljcec, const int& sljcec,
                 const cbool& itxtype, const cbool iincludes_novel_ss)

    cdef cppclass MXE_key:
        long mxe_first
        long mxe_second
        long first
        long second
        long third
        long fourth
        string chrom

    cdef cppclass MXE_info:
        int iid
        string gID
        SupInfo supInfo
        long ts
        long te
        long ss
        long se
        long us
        long ue
        long ds
        long de
        size_t tidx
        size_t sidx
        size_t uidx
        size_t didx
        int inc_len
        int skp_len
        int inc_len_jcec
        int skp_len_jcec
        cbool txtype
        cbool includes_novel_ss
        void set(const int& iiid, const string& igID, const SupInfo& isupInfo,
                 const long& its, const long& ite, const long& iss, const long& ise,
                 const long& ius, const long& iue, const long& ids, const long& ide,
                 const size_t& itidx, const size_t& isidx, const size_t& iuidx,
                 const size_t& ididx, const int& il, const int& sl,
                 const int& iljcec, const int& sljcec, const cbool& itxtype,
                 const cbool iincludes_novel_ss)

    cdef cppclass ALT35_key:
        long first
        long second
        long third
        string chrom

    cdef cppclass ALT35_info:
        int iid
        string gID
        SupInfo supInfo
        long ls
        long le
        long ss
        long se
        long fs
        long fe
        size_t lidx
        size_t sidx
        size_t fidx
        int inc_len
        int skp_len
        int inc_len_jcec
        int skp_len_jcec
        cbool txtype
        cbool includes_novel_ss
        void set(const int& iiid, const string& igID, const SupInfo& isupInfo,
                 const long& ils, const long& ile, const long& iss, const long& ise,
                 const long& ifs, const long& ife, const size_t& ilidx,
                 const size_t& isidx, const size_t& ifidx, const int& il,
                 const int& sl, const int& iljcec, const int& sljcec,
                 const cbool& itxtype, const cbool iincludes_novel_ss)

    cdef cppclass RI_key:
        long first
        long second
        string chrom

    cdef cppclass RI_info:
        int iid
        string gID
        SupInfo supInfo
        long rs
        long re
        long us
        long ue
        long ds
        long de
        size_t ridx
        size_t uidx
        size_t didx
        int inc_len
        int skp_len
        int inc_len_jcec
        int skp_len_jcec
        cbool txtype
        cbool includes_novel_ss
        void set(const int& iiid, const string& igID, const SupInfo& isupInfo,
                 const long& irs, const long& ire, const long& ius, const long& iue,
                 const long& ids, const long& ide, const size_t& iridx,
                 const size_t& iuidx, const size_t& ididx, const int& il,
                 const int& sl, const int& iljcec, const int& sljcec,
                 const cbool& itxtype, const cbool iincludes_novel_ss)

    cdef cppclass Tuple:
        Tuple()
        int first
        int second
        cbool operator<(Tuple&)

    cdef cppclass Triad:
        Triad()
        long left
        long mid
        long right
        cbool operator<(Triad&)
        void set(const long& ileft, const long& imid, const long& iright)

    cdef cppclass Tetrad:
        Tetrad()
        long first
        long second
        long third
        long fourth
        cbool operator<(Tetrad&)
        void set(const long& ifirst, const long& isecond, const long& ithird, const long& ifourth)

    cdef cppclass MNMNM:
        long coord1
        long coord2
        long coord3
        long coord4
        long coord5
        long coord6
        cbool operator<(MNMNM& t)
        void set(const long& i1, const long& i2, const long& i3, const long& i4, const long& i5, const long& i6)

    cdef cppclass SupInfo:
        string g_name
        string chrom
        char strand
        void set_info(string, string, string)

    cdef cppclass Transcript:
        vector[pair[long,long]] exons
        unordered_map[long,size_t] first
        unordered_map[long,size_t] second

    cdef cppclass Gene:
        Gene()
        cmap[pair[long,long],size_t] exon_idx
        vector[pair[long,long]] idx_exon
        vector[vector[cset[pair[size_t,cbool]]]] sg
        unordered_map[string,Transcript] trans

    cdef cppclass Exon:
        char strand
        Tetrad tetrad
        void set(const char& istrand, const long& ifirst, const long& isecond, const long& ithird, const long& ifourth)

    cdef cppclass Inc_skp_len:
        Inc_skp_len()
        int inc_len
        int skp_len

    cdef cppclass Inc_skp_count:
        Inc_skp_count()
        int inc_count
        int skp_count

    cdef cppclass String_and_stream:
        const char* c_str()

    cdef cppclass SE_counts_for_event:
        SE_counts_for_event()
        Inc_skp_count jc_counts
        Inc_skp_count jcec_counts
        int upstream_to_target_count
        int target_to_downstream_count
        int target_count
        int upstream_to_downstream_count

    cdef cppclass SE_joined_count_strings:
        String_and_stream jc_inc_1
        String_and_stream jc_skp_1
        String_and_stream jc_inc_2
        String_and_stream jc_skp_2
        String_and_stream jcec_inc_1
        String_and_stream jcec_skp_1
        String_and_stream jcec_inc_2
        String_and_stream jcec_skp_2
        String_and_stream upstream_to_target
        String_and_stream target_to_downstream
        String_and_stream target
        String_and_stream upstream_to_downstream

        void clear()

    cdef cppclass SE_counts_for_event_by_bam:
      char strand
      Inc_skp_len jc_lengths
      Inc_skp_len jcec_lengths
      vector[SE_counts_for_event] counts

      void join_counts_across_bams(
          int sam1len,
          SE_joined_count_strings* joined_strings) const

    cdef cppclass MXE_counts_for_event:
        MXE_counts_for_event()
        Inc_skp_count jc_counts
        Inc_skp_count jcec_counts
        int upstream_to_first_count
        int first_to_downstream_count
        int first_count
        int upstream_to_second_count
        int second_to_downstream_count
        int second_count

    cdef cppclass MXE_joined_count_strings:
        String_and_stream jc_inc_1
        String_and_stream jc_skp_1
        String_and_stream jc_inc_2
        String_and_stream jc_skp_2
        String_and_stream jcec_inc_1
        String_and_stream jcec_skp_1
        String_and_stream jcec_inc_2
        String_and_stream jcec_skp_2
        String_and_stream upstream_to_first
        String_and_stream first_to_downstream
        String_and_stream first
        String_and_stream upstream_to_second
        String_and_stream second_to_downstream
        String_and_stream second

        void clear()

    cdef cppclass MXE_counts_for_event_by_bam:
      char strand
      Inc_skp_len jc_lengths
      Inc_skp_len jcec_lengths
      vector[MXE_counts_for_event] counts

      void join_counts_across_bams(
          int sam1len,
          MXE_joined_count_strings* joined_strings) const

    cdef cppclass ALT35_counts_for_event:
        ALT35_counts_for_event()
        Inc_skp_count jc_counts
        Inc_skp_count jcec_counts
        int across_short_boundary_count
        int long_to_flanking_count
        int exclusive_to_long_count
        int short_to_flanking_count

    cdef cppclass ALT35_joined_count_strings:
        String_and_stream jc_inc_1
        String_and_stream jc_skp_1
        String_and_stream jc_inc_2
        String_and_stream jc_skp_2
        String_and_stream jcec_inc_1
        String_and_stream jcec_skp_1
        String_and_stream jcec_inc_2
        String_and_stream jcec_skp_2
        String_and_stream across_short_boundary
        String_and_stream long_to_flanking
        String_and_stream exclusive_to_long
        String_and_stream short_to_flanking

        void clear()

    cdef cppclass ALT35_counts_for_event_by_bam:
      char strand
      Inc_skp_len jc_lengths
      Inc_skp_len jcec_lengths
      vector[ALT35_counts_for_event] counts

      void join_counts_across_bams(
          int sam1len,
          ALT35_joined_count_strings* joined_strings) const

    cdef cppclass RI_counts_for_event:
        RI_counts_for_event()
        Inc_skp_count jc_counts
        Inc_skp_count jcec_counts
        int upstream_to_intron_count
        int intron_to_downstream_count
        int intron_count
        int upstream_to_downstream_count

    cdef cppclass RI_joined_count_strings:
        String_and_stream jc_inc_1
        String_and_stream jc_skp_1
        String_and_stream jc_inc_2
        String_and_stream jc_skp_2
        String_and_stream jcec_inc_1
        String_and_stream jcec_skp_1
        String_and_stream jcec_inc_2
        String_and_stream jcec_skp_2
        String_and_stream upstream_to_intron
        String_and_stream intron_to_downstream
        String_and_stream intron
        String_and_stream upstream_to_downstream

        void clear()

    cdef cppclass RI_counts_for_event_by_bam:
      char strand
      Inc_skp_len jc_lengths
      Inc_skp_len jcec_lengths
      vector[RI_counts_for_event] counts

      void join_counts_across_bams(
          int sam1len,
          RI_joined_count_strings* joined_strings) const

    cdef cppclass Str_ptr:
        const string *p
        Str_ptr(const string& s)
        const string& get() const

    void insert_str_ptr(cset[Str_ptr]& iset, cset[string].iterator first, cset[string].iterator last)

    void cprint[T](T)
    T3 at_map[T1, T2, T3](cmap[T1, T2] imap, int idx)
    vector[string]& split(const string, char, vector[string])
    vector[string] split(const string, char)
    string cjoin[T](T, T, const char&)
    string join_cigar[T](long& mc, vector[T]& CigarData, char& sep)
    char* join_tetrad(char* numstr, Tetrad& tetrad, char& sep)
    char* join_pair(char* numstr, long left, long right, char& sep)
    string join_pair_vector(vector[pair[long,long]].iterator,
                            vector[pair[long,long]].iterator, char&)
    string to_novel_txname(string, string)
    string num2str[T](T Number)
    void insert_set[InputType](cset[InputType]& iset, cset[InputType].iterator first, cset[InputType].iterator last)
    void insert_unset[InputType](unordered_set[InputType]& iset, cset[InputType].iterator first, cset[InputType].iterator last)
