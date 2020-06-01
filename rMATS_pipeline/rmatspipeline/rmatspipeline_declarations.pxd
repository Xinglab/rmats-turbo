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

    cdef cppclass Read_count_table:
        vector[int] incv
        vector[int] skpv
        int inc_len
        int skp_len
        char strand

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
