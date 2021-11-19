from libcpp cimport bool as cbool
from libcpp.set cimport set as cset
from libcpp.map cimport map as cmap
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.algorithm cimport sort, lower_bound
from libc.stdio cimport FILE, fprintf, fopen, fclose


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
        int inc_len
        int skp_len
        int inc_len_jcec
        int skp_len_jcec
        cbool txtype
        cbool includes_novel_ss
        void set(const int& iiid, const string& igID, const SupInfo& isupInfo,
                 const long& its, const long& ite, const long& ius, const long& iue,
                 const long& ids, const long& ide, const int& il,
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
        int inc_len
        int skp_len
        int inc_len_jcec
        int skp_len_jcec
        cbool txtype
        cbool includes_novel_ss
        void set(const int& iiid, const string& igID, const SupInfo& isupInfo,
                 const long& its, const long& ite, const long& iss, const long& ise,
                 const long& ius, const long& iue, const long& ids, const long& ide,
                 const int& il, const int& sl,
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
        int inc_len
        int skp_len
        int inc_len_jcec
        int skp_len_jcec
        cbool txtype
        cbool includes_novel_ss
        void set(const int& iiid, const string& igID, const SupInfo& isupInfo,
                 const long& ils, const long& ile, const long& iss, const long& ise,
                 const long& ifs, const long& ife, const int& il,
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
        int inc_len
        int skp_len
        int inc_len_jcec
        int skp_len_jcec
        cbool txtype
        cbool includes_novel_ss
        void set(const int& iiid, const string& igID, const SupInfo& isupInfo,
                 const long& irs, const long& ire, const long& ius, const long& iue,
                 const long& ids, const long& ide, const int& il,
                 const int& sl, const int& iljcec, const int& sljcec,
                 const cbool& itxtype, const cbool iincludes_novel_ss)

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

    string cjoin[T](T, T, const char&)
    char* join_pair(char* numstr, long left, long right, char& sep)
    string num2str[T](T Number)
    void insert_set[InputType](cset[InputType]& iset, cset[InputType].iterator first, cset[InputType].iterator last)
