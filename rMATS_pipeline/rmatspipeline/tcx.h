#ifndef RMATS_TCX
#define RMATS_TCX


#include <set>
#include <map>
#include <tuple>
#include <vector>
#include <string>
#include <utility>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <iostream>


namespace rmats {
    struct SupInfo
    {
        std::string g_name;
        std::string chrom;
        char strand;
        void set_info(std::string& iname, std::string& ichrom, std::string& istrand) {
            this->g_name = iname;
            this->chrom = ichrom;
            this->strand = istrand.c_str()[0];
        }
    };

    struct SE_key
    {
        long first;
        long second;
        long third;
        long fourth;
        std::string chrom;
        bool operator<(const SE_key& t) const {
            return std::tie(this->first, this->second, this->third, this->fourth, this->chrom) < std::tie(t.first, t.second, t.third, t.fourth, t.chrom);
        }
    };

    struct SE_info
    {
        int iid;
        std::string gID;
        SupInfo supInfo;
        long ts;
        long te;
        long us;
        long ue;
        long ds;
        long de;
        size_t tidx;
        size_t uidx;
        size_t didx;
        int inc_len;
        int skp_len;
        int inc_len_jcec;
        int skp_len_jcec;
        bool txtype;
        bool includes_novel_ss;
        void set(const int& iiid, const std::string& igID, const SupInfo& isupInfo,
                 const long& its, const long& ite, const long& ius, const long& iue,
                 const long& ids, const long& ide, const size_t& itidx,
                 const size_t& iuidx, const size_t& ididx, const int& il,
                 const int& sl, const int& iljcec, const int& sljcec,
                 const bool& itxtype, const bool iincludes_novel_ss) {
            if (iiid >= 0) {
                this->iid = iiid;
            }
            this->gID = igID;
            this->supInfo = isupInfo;
            this->ts = its;
            this->te = ite;
            this->us = ius;
            this->ue = iue;
            this->ds = ids;
            this->de = ide;
            this->tidx = itidx;
            this->uidx = iuidx;
            this->didx = ididx;
            this->inc_len = il;
            this->skp_len = sl;
            this->inc_len_jcec = iljcec;
            this->skp_len_jcec = sljcec;
            this->txtype = itxtype;
            this->includes_novel_ss = iincludes_novel_ss;
        }
        bool operator<(const SE_info& t) const {
            return std::tie(this->iid) < std::tie(t.iid);
        }
    };

    struct MXE_key
    {
        long mxe_first;
        long mxe_second;
        long first;
        long second;
        long third;
        long fourth;
        std::string chrom;
        bool operator<(const MXE_key& t) const {
            return std::tie(this->mxe_first, this->mxe_second, this->first, this->second, this->third, this->fourth, this->chrom) < std::tie(t.mxe_first, t.mxe_second, t.first, t.second, t.third, t.fourth, t.chrom);
        }
    };

    struct MXE_info
    {
        int iid;
        std::string gID;
        SupInfo supInfo;
        long ts;
        long te;
        long ss;
        long se;
        long us;
        long ue;
        long ds;
        long de;
        size_t tidx;
        size_t sidx;
        size_t uidx;
        size_t didx;
        int inc_len;
        int skp_len;
        int inc_len_jcec;
        int skp_len_jcec;
        bool txtype;
        bool includes_novel_ss;
        void set(const int& iiid, const std::string& igID, const SupInfo& isupInfo,
                 const long& its, const long& ite, const long& iss, const long& ise,
                 const long& ius, const long& iue, const long& ids, const long& ide,
                 const size_t& itidx, const size_t& isidx, const size_t& iuidx,
                 const size_t& ididx, const int& il, const int& sl,
                 const int& iljcec, const int& sljcec, const bool& itxtype,
                 const bool iincludes_novel_ss) {
            if (iiid >= 0) {
                this->iid = iiid;
            }
            this->gID = igID;
            this->supInfo = isupInfo;
            this->ts = its;
            this->te = ite;
            this->ss = iss;
            this->se = ise;
            this->us = ius;
            this->ue = iue;
            this->ds = ids;
            this->de = ide;
            this->tidx = itidx;
            this->sidx = isidx;
            this->uidx = iuidx;
            this->didx = ididx;
            this->inc_len = il;
            this->skp_len = sl;
            this->inc_len_jcec = iljcec;
            this->skp_len_jcec = sljcec;
            this->txtype = itxtype;
            this->includes_novel_ss = iincludes_novel_ss;
        }
        bool operator<(const MXE_info& t) const {
            return std::tie(this->iid) < std::tie(t.iid);
        }
    };

    struct ALT35_key
    {
        long first;
        long second;
        long third;
        std::string chrom;
        bool operator<(const ALT35_key& t) const {
            return std::tie(this->first, this->second, this->third, this->chrom) < std::tie(t.first, t.second, t.third, t.chrom);
        }
    };

    struct ALT35_info
    {
        int iid;
        std::string gID;
        SupInfo supInfo;
        long ls;
        long le;
        long ss;
        long se;
        long fs;
        long fe;
        size_t lidx;
        size_t sidx;
        size_t fidx;
        int inc_len;
        int skp_len;
        int inc_len_jcec;
        int skp_len_jcec;
        bool txtype;
        bool includes_novel_ss;
        void set(const int& iiid, const std::string& igID, const SupInfo& isupInfo,
                 const long& ils, const long& ile, const long& iss, const long& ise,
                 const long& ifs, const long& ife, const size_t& ilidx,
                 const size_t& isidx, const size_t& ifidx, const int& il,
                 const int& sl, const int& iljcec, const int& sljcec,
                 const bool& itxtype, const bool iincludes_novel_ss) {
            if (iiid >= 0) {
                this->iid = iiid;
            }
            this->gID = igID;
            this->supInfo = isupInfo;
            this->ls = ils;
            this->le = ile;
            this->ss = iss;
            this->se = ise;
            this->fs = ifs;
            this->fe = ife;
            this->lidx = ilidx;
            this->sidx = isidx;
            this->fidx = ifidx;
            this->inc_len = il;
            this->skp_len = sl;
            this->inc_len_jcec = iljcec;
            this->skp_len_jcec = sljcec;
            this->txtype = itxtype;
            this->includes_novel_ss = iincludes_novel_ss;
        }
        bool operator<(const ALT35_info& t) const {
            return std::tie(this->iid) < std::tie(t.iid);
        }
    };

    struct RI_key
    {
        long first;
        long second;
        std::string chrom;
        bool operator<(const RI_key& t) const {
            return std::tie(this->first, this->second, this->chrom) < std::tie(t.first, t.second, t.chrom);
        }
    };

    struct RI_info
    {
        int iid;
        std::string gID;
        SupInfo supInfo;
        long rs;
        long re;
        long us;
        long ue;
        long ds;
        long de;
        size_t ridx;
        size_t uidx;
        size_t didx;
        int inc_len;
        int skp_len;
        int inc_len_jcec;
        int skp_len_jcec;
        bool txtype;
        bool includes_novel_ss;
        void set(const int& iiid, const std::string& igID, const SupInfo& isupInfo,
                 const long& irs, const long& ire, const long& ius, const long& iue,
                 const long& ids, const long& ide, const size_t& iridx,
                 const size_t& iuidx, const size_t& ididx, const int& il,
                 const int& sl, const int& iljcec, const int& sljcec,
                 const bool& itxtype, const bool iincludes_novel_ss) {
            if (iiid >= 0) {
                this->iid = iiid;
            }
            this->gID = igID;
            this->supInfo = isupInfo;
            this->rs = irs;
            this->re = ire;
            this->us = ius;
            this->ue = iue;
            this->ds = ids;
            this->de = ide;
            this->ridx = iridx;
            this->uidx = iuidx;
            this->didx = ididx;
            this->inc_len = il;
            this->skp_len = sl;
            this->inc_len_jcec = iljcec;
            this->skp_len_jcec = sljcec;
            this->txtype = itxtype;
            this->includes_novel_ss = iincludes_novel_ss;
        }
        bool operator<(const RI_info& t) const {
            return std::tie(this->iid) < std::tie(t.iid);
        }
    };

    struct Tuple
    {
        int first;
        int second;
        bool operator<(const Tuple& t) const {
            return std::tie(this->first, this->second) < std::tie(t.first, t.second);
        }
    };

    struct Triad
    {
        long left;
        long mid;
        long right;
        bool operator<(const Triad& t) const {
            return std::tie(this->left, this->mid, this->right) < std::tie(t.left, t.mid, t.right);
        }
        void set(const long& ileft, const long& imid, const long& iright) {
            this->left = ileft;
            this->mid = imid;
            this->right = iright;
            return;
        }
    };

    struct Tetrad
    {
        long first;
        long second;
        long third;
        long fourth;
        bool operator<(const Tetrad& t) const {
            return std::tie(this->first, this->second, this->third, this->fourth) < std::tie(t.first, t.second, t.third, t.fourth);
        }
        void set(const long& ifirst, const long& isecond, const long& ithird, const long& ifourth) {
            this->first = ifirst;
            this->second = isecond;
            this->third = ithird;
            this->fourth = ifourth;
            return;
        }
    };

    struct MNMNM
    {
        long coord1;
        long coord2;
        long coord3;
        long coord4;
        long coord5;
        long coord6;
        bool operator<(const MNMNM& t) const {
            return std::tie(this->coord1, this->coord2, this->coord3, this->coord4, this->coord5, this->coord6) < std::tie(t.coord1, t.coord2, t.coord3, t.coord4, t.coord5, t.coord6);
        }
        void set(const long& i1, const long& i2, const long& i3, const long& i4, const long& i5, const long& i6) {
            this->coord1 = i1;
            this->coord2 = i2;
            this->coord3 = i3;
            this->coord4 = i4;
            this->coord5 = i5;
            this->coord6 = i6;
            return;
        }
    };

    struct Transcript
    {
        std::vector<std::pair<long,long> > exons;
        std::unordered_map<long,size_t> first;
        std::unordered_map<long,size_t> second;
    };

    struct Gene
    {
        std::map<std::pair<long,long>,size_t> exon_idx;
        std::vector<std::pair<long,long> > idx_exon;
        std::vector<std::vector<std::set<std::pair<size_t,bool> > > > sg;
        std::unordered_map<std::string,Transcript> trans;
    };

    struct Exon
    {
        char strand;
        Tetrad tetrad;
        bool operator<(const Exon& t) const {
            return std::tie(this->tetrad, this->strand) < std::tie(t.tetrad, t.strand);
        }
        void set(const char& istrand, const long& ifirst, const long& isecond, const long& ithird, const long& ifourth) {
            this->strand = istrand;
            this->tetrad.set(ifirst, isecond, ithird, ifourth);
            return;
        }
    };

    struct Read_count_table
    {
        std::vector<int> incv;
        std::vector<int> skpv;
        int inc_len;
        int skp_len;
        char strand;
    };

    struct Str_ptr
    {
        const std::string* p;
        Str_ptr(const std::string& s): p(&s) {}
        const std::string& get() const {
            return *(this->p);
        }
        bool operator<(const Str_ptr& t) const {
            return std::tie(*(this->p)) < std::tie(*(t.p));
        }
        bool operator<(const std::string& t) const {
            return std::tie(*(this->p)) < std::tie(t);
        }
    };

    void insert_str_ptr(std::set<Str_ptr>& iset, std::set<std::string>::iterator first, std::set<std::string>::iterator last) {
        std::set<std::string>::iterator tmp = first;
        while (tmp != last) {
            iset.insert(Str_ptr(*tmp));
            ++tmp;
        }
    }

    template <typename T>
    void cprint(T t) {
        std::cout << t << std::endl;
    }


    template <typename T1, typename T2, typename T3>
    T3 at_map(std::map<T1, T2> imap, int idx) {
        T3 iter = imap.begin();
        return iter+idx;
    }


    std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        return elems;
    }


    std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split(s, delim, elems);
        return elems;
    }


    template <typename Iter>
    std::string cjoin(Iter begin, Iter end, const char& separator)
    {
        std::ostringstream result;
        if (begin != end)
            result << *begin++;
        while (begin != end)
            result << separator << *begin++;
        return result.str();
    }


    template <typename T>
    std::string join_cigar(long& mc, std::vector<T>& CigarData, char& sep) {
        int i;
        std::ostringstream result;

        if (CigarData.size() > 0) {
            result << mc;
        }
        for (i = 0; i < CigarData.size(); ++i) {
            result << sep << CigarData[i].Length;
        }
        return result.str();
    }


    // is brittle (you must supply a large enough buffer), fast, and verbose; requires nothing (is standard C++); all platforms
    char* join_tetrad(char* numstr, Tetrad& tetrad, char& sep) {
        sprintf(numstr, "%ld%c%ld%c%ld%c%ld", tetrad.first, sep, tetrad.second, sep, tetrad.third, sep, tetrad.fourth);
        return numstr;
    }


    char* join_pair(char* numstr, long left, long right, char& sep) {
        sprintf(numstr, "%ld%c%ld", left, sep, right);
        return numstr;
    }


    std::string join_pair_vector(std::vector<std::pair<long,long> >::iterator begin,
                                 std::vector<std::pair<long,long> >::iterator end,
                                 char& separator) {
        std::ostringstream result;
        if (begin != end) {
            result << (*begin).first;
            result << separator << (*begin++).second;
        }
        while (begin != end) {
            result << separator << (*begin).first;
            result << separator << (*begin++).second;
        }
        return result.str();
    }


    std::string to_novel_txname(std::string cg, std::string novel) {
        return cg + ".n." + novel;
    }


    template<typename T> std::string num2str(T Number) {
        std::ostringstream ss;
        ss << Number;
        return ss.str();
    }

    template<class InputType, class InputIt>
    void insert_set(std::set<InputType>& iset, InputIt first,InputIt last) {
        iset.insert(first, last);
    }

    template<class InputType, class InputIt>
    void insert_unset(std::unordered_set<InputType>& iset, InputIt first,InputIt last) {
        iset.insert(first, last);
    }
}


#endif
