#ifndef RMATS_TCX
#define RMATS_TCX


#include <set>
#include <map>
#include <tuple>
#include <vector>
#include <string>
#include <utility>
#include <sstream>
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
        int inc_len;
        int skp_len;
        int inc_len_jcec;
        int skp_len_jcec;
        bool txtype;
        bool includes_novel_ss;
        void set(const int& iiid, const std::string& igID, const SupInfo& isupInfo,
                 const long& its, const long& ite, const long& ius, const long& iue,
                 const long& ids, const long& ide, const int& il,
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
        int inc_len;
        int skp_len;
        int inc_len_jcec;
        int skp_len_jcec;
        bool txtype;
        bool includes_novel_ss;
        void set(const int& iiid, const std::string& igID, const SupInfo& isupInfo,
                 const long& its, const long& ite, const long& iss, const long& ise,
                 const long& ius, const long& iue, const long& ids, const long& ide,
                 const int& il, const int& sl,
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
        int inc_len;
        int skp_len;
        int inc_len_jcec;
        int skp_len_jcec;
        bool txtype;
        bool includes_novel_ss;
        void set(const int& iiid, const std::string& igID, const SupInfo& isupInfo,
                 const long& ils, const long& ile, const long& iss, const long& ise,
                 const long& ifs, const long& ife, const int& il,
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
        int inc_len;
        int skp_len;
        int inc_len_jcec;
        int skp_len_jcec;
        bool txtype;
        bool includes_novel_ss;
        void set(const int& iiid, const std::string& igID, const SupInfo& isupInfo,
                 const long& irs, const long& ire, const long& ius, const long& iue,
                 const long& ids, const long& ide, const int& il,
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

    char* join_pair(char* numstr, long left, long right, char& sep) {
        sprintf(numstr, "%ld%c%ld", left, sep, right);
        return numstr;
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
}


#endif
