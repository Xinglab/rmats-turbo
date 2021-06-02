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

    struct Inc_skp_len
    {
        Inc_skp_len() : inc_len(0), skp_len(0) {}
        int inc_len;
        int skp_len;
    };

    struct Inc_skp_count
    {
        Inc_skp_count() : inc_count(0), skp_count(0) {}
        int inc_count;
        int skp_count;
    };

    // oss.str() needs to be stored in order to use the pointer  from .c_str()
    struct String_and_stream
    {
        std::ostringstream oss;
        std::string s;

        void clear()
        {
            oss.str("");
            s.clear();
        }

        const char* c_str()
        {
            s = oss.str();
            return s.c_str();
        }
    };

    struct SE_counts_for_event
    {
        SE_counts_for_event()
            : upstream_to_target_count(0),
              target_to_downstream_count(0),
              target_count(0),
              upstream_to_downstream_count(0) {}
        Inc_skp_count jc_counts;
        Inc_skp_count jcec_counts;
        int upstream_to_target_count;
        int target_to_downstream_count;
        int target_count;
        int upstream_to_downstream_count;
    };

    struct SE_joined_count_strings
    {
        String_and_stream jc_inc_1;
        String_and_stream jc_skp_1;
        String_and_stream jc_inc_2;
        String_and_stream jc_skp_2;
        String_and_stream jcec_inc_1;
        String_and_stream jcec_skp_1;
        String_and_stream jcec_inc_2;
        String_and_stream jcec_skp_2;
        String_and_stream upstream_to_target;
        String_and_stream target_to_downstream;
        String_and_stream target;
        String_and_stream upstream_to_downstream;

        void clear()
        {
            jc_inc_1.clear();
            jc_skp_1.clear();
            jc_inc_2.clear();
            jc_skp_2.clear();
            jcec_inc_1.clear();
            jcec_skp_1.clear();
            jcec_inc_2.clear();
            jcec_skp_2.clear();
            upstream_to_target.clear();
            target_to_downstream.clear();
            target.clear();
            upstream_to_downstream.clear();
        }
    };

    struct SE_counts_for_event_by_bam
    {
        char strand;
        Inc_skp_len jc_lengths;
        Inc_skp_len jcec_lengths;
        std::vector<SE_counts_for_event> counts;

        void join_counts_across_bams(
            int sam1len,
            SE_joined_count_strings* joined_strings) const
        {
            joined_strings->clear();
            std::string maybe_comma = "";
            std::string sam2_maybe_comma = "";
            for (int i = 0; i < counts.size(); ++i)
            {
                const SE_counts_for_event& for_event = counts[i];
                if (i < sam1len)
                {
                    joined_strings->jc_inc_1.oss
                        << maybe_comma << for_event.jc_counts.inc_count;
                    joined_strings->jc_skp_1.oss
                        << maybe_comma << for_event.jc_counts.skp_count;
                    joined_strings->jcec_inc_1.oss
                        << maybe_comma << for_event.jcec_counts.inc_count;
                    joined_strings->jcec_skp_1.oss
                        << maybe_comma << for_event.jcec_counts.skp_count;
                }
                else
                {
                    joined_strings->jc_inc_2.oss
                        << sam2_maybe_comma << for_event.jc_counts.inc_count;
                    joined_strings->jc_skp_2.oss
                        << sam2_maybe_comma << for_event.jc_counts.skp_count;
                    joined_strings->jcec_inc_2.oss
                        << sam2_maybe_comma << for_event.jcec_counts.inc_count;
                    joined_strings->jcec_skp_2.oss
                        << sam2_maybe_comma << for_event.jcec_counts.skp_count;
                    sam2_maybe_comma = ",";
                }

                joined_strings->upstream_to_target.oss
                    << maybe_comma << for_event.upstream_to_target_count;
                joined_strings->target_to_downstream.oss
                    << maybe_comma << for_event.target_to_downstream_count;
                joined_strings->target.oss
                    << maybe_comma << for_event.target_count;
                joined_strings->upstream_to_downstream.oss
                    << maybe_comma << for_event.upstream_to_downstream_count;
                maybe_comma = ",";
            }
        }
    };

    struct MXE_counts_for_event
    {
        MXE_counts_for_event()
            : upstream_to_first_count(0),
              first_to_downstream_count(0),
              first_count(0),
              upstream_to_second_count(0),
              second_to_downstream_count(0),
              second_count(0) {}
        Inc_skp_count jc_counts;
        Inc_skp_count jcec_counts;
        int upstream_to_first_count;
        int first_to_downstream_count;
        int first_count;
        int upstream_to_second_count;
        int second_to_downstream_count;
        int second_count;
    };

    struct MXE_joined_count_strings
    {
        String_and_stream jc_inc_1;
        String_and_stream jc_skp_1;
        String_and_stream jc_inc_2;
        String_and_stream jc_skp_2;
        String_and_stream jcec_inc_1;
        String_and_stream jcec_skp_1;
        String_and_stream jcec_inc_2;
        String_and_stream jcec_skp_2;
        String_and_stream upstream_to_first;
        String_and_stream first_to_downstream;
        String_and_stream first;
        String_and_stream upstream_to_second;
        String_and_stream second_to_downstream;
        String_and_stream second;

        void clear()
        {
            jc_inc_1.clear();
            jc_skp_1.clear();
            jc_inc_2.clear();
            jc_skp_2.clear();
            jcec_inc_1.clear();
            jcec_skp_1.clear();
            jcec_inc_2.clear();
            jcec_skp_2.clear();
            upstream_to_first.clear();
            first_to_downstream.clear();
            first.clear();
            upstream_to_second.clear();
            second_to_downstream.clear();
            second.clear();
        }
    };

    struct MXE_counts_for_event_by_bam
    {
        char strand;
        Inc_skp_len jc_lengths;
        Inc_skp_len jcec_lengths;
        std::vector<MXE_counts_for_event> counts;

        void join_counts_across_bams(
            int sam1len,
            MXE_joined_count_strings* joined_strings) const
        {
            joined_strings->clear();
            std::string maybe_comma = "";
            std::string sam2_maybe_comma = "";
            for (int i = 0; i < counts.size(); ++i)
            {
                const MXE_counts_for_event& for_event = counts[i];
                if (i < sam1len)
                {
                    joined_strings->jc_inc_1.oss
                        << maybe_comma << for_event.jc_counts.inc_count;
                    joined_strings->jc_skp_1.oss
                        << maybe_comma << for_event.jc_counts.skp_count;
                    joined_strings->jcec_inc_1.oss
                        << maybe_comma << for_event.jcec_counts.inc_count;
                    joined_strings->jcec_skp_1.oss
                        << maybe_comma << for_event.jcec_counts.skp_count;
                }
                else
                {
                    joined_strings->jc_inc_2.oss
                        << sam2_maybe_comma << for_event.jc_counts.inc_count;
                    joined_strings->jc_skp_2.oss
                        << sam2_maybe_comma << for_event.jc_counts.skp_count;
                    joined_strings->jcec_inc_2.oss
                        << sam2_maybe_comma << for_event.jcec_counts.inc_count;
                    joined_strings->jcec_skp_2.oss
                        << sam2_maybe_comma << for_event.jcec_counts.skp_count;
                    sam2_maybe_comma = ",";
                }

                joined_strings->upstream_to_first.oss
                    << maybe_comma << for_event.upstream_to_first_count;
                joined_strings->first_to_downstream.oss
                    << maybe_comma << for_event.first_to_downstream_count;
                joined_strings->first.oss
                    << maybe_comma << for_event.first_count;
                joined_strings->upstream_to_second.oss
                    << maybe_comma << for_event.upstream_to_second_count;
                joined_strings->second_to_downstream.oss
                    << maybe_comma << for_event.second_to_downstream_count;
                joined_strings->second.oss
                    << maybe_comma << for_event.second_count;
                maybe_comma = ",";
            }
        }
    };

    struct ALT35_counts_for_event
    {
        ALT35_counts_for_event()
            : across_short_boundary_count(0),
              long_to_flanking_count(0),
              exclusive_to_long_count(0),
              short_to_flanking_count(0) {}
        Inc_skp_count jc_counts;
        Inc_skp_count jcec_counts;
        int across_short_boundary_count;
        int long_to_flanking_count;
        int exclusive_to_long_count;
        int short_to_flanking_count;
    };

    struct ALT35_joined_count_strings
    {
        String_and_stream jc_inc_1;
        String_and_stream jc_skp_1;
        String_and_stream jc_inc_2;
        String_and_stream jc_skp_2;
        String_and_stream jcec_inc_1;
        String_and_stream jcec_skp_1;
        String_and_stream jcec_inc_2;
        String_and_stream jcec_skp_2;
        String_and_stream across_short_boundary;
        String_and_stream long_to_flanking;
        String_and_stream exclusive_to_long;
        String_and_stream short_to_flanking;

        void clear()
        {
            jc_inc_1.clear();
            jc_skp_1.clear();
            jc_inc_2.clear();
            jc_skp_2.clear();
            jcec_inc_1.clear();
            jcec_skp_1.clear();
            jcec_inc_2.clear();
            jcec_skp_2.clear();
            across_short_boundary.clear();
            long_to_flanking.clear();
            exclusive_to_long.clear();
            short_to_flanking.clear();
        }
    };

    struct ALT35_counts_for_event_by_bam
    {
        char strand;
        Inc_skp_len jc_lengths;
        Inc_skp_len jcec_lengths;
        std::vector<ALT35_counts_for_event> counts;

        void join_counts_across_bams(
            int sam1len,
            ALT35_joined_count_strings* joined_strings) const
        {
            joined_strings->clear();
            std::string maybe_comma = "";
            std::string sam2_maybe_comma = "";
            for (int i = 0; i < counts.size(); ++i)
            {
                const ALT35_counts_for_event& for_event = counts[i];
                if (i < sam1len)
                {
                    joined_strings->jc_inc_1.oss
                        << maybe_comma << for_event.jc_counts.inc_count;
                    joined_strings->jc_skp_1.oss
                        << maybe_comma << for_event.jc_counts.skp_count;
                    joined_strings->jcec_inc_1.oss
                        << maybe_comma << for_event.jcec_counts.inc_count;
                    joined_strings->jcec_skp_1.oss
                        << maybe_comma << for_event.jcec_counts.skp_count;
                }
                else
                {
                    joined_strings->jc_inc_2.oss
                        << sam2_maybe_comma << for_event.jc_counts.inc_count;
                    joined_strings->jc_skp_2.oss
                        << sam2_maybe_comma << for_event.jc_counts.skp_count;
                    joined_strings->jcec_inc_2.oss
                        << sam2_maybe_comma << for_event.jcec_counts.inc_count;
                    joined_strings->jcec_skp_2.oss
                        << sam2_maybe_comma << for_event.jcec_counts.skp_count;
                    sam2_maybe_comma = ",";
                }

                joined_strings->across_short_boundary.oss
                    << maybe_comma << for_event.across_short_boundary_count;
                joined_strings->long_to_flanking.oss
                    << maybe_comma << for_event.long_to_flanking_count;
                joined_strings->exclusive_to_long.oss
                    << maybe_comma << for_event.exclusive_to_long_count;
                joined_strings->short_to_flanking.oss
                    << maybe_comma << for_event.short_to_flanking_count;
                maybe_comma = ",";
            }
        }
    };

    struct RI_counts_for_event
    {
        RI_counts_for_event()
            : upstream_to_intron_count(0),
              intron_to_downstream_count(0),
              intron_count(0),
              upstream_to_downstream_count(0) {}
        Inc_skp_count jc_counts;
        Inc_skp_count jcec_counts;
        int upstream_to_intron_count;
        int intron_to_downstream_count;
        int intron_count;
        int upstream_to_downstream_count;
    };

    struct RI_joined_count_strings
    {
        String_and_stream jc_inc_1;
        String_and_stream jc_skp_1;
        String_and_stream jc_inc_2;
        String_and_stream jc_skp_2;
        String_and_stream jcec_inc_1;
        String_and_stream jcec_skp_1;
        String_and_stream jcec_inc_2;
        String_and_stream jcec_skp_2;
        String_and_stream upstream_to_intron;
        String_and_stream intron_to_downstream;
        String_and_stream intron;
        String_and_stream upstream_to_downstream;

        void clear()
        {
            jc_inc_1.clear();
            jc_skp_1.clear();
            jc_inc_2.clear();
            jc_skp_2.clear();
            jcec_inc_1.clear();
            jcec_skp_1.clear();
            jcec_inc_2.clear();
            jcec_skp_2.clear();
            upstream_to_intron.clear();
            intron_to_downstream.clear();
            intron.clear();
            upstream_to_downstream.clear();
        }
    };

    struct RI_counts_for_event_by_bam
    {
        char strand;
        Inc_skp_len jc_lengths;
        Inc_skp_len jcec_lengths;
        std::vector<RI_counts_for_event> counts;

        void join_counts_across_bams(
            int sam1len,
            RI_joined_count_strings* joined_strings) const
        {
            joined_strings->clear();
            std::string maybe_comma = "";
            std::string sam2_maybe_comma = "";
            for (int i = 0; i < counts.size(); ++i)
            {
                const RI_counts_for_event& for_event = counts[i];
                if (i < sam1len)
                {
                    joined_strings->jc_inc_1.oss
                        << maybe_comma << for_event.jc_counts.inc_count;
                    joined_strings->jc_skp_1.oss
                        << maybe_comma << for_event.jc_counts.skp_count;
                    joined_strings->jcec_inc_1.oss
                        << maybe_comma << for_event.jcec_counts.inc_count;
                    joined_strings->jcec_skp_1.oss
                        << maybe_comma << for_event.jcec_counts.skp_count;
                }
                else
                {
                    joined_strings->jc_inc_2.oss
                        << sam2_maybe_comma << for_event.jc_counts.inc_count;
                    joined_strings->jc_skp_2.oss
                        << sam2_maybe_comma << for_event.jc_counts.skp_count;
                    joined_strings->jcec_inc_2.oss
                        << sam2_maybe_comma << for_event.jcec_counts.inc_count;
                    joined_strings->jcec_skp_2.oss
                        << sam2_maybe_comma << for_event.jcec_counts.skp_count;
                    sam2_maybe_comma = ",";
                }

                joined_strings->upstream_to_intron.oss
                    << maybe_comma << for_event.upstream_to_intron_count;
                joined_strings->intron_to_downstream.oss
                    << maybe_comma << for_event.intron_to_downstream_count;
                joined_strings->intron.oss
                    << maybe_comma << for_event.intron_count;
                joined_strings->upstream_to_downstream.oss
                    << maybe_comma << for_event.upstream_to_downstream_count;
                maybe_comma = ",";
            }
        }
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
