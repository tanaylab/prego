#ifndef seqpack_KMerMultiStat_h
#define seqpack_KMerMultiStat_h 1

/*=================================================
=================================================*/

#include "BitVec.h"
#include "options.h"
#include <Rcpp.h>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

class KMerMultiStat {

  protected:
    int m_k;

    int m_min_gap_size;
    int m_max_gap_size;
    int m_gap_margin;

    vector<string> *m_sequences;
    vector<int> *m_is_train;

    int m_resp_dim;
    vector<float> m_interv_flat_stat;

    int m_max_degen;

    int m_norm_size;
    float m_norm_factor;

    string m_cur_pat;

    int m_cur_locid;

    ds_bitvec m_tmp_bv;

    int m_max_multi;
    unordered_map<string, int> m_kmer_counter;
    map<const string, vector<pair<int, vector<float>>>> m_kmer_multi_stat;

    ds_bitvec m_null_fp;

    bool m_should_filter;

    const set<string> *m_pat_filter;

  public:
    bool m_logit;

    void set_gap_margin(int mg) { m_gap_margin = mg; }

    pair<int, vector<float>> get_kmer_multi_stat(const string &kmer, int multi) {
        if (multi >= m_max_multi || m_kmer_multi_stat.count(kmer) == 0) {
            return (pair<int, vector<float>>(0, m_resp_dim));
        }
        return (m_kmer_multi_stat[kmer][multi]);
    }

    map<const string, vector<pair<int, vector<float>>>>::const_iterator get_pat_begin() const {
        return (m_kmer_multi_stat.begin());
    }
    map<const string, vector<pair<int, vector<float>>>>::const_iterator get_pat_end() const {
        return (m_kmer_multi_stat.end());
    }

    int get_pat_size() const { return (m_kmer_multi_stat.size()); }

    KMerMultiStat() {}

    KMerMultiStat(int k, int degen, int min_gap, int max_gap, vector<string> *focus_loci,
                  vector<int> *is_train, int bin_num, int norm, float norm_factor,
                  const vector<vector<float>> &interv_stat, int range_min, int range_max,
                  int gap_marg = 2, const bool &logit = true, const set<string> *pat_filter = 0);

    void init(int k, int degen, int min_gap, int max_gap, vector<string> *sequences,
              vector<int> *is_train, int norm, float norm_factor,
              const vector<vector<float>> &interv_stat, int range_min, int range_max,
              int gap_marg = 2);

    void init_flat_stat(const vector<vector<float>> &interv_stat);

    void build_kmers(int range_min, int range_max);

    void fill_pat(string::const_iterator seq, int k1, int k2, int g);

    bool have_kmer(const string &pat) { return (m_kmer_multi_stat.count(pat)); }
};

#endif // KMerMultiStat
