#include "port.h"

#include "KMerMultiStat.h"

KMerMultiStat::KMerMultiStat(int k, int max_degen, int gap_min_size, int gap_max_size,
                             vector<string> *sequences, vector<int> *is_train, int bin_num,
                             int norm_size, float norm_factor,
                             const vector<vector<float>> &interv_stat, int range_min, int range_max,
                             int gap_marg, const bool &logit, const set<string> *pat_filter)
    : m_k(k), m_min_gap_size(gap_min_size), m_max_gap_size(gap_max_size), m_sequences(sequences),
      m_is_train(is_train), m_resp_dim(interv_stat.size()), m_max_degen(max_degen),
      m_norm_size(norm_size), m_norm_factor(norm_factor), m_should_filter(pat_filter != 0),
      m_pat_filter(pat_filter), m_logit(logit) {
    m_tmp_bv.set_size(sequences->size());
    m_max_multi = bin_num;
    m_gap_margin = gap_marg;
    init_flat_stat(interv_stat);
    build_kmers(range_min, range_max);
}

void KMerMultiStat::init(int k, int max_degen, int gap_min_size, int gap_max_size,
                         vector<string> *sequences, vector<int> *is_train, int norm_size,
                         float norm_factor, const vector<vector<float>> &interv_stat, int range_min,
                         int range_max, int gap_marg) {
    m_kmer_multi_stat.clear();
    m_k = k;
    m_max_degen = max_degen;
    m_min_gap_size = gap_min_size;
    m_max_gap_size = gap_max_size;
    m_sequences = sequences;
    m_is_train = is_train;
    m_tmp_bv.set_size(m_sequences->size());
    m_max_multi = 10;
    m_norm_size = norm_size;
    m_norm_factor = norm_factor;
    m_gap_margin = gap_marg;
    m_resp_dim = interv_stat.size();
    init_flat_stat(interv_stat);
    build_kmers(range_min, range_max);
}

void KMerMultiStat::init_flat_stat(const vector<vector<float>> &stats) {
    m_interv_flat_stat.resize(stats.size() * stats[0].size());

    if (m_is_train->size() != stats[0].size()) {
        Rcpp::Rcerr
            << "ERROR: mismatch sizes between train mask and response stat when adding response "
               "to PWML regression"
            << endl;
        return;
    }
    vector<float>::iterator i_multi = m_interv_flat_stat.begin();

    for (size_t seq_i = 0; seq_i < m_is_train->size(); seq_i++) {
        for (int rd = 0; rd < m_resp_dim; rd++) {
            *i_multi = stats[rd][seq_i];
            ++i_multi;
        }
    }
}

void KMerMultiStat::build_kmers(int range_min, int range_max) {

    if (range_max == -1) {
        range_max = (*m_sequences)[0].size();
    }

    if (m_logit) {
        Rcpp::Rcerr << "Will build kmer multi stat"
                    << " interv stats " << m_interv_flat_stat.size() / m_resp_dim << " max gap "
                    << m_max_gap_size << endl;
    }

    vector<float>::iterator stat = m_interv_flat_stat.begin();

    for (uint locid = 0; locid < m_sequences->size(); locid++) {
        if (!(*m_is_train)[locid]) {
            stat += m_resp_dim;
            continue;
        }
        m_cur_locid = locid;
        if (m_logit && locid % 20 == 0) {
            Rcpp::Rcerr << "Processing id " << locid << "\r";
        }
        string::const_iterator seq = (*m_sequences)[locid].begin() + range_min;
        string::const_iterator seq_end =
            (*m_sequences)[locid].begin() + range_max - m_k - m_max_gap_size;
        int cur_length = seq_end - seq;
        m_kmer_counter.clear();
        for (string::const_iterator j = seq; j <= seq_end; ++j) {
            m_cur_pat.resize(m_k + 1);
            for (int g = m_min_gap_size; g <= m_max_gap_size; g++) {
                if (g == 0) {
                    m_cur_pat.resize(m_k);
                    fill_pat(j, m_k, 0, 0);
                    if (m_should_filter && m_pat_filter->count(m_cur_pat) == 0) {
                        continue;
                    }
                    if (!m_kmer_counter.count(m_cur_pat)) {
                        m_kmer_counter[m_cur_pat] = 1;
                    } else {
                        m_kmer_counter[m_cur_pat]++;
                    }
                    m_cur_pat.resize(m_k + 1);
                } else {
                    for (int k1 = m_gap_margin; k1 <= m_k - m_gap_margin; k1++) {
                        fill_pat(j, k1, m_k - k1, g);
                        if (!m_kmer_counter.count(m_cur_pat)) {
                            m_kmer_counter[m_cur_pat] = 1;
                        } else {
                            m_kmer_counter[m_cur_pat]++;
                        }
                    }
                }
            }
        }
        // update stat
        for (unordered_map<string, int>::iterator kmer = m_kmer_counter.begin();
             kmer != m_kmer_counter.end(); kmer++) {
            int multi;
            if (m_norm_size) {
                multi = int(kmer->second * m_norm_factor / float(cur_length));
            } else {
                multi = kmer->second;
            }
            if (multi >= m_max_multi) {
                multi = m_max_multi - 1;
            }
            if (m_kmer_multi_stat.count(kmer->first) == 0) {
                m_kmer_multi_stat[kmer->first].resize(m_max_multi,
                                                      pair<int, vector<float>>(0, m_resp_dim));
            }
            pair<int, vector<float>> &p = m_kmer_multi_stat[kmer->first][multi];
            p.first++;
            for (int r = 0; r < m_resp_dim; r++) {
                p.second[r] += *(stat + r);
            }
        }
        stat += m_resp_dim;
    }
}

void KMerMultiStat::fill_pat(string::const_iterator seq, int k1, int k2, int g) {
    int i = 0;
    string::iterator cp = m_cur_pat.begin();
    for (; i < k1; i++) {
        *cp = *(seq + i);
        cp++;
    }
    if (g == 0) {
        return;
    }
    i += g;
    *cp = '0' + g;
    cp++;
    int max_i = k1 + g + k2;
    for (; i < max_i; i++) {
        *cp = *(seq + i);
        cp++;
    }
}
