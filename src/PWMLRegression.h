#ifndef seqpack_PWMLRegression_h
#define seqpack_PWMLRegression_h 1

/*=================================================
=================================================*/

#include "BitVec.h"
#include "DnaPWML.h"
#include <Rcpp.h>
#include <string.h>
#include <vector>

class PWMLRegression {

  private:
    struct UpdAux {
        vector<float>::iterator p;
        float factor;
    };

    struct NeighStep {
        char nuc;
        float diff;
        NeighStep(char n, float d) : nuc(n), diff(d) {}
    };

  protected:
    // Relevant chromosomes and probes (if in tiling mode)
    const vector<string> &m_sequences;
    const vector<int> &m_train_mask;
    int m_train_n = 0;

    // response dimension
    int m_rdim;

    // response statistics in case of categorical response
    int m_ncat = 0;

    // epigenomic readouts to fit to (blocks of m_rdim in the vector)
    vector<float> m_interv_stat;
    vector<double> m_data_avg;
    vector<double> m_data_var;

    int m_min_range;
    int m_max_range;

    vector<bool> m_is_wildcard;

    // m_bidirect[i] is mixture pssm i bidirectional?
    bool m_bidirect;

    // m_chars[mix_id][pos]['A'] = weight for A in position pos
    vector<vector<float>> m_nuc_factors;
    vector<vector<vector<float>>> m_derivs; // deriv[vid][pos][char]

    float m_min_prob;

    // step sizes at each phase
    vector<float> m_resolutions;
    vector<float> m_spat_resolutions;

    vector<float> m_spat_factors;
    vector<vector<float>> m_spat_derivs; // deriv[vid][pos][char]
    int m_spat_bin_size;
    float m_spat_factor_step;

    float m_cur_score;

    // cur best fit line m_a + m_b*pred(seq)
    vector<float> m_a;
    vector<float> m_b;

    float m_min_rms_for_star;

    float m_unif_prior;

    vector<UpdAux> m_aux_upds;

    // defining the local search policty
    vector<vector<NeighStep>> m_cur_neigh;

    // imporvement epsilon
    float m_imporve_epsilon;

    // score metric
    string m_score_metric;

    // random number for each sequence in order to break ties
    vector<float> m_data_epsilon;

    // aux for storing temp energy predictions
    vector<pair<float, int>> m_aux_preds;

  public:
    bool m_logit;

    float get_cur_score() const { return (m_cur_score); }

    // init in tiling mode
    PWMLRegression(const vector<string> &loci, const vector<int> &train_mask, int min_range,
                   int max_range, float min_prob, int spat_bin_size,
                   const vector<float> &resolutions, const vector<float> &s_resolutions, float eps,
                   float min_improv_for_star, float unif_prior, const string &score_metric);

    void add_responses(const vector<vector<float>> &stats);

    void init_seed(const string &init_mot, int isbid);
    void init_pwm(DnaPSSM &pwm);
    void init_pwm_spat(DnaPSSM &pwm, const vector<float> &spat_factors);

    void optimize();

    void output_pssm(ostream &out, ostream &spat, int psid);
    Rcpp::DataFrame output_pssm_df(int psid);
    Rcpp::DataFrame output_spat_df(int psid);

    void report_cur_lpwm();

    void init_energies();
    void fill_predictions(vector<float> &pred);
    void get_diff_vector(vector<float> &to_diff, float &avg, float &var);
    float get_r2();
    void get_model(DnaPWML &pwml);

  private:
    void init_neighborhood(float resolution);

    void update_seq_interval(int seq_id, string::const_iterator min_i, string::const_iterator max_i,
                             int sign, int pos);

    void take_best_step();

    float compute_cur_spat_score();
    float compute_cur_score(const int &pos, const vector<float> &probs);

    float compute_cur_r2(const int &pos, const vector<float> &probs);
    float compute_cur_r2_spat();

    float compute_cur_ks(const int &pos, const vector<float> &probs);
    float compute_cur_ks_spat();
};

#endif // PWMLRegression_h
