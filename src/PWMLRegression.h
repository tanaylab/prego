#ifndef seqpack_PWMLRegression_h
#define seqpack_PWMLRegression_h 1

/*=================================================
=================================================*/

#include "FunctionInterpolator.h"
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
    int m_step_num = 0;

    // response dimension
    int m_rdim;

    // response statistics in case of categorical response
    int m_ncat = 0;

    // epigenomic readouts to fit to (blocks of m_rdim in the vector)
    vector<float> m_interv_stat;
    vector<float> m_data_avg;
    vector<float> m_data_var;
    vector<vector<double> > m_data_avg_fold;
    vector<vector<double> > m_data_var_fold;

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
    int m_spat_bins_num;

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

    // number of folds for cross validation
    int m_num_folds;
    vector<int> m_folds;
    vector<int> m_fold_sizes;

    // random number for each sequence in order to break ties
    vector<float> m_data_epsilon;

    // aux for storing temp energy predictions
    vector<pair<float, int>> m_aux_preds;

    bool m_log_energy;

    // energy epsilon
    float m_energy_epsilon;

    // energy function callback
    FunctionInterpolator m_energy_func;

    bool m_optimize_pwm;
    bool m_optimize_spat;

  public:
    bool m_logit;

    float get_cur_score() const { return (m_cur_score); }

    // init in tiling mode
    PWMLRegression(const vector<string> &loci, const vector<int> &train_mask, int min_range,
                   int max_range, float min_prob, int spat_bin_size,
                   const vector<float> &resolutions, const vector<float> &s_resolutions, float eps,
                   float min_improv_for_star, float unif_prior, const string &score_metric,
                   const int &num_folds, const bool &log_energy, const float &energy_epsilon,
                   Rcpp::Nullable<Rcpp::Function> m_energy_func = R_NilValue, const float &xmin = -100,
                   const float &xmax = 100, const int &npts = 1000, const bool &optimize_pwm = true, const bool &optimize_spat = true);

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

    void compute_step_probs(const int &pos, const int& step, vector<float> &probs);
    tuple<int, int, float> choose_best_move();
    void take_best_step();
    void apply_move(const int &best_pos, const int &best_step, const float &best_score);
    void optimize_spatial_factors();
    pair<float, float> check_spat_bin(float best_spat_score, int spat_bin);
    void normalize_spat_factors(float best_spat_diff);

    float compute_cur_spat_score();
    float compute_cur_score(const int &pos, const vector<float> &probs);
    float compute_cur_fold_score(const int &pos, const vector<float> &probs, const int &fold);

    float compute_cur_r2(const int &pos, const vector<float> &probs);
    float compute_cur_r2_fold(const int &pos, const vector<float> &probs, const int &fold);
    float compute_cur_r2_spat();

    float compute_cur_ks(const int &pos, const vector<float> &probs);
    float compute_cur_ks_fold(const int &pos, const vector<float> &probs, const int &fold);
    float compute_cur_ks_spat();

    int pos_to_spat_bin(const int &pos);
    void symmetrize_spat_factors();
};

#endif // PWMLRegression_h
