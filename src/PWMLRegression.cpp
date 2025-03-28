#include "port.h"
#include <chrono>
#include <cmath>
#include <random>
#include <algorithm>

#include "LeastSquare.h"
#include "PWMLRegression.h"

PWMLRegression::PWMLRegression(const vector<string> &seqs, const vector<int> &train_mask,
                               int min_range, int max_range, float min_prob, int spat_bin_size,
                               const vector<float> &resolutions, const vector<float> &s_resolutions,
                               float eps, float min_improv_for_star, float unif_prior,
                               const string &score_metric, const int &num_folds, const bool &log_energy, 
                               const float &energy_epsilon, Rcpp::Nullable<Rcpp::Function> energy_func, 
                               const float &xmin, const float &xmax, const int &npts, const bool &optimize_pwm, 
                               const bool &optimize_spat, const bool &symmetrize_spat)
    : m_sequences(seqs), m_train_mask(train_mask), m_min_range(min_range), m_max_range(max_range),
      m_min_prob(min_prob), m_resolutions(resolutions), m_spat_resolutions(s_resolutions),
      m_spat_bin_size(spat_bin_size),
      m_unif_prior(unif_prior), m_imporve_epsilon(eps), m_score_metric(score_metric), m_num_folds(num_folds),  
      m_log_energy(log_energy), m_energy_epsilon(energy_epsilon), m_optimize_pwm(optimize_pwm), 
      m_optimize_spat(optimize_spat), m_symmetrize_spat(symmetrize_spat) {
    if (m_num_folds < 1) {
        Rcpp::stop("number of folds must be at least 1");
    } 
    // generate an indicator vector the length of the number of sequences
    m_folds.resize(m_sequences.size());
    m_fold_sizes.resize(m_num_folds, 0);
    for (size_t i = 0; i < m_sequences.size(); i++) {
        m_folds[i] = i % m_num_folds;
        if (m_train_mask[i]){
            m_fold_sizes[m_folds[i]]++;
        }
    }
    if (m_num_folds > 1){
        Rcpp::RNGScope rngScope;
        

        std::random_device rd;
        std::mt19937 gen(rd());        
        std::shuffle(m_folds.begin(), m_folds.end(), gen);
    }

    if (energy_func.isNotNull()) {
        m_energy_func.init(Rcpp::as<Rcpp::Function>(energy_func), xmin, xmax, npts);
    }
}

void PWMLRegression::add_responses(const vector<vector<float>> &stats) {
    m_rdim = stats.size();
    if (m_score_metric == "ks" && m_rdim > 1) {
        Rcpp::stop("Warning: KS test is only for binary response");
    }
    if (m_logit) {
        Rcpp::Rcerr << "will init response "
                    << " m_rdim " << m_rdim << " size of vec " << stats[0].size() << endl;
    }

    m_data_avg.resize(m_rdim, 0);
    m_data_var.resize(m_rdim, 0);
    
    m_data_avg_fold.resize(m_num_folds);
    m_data_var_fold.resize(m_num_folds);
    for (int i = 0; i < m_num_folds; i++) {
        m_data_avg_fold[i].resize(m_rdim, 0);
        m_data_var_fold[i].resize(m_rdim, 0);
    }    

    m_interv_stat.resize(stats.size() * stats[0].size());

    if (m_train_mask.size() != stats[0].size()) {
        Rcpp::Rcerr
            << "mismatch sizes between train mask and response stat when adding response to PWML "
               "regression"
            << endl;
        return;
    }
    m_train_n = 0;
    vector<float>::iterator i_multi = m_interv_stat.begin();
    int seq_i = 0;
    int fold = 0;
    for (vector<int>::const_iterator mask = m_train_mask.begin(); mask != m_train_mask.end();
         mask++) {
        for (int rd = 0; rd < m_rdim; rd++) {
            *i_multi = stats[rd][seq_i];
            fold = m_folds[seq_i];
            if (*mask) {
                m_data_avg[rd] += *i_multi;
                m_data_var[rd] += (*i_multi) * (*i_multi);
                m_data_avg_fold[fold][rd] += *i_multi;
                m_data_var_fold[fold][rd] += (*i_multi) * (*i_multi);
            }
            ++i_multi;
        }
        if (*mask) {
            m_train_n++;
            if (m_score_metric == "ks") {
                if (stats[0][seq_i] == 1) {
                    m_ncat++;
                }
            }
        }
        seq_i++;
    }

    if (m_score_metric == "ks") {
        Rcpp::RNGScope scope; 
        Rcpp::NumericVector x = Rcpp::runif(m_sequences.size());
        m_data_epsilon.resize(m_sequences.size());
        for (size_t i = 0; i < m_data_epsilon.size(); i++) {
            m_data_epsilon[i] = x[i] * 1e-5;
        }
        m_aux_preds.reserve(m_train_n);
    }

    for (int rd = 0; rd < m_rdim; rd++) {
        m_data_avg[rd] /= m_train_n;
        m_data_var[rd] /= m_train_n;
        m_data_var[rd] -= m_data_avg[rd] * m_data_avg[rd];
    }

    for (int f = 0; f < m_num_folds; f++){
        for (int rd = 0; rd < m_rdim; rd++) {
            m_data_avg_fold[f][rd] /= m_fold_sizes[f];
            m_data_var_fold[f][rd] /= m_fold_sizes[f];
            m_data_var_fold[f][rd] -= m_data_avg_fold[f][rd] * m_data_avg_fold[f][rd];
        }
    }
    if (m_logit) {
        Rcpp::Rcerr << "response avg " << m_data_avg[0] << " response var " << m_data_var[0]
                    << endl;
        Rcpp::Rcerr << "ranges " << m_min_range << " " << m_max_range << " bin " << m_spat_bin_size
                    << endl;
    }
}

void PWMLRegression::init_seed(const string &init_mot, int isbid) {
    m_nuc_factors.resize(init_mot.size(), vector<float>('T' + 1));
    m_is_wildcard.resize(init_mot.size(), false);
    m_bidirect = isbid;

    m_spat_bins_num = (m_max_range - m_min_range) / m_spat_bin_size;    

    // abort if the number of bins is even and bidirect
    if (m_bidirect && m_spat_bins_num % 2 == 0) {
        Rcpp::stop("number of spatial bins must be odd when bidirect is true");
    }

    m_spat_factors.resize(m_spat_bins_num);

    int pos = 0;
    for (string::const_iterator i = init_mot.begin(); i != init_mot.end(); i++) {
        if (*i == '*') {
            fill(m_nuc_factors[pos].begin(), m_nuc_factors[pos].end(), 0.25);
            m_is_wildcard[pos] = true;
        } else {
            fill(m_nuc_factors[pos].begin(), m_nuc_factors[pos].end(), m_unif_prior);
            m_nuc_factors[pos][*i] = 1 - m_unif_prior * 3;
            m_is_wildcard[pos] = false;
        }
        pos++;
    }

    m_aux_upds.resize(pos);

    m_derivs.resize(m_sequences.size());
    for (size_t seq_id = 0; seq_id < m_sequences.size(); seq_id++) {
        m_derivs[seq_id].resize(pos, vector<float>('T' + 1));
    }
    int max_spat_bin = m_spat_factors.size();
    fill(m_spat_factors.begin(), m_spat_factors.end(), 1.0 / max_spat_bin);
    m_spat_derivs.resize(m_sequences.size(), vector<float>(max_spat_bin));
}

void PWMLRegression::init_seed_spat(const string &init_mot, const vector<float>& spat_factors, int isbid){
    init_seed(init_mot, isbid);
    m_spat_factors.resize(spat_factors.size());
    copy(spat_factors.begin(), spat_factors.end(), m_spat_factors.begin());
    int max_spat_bin = m_spat_factors.size();    
    m_spat_derivs.resize(m_sequences.size(), vector<float>(max_spat_bin));
    if (m_logit) {
        Rcpp::Rcerr << "init pwm spat: " << endl;    
        for (int bin = 0; bin < max_spat_bin; bin++) {
            Rcpp::Rcerr << m_spat_bin_size * bin << "\t" << m_spat_factors[bin] << endl;
        }
    }
}

void PWMLRegression::init_pwm_spat(DnaPSSM &pwm, const vector<float>& spat_factors){
    init_pwm(pwm);
    m_spat_factors.resize(spat_factors.size());
    copy(spat_factors.begin(), spat_factors.end(), m_spat_factors.begin());
    int max_spat_bin = m_spat_factors.size();    
    m_spat_derivs.resize(m_sequences.size(), vector<float>(max_spat_bin));
    if (m_logit) {
        Rcpp::Rcerr << "init pwm spat: " << endl;    
        for (int bin = 0; bin < max_spat_bin; bin++) {
            Rcpp::Rcerr << m_spat_bin_size * bin << "\t" << m_spat_factors[bin] << endl;
        }
    }
}

void PWMLRegression::init_pwm(DnaPSSM &pwm) {    
    m_nuc_factors.resize(pwm.size(), vector<float>('T' + 1));    

    m_spat_bins_num = (m_max_range - m_min_range) / m_spat_bin_size;    
    // abort if the number of bins is even and bidirect
    if (m_bidirect && m_spat_bins_num % 2 == 0) {
        Rcpp::stop("number of spatial bins must be odd when bidirect is true");
    }

    m_spat_factors.resize(m_spat_bins_num);

    m_is_wildcard.resize(pwm.size(), false);    
    m_bidirect = pwm.is_bidirect();    

    for (size_t i = 0; i < pwm.size(); i++) {
        m_is_wildcard[i] = false;
        m_nuc_factors[i]['A'] = pwm[i].get_prob('A');
        m_nuc_factors[i]['C'] = pwm[i].get_prob('C');
        m_nuc_factors[i]['G'] = pwm[i].get_prob('G');
        m_nuc_factors[i]['T'] = pwm[i].get_prob('T');
        // Rcpp::Rcerr << "set pos " << i + 1 << " to " << m_nuc_factors[i]['A'] << " "
        //             << m_nuc_factors[i]['C'] << " " << m_nuc_factors[i]['G'] << " "
        //             << m_nuc_factors[i]['T'] << endl;
    }

    m_aux_upds.resize(pwm.size());

    m_derivs.resize(m_sequences.size(), vector<vector<float>>(pwm.size(), vector<float>('T' + 1)));
    m_derivs.resize(m_sequences.size());
    for (size_t seq_id = 0; seq_id < m_sequences.size(); seq_id++) {
        m_derivs[seq_id].resize(pwm.size());
        for (size_t pos = 0; pos < pwm.size(); pos++) {
            m_derivs[seq_id][pos].resize('T' + 1);
        }
    }
    int max_spat_bin = m_spat_factors.size();
    fill(m_spat_factors.begin(), m_spat_factors.end(), 1.0 / max_spat_bin);
    m_spat_derivs.resize(m_sequences.size(), vector<float>(max_spat_bin));
}

void PWMLRegression::init_neighborhood(float resolution) {
    // ugly, but faster to program the the full recurrence

    m_cur_neigh.resize(0);
    m_cur_neigh.resize(20);

    m_cur_neigh[0].push_back(NeighStep('A', resolution));
    m_cur_neigh[1].push_back(NeighStep('A', -resolution));
    m_cur_neigh[2].push_back(NeighStep('C', resolution));
    m_cur_neigh[3].push_back(NeighStep('C', -resolution));
    m_cur_neigh[4].push_back(NeighStep('G', resolution));
    m_cur_neigh[5].push_back(NeighStep('G', -resolution));
    m_cur_neigh[6].push_back(NeighStep('T', resolution));
    m_cur_neigh[7].push_back(NeighStep('T', -resolution));

    m_cur_neigh[8].push_back(NeighStep('A', resolution));
    m_cur_neigh[8].push_back(NeighStep('C', resolution));

    m_cur_neigh[9].push_back(NeighStep('A', resolution));
    m_cur_neigh[9].push_back(NeighStep('G', resolution));

    m_cur_neigh[10].push_back(NeighStep('A', resolution));
    m_cur_neigh[10].push_back(NeighStep('T', resolution));

    m_cur_neigh[11].push_back(NeighStep('A', -resolution));
    m_cur_neigh[11].push_back(NeighStep('C', -resolution));

    m_cur_neigh[12].push_back(NeighStep('A', -resolution));
    m_cur_neigh[12].push_back(NeighStep('G', -resolution));

    m_cur_neigh[13].push_back(NeighStep('A', -resolution));
    m_cur_neigh[13].push_back(NeighStep('T', -resolution));

    m_cur_neigh[14].push_back(NeighStep('C', resolution));
    m_cur_neigh[14].push_back(NeighStep('G', resolution));

    m_cur_neigh[15].push_back(NeighStep('C', resolution));
    m_cur_neigh[15].push_back(NeighStep('T', resolution));

    m_cur_neigh[16].push_back(NeighStep('C', -resolution));
    m_cur_neigh[16].push_back(NeighStep('G', -resolution));

    m_cur_neigh[17].push_back(NeighStep('C', -resolution));
    m_cur_neigh[17].push_back(NeighStep('T', -resolution));

    m_cur_neigh[18].push_back(NeighStep('G', resolution));
    m_cur_neigh[18].push_back(NeighStep('T', resolution));

    m_cur_neigh[19].push_back(NeighStep('G', -resolution));
    m_cur_neigh[19].push_back(NeighStep('T', -resolution));
}

void PWMLRegression::optimize() {
    // init
    float prev_score = 0;
    m_cur_score = 0;

    for (size_t phase = 0; phase < m_resolutions.size(); phase++) {
        if (m_logit) {
            Rcpp::Rcerr << "Start phase " << phase << " resol " << m_resolutions[phase] << endl;
        }
        init_neighborhood(m_resolutions[phase]); // prepare all possible steps
        m_spat_factor_step = m_spat_resolutions[phase];
        do {
            prev_score = m_cur_score;

            Rcpp::checkUserInterrupt();
            // auto start = chrono::high_resolution_clock::now();
            init_energies();
            // auto stop = chrono::high_resolution_clock::now();
            // auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
            // Rcpp::Rcout << "init_energies took " << duration.count() << " milliseconds" << endl;
            Rcpp::checkUserInterrupt();

            // start = chrono::high_resolution_clock::now();
            take_best_step();
            // stop = chrono::high_resolution_clock::now();
            // duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
            // Rcpp::Rcout << "take_best_step took " << duration.count() << " milliseconds" << endl;

            if (m_logit) {
                Rcpp::Rcerr << "S -lrtEP prev " << prev_score << " " << m_cur_score << endl;
            }
            Rcpp::checkUserInterrupt();
            m_step_num++;
        } while (m_cur_score > prev_score + m_imporve_epsilon);

        if (m_symmetrize_spat){
            symmetrize_spat_factors();
        }        
    }
}

void PWMLRegression::init_energies() {
    // iter start and end
    // each time increase,decrease derivs

    if (m_logit) {
        Rcpp::Rcerr << "Will init til energ "
                    << " num of sequences " << m_sequences.size() << endl;
    }

    int max_pos = m_derivs[0].size();
    int max_seq_id = m_sequences.size();

    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (!m_train_mask[seq_id]) {
            continue;
        }
        string::const_iterator new_start = m_sequences[seq_id].begin();
        string::const_iterator new_end = m_sequences[seq_id].end();
        for (int pos = 0; pos < max_pos; pos++) {
            fill(m_derivs[seq_id][pos].begin(), m_derivs[seq_id][pos].end(), 0);
        }
        fill(m_spat_derivs[seq_id].begin(), m_spat_derivs[seq_id].end(), 0);
        update_seq_interval(seq_id, new_start, new_end, +1, 0);
    }

    if (m_logit) {
        Rcpp::Rcerr << "done init energies " << endl;
    }
}

int PWMLRegression::pos_to_spat_bin(const int& pos) {
    int bin = int(pos / m_spat_bin_size);
    if (m_symmetrize_spat && m_bidirect){        
        int center_bin = int(m_spat_bins_num / 2);

        // if bin is in the right side of the center symmetrize it
        if (bin > center_bin){
            bin = bin - (bin - center_bin) * 2;
        }
    }
    return bin;
}

void PWMLRegression::update_seq_interval(int seq_id, string::const_iterator min_i,
                                         string::const_iterator max_i, int sign, int pos) {
    for (string::const_iterator i = min_i; i < max_i; i++) {
        int spat_bin = pos_to_spat_bin(pos);        
        pos++;
        string::const_iterator j = i;
        double prod = sign;
        vector<UpdAux>::iterator upds = m_aux_upds.begin();
        vector<vector<float>>::iterator deriv = m_derivs[seq_id].begin();

        for (vector<vector<float>>::const_iterator p = m_nuc_factors.begin();
             p != m_nuc_factors.end(); p++) {
            if (!(*j) || *j == 'N' || *j == '*') {
                prod = 0;
                break;
            }
            float factor = (*p)[*j];
            upds->p = deriv->begin() + (*j);
            upds->factor = factor;
            prod *= factor;
            upds++;
            deriv++;
            j++;
        }
        if (!prod) {
            continue;
        }
        m_spat_derivs[seq_id][spat_bin] += prod;
        prod *= m_spat_factors[spat_bin];
        for (upds = m_aux_upds.begin(); upds != m_aux_upds.end(); upds++) {
            *(upds->p) += prod / upds->factor;
        }
        // update all derivs -

        if (m_bidirect) {
            float rprod = sign;
            j = i;
            vector<UpdAux>::iterator upds = m_aux_upds.begin();
            vector<vector<float>>::reverse_iterator deriv = m_derivs[seq_id].rbegin();
            for (vector<vector<float>>::reverse_iterator p = m_nuc_factors.rbegin();
                 p != m_nuc_factors.rend(); p++) {
                if (!(*j) || *j == 'N' || *j == '*') {
                    prod = 0;
                    break;
                }
                float factor = 0.25;
                switch (*j) {
                case 'A':
                    factor = (*p)['T'];
                    upds->p = deriv->begin() + 'T';
                    break;
                case 'T':
                    factor = (*p)['A'];
                    upds->p = deriv->begin() + 'A';
                    break;
                case 'C':
                    factor = (*p)['G'];
                    upds->p = deriv->begin() + 'G';
                    break;
                case 'G':
                    factor = (*p)['C'];
                    upds->p = deriv->begin() + 'C';
                    break;
                default:
                    break;
                }
                upds->factor = factor;
                rprod *= factor;
                upds++;
                deriv++;
                j++;
            }
            m_spat_derivs[seq_id][spat_bin] += rprod;
            rprod *= m_spat_factors[spat_bin];
            for (upds = m_aux_upds.begin(); upds != m_aux_upds.end(); upds++) {
                *(upds->p) += rprod / upds->factor;
            }
        }
    }
    if (m_symmetrize_spat){
        symmetrize_spat_factors();
    }
    
}

void PWMLRegression::symmetrize_spat_factors(){
    if (m_bidirect){
        int center_bin = int(m_spat_bins_num / 2);
        for (int bin = center_bin + 1; bin < m_spat_bins_num; bin++){
            m_spat_factors[bin] = m_spat_factors[bin - (bin - center_bin) * 2];
        }
    }
}

void PWMLRegression::compute_step_probs(const int &pos, const int &step, vector<float> &probs) {
    probs = m_nuc_factors[pos];
    for (auto delta = m_cur_neigh[step].begin(); delta != m_cur_neigh[step].end(); delta++) {
        probs[delta->nuc] += delta->diff;
        if (probs[delta->nuc] <= 0) {
            probs[delta->nuc] = m_min_prob;
        }
    }
}

tuple<int, int, float> PWMLRegression::choose_best_move() {
    int best_pos = 0;
    int best_step = 0;
    float best_score = m_cur_score;

    int max_pos = m_nuc_factors.size();
    int neigh_size = m_cur_neigh.size();
    vector<pair<int, int>> steps;
    vector<vector<float>> scores(m_num_folds, vector<float>(max_pos * neigh_size));
    vector<vector<int>> ranks(m_num_folds);
    vector<int> avg_ranks(max_pos * neigh_size, 0);
    
    int cur_step = 0;

    for (int pos = 0; pos < max_pos; pos++) {
        for (int step = 0; step < neigh_size; step++) {
            steps.emplace_back(pos, step);
            vector<float> probs;
            compute_step_probs(pos, step, probs);
            for (int i = 0; i < m_num_folds; i++) {
                scores[i][cur_step] = compute_cur_fold_score(pos, probs, i);
            }
            cur_step++;
        }
    }

    for (int i = 0; i < m_num_folds; i++) {
        rank_vector(scores[i], ranks[i]);        
    }

    int best_avg_rank = 0;
    for (int i = 0; i < max_pos * neigh_size; i++) {
        for (int j = 0; j < m_num_folds; j++) {
            avg_ranks[i] += ranks[j][i];            
        }
        avg_ranks[i] /= m_num_folds;
        if (avg_ranks[i] > best_avg_rank) {
            best_avg_rank = avg_ranks[i];
            best_pos = steps[i].first;
            best_step = steps[i].second;
            vector<float> probs;
            compute_step_probs(best_pos, best_step, probs);
            best_score = compute_cur_score(best_pos, probs);
        }
    }

    if (m_logit) {
        Rcpp::Rcerr << "best step was " << best_step << " pos " << best_pos << " score " << best_score << " best rank " << best_avg_rank << endl;
    }

    return make_tuple(best_pos, best_step, best_score);
}

void PWMLRegression::apply_move(const int &best_pos, const int &best_step, const float &best_score) {
    auto &pos_probs = m_nuc_factors[best_pos];
    for (auto &delta : m_cur_neigh[best_step]) {
        pos_probs[delta.nuc] += delta.diff;
        pos_probs[delta.nuc] = max(pos_probs[delta.nuc], m_min_prob);
    }

    float tot = pos_probs['A'] + pos_probs['C'] + pos_probs['G'] + pos_probs['T'];
    if (m_logit) {
        Rcpp::Rcerr << "update pos " << best_pos << " tot = " << tot << endl;
    }

    pos_probs['A'] /= tot;
    pos_probs['C'] /= tot;
    pos_probs['G'] /= tot;
    pos_probs['T'] /= tot;
    m_cur_score = best_score;
}

void PWMLRegression::take_best_step() {
    if (m_optimize_pwm) {
        auto [best_pos, best_step, best_score] = choose_best_move();
        if (best_score == m_cur_score) {
            if (m_logit) {
                Rcpp::Rcerr << "no improvement" << endl;
            }
            return;
        }
        apply_move(best_pos, best_step, best_score);
    }

    Rcpp::checkUserInterrupt();

    if (m_optimize_spat) {
        optimize_spatial_factors();
    }
}

void PWMLRegression::optimize_spatial_factors() {
    int max_spat_bin = m_spat_factors.size();
    float best_spat_score = m_cur_score;
    int best_spat_bin = -1;
    float best_spat_diff = 0;
    float spat_score = compute_cur_spat_score();
    if (m_logit) {
        Rcpp::Rcerr << "spat score without change = " << spat_score << " cur is " << m_cur_score << endl;
    }
    for (int spat_bin = 0; spat_bin < max_spat_bin; spat_bin++) {
        auto [new_best_score, new_best_diff] = check_spat_bin(best_spat_score, spat_bin);
        if(new_best_score > best_spat_score) {
            best_spat_score = new_best_score;
            best_spat_bin = spat_bin;
            best_spat_diff = new_best_diff;
        }
    }

    if (m_logit) {
        Rcpp::Rcerr << "best spat step was bin " << best_spat_bin << " diff " << best_spat_diff << " score " << best_spat_score << endl;
    }

    if (best_spat_score > m_cur_score) {
        if (m_logit) {
            Rcpp::Rcerr << "update spat bin " << best_spat_bin << " diff " << best_spat_diff << endl;
        }
        m_spat_factors[best_spat_bin] += best_spat_diff;
        normalize_spat_factors(best_spat_diff);
        m_cur_score = best_spat_score;
    }
}

pair<float, float> PWMLRegression::check_spat_bin(float best_spat_score, int spat_bin) {
    m_spat_factors[spat_bin] += m_spat_factor_step;
    float spat_score = compute_cur_spat_score();
    float best_spat_diff = m_spat_factor_step;
    if(spat_score > best_spat_score) {
        best_spat_score = spat_score;
    } else {
        m_spat_factors[spat_bin] -= 2 * m_spat_factor_step;
        if (m_spat_factors[spat_bin] >= 0) {
            spat_score = compute_cur_spat_score();
            if(spat_score > best_spat_score) {
                best_spat_score = spat_score;
                best_spat_diff = -m_spat_factor_step;
            }
        }
        m_spat_factors[spat_bin] += m_spat_factor_step;
    }
    return make_pair(best_spat_score, best_spat_diff);
}

void PWMLRegression::normalize_spat_factors(float best_spat_diff) {
    for(auto &bin: m_spat_factors) {
        bin /= (1 + best_spat_diff);
    }
}



float PWMLRegression::compute_cur_fold_score(const int &pos, const vector<float> &probs, const int &fold) {
    float score;
    if (m_score_metric == "r2") {
        score = compute_cur_r2_fold(pos, probs, fold);
    } else if (m_score_metric == "ks") {
        score = compute_cur_ks_fold(pos, probs, fold);
    } else {
        Rcpp::stop("Unknown score metric (can be either 'r2' or 'ks')");
    }
    return score;
}

float PWMLRegression::compute_cur_score(const int &pos, const vector<float> &probs) {
    float score;
    if (m_score_metric == "r2") {
        score = compute_cur_r2(pos, probs);
    } else if (m_score_metric == "ks") {
        score = compute_cur_ks(pos, probs);
    } else {
        Rcpp::stop("Unknown score metric (can be either 'r2' or 'ks')");
    }
    return score;
}

float PWMLRegression::compute_cur_ks(const int &pos, const vector<float> &probs) {
    int max_seq_id = m_sequences.size();
    m_aux_preds.resize(0);

    vector<vector<vector<float>>>::const_iterator seq_deriv = m_derivs.begin();
    vector<float>::const_iterator resp = m_interv_stat.begin();
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {
            const vector<float> &deriv = (*seq_deriv)[pos];
            float v = probs['A'] * deriv['A'] + probs['C'] * deriv['C'] + probs['G'] * deriv['G'] +
                      probs['T'] * deriv['T'];

            // push predictions according to category
            
            float epsilon = m_data_epsilon[seq_id];
            m_aux_preds.push_back(make_pair(-v * (1 + epsilon), *resp));
        }
        resp++;
        seq_deriv++;
    }

    sort(m_aux_preds.begin(), m_aux_preds.end());

    float max_diff = 0;
    float cur_diff = 0;
    float n_0 = m_train_n - m_ncat;
    float n_1 = m_ncat;
    for (size_t i = 0; i < m_aux_preds.size(); i++) {
        if (m_aux_preds[i].second == 0) {
            cur_diff -= 1 / n_0;
        } else {
            cur_diff += 1 / n_1;
        }
        if (cur_diff > max_diff) {
            max_diff = cur_diff;
        }
    }

    return max_diff;
}

float PWMLRegression::compute_cur_ks_fold(const int &pos, const vector<float> &probs, const int &fold) {
    int max_seq_id = m_sequences.size();
    m_aux_preds.resize(0);

    vector<vector<vector<float>>>::const_iterator seq_deriv = m_derivs.begin();
    vector<float>::const_iterator resp = m_interv_stat.begin();
    float n_0 = 0;
    float n_1 = 0;
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id] && m_folds[seq_id] == fold) {
            const vector<float> &deriv = (*seq_deriv)[pos];
            float v = probs['A'] * deriv['A'] + probs['C'] * deriv['C'] + probs['G'] * deriv['G'] +
                      probs['T'] * deriv['T'];

            // push predictions according to category

            float epsilon = m_data_epsilon[seq_id];
            m_aux_preds.push_back(make_pair(-v * (1 + epsilon), *resp));
            if (*resp == 0) {
                n_0++;
            } else {
                n_1++;
            }
        }
        resp++;
        seq_deriv++;
    }

    sort(m_aux_preds.begin(), m_aux_preds.end());

    float max_diff = 0;
    float cur_diff = 0;
    
    for (size_t i = 0; i < m_aux_preds.size(); i++) {
        if (m_aux_preds[i].second == 0) {
            cur_diff -= 1 / n_0;
        } else {
            cur_diff += 1 / n_1;
        }
        if (cur_diff > max_diff) {
            max_diff = cur_diff;
        }
    }

    return max_diff;
}

float PWMLRegression::compute_cur_r2(const int &pos, const vector<float> &probs) {
    int max_seq_id = m_sequences.size();

    // calculate energies
    vector<vector<vector<float>>>::const_iterator seq_deriv = m_derivs.begin();
    
    vector<float> energies(max_seq_id, 0);   
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {
            const vector<float> &deriv = (*seq_deriv)[pos];
            float v = probs['A'] * deriv['A'] + probs['C'] * deriv['C'] + probs['G'] * deriv['G'] +
                      probs['T'] * deriv['T'];
            if (m_log_energy) {
                v = log(v + m_energy_epsilon);
            }  
            energies[seq_id] = v;            
        }
        seq_deriv++;
    }

    if (m_energy_func.is_initialized()){           
        if (m_logit) {
            Rcpp::Rcerr << "energies before " << energies[0] << " " << energies[1] << " " << energies[2] << " " << energies[3] << " " << energies[4] << " " << energies[5] << endl;
        }
        
        energies = m_energy_func.interpolate(energies);
        // energies = Rcpp::as<vector<float>>(Rcpp::as<Rcpp::Function>(m_energy_func)(energies));
        if (m_logit) {
            Rcpp::Rcerr << "energies after " << energies[0] << " " << energies[1] << " " << energies[2] << " " << energies[3] << " " << energies[4] << " " << energies[5] << endl;
        }        
        

        if (energies.size() != (size_t)max_seq_id){
            Rcpp::stop("Energy function must return a vector of the same length as the number of sequences");
        }
    }


    // calculate statistics
    vector<double> xy(m_rdim, 0);
    double ex = 0;
    double ex2 = 0;
    seq_deriv = m_derivs.begin();
    vector<float>::const_iterator resp = m_interv_stat.begin();
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {            
            float v = energies[seq_id];
          
            ex += v;
            ex2 += v * v;
            for (int rd = 0; rd < m_rdim; rd++) {
                xy[rd] += v * *resp;
                resp++;
            }
        }
        seq_deriv++;
    }
    // if (m_logit) {
    //     Rcpp::Rcerr << "done initing xy for all sequences, xy 0 is " << xy[0] << endl;
    // }
    ex /= m_train_n;
    ex2 /= m_train_n;
    double pred_var = ex2 - ex * ex;
    float tot_r2 = 0;
    m_a.resize(m_rdim);
    m_b.resize(m_rdim);
    for (int rd = 0; rd < m_rdim; rd++) {
        xy[rd] /= m_train_n;
        float cov = xy[rd] - ex * m_data_avg[rd];
        float r2 = cov * cov / (pred_var * m_data_var[rd]);
        m_a[rd] = (ex2 * m_data_avg[rd] - ex * xy[rd]) / pred_var;
        m_b[rd] = cov / pred_var;
        tot_r2 += r2;
    }

    if (m_logit && std::isnan(tot_r2)) {
        Rcpp::Rcerr << "Nan at at r2 var " << pred_var << " " << ex << " " << ex2 << " "
                    << m_data_avg[0] << endl;
    }
    return (tot_r2);
}

float PWMLRegression::compute_cur_r2_fold(const int &pos, const vector<float> &probs, const int &fold) {
    int max_seq_id = m_sequences.size();

    // calculate energies
    vector<vector<vector<float>>>::const_iterator seq_deriv = m_derivs.begin();    
    vector<float> energies(max_seq_id, 0);   
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {
            const vector<float> &deriv = (*seq_deriv)[pos];
            float v = probs['A'] * deriv['A'] + probs['C'] * deriv['C'] + probs['G'] * deriv['G'] +
                      probs['T'] * deriv['T'];
            if (m_log_energy) {
                v = log(v + m_energy_epsilon);
            }  
            energies[seq_id] = v;            
        }
        seq_deriv++;
    }

    if (m_energy_func.is_initialized()) {
        energies = m_energy_func.interpolate(energies);
        // energies = Rcpp::as<vector<float>>(Rcpp::as<Rcpp::Function>(m_energy_func)(energies));
    }

    // calculate statistics
    vector<double> xy(m_rdim, 0);
    double ex = 0;
    double ex2 = 0;
    seq_deriv = m_derivs.begin();
    vector<float>::const_iterator resp = m_interv_stat.begin();
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {            
            float v = energies[seq_id];
            ex += v;
            ex2 += v * v;
            for (int rd = 0; rd < m_rdim; rd++) {                
                if (m_folds[seq_id] == fold) {
                    xy[rd] += v * *resp;
                }
                resp++;
            }
        }
        seq_deriv++;
    }
    // if (m_logit) {
    //     Rcpp::Rcerr << "done initing xy for all sequences, xy 0 is " << xy[0] << endl;
    // }
    ex /= m_fold_sizes[fold];
    ex2 /= m_fold_sizes[fold];
    double pred_var = ex2 - ex * ex;
    float tot_r2 = 0;
    m_a.resize(m_rdim);
    m_b.resize(m_rdim);
    for (int rd = 0; rd < m_rdim; rd++) {
        xy[rd] /= m_fold_sizes[fold];
        float cov = xy[rd] - ex * m_data_avg_fold[fold][rd];
        float r2 = cov * cov / (pred_var * m_data_var_fold[fold][rd]);
        m_a[rd] = (ex2 * m_data_avg_fold[fold][rd] - ex * xy[rd]) / pred_var;
        m_b[rd] = cov / pred_var;
        tot_r2 += r2;
    }

    if (m_logit && std::isnan(tot_r2)) {
        Rcpp::Rcerr << "Nan at at r2 var " << pred_var << " " << ex << " " << ex2 << " " << fold << " "  
                    << m_data_avg_fold[fold][0] << endl;
    }
    return (tot_r2);
}

float PWMLRegression::compute_cur_spat_score() {
    float score;
    if (m_score_metric == "r2") {
        score = compute_cur_r2_spat();
    } else if (m_score_metric == "ks") {
        score = compute_cur_ks_spat();
    } else {
        Rcpp::stop("Unknown score metric (can be either 'r2' or 'ks')");
    }
    return (score);
}

float PWMLRegression::compute_cur_ks_spat() {
    int max_seq_id = m_sequences.size();
    m_aux_preds.resize(0);

    vector<vector<float>>::const_iterator seq_derivs = m_spat_derivs.begin();
    vector<float>::const_iterator resp = m_interv_stat.begin();

    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {
            float v = 0;
            vector<float>::const_iterator fact = m_spat_factors.begin();
            for (vector<float>::const_iterator bin = seq_derivs->begin(); bin != seq_derivs->end();
                 bin++) {
                v += *bin * *fact;
                fact++;
            }

            // push predictions according to category
            
            float epsilon = m_data_epsilon[seq_id];
            m_aux_preds.push_back(make_pair(-v * (1 + epsilon), *resp));
            resp++;
        }
        seq_derivs++;
    }

    sort(m_aux_preds.begin(), m_aux_preds.end());

    float max_diff = 0;
    float cur_diff = 0;
    float n_0 = m_train_n - m_ncat;
    float n_1 = m_ncat;
    for (size_t i = 0; i < m_aux_preds.size(); i++) {
        if (m_aux_preds[i].second == 0) {
            cur_diff -= 1 / n_0;
        } else {
            cur_diff += 1 / n_1;
        }
        if (cur_diff > max_diff) {
            max_diff = cur_diff;
        }
    }

    return max_diff;
}

float PWMLRegression::compute_cur_r2_spat() {
    int max_seq_id = m_sequences.size();

    vector<vector<float>>::const_iterator seq_derivs = m_spat_derivs.begin();    
    vector<float> energies(max_seq_id, 0);   
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {
            float v = 0;
            vector<float>::const_iterator fact = m_spat_factors.begin();
            for (vector<float>::const_iterator bin = seq_derivs->begin(); bin != seq_derivs->end();
                 bin++) {
                v += *bin * *fact;
                fact++;
            }
            if (m_log_energy) {
                v = log(v + m_energy_epsilon);
            }
            energies[seq_id] = v;            
        }
        seq_derivs++;
    }

    if (m_energy_func.is_initialized()) {
        energies = m_energy_func.interpolate(energies);
        // energies = Rcpp::as<vector<float>>(Rcpp::as<Rcpp::Function>(m_energy_func)(energies));
    }

    vector<double> xy(m_rdim, 0);
    float ex = 0;
    float ex2 = 0;
    seq_derivs = m_spat_derivs.begin();
    vector<float>::const_iterator resp = m_interv_stat.begin();
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {
            float v = energies[seq_id];            
            ex += v;
            ex2 += v * v;
            for (int rd = 0; rd < m_rdim; rd++) {
                xy[rd] += v * *resp;
                resp++;
            }
        }
        seq_derivs++;
    }

    ex /= m_train_n;
    ex2 /= m_train_n;
    double pred_var = ex2 - ex * ex;
    float tot_r2 = 0;
    m_a.resize(m_rdim);
    m_b.resize(m_rdim);
    for (int rd = 0; rd < m_rdim; rd++) {
        xy[rd] /= m_train_n;
        float cov = xy[rd] - ex * m_data_avg[rd];
        float r2 = cov * cov / (pred_var * m_data_var[rd]);
        m_a[rd] = (ex2 * m_data_avg[rd] - ex * xy[rd]) / pred_var;
        m_b[rd] = cov / pred_var;
        tot_r2 += r2;
    }

    return (tot_r2);
}

void PWMLRegression::report_cur_lpwm() {
    int max_pos = m_nuc_factors.size();
    for (int pos = 0; pos < max_pos; pos++) {
        REprintf("(%.2f %.2f %.2f %.2f) ", m_nuc_factors[pos]['A'], m_nuc_factors[pos]['G'],
                 m_nuc_factors[pos]['C'], m_nuc_factors[pos]['T']);
    }
    Rcpp::Rcerr << " ";
    for (size_t sbin = 0; sbin < m_spat_factors.size(); sbin++) {
        REprintf("%.2f ", m_spat_factors[sbin]);
    }
}

Rcpp::DataFrame PWMLRegression::output_pssm_df(int psid) {
    int max_pos = m_nuc_factors.size();

    vector<int> poss(max_pos);
    vector<float> pssm_A(max_pos, 0);
    vector<float> pssm_C(max_pos, 0);
    vector<float> pssm_G(max_pos, 0);
    vector<float> pssm_T(max_pos, 0);

    for (int pos = 0; pos < max_pos; pos++) {
        poss[pos] = pos;
        pssm_A[pos] = m_nuc_factors[pos]['A'];
        pssm_C[pos] = m_nuc_factors[pos]['C'];
        pssm_G[pos] = m_nuc_factors[pos]['G'];
        pssm_T[pos] = m_nuc_factors[pos]['T'];
    }

    Rcpp::DataFrame pssm = Rcpp::DataFrame::create(
        Rcpp::Named("pos") = poss, Rcpp::Named("A") = pssm_A, Rcpp::Named("C") = pssm_C,
        Rcpp::Named("G") = pssm_G, Rcpp::Named("T") = pssm_T);

    return pssm;
}

Rcpp::DataFrame PWMLRegression::output_spat_df(int psid) {
    int max_spat_bin = m_spat_factors.size();
    vector<int> spat_bins(max_spat_bin);
    vector<float> spat_factors(max_spat_bin, 0);

    for (int bin = 0; bin < max_spat_bin; bin++) {
        spat_bins[bin] = m_spat_bin_size * bin;
        spat_factors[bin] = m_spat_factors[bin];
    }

    Rcpp::DataFrame spat_df = Rcpp::DataFrame::create(Rcpp::Named("bin") = spat_bins,
                                                      Rcpp::Named("spat_factor") = spat_factors);

    return (spat_df);
}

void PWMLRegression::output_pssm(ostream &out, ostream &spat, int psid) {
    Rcpp::Rcerr << "will out pwm" << endl;
    int max_pos = m_nuc_factors.size();
    for (int pos = 0; pos < max_pos; pos++) {
        out << psid << "\t" << pos;
        out << "\t" << m_nuc_factors[pos]['A'] << "\t" << m_nuc_factors[pos]['C'] << "\t"
            << m_nuc_factors[pos]['G'] << "\t" << m_nuc_factors[pos]['T'] << endl;
    }
    int max_spat_bin = m_spat_factors.size();
    for (int bin = 0; bin < max_spat_bin; bin++) {
        spat << psid << "\t" << m_spat_bin_size * bin << "\t" << m_spat_factors[bin] << endl;
    }
}

float PWMLRegression::get_r2() {
    init_energies();
    float r2 = compute_cur_r2_spat();
    return (r2);
}

void PWMLRegression::get_model(DnaPWML &model) {
    model.get_pssm().resize(m_nuc_factors.size());
    model.get_pssm().set_range(m_min_range, m_max_range);
    for (size_t i = 0; i < m_nuc_factors.size(); i++) {
        model.get_pssm()[i].set_weight('A', m_nuc_factors[i]['A']);
        model.get_pssm()[i].set_weight('C', m_nuc_factors[i]['C']);
        model.get_pssm()[i].set_weight('G', m_nuc_factors[i]['G']);
        model.get_pssm()[i].set_weight('T', m_nuc_factors[i]['T']);
    }
    model.get_pssm().normalize();
    model.set_spat_bin_size(m_spat_bin_size);
    for (size_t bin = 0; bin < m_spat_factors.size(); bin++) {
        model.set_spat_factor(bin, m_spat_factors[bin]);
    }
    model.get_pssm().set_bidirect(m_bidirect);
}

void PWMLRegression::fill_predictions(vector<float> &preds) {
    // updating the m_a, m_b
    init_energies();
    Rcpp::Rcerr << "done init energ" << endl;
    // float r2 = compute_cur_r2_spat();
    // Rcpp::Rcerr << "done comp r2 = " << r2 << endl;

    preds.resize(m_interv_stat.size());
    vector<float>::iterator prs = preds.begin();
    vector<vector<vector<float>>>::iterator seq_deriv = m_derivs.begin();
    vector<float>::iterator resp = m_interv_stat.begin();
    vector<float> &factors = m_nuc_factors[0];
    for (size_t seq_id = 0; seq_id < m_interv_stat.size(); seq_id++) {
        if (m_train_mask[seq_id]) {
            vector<float> &deriv = (*seq_deriv)[0];

            float predict = factors['A'] * deriv['A'] + factors['C'] * deriv['C'] +
                            factors['G'] * deriv['G'] + factors['T'] * deriv['T'];

            *prs = predict;
        } else {
            *prs = -1e+30;
        }
        prs++;
        seq_deriv++;
        resp++;
    }
}
