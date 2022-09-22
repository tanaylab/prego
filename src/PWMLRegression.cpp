#include "port.h"
#include <chrono>
#include <cmath>

BASE_CC_FILE
#include "LeastSquare.h"
#include "PWMLRegression.h"

PWMLRegression::PWMLRegression(const vector<string> &seqs, const vector<int> &train_mask,
                               int min_range, int max_range, float min_prob, int spat_bin_size,
                               const vector<float> &resolutions, const vector<float> &s_resolutions,
                               float eps, float min_improv_for_star, float unif_prior,
                               const string &score_metric)
    : m_sequences(seqs), m_train_mask(train_mask), m_min_range(min_range), m_max_range(max_range),
      m_min_prob(min_prob), m_resolutions(resolutions), m_spat_resolutions(s_resolutions),
      m_spat_bin_size(spat_bin_size), // no spat bin for tiling,
      m_unif_prior(unif_prior), m_imporve_epsilon(eps), m_score_metric(score_metric) {}

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
    for (vector<int>::const_iterator mask = m_train_mask.begin(); mask != m_train_mask.end();
         mask++) {
        for (int rd = 0; rd < m_rdim; rd++) {
            *i_multi = stats[rd][seq_i];
            if (*mask) {
                m_data_avg[rd] += *i_multi;
                m_data_var[rd] += (*i_multi) * (*i_multi);
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
    if (m_logit) {
        Rcpp::Rcerr << "response avg " << m_data_avg[0] << " response var " << m_data_var[0]
                    << endl;
        Rcpp::Rcerr << "ranges " << m_min_range << " " << m_max_range << " bin " << m_spat_bin_size
                    << endl;
    }
}

void PWMLRegression::init_seed(const string &init_mot, int isbid) {
    m_nuc_factors.resize(init_mot.size(), vector<float>('T' + 1));
    m_spat_factors.resize((m_max_range - m_min_range) / m_spat_bin_size + 1);

    m_is_wildcard.resize(init_mot.size(), false);
    m_bidirect = isbid;

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
    for (int seq_id = 0; seq_id < m_sequences.size(); seq_id++) {
        m_derivs[seq_id].resize(pos, vector<float>('T' + 1));
    }
    int max_spat_bin = m_spat_factors.size();
    fill(m_spat_factors.begin(), m_spat_factors.end(), 1.0 / max_spat_bin);
    m_spat_derivs.resize(m_sequences.size(), vector<float>(max_spat_bin));
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
    Rcpp::Rcerr << "init pwm with size " << pwm.size() << endl;
    m_nuc_factors.resize(pwm.size(), vector<float>('T' + 1));
    m_spat_factors.resize((m_max_range - m_min_range) / m_spat_bin_size + 1);
    m_is_wildcard.resize(pwm.size(), false);
    m_bidirect = pwm.is_bidirect();

    for (int i = 0; i < pwm.size(); i++) {
        m_is_wildcard[i] = false;
        m_nuc_factors[i]['A'] = pwm[i].get_prob('A');
        m_nuc_factors[i]['C'] = pwm[i].get_prob('C');
        m_nuc_factors[i]['G'] = pwm[i].get_prob('G');
        m_nuc_factors[i]['T'] = pwm[i].get_prob('T');
        Rcpp::Rcerr << "set pos " << i + 1 << " to " << m_nuc_factors[i]['A'] << " "
                    << m_nuc_factors[i]['C'] << " " << m_nuc_factors[i]['G'] << " "
                    << m_nuc_factors[i]['T'] << endl;
    }

    m_aux_upds.resize(pwm.size());

    m_derivs.resize(m_sequences.size(), vector<vector<float>>(pwm.size(), vector<float>('T' + 1)));
    m_derivs.resize(m_sequences.size());
    for (int seq_id = 0; seq_id < m_sequences.size(); seq_id++) {
        m_derivs[seq_id].resize(pwm.size());
        for (int pos = 0; pos < pwm.size(); pos++) {
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

    for (int phase = 0; phase < m_resolutions.size(); phase++) {
        if (m_logit) {
            Rcpp::Rcerr << "Start phase " << phase << " resol " << m_resolutions[phase] << endl;
        }
        init_neighborhood(m_resolutions[phase]); // prepare all possible steps
        m_spat_factor_step = m_spat_resolutions[phase];
        do {
            prev_score = m_cur_score;
            // auto start = chrono::high_resolution_clock::now();
            init_energies();
            // auto stop = chrono::high_resolution_clock::now();
            // auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
            // Rcpp::Rcout << "init_energies took " << duration.count() << " milliseconds" << endl;

            // start = chrono::high_resolution_clock::now();
            take_best_step();
            // stop = chrono::high_resolution_clock::now();
            // duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
            // Rcpp::Rcout << "take_best_step took " << duration.count() << " milliseconds" << endl;

            if (m_logit) {
                Rcpp::Rcerr << "S -lrtEP prev " << prev_score << " " << m_cur_score << endl;
            }
        } while (m_cur_score > prev_score + m_imporve_epsilon);
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

void PWMLRegression::update_seq_interval(int seq_id, string::const_iterator min_i,
                                         string::const_iterator max_i, int sign, int pos) {
    for (string::const_iterator i = min_i; i < max_i; i++) {
        int spat_bin = int(pos / m_spat_bin_size);
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
}

void PWMLRegression::take_best_step() {

    int best_pos = -1;
    int best_step = -1;
    float best_score = m_cur_score;

    int max_pos = m_nuc_factors.size();
    int neigh_size = m_cur_neigh.size();
    vector<tuple<float, int, int>> pos_output(max_pos * neigh_size);

    // create a vector of positions and steps
    vector<pair<int, int>> pos_steps;
    pos_steps.reserve(max_pos * neigh_size);
    for (int pos = 0; pos < max_pos; pos++) {
        for (int step = 0; step < neigh_size; step++) {
            pos_steps.push_back(make_pair(pos, step));
        }
    }

    // score a single step
    auto score_step = [&](const pair<int, int> &pos_step) {
        int pos = pos_step.first;
        int step = pos_step.second;
        vector<float> probs = m_nuc_factors[pos];
        for (auto delta = m_cur_neigh[step].begin(); delta != m_cur_neigh[step].end(); delta++) {
            probs[delta->nuc] += delta->diff;
            if (probs[delta->nuc] <= 0) {
                probs[delta->nuc] = m_min_prob;
            }
        }
        float score = compute_cur_score(pos, probs);
        if (m_cur_score > score) {
            score = m_cur_score;
        }
        return make_tuple(score, pos, step);
    };

    // choose best step
    transform(pos_steps.begin(), pos_steps.end(), pos_output.begin(), score_step);

    sort(pos_output.begin(), pos_output.end(), greater<tuple<float, int, int>>());
    best_score = get<0>(pos_output[0]);
    best_pos = get<1>(pos_output[0]);
    best_step = get<2>(pos_output[0]);
    Rcpp::checkUserInterrupt();

    if (m_logit) {
        Rcpp::Rcerr << "best step was " << best_step << " pos " << best_pos << " score "
                    << best_score << endl;
    }

    int max_spat_bin = m_spat_factors.size();
    float best_spat_score = m_cur_score;
    int best_spat_bin = -1;
    float best_spat_diff = 0;
    float spat_score = compute_cur_spat_score();
    if (m_logit) {
        Rcpp::Rcerr << "spat score without change = " << spat_score << " cur is " << m_cur_score
                    << endl;
    }
    for (int spat_bin = 0; spat_bin < max_spat_bin; spat_bin++) {
        m_spat_factors[spat_bin] += m_spat_factor_step;
        float spat_score = compute_cur_spat_score();
        // Rcpp::Rcerr << "bin " << spat_bin << " add spat score " << spat_score << endl;
        if (spat_score > best_spat_score) {
            best_spat_score = spat_score;
            best_spat_bin = spat_bin;
            best_spat_diff = m_spat_factor_step;
        }
        m_spat_factors[spat_bin] -= 2 * m_spat_factor_step;
        if (m_spat_factors[spat_bin] >= 0) {
            spat_score = compute_cur_spat_score();
            // Rcpp::Rcerr << "bin " << spat_bin << " del spat score " << spat_score << endl;
            if (spat_score > best_spat_score) {
                best_spat_score = spat_score;
                best_spat_bin = spat_bin;
                best_spat_diff = -m_spat_factor_step;
            }
        }
        m_spat_factors[spat_bin] += m_spat_factor_step;
    }
    if (m_logit) {
        Rcpp::Rcerr << "best spat step was bin " << best_spat_bin << " diff " << best_spat_diff
                    << " score " << best_spat_score << endl;
    }
    if (best_score == m_cur_score) {
        if (m_logit) {
            Rcpp::Rcerr << "no improvement" << endl;
        }
        return;
    }

    if (best_score > best_spat_score) {
        vector<float> &pos_probs = m_nuc_factors[best_pos];
        for (vector<NeighStep>::iterator delta = m_cur_neigh[best_step].begin();
             delta != m_cur_neigh[best_step].end(); delta++) {
            pos_probs[delta->nuc] += delta->diff;
            if (pos_probs[delta->nuc] <= 0) {
                pos_probs[delta->nuc] = m_min_prob;
                // normalize
            }
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
    } else {
        if (m_logit) {
            Rcpp::Rcerr << "update spat bin " << best_spat_bin << " diff " << best_spat_diff
                        << endl;
        }

        m_spat_factors[best_spat_bin] += best_spat_diff;

        // normalize
        for (int bin = 0; bin < m_spat_factors.size(); bin++) {
            m_spat_factors[bin] /= (1 + best_spat_diff);
        }
        m_cur_score = best_spat_score;
    }
    // if (m_logit) {
    //     Rcpp::Rcerr << "after step";
    //     report_cur_lpwm();
    //     Rcpp::Rcerr << endl;
    // }
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
            assert(*resp == 0 || *resp == 1);
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
    for (int i = 0; i < m_aux_preds.size(); i++) {
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
    vector<double> xy(m_rdim, 0);

    double ex = 0;
    double ex2 = 0;
    int max_seq_id = m_sequences.size();

    vector<vector<vector<float>>>::const_iterator seq_deriv = m_derivs.begin();
    vector<float>::const_iterator resp = m_interv_stat.begin();
    for (int seq_id = 0; seq_id < max_seq_id; seq_id++) {
        if (m_train_mask[seq_id]) {
            const vector<float> &deriv = (*seq_deriv)[pos];
            float v = probs['A'] * deriv['A'] + probs['C'] * deriv['C'] + probs['G'] * deriv['G'] +
                      probs['T'] * deriv['T'];
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
            assert(*resp == 0 || *resp == 1);
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
    for (int i = 0; i < m_aux_preds.size(); i++) {
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
    vector<double> xy(m_rdim, 0);
    float ex = 0;
    float ex2 = 0;
    int max_seq_id = m_sequences.size();

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
    for (int sbin = 0; sbin < m_spat_factors.size(); sbin++) {
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
    for (int i = 0; i < m_nuc_factors.size(); i++) {
        model.get_pssm()[i].set_weight('A', m_nuc_factors[i]['A']);
        model.get_pssm()[i].set_weight('C', m_nuc_factors[i]['C']);
        model.get_pssm()[i].set_weight('G', m_nuc_factors[i]['G']);
        model.get_pssm()[i].set_weight('T', m_nuc_factors[i]['T']);
    }
    model.get_pssm().normalize();
    model.set_spat_bin_size(m_spat_bin_size);
    for (int bin = 0; bin < m_spat_factors.size(); bin++) {
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
    for (int seq_id = 0; seq_id < m_interv_stat.size(); seq_id++) {
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
