#include "port.h"
#include <cmath>
BASE_CC_FILE
#include "LeastSquare.h"
#include "PWMLRegression.h"

PWMLRegression::PWMLRegression(const vector<string> &seqs, const vector<int> &train_mask,
                               int min_range, int max_range, float min_prob, int spat_bin_size,
                               const vector<float> &resolutions, const vector<float> &s_resolutions,
                               float eps, float min_improv_for_star, float unif_prior)
    : m_sequences(seqs), m_train_mask(train_mask), m_min_range(min_range), m_max_range(max_range),
      m_min_prob(min_prob), m_spat_bin_size(spat_bin_size), // no spat bin for tiling
      m_resolutions(resolutions), m_spat_resolutions(s_resolutions), m_unif_prior(unif_prior), 
      m_imporve_epsilon(eps) {}

void PWMLRegression::add_responses(const vector<vector<float>> &stats) {
    m_rdim = stats.size();
    Rcpp::Rcerr << "will init response "
         << " m_rdim " << m_rdim << " size of vec " << stats[0].size() << endl;
    m_data_avg.resize(m_rdim, 0);
    m_data_var.resize(m_rdim, 0);
    m_interv_stat.resize(stats.size() * stats[0].size());

    if (m_train_mask.size() != stats[0].size()) {
        Rcpp::Rcerr << "mismatch sizes between train mask and response stat when adding response to PWML "
                "regression"
             << endl;
        return;
    }
    m_train_n = 0;
    vector<float>::iterator i_multi = m_interv_stat.begin();
    int prb_i = 0;
    for (vector<int>::const_iterator mask = m_train_mask.begin(); mask != m_train_mask.end();
         mask++) {
        for (int rd = 0; rd < m_rdim; rd++) {
            *i_multi = stats[rd][prb_i];
            if (*mask) {
                m_data_avg[rd] += *i_multi;
                m_data_var[rd] += (*i_multi) * (*i_multi);
            }
            ++i_multi;
        }
        prb_i++;
        if (*mask) {
            m_train_n++;
        }
    }

    for (int rd = 0; rd < m_rdim; rd++) {
        m_data_avg[rd] /= m_train_n;
        m_data_var[rd] /= m_train_n;
        m_data_var[rd] -= m_data_avg[rd] * m_data_avg[rd];
    }
    Rcpp::Rcerr << "response avg " << m_data_avg[0] << " response var " << m_data_var[0] << endl;
    Rcpp::Rcerr << "ranges " << m_min_range << " " << m_max_range << " bin " << m_spat_bin_size << endl;
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
    for (int prb_id = 0; prb_id < m_sequences.size(); prb_id++) {
        m_derivs[prb_id].resize(pos, vector<float>('T' + 1));
    }
    int max_spat_bin = m_spat_factors.size();
    fill(m_spat_factors.begin(), m_spat_factors.end(), 1.0 / max_spat_bin);
    m_spat_derivs.resize(m_sequences.size(), vector<float>(max_spat_bin));
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
        Rcpp::Rcerr << "set pos i to " << m_nuc_factors[i]['A'] << " " << m_nuc_factors[i]['C'] << " "
             << m_nuc_factors[i]['G'] << " " << m_nuc_factors[i]['T'] << endl;
    }

    m_aux_upds.resize(pwm.size());

    m_derivs.resize(m_sequences.size(), vector<vector<float>>(pwm.size(), vector<float>('T' + 1)));
    m_derivs.resize(m_sequences.size());
    for (int prb_id = 0; prb_id < m_sequences.size(); prb_id++) {
        m_derivs[prb_id].resize(pwm.size());
        for (int pos = 0; pos < pwm.size(); pos++) {
            m_derivs[prb_id][pos].resize('T' + 1);
        }
    }
    int max_spat_bin = m_spat_factors.size();
    fill(m_spat_factors.begin(), m_spat_factors.end(), 1.0 / max_spat_bin);
    m_spat_derivs.resize(m_sequences.size(), vector<float>(max_spat_bin));
}

void PWMLRegression::init_neighborhood(float resolution) {
    // ugly, but faster to program the the full recurence

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
    float prev_r2 = 0;
    m_cur_r2 = 0;    

    for (int phase = 0; phase < m_resolutions.size(); phase++) {
        Rcpp::Rcerr << "Start phase " << phase << " resol " << m_resolutions[phase] << endl;
        init_neighborhood(m_resolutions[phase]);
        m_spat_factor_step = m_spat_resolutions[phase];
        do {
            prev_r2 = m_cur_r2;
            init_energies();
            take_best_step();
            Rcpp::Rcerr << "S -lrtEP prev " << prev_r2 << " " << m_cur_r2 << endl;
			Rcpp::checkUserInterrupt();
        } while (m_cur_r2 > prev_r2 + m_imporve_epsilon);
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
    int max_prb_id = m_sequences.size();
    for (int prb_id = 0; prb_id < max_prb_id; prb_id++) {
        if (!m_train_mask[prb_id]) {
            continue;
        }
        string::const_iterator new_start = m_sequences[prb_id].begin();
        string::const_iterator new_end = m_sequences[prb_id].end();
        for (int pos = 0; pos < max_pos; pos++) {
            fill(m_derivs[prb_id][pos].begin(), m_derivs[prb_id][pos].end(), 0);
        }
        fill(m_spat_derivs[prb_id].begin(), m_spat_derivs[prb_id].end(), 0);
        update_seq_interval(prb_id, new_start, new_end, +1, 0);
    }
	
	if (m_logit) {
    	Rcpp::Rcerr << "done init energies " << endl;
	}
}

void PWMLRegression::update_seq_interval(int prb_id, string::const_iterator min_i,
                                         string::const_iterator max_i, int sign, int pos) {
    for (string::const_iterator i = min_i; i < max_i; i++) {
        int spat_bin = int(pos / m_spat_bin_size);
        pos++;
        string::const_iterator j = i;
        double prod = sign;
        vector<UpdAux>::iterator upds = m_aux_upds.begin();
        vector<vector<float>>::iterator deriv = m_derivs[prb_id].begin();
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
        m_spat_derivs[prb_id][spat_bin] += prod;
        prod *= m_spat_factors[spat_bin];
        for (upds = m_aux_upds.begin(); upds != m_aux_upds.end(); upds++) {
            *(upds->p) += prod / upds->factor;
        }
        // update all derivs -

        if (m_bidirect) {
            float rprod = sign;
            j = i;
            vector<UpdAux>::iterator upds = m_aux_upds.begin();
            vector<vector<float>>::reverse_iterator deriv = m_derivs[prb_id].rbegin();
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
            m_spat_derivs[prb_id][spat_bin] += rprod;
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
    float best_r2 = m_cur_r2;

    int max_pos = m_nuc_factors.size();

	// choose best step
    vector<float> probs(256, 0);
    for (int pos = 0; pos < max_pos; pos++) {
        // iter on neighborhood
        int neigh_size = m_cur_neigh.size();
        float best_pos_r2 = 0;
        for (int step_i = 0; step_i < neigh_size; step_i++) {
            // update probs
            probs = m_nuc_factors[pos];
            for (vector<NeighStep>::iterator delta = m_cur_neigh[step_i].begin();
                 delta != m_cur_neigh[step_i].end(); delta++) {
                probs[delta->nuc] += delta->diff;
                if (probs[delta->nuc] <= 0) {
                    probs[delta->nuc] = m_min_prob;
                    // normalize
                }
            }
            float cur_step_r2 = compute_cur_r2(pos, probs);
            // Rcpp::Rcerr << "r2 at step " << step_i << " was " << cur_step_r2 << endl;
            if (best_pos_r2 < cur_step_r2) {
                best_pos_r2 = cur_step_r2;
            }
            if (cur_step_r2 > best_r2) {
                best_step = step_i;
                best_pos = pos;
                best_r2 = cur_step_r2;
            }
        }
        if (m_logit) {
            Rcpp::Rcerr << "best r2 at pos " << pos << " = " << best_pos_r2 << endl;
        }
        // iter on all changes: single nuc and double nucs (regression?)
    }
    Rcpp::Rcerr << "best step was " << best_step << " pos " << best_pos << " r2 " << best_r2 << endl;
    int max_spat_bin = m_spat_factors.size();
    float best_spat_r2 = m_cur_r2;
    int best_spat_bin = -1;
    float best_spat_diff = 0;    
    float spat_r2 = compute_cur_r2_spat();
    Rcpp::Rcerr << "spat r2 without change = " << spat_r2 << " cur is " << m_cur_r2 << endl;
    for (int spat_bin = 0; spat_bin < max_spat_bin; spat_bin++) {
        m_spat_factors[spat_bin] += m_spat_factor_step;
        float spat_r2 = compute_cur_r2_spat();
        // Rcpp::Rcerr << "bin " << spat_bin << " add spat r2 " << spat_r2 << endl;
        if (spat_r2 > best_spat_r2) {
            best_spat_r2 = spat_r2;
            best_spat_bin = spat_bin;
            best_spat_diff = m_spat_factor_step;
        }
        m_spat_factors[spat_bin] -= 2 * m_spat_factor_step;
        if (m_spat_factors[spat_bin] >= 0) {
            spat_r2 = compute_cur_r2_spat();
            // Rcpp::Rcerr << "bin " << spat_bin << " del spat r2 " << spat_r2 << endl;
            if (spat_r2 > best_spat_r2) {
                best_spat_r2 = spat_r2;
                best_spat_bin = spat_bin;
                best_spat_diff = -m_spat_factor_step;
            }
        }
        m_spat_factors[spat_bin] += m_spat_factor_step;
    }
    Rcpp::Rcerr << "best spat step was bin " << best_spat_bin << " diff " << best_spat_diff << " r2 "
         << best_spat_r2 << endl;
    if (best_r2 == m_cur_r2) {
        Rcpp::Rcerr << "no improvement" << endl;
        return;
    }
    if (best_r2 > best_spat_r2) {
        vector<float> &pos_probs = m_nuc_factors[best_pos];
        for (vector<NeighStep>::iterator delta = m_cur_neigh[best_step].begin();
             delta != m_cur_neigh[best_step].end(); delta++) {
            pos_probs[delta->nuc] += delta->diff;
            if (pos_probs[delta->nuc] <= 0) {
                pos_probs[delta->nuc] = m_min_prob;
                // noramalize
            }
        }
        float tot = pos_probs['A'] + pos_probs['C'] + pos_probs['G'] + pos_probs['T'];
        Rcpp::Rcerr << "update pos " << best_pos << " tot = " << tot << endl;
        pos_probs['A'] /= tot;
        pos_probs['C'] /= tot;
        pos_probs['G'] /= tot;
        pos_probs['T'] /= tot;
        m_cur_r2 = best_r2;
    } else {
        Rcpp::Rcerr << "update spat bin " << best_spat_bin << " diff " << best_spat_diff << endl;
        m_spat_factors[best_spat_bin] += best_spat_diff;
        // normalize
        for (int bin = 0; bin < m_spat_factors.size(); bin++) {
            m_spat_factors[bin] /= (1 + best_spat_diff);
        }
        m_cur_r2 = best_spat_r2;
    }
    if (m_logit) {
        Rcpp::Rcerr << "after step";
        report_cur_lpwm();
        Rcpp::Rcerr << endl;
    }
}

// float PWMLRegression::compute_cur_wilcox() {}

float PWMLRegression::compute_cur_r2(int pos, vector<float> &probs) {
    vector<double> xy(m_rdim, 0);

    double ex = 0;
    double ex2 = 0;
    int max_prb_id = m_sequences.size();
	
    vector<vector<vector<float>>>::iterator prb_deriv = m_derivs.begin();
    vector<float>::iterator resp = m_interv_stat.begin();
    for (int prb_id = 0; prb_id < max_prb_id; prb_id++) {
        if (m_train_mask[prb_id]) {
            vector<float> &deriv = (*prb_deriv)[pos];
            float v = probs['A'] * deriv['A'] + probs['C'] * deriv['C'] + probs['G'] * deriv['G'] +
                      probs['T'] * deriv['T'];
            // push v,response to priority queue
            // compute wilcoxon
            ex += v;
            ex2 += v * v;
            for (int rd = 0; rd < m_rdim; rd++) {
                xy[rd] += v * *resp;
                resp++;
            }
        }
        prb_deriv++;
    }
	if (m_logit) {
    	Rcpp::Rcerr << "done initing xy for all sequences, xy 0 is " << xy[0] << endl;
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
    //	Rcpp::Rcerr << "n = " << n << " ex " << ex << " ex2 " << ex2 << " xy " << xy << " pred var " <<
    //pred_var << " cov = " << cov << " resp  avg " << m_data_avg << " resp _var " << m_data_var <<
    //endl;

    if (std::isnan(tot_r2)) {
        Rcpp::Rcerr << "Nan at at r2 var " << pred_var << " " << ex << " " << ex2 << " " << m_data_avg[0]
             << endl;
    }
    return (tot_r2);
}

float PWMLRegression::compute_cur_r2_spat() {	
    vector<double> xy(m_rdim, 0);
    float ex = 0;
    float ex2 = 0;    
    int max_prb_id = m_sequences.size();	

    vector<vector<float>>::iterator prb_derivs = m_spat_derivs.begin();    
    vector<float>::iterator resp = m_interv_stat.begin();   	
    for (int prb_id = 0; prb_id < max_prb_id; prb_id++) {        
        if (m_train_mask[prb_id]) {
            float v = 0;            
            vector<float>::iterator fact = m_spat_factors.begin();
            for (vector<float>::iterator bin = prb_derivs->begin(); bin != prb_derivs->end();
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
        prb_derivs++;
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

Rcpp::DataFrame PWMLRegression::output_pssm_df(int psid){
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
		Rcpp::Named("psid") = psid,
		Rcpp::Named("pos") = poss,
		Rcpp::Named("A") = pssm_A, 
		Rcpp::Named("C") = pssm_C,
		Rcpp::Named("G") = pssm_G,
		Rcpp::Named("T") = pssm_T
	);

	return pssm;
}

Rcpp::DataFrame PWMLRegression::output_spat_df(int psid){
	int max_spat_bin = m_spat_factors.size();
	vector<int> spat_bins(max_spat_bin);
	vector<float> spat_factors(max_spat_bin, 0);

	for (int bin = 0; bin < max_spat_bin; bin++) {
		spat_bins[bin] = m_spat_bin_size * bin;
		spat_factors[bin] = m_spat_factors[bin];
	}

	Rcpp::DataFrame spat_df = Rcpp::DataFrame::create(
		Rcpp::Named("psid") = psid,
		Rcpp::Named("bin") = spat_bins,
		Rcpp::Named("spat_factor") = spat_factors		
	);

	return(spat_df);
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
}

void PWMLRegression::fill_predictions(vector<float> &preds) {
    // updateting the m_a, m_b
    init_energies();
    Rcpp::Rcerr << "done init energ" << endl;
    float r2 = compute_cur_r2_spat();
    Rcpp::Rcerr << "done comp r2 = " << r2 << endl;

    preds.resize(m_interv_stat.size());
    vector<float>::iterator prs = preds.begin();
    vector<vector<vector<float>>>::iterator prb_deriv = m_derivs.begin();
    vector<float>::iterator resp = m_interv_stat.begin();
    vector<float> &factors = m_nuc_factors[0];
    for (int prb_id = 0; prb_id < m_interv_stat.size(); prb_id++) {
        if (m_train_mask[prb_id]) {
            vector<float> &deriv = (*prb_deriv)[0];

            float predict = factors['A'] * deriv['A'] + factors['C'] * deriv['C'] +
                            factors['G'] * deriv['G'] + factors['T'] * deriv['T'];

            //		Rcpp::Rcerr << "prb " << prb_id << " deriv " << deriv['A'] << " " << deriv['C'] << " " <<
            //deriv['G'] << " " << deriv['T'] << " factors " << factors['A'] << " " << factors['C']
            //<< " " << factors['G'] << " " << factors['T'] << endl;

            *prs = predict;
        } else {
            *prs = -1e+30;
        }
        prs++;
        prb_deriv++;
        resp++;
    }
}
