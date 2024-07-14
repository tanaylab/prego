#include "port.h"

#include "BitVecIter.h"
#include "LeastSquare.h"
#include "ProgressReporter.h"
#include "Random.h"
#include "SpecialFunc.h"
#include "dnastrutil.h"
#include "FunctionInterpolator.h"



#include "KMerMultiStat.h"
#include "PWMLRegression.h"
#include "PssmRegression.h"
#include <Rcpp.h>
#include <string.h>
#include <vector>
#include <unordered_map>
using namespace std;

// [[Rcpp::plugins("cpp17")]]

// [[Rcpp::export]]
std::string get_consensus_cpp(const Rcpp::NumericMatrix &pssm_mat, const float &single_thresh,
                              const float &double_thresh) {

    DnaPSSM pssm;

    pssm.resize(pssm_mat.nrow());
    for (int i = 0; i < pssm_mat.nrow(); i++) {
        pssm[i].set_weight('A', pssm_mat(i, 0));
        pssm[i].set_weight('C', pssm_mat(i, 1));
        pssm[i].set_weight('G', pssm_mat(i, 2));
        pssm[i].set_weight('T', pssm_mat(i, 3));
    }
    pssm.normalize();

    return (pssm.get_consensus(single_thresh, double_thresh));
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_pwm_cpp(const Rcpp::StringVector &sequences,
                                    const Rcpp::NumericMatrix &pssm_mat, const bool &is_bidirect,
                                    const int &spat_min, const int &spat_max,
                                    const Rcpp::NumericVector &spat_factor, const int &bin_size, const bool &use_max = false) {
    vector<string> seqs = Rcpp::as<vector<string>>(sequences);
    vector<float> spat_fac = Rcpp::as<vector<float>>(spat_factor);
    int smin = spat_min;
    int smax = spat_max;

    DnaPSSM pssm;

    pssm.set_bidirect(is_bidirect);
    pssm.resize(pssm_mat.nrow());
    pssm.set_range(smin, smax);
    for (int i = 0; i < pssm_mat.nrow(); i++) {
        pssm[i].set_weight('A', pssm_mat(i, 0));
        pssm[i].set_weight('C', pssm_mat(i, 1));
        pssm[i].set_weight('G', pssm_mat(i, 2));
        pssm[i].set_weight('T', pssm_mat(i, 3));
    }
    pssm.normalize();

    DnaPWML pwml(pssm, spat_fac, bin_size);
    Rcpp::NumericVector preds(seqs.size());

    for (size_t i = 0; i < seqs.size(); i++) {
        float energy;
        if (use_max){
            pwml.integrate_energy_max(seqs[i], energy);
        } else {
            pwml.integrate_energy(seqs[i], energy);
        }
        
        preds[i] = energy;
    }

    return (preds);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_local_pwm_cpp(const Rcpp::StringVector &sequences,
                                    const Rcpp::NumericMatrix &pssm_mat, const bool &is_bidirect,
                                    const int &spat_min, const int &spat_max,
                                    const Rcpp::NumericVector &spat_factor, const int &bin_size){
    vector<string> seqs = Rcpp::as<vector<string>>(sequences);
    vector<float> spat_fac = Rcpp::as<vector<float>>(spat_factor);
    int smin = spat_min;
    int smax = spat_max;
    int motif_len = pssm_mat.nrow();

    DnaPSSM pssm;

    pssm.set_bidirect(is_bidirect);
    pssm.resize(motif_len);
    pssm.set_range(smin, smax);
    for (int i = 0; i < motif_len; i++) {
        pssm[i].set_weight('A', pssm_mat(i, 0));
        pssm[i].set_weight('C', pssm_mat(i, 1));
        pssm[i].set_weight('G', pssm_mat(i, 2));
        pssm[i].set_weight('T', pssm_mat(i, 3));
    }
    pssm.normalize();

    DnaPWML pwml(pssm, spat_fac, bin_size);
    Rcpp::NumericVector preds(seqs.size());

    Rcpp::NumericMatrix local_preds(seqs.size(), seqs[0].length());

    for (size_t i = 0; i < seqs.size(); i++) {
        // go over windows of motif_len and compute local scores
        for (size_t j = 0; j < seqs[i].length() - motif_len + 1; j++) {
            float energy;
            pwml.integrate_energy(seqs[i].substr(j, motif_len), energy);
            local_preds(i, j) = energy;
        }

        // fill in the rest with NA
        for (size_t j = seqs[i].length() - motif_len + 1; j < seqs[i].length(); j++) {
            local_preds(i, j) = Rcpp::NumericVector::get_na();
        }
    }

    return (local_preds);
}

// [[Rcpp::export]]
Rcpp::StringVector mask_sequences_cpp(const Rcpp::StringVector &sequences,
                                      const Rcpp::NumericMatrix &pssm_mat, const bool &is_bidirect,
                                      const int &spat_min, const int &spat_max,
                                      const Rcpp::NumericVector &spat_factor, const int &bin_size,
                                      const float &mask_thresh, const Rcpp::LogicalVector &pos_mask) {
    vector<string> seqs = Rcpp::as<vector<string>>(sequences);
    vector<float> spat_fac = Rcpp::as<vector<float>>(spat_factor);
    int smin = spat_min;
    int smax = spat_max;
    int motif_len = pssm_mat.nrow();    

    // clone sequences
    vector<string> output = seqs;

    DnaPSSM pssm;

    pssm.set_bidirect(is_bidirect);
    pssm.resize(pssm_mat.nrow());
    pssm.set_range(smin, smax);
    for (int i = 0; i < pssm_mat.nrow(); i++) {
        pssm[i].set_weight('A', pssm_mat(i, 0));
        pssm[i].set_weight('C', pssm_mat(i, 1));
        pssm[i].set_weight('G', pssm_mat(i, 2));
        pssm[i].set_weight('T', pssm_mat(i, 3));
    }
    pssm.normalize();    

    DnaPWML pwml(pssm, spat_fac, bin_size);
    Rcpp::StringVector preds(seqs.size());

    for (size_t i = 0; i < seqs.size(); i++) {
        // go over windows of motif_len and mask if the score is above the threshold
        for (int j = 0; j < (int)seqs[i].size() - motif_len; j++) {
            float energy;
            pwml.integrate_energy(seqs[i].substr(j, motif_len), energy);
            if (energy > mask_thresh) {
                // replace only the positions that are not masked (i.e. are informative)
                for (int k = 0; k < motif_len; k++) {
                    if (pos_mask[k]) {
                        output[i][j + k] = 'N';
                    }
                }                
            }
        }    

        // mask the few bases that were shorter than motif_len
        for (int j = (int)seqs[i].size() - motif_len; j < (int)seqs[i].size(); j++) {
            output[i][j] = 'N';
        }
    }

    return (Rcpp::wrap(output));
}

// [[Rcpp::export]]
Rcpp::List regress_pwm_cpp(
    const Rcpp::StringVector &sequences, const Rcpp::DataFrame &response,
    const Rcpp::LogicalVector &is_train_logical, const std::string &motif, const int &spat_min,
    const int &spat_max, const float &min_nuc_prob, const int &spat_bin,
    const float &improve_epsilon, const bool &is_bidirect, const float &unif_prior,
    const std::string &score_metric, const int &verbose, const int &seed,
    const Rcpp::NumericMatrix &pssm_mat, const Rcpp::Nullable<Rcpp::NumericVector> &spat_factor,
    const float &consensus_single_thresh, const float &consensus_double_thresh,
    const int &num_folds = 1, const float &energy_epsilon = 0, const bool &log_energy = false,
    Rcpp::Nullable<Rcpp::Function> energy_func = R_NilValue, const float &xmin = -100,
    const float &xmax = 100, const int &npts = 1000, const bool &optimize_pwm = true,
    const bool &optimize_spat = true, const bool &symmetrize_spat = true) {
    Random::reset(seed);
    vector<vector<float>> response_stat = Rcpp::as<vector<vector<float>>>(response);    

    vector<string> seqs = Rcpp::as<vector<string>>(sequences);

    vector<int> is_train = Rcpp::as<vector<int>>(is_train_logical);

    // initialize resolution objects
    vector<float> res(4);
    res[0] = 0.05;
    res[1] = 0.02;
    res[2] = 0.01;
    res[3] = 0.005;

    vector<float> spres(4);
    spres[0] = 0.01;
    spres[1] = 0.01;
    spres[2] = 0.01;
    spres[3] = 0.005;

    int smin = spat_min;
    int smax = spat_max;

    if (verbose) {
        Rcpp::Rcerr << "into pwmlreg" << endl;
    }

    if (float(int((smax - smin) / spat_bin)) != float((smax - smin)) / float(spat_bin)) {
        Rcpp::stop("spat_bin must divide spat_max - spat_min");
    }

    PWMLRegression pwmlreg(seqs, is_train, smin, smax, min_nuc_prob, spat_bin, res, spres,
                           improve_epsilon, 0.001, unif_prior, score_metric, num_folds, log_energy,
                           energy_epsilon, energy_func, xmin, xmax, npts, optimize_pwm, optimize_spat, 
                           symmetrize_spat);

    pwmlreg.m_logit = verbose;

    // add responses (and compute their statistics - avg and variance)
    pwmlreg.add_responses(response_stat);

    string seedmot(motif);
    if (motif.empty()) { // initialize using pssm
        if (verbose) {
            Rcpp::Rcerr << "using pre-computed pssm" << std::endl;
        }

        DnaPSSM pssm;
        pssm.set_bidirect(is_bidirect);
        pssm.resize(pssm_mat.nrow());
        pssm.set_range(smin, smax);
        for (int i = 0; i < pssm_mat.nrow(); i++) {
            pssm[i].set_weight('A', pssm_mat(i, 0));
            pssm[i].set_weight('C', pssm_mat(i, 1));
            pssm[i].set_weight('G', pssm_mat(i, 2));
            pssm[i].set_weight('T', pssm_mat(i, 3));
        }
        pssm.normalize();        
        if (spat_factor.isNotNull()) {
            vector<float> spat_fac = Rcpp::as<vector<float>>(spat_factor);            
            pwmlreg.init_pwm_spat(pssm, spat_fac);            
        } else {
            pwmlreg.init_pwm(pssm);
        }
    } else { // initialize using a seed motif
        if (spat_factor.isNotNull()) {
            vector<float> spat_fac = Rcpp::as<vector<float>>(spat_factor);
            pwmlreg.init_seed_spat(seedmot, spat_fac, is_bidirect);
        } else {
            pwmlreg.init_seed(seedmot, is_bidirect);
        }
        
        if (verbose) {
            Rcpp::Rcerr << "done init seed " << seedmot << endl;
        }
    }

    // main loop
    pwmlreg.optimize();

    // get predictions
    DnaPWML pwml;
    pwmlreg.get_model(pwml);

    vector<float> preds(seqs.size());

    for (size_t i = 0; i < seqs.size(); i++) {
        float energy;
        pwml.integrate_energy(seqs[i], energy);
        preds[i] = energy;
    }

    if (energy_func.isNotNull()) {
        preds = Rcpp::as<vector<float>>(Rcpp::as<Rcpp::Function>(energy_func)(preds));

        if (preds.size() != seqs.size()) {
            Rcpp::stop("Energy function must return a vector of the same length as the number of "
                       "sequences");
        }
    }

    // prepare output
    Rcpp::List res_list = Rcpp::List::create(
        Rcpp::Named("pssm") = pwmlreg.output_pssm_df(0),
        Rcpp::Named("spat") = pwmlreg.output_spat_df(0), Rcpp::Named("pred") = preds,
        Rcpp::Named("consensus") =
            pwml.get_pssm().get_consensus(consensus_single_thresh, consensus_double_thresh));
    return (res_list);
}

// [[Rcpp::export]]
Rcpp::DataFrame screen_kmers_cpp(const Rcpp::StringVector &sequences,
                                 const Rcpp::DataFrame &response,
                                 const Rcpp::LogicalVector &is_train_logical, const int &L,
                                 const int &from_range, const int &to_range, const float &min_cor,
                                 const int &min_gap, const int &max_gap,
                                 const int &n_in_train, const int &seed, const bool &verbose) {

    Random::reset(seed);
    vector<vector<float>> response_stat = Rcpp::as<vector<vector<float>>>(response);
    int resp_dim = response_stat.size();

    vector<string> seqs = Rcpp::as<vector<string>>(sequences);

    vector<int> is_train = Rcpp::as<vector<int>>(is_train_logical);

    float r2_thresh = min_cor * min_cor;
    int bin_num = 20;
    int norm = 0;
    float norm_factor = 1;

    KMerMultiStat multi(L, 0, min_gap, max_gap, &seqs, &is_train, bin_num, norm, norm_factor,
                        response_stat, from_range, to_range, 2, verbose);

    string best_mot = "";
    float best_r2 = 0;
    vector<string> foc_mots;
    vector<float> foc_scores;
    vector<int> foc_ids;
    vector<float> response_avg(resp_dim, 0);
    vector<float> response_var(resp_dim, 0);

    // normalize response variables
    for (int ri = 0; ri < resp_dim; ri++) {
        vector<int>::iterator train = is_train.begin();
        for (vector<float>::iterator i = response_stat[ri].begin(); i != response_stat[ri].end();
             i++) {
            if (*train) {
                response_avg[ri] += *i;
                response_var[ri] += (*i) * (*i);
            }
            train++;
        }
        response_avg[ri] /= n_in_train;
        response_var[ri] /= n_in_train;
        response_var[ri] -= response_avg[ri] * response_avg[ri];
    }

    if (verbose) {
        Rcpp::Rcerr << "done normalizing response " << endl;
    }

    // results vectors
    vector<string> res_kmer;
    vector<float> res_max_r2;
    vector<float> res_avg_multi;
    vector<float> res_multi_var;
    vector<vector<float>> res_cors(resp_dim);

    // iterate over all kmers
    ProgressReporter progress;
    if (verbose) {
        progress.init(multi.get_pat_size(), 1);
    }

    for (auto k = multi.get_pat_begin(); k != multi.get_pat_end(); k++) {
        // if pattern has 'N' skip
        if (k->first.find('N') != string::npos) {
            continue;
        }
        vector<float> cov(resp_dim, 0);
        vector<float> corr(resp_dim, 0);
        float avg_multi = 0;
        float tot_multi2 = 0;
        // compute cov over the rdim
        const vector<pair<int, vector<float>>> &multi = k->second;
        for (size_t m = 1; m < multi.size(); m++) {
            avg_multi += multi[m].first * m;
            tot_multi2 += multi[m].first * m * m;
            for (int ri = 0; ri < resp_dim; ri++) {
                cov[ri] += m * multi[m].second[ri];
            }
        }
        avg_multi /= n_in_train;
        float multi_var = tot_multi2 / n_in_train - avg_multi * avg_multi;
        float max_r2 = 0;
        for (int ri = 0; ri < resp_dim; ri++) {
            cov[ri] /= n_in_train;
            cov[ri] -= avg_multi * response_avg[ri];
            corr[ri] = cov[ri] / sqrt(multi_var * response_var[ri]);
            if (max_r2 < corr[ri] * corr[ri]) {
                max_r2 = corr[ri] * corr[ri];
            }
        }
        if (max_r2 > r2_thresh) {
            res_kmer.push_back(k->first);
            res_max_r2.push_back(max_r2);
            res_avg_multi.push_back(avg_multi);
            res_multi_var.push_back(multi_var);
            for (int ri = 0; ri < resp_dim; ri++) {
                res_cors[ri].push_back(corr[ri]);
            }

            foc_mots.push_back(k->first);
            foc_scores.push_back(max_r2);
            foc_ids.push_back(foc_ids.size());
            if (max_r2 > best_r2) {
                best_r2 = max_r2;
                best_mot = k->first;
                if (verbose) {
                    Rcpp::Rcerr << "new best motif: " << best_mot << " r2: " << best_r2 << endl;
                }
            }
        }
        if (verbose) {
            progress.report(1);
        }
    }

    if (verbose) {
        progress.report_last();
    }

    if (verbose) {
        Rcpp::Rcerr << "done screening " << endl;
    }

    // assemble results
    Rcpp::DataFrame res = Rcpp::DataFrame::create(
        Rcpp::Named("kmer") = Rcpp::wrap(res_kmer), Rcpp::Named("max_r2") = Rcpp::wrap(res_max_r2),
        Rcpp::Named("avg_n") = Rcpp::wrap(res_avg_multi),
        Rcpp::Named("avg_var") = Rcpp::wrap(res_multi_var));

    for (int ri = 0; ri < resp_dim; ri++) {
        res.push_back(Rcpp::wrap(res_cors[ri]),
                      Rcpp::as<string>(Rcpp::as<Rcpp::CharacterVector>(response.names())[ri]));
    }
    return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector interpolateFunction(Rcpp::Function func, float xmin, float xmax, int npts, Rcpp::NumericVector x) {
    FunctionInterpolator interp(func, xmin, xmax, npts);
    vector<float> xvals = Rcpp::as<vector<float>>(x);
    
    return Rcpp::wrap(interp.interpolate(xvals));
}

// [[Rcpp::export]]
Rcpp::DataFrame dinuc_distribution(Rcpp::StringVector sequences, int size = 1000) {
    if(size < 1) {
        Rcpp::stop("Size must be at least 1");
    }
    
    // Prepare a vector of unordered_maps to hold dinucleotide frequencies at each position
    vector<unordered_map<string, int>> nucleotide_counts(size);

    // Convert StringVector to vector of strings
    vector<std::string> seqs = Rcpp::as<vector<std::string>>(sequences);
    
    // Iterate over each sequence
    for(const auto& sequence : seqs) {
        // Validate that the sequence length is at least 'size'
        if((int)sequence.size() < size) {
            Rcpp::stop("All sequences must have a length of at least 'size'");
        }
        
        // In each sequence, iterate over each position and extract the dinucleotide
        for(int pos = 0; pos < size - 1; ++pos) {
            std::string dinucleotide = sequence.substr(pos, 2);
            
            // Update the counts in the unordered_map for the corresponding dinucleotide and position
            nucleotide_counts[pos][dinucleotide]++;
        }
    }
    
    // After iterating over all sequences and positions, calculate the percentages
    int total_sequences = sequences.size();

    // Create an Rcpp data frame with the required structure and fill it with the calculated percentages
    Rcpp::NumericMatrix percentage_matrix(size, 16);
    Rcpp::CharacterVector dinucleotides = {"AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"};
    
    for(int pos = 0; pos < size; ++pos) {
        for(int dinuc_i = 0; dinuc_i < 16; ++dinuc_i) {
            string dinucleotide = Rcpp::as<string>(dinucleotides[dinuc_i]);
            percentage_matrix(pos, dinuc_i) = (double)nucleotide_counts[pos][dinucleotide] / total_sequences;
        }
    }
    
    return Rcpp::DataFrame::create(
        Rcpp::Named("pos") = Rcpp::seq(1, size), 
        Rcpp::Named("AA") = percentage_matrix(Rcpp::_, 0), 
        Rcpp::Named("AC") = percentage_matrix(Rcpp::_, 1),
        Rcpp::Named("AG") = percentage_matrix(Rcpp::_, 2), 
        Rcpp::Named("AT") = percentage_matrix(Rcpp::_, 3),
        Rcpp::Named("CA") = percentage_matrix(Rcpp::_, 4), 
        Rcpp::Named("CC") = percentage_matrix(Rcpp::_, 5),
        Rcpp::Named("CG") = percentage_matrix(Rcpp::_, 6), 
        Rcpp::Named("CT") = percentage_matrix(Rcpp::_, 7),
        Rcpp::Named("GA") = percentage_matrix(Rcpp::_, 8), 
        Rcpp::Named("GC") = percentage_matrix(Rcpp::_, 9),
        Rcpp::Named("GG") = percentage_matrix(Rcpp::_, 10), 
        Rcpp::Named("GT") = percentage_matrix(Rcpp::_, 11),
        Rcpp::Named("TA") = percentage_matrix(Rcpp::_, 12), 
        Rcpp::Named("TC") = percentage_matrix(Rcpp::_, 13),
        Rcpp::Named("TG") = percentage_matrix(Rcpp::_, 14), 
        Rcpp::Named("TT") = percentage_matrix(Rcpp::_, 15)
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame trinuc_distribution(Rcpp::StringVector sequences, int size = 1000) {
    if (size < 1) {
        Rcpp::stop("Size must be at least 1");
    }

    // Prepare a vector of unordered_maps to hold trinucleotide frequencies at each position
    std::vector<std::unordered_map<std::string, int>> nucleotide_counts(size);

    // Convert StringVector to vector of strings
    std::vector<std::string> seqs = Rcpp::as<std::vector<std::string>>(sequences);

    // Iterate over each sequence
    for (const auto &sequence : seqs) {
        // Validate that the sequence length is at least 'size'
        if ((int)sequence.size() < size) {
            Rcpp::stop("All sequences must have a length of at least 'size'");
        }

        // In each sequence, iterate over each position and extract the trinucleotide
        for (int pos = 0; pos < size - 2; ++pos) {
            std::string trinucleotide = sequence.substr(pos, 3);

            // Update the counts in the unordered_map for the corresponding trinucleotide and
            // position
            nucleotide_counts[pos][trinucleotide]++;
        }
    }

    // After iterating over all sequences and positions, calculate the percentages
    int total_sequences = sequences.size();

    // Create an Rcpp data frame with the required structure and fill it with the calculated
    // percentages
    Rcpp::NumericMatrix percentage_matrix(size, 64);
    Rcpp::CharacterVector trinucleotides = {
        "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA",
        "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC",
        "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG",
        "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT",
        "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

    for (int pos = 0; pos < size; ++pos) {
        for (int trinuc_i = 0; trinuc_i < 64; ++trinuc_i) {
            std::string trinucleotide = Rcpp::as<std::string>(trinucleotides[trinuc_i]);
            percentage_matrix(pos, trinuc_i) =
                (double)nucleotide_counts[pos][trinucleotide] / total_sequences;
        }
    }

    // Create the DataFrame
    Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::Named("pos") = Rcpp::seq(1, size));

    // Add columns for each trinucleotide
    for (int i = 0; i < 64; ++i) {
        df.push_back(percentage_matrix(Rcpp::_, i), Rcpp::as<std::string>(trinucleotides[i]));
    }

    return df;
}