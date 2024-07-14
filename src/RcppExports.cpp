// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// kmer_matrix_cpp
Rcpp::IntegerMatrix kmer_matrix_cpp(Rcpp::CharacterVector sequences, int kmer_length, int from_range, Rcpp::Nullable<int> to_range, Rcpp::Nullable<Rcpp::CharacterVector> mask, bool add_mask, size_t max_gap);
RcppExport SEXP _prego_kmer_matrix_cpp(SEXP sequencesSEXP, SEXP kmer_lengthSEXP, SEXP from_rangeSEXP, SEXP to_rangeSEXP, SEXP maskSEXP, SEXP add_maskSEXP, SEXP max_gapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< int >::type kmer_length(kmer_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type from_range(from_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type to_range(to_rangeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< bool >::type add_mask(add_maskSEXP);
    Rcpp::traits::input_parameter< size_t >::type max_gap(max_gapSEXP);
    rcpp_result_gen = Rcpp::wrap(kmer_matrix_cpp(sequences, kmer_length, from_range, to_range, mask, add_mask, max_gap));
    return rcpp_result_gen;
END_RCPP
}
// get_consensus_cpp
std::string get_consensus_cpp(const Rcpp::NumericMatrix& pssm_mat, const float& single_thresh, const float& double_thresh);
RcppExport SEXP _prego_get_consensus_cpp(SEXP pssm_matSEXP, SEXP single_threshSEXP, SEXP double_threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type pssm_mat(pssm_matSEXP);
    Rcpp::traits::input_parameter< const float& >::type single_thresh(single_threshSEXP);
    Rcpp::traits::input_parameter< const float& >::type double_thresh(double_threshSEXP);
    rcpp_result_gen = Rcpp::wrap(get_consensus_cpp(pssm_mat, single_thresh, double_thresh));
    return rcpp_result_gen;
END_RCPP
}
// compute_pwm_cpp
Rcpp::NumericVector compute_pwm_cpp(const Rcpp::StringVector& sequences, const Rcpp::NumericMatrix& pssm_mat, const bool& is_bidirect, const int& spat_min, const int& spat_max, const Rcpp::NumericVector& spat_factor, const int& bin_size, const bool& use_max);
RcppExport SEXP _prego_compute_pwm_cpp(SEXP sequencesSEXP, SEXP pssm_matSEXP, SEXP is_bidirectSEXP, SEXP spat_minSEXP, SEXP spat_maxSEXP, SEXP spat_factorSEXP, SEXP bin_sizeSEXP, SEXP use_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type pssm_mat(pssm_matSEXP);
    Rcpp::traits::input_parameter< const bool& >::type is_bidirect(is_bidirectSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_min(spat_minSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_max(spat_maxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type spat_factor(spat_factorSEXP);
    Rcpp::traits::input_parameter< const int& >::type bin_size(bin_sizeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type use_max(use_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_pwm_cpp(sequences, pssm_mat, is_bidirect, spat_min, spat_max, spat_factor, bin_size, use_max));
    return rcpp_result_gen;
END_RCPP
}
// compute_local_pwm_cpp
Rcpp::NumericMatrix compute_local_pwm_cpp(const Rcpp::StringVector& sequences, const Rcpp::NumericMatrix& pssm_mat, const bool& is_bidirect, const int& spat_min, const int& spat_max, const Rcpp::NumericVector& spat_factor, const int& bin_size);
RcppExport SEXP _prego_compute_local_pwm_cpp(SEXP sequencesSEXP, SEXP pssm_matSEXP, SEXP is_bidirectSEXP, SEXP spat_minSEXP, SEXP spat_maxSEXP, SEXP spat_factorSEXP, SEXP bin_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type pssm_mat(pssm_matSEXP);
    Rcpp::traits::input_parameter< const bool& >::type is_bidirect(is_bidirectSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_min(spat_minSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_max(spat_maxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type spat_factor(spat_factorSEXP);
    Rcpp::traits::input_parameter< const int& >::type bin_size(bin_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_local_pwm_cpp(sequences, pssm_mat, is_bidirect, spat_min, spat_max, spat_factor, bin_size));
    return rcpp_result_gen;
END_RCPP
}
// mask_sequences_cpp
Rcpp::StringVector mask_sequences_cpp(const Rcpp::StringVector& sequences, const Rcpp::NumericMatrix& pssm_mat, const bool& is_bidirect, const int& spat_min, const int& spat_max, const Rcpp::NumericVector& spat_factor, const int& bin_size, const float& mask_thresh, const Rcpp::LogicalVector& pos_mask);
RcppExport SEXP _prego_mask_sequences_cpp(SEXP sequencesSEXP, SEXP pssm_matSEXP, SEXP is_bidirectSEXP, SEXP spat_minSEXP, SEXP spat_maxSEXP, SEXP spat_factorSEXP, SEXP bin_sizeSEXP, SEXP mask_threshSEXP, SEXP pos_maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type pssm_mat(pssm_matSEXP);
    Rcpp::traits::input_parameter< const bool& >::type is_bidirect(is_bidirectSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_min(spat_minSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_max(spat_maxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type spat_factor(spat_factorSEXP);
    Rcpp::traits::input_parameter< const int& >::type bin_size(bin_sizeSEXP);
    Rcpp::traits::input_parameter< const float& >::type mask_thresh(mask_threshSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type pos_mask(pos_maskSEXP);
    rcpp_result_gen = Rcpp::wrap(mask_sequences_cpp(sequences, pssm_mat, is_bidirect, spat_min, spat_max, spat_factor, bin_size, mask_thresh, pos_mask));
    return rcpp_result_gen;
END_RCPP
}
// regress_pwm_cpp
Rcpp::List regress_pwm_cpp(const Rcpp::StringVector& sequences, const Rcpp::DataFrame& response, const Rcpp::LogicalVector& is_train_logical, const std::string& motif, const int& spat_min, const int& spat_max, const float& min_nuc_prob, const int& spat_bin, const float& improve_epsilon, const bool& is_bidirect, const float& unif_prior, const std::string& score_metric, const int& verbose, const int& seed, const Rcpp::NumericMatrix& pssm_mat, const Rcpp::Nullable<Rcpp::NumericVector>& spat_factor, const float& consensus_single_thresh, const float& consensus_double_thresh, const int& num_folds, const float& energy_epsilon, const bool& log_energy, Rcpp::Nullable<Rcpp::Function> energy_func, const float& xmin, const float& xmax, const int& npts, const bool& optimize_pwm, const bool& optimize_spat, const bool& symmetrize_spat);
RcppExport SEXP _prego_regress_pwm_cpp(SEXP sequencesSEXP, SEXP responseSEXP, SEXP is_train_logicalSEXP, SEXP motifSEXP, SEXP spat_minSEXP, SEXP spat_maxSEXP, SEXP min_nuc_probSEXP, SEXP spat_binSEXP, SEXP improve_epsilonSEXP, SEXP is_bidirectSEXP, SEXP unif_priorSEXP, SEXP score_metricSEXP, SEXP verboseSEXP, SEXP seedSEXP, SEXP pssm_matSEXP, SEXP spat_factorSEXP, SEXP consensus_single_threshSEXP, SEXP consensus_double_threshSEXP, SEXP num_foldsSEXP, SEXP energy_epsilonSEXP, SEXP log_energySEXP, SEXP energy_funcSEXP, SEXP xminSEXP, SEXP xmaxSEXP, SEXP nptsSEXP, SEXP optimize_pwmSEXP, SEXP optimize_spatSEXP, SEXP symmetrize_spatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type is_train_logical(is_train_logicalSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type motif(motifSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_min(spat_minSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_max(spat_maxSEXP);
    Rcpp::traits::input_parameter< const float& >::type min_nuc_prob(min_nuc_probSEXP);
    Rcpp::traits::input_parameter< const int& >::type spat_bin(spat_binSEXP);
    Rcpp::traits::input_parameter< const float& >::type improve_epsilon(improve_epsilonSEXP);
    Rcpp::traits::input_parameter< const bool& >::type is_bidirect(is_bidirectSEXP);
    Rcpp::traits::input_parameter< const float& >::type unif_prior(unif_priorSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type score_metric(score_metricSEXP);
    Rcpp::traits::input_parameter< const int& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type pssm_mat(pssm_matSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericVector>& >::type spat_factor(spat_factorSEXP);
    Rcpp::traits::input_parameter< const float& >::type consensus_single_thresh(consensus_single_threshSEXP);
    Rcpp::traits::input_parameter< const float& >::type consensus_double_thresh(consensus_double_threshSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_folds(num_foldsSEXP);
    Rcpp::traits::input_parameter< const float& >::type energy_epsilon(energy_epsilonSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_energy(log_energySEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::Function> >::type energy_func(energy_funcSEXP);
    Rcpp::traits::input_parameter< const float& >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< const float& >::type xmax(xmaxSEXP);
    Rcpp::traits::input_parameter< const int& >::type npts(nptsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type optimize_pwm(optimize_pwmSEXP);
    Rcpp::traits::input_parameter< const bool& >::type optimize_spat(optimize_spatSEXP);
    Rcpp::traits::input_parameter< const bool& >::type symmetrize_spat(symmetrize_spatSEXP);
    rcpp_result_gen = Rcpp::wrap(regress_pwm_cpp(sequences, response, is_train_logical, motif, spat_min, spat_max, min_nuc_prob, spat_bin, improve_epsilon, is_bidirect, unif_prior, score_metric, verbose, seed, pssm_mat, spat_factor, consensus_single_thresh, consensus_double_thresh, num_folds, energy_epsilon, log_energy, energy_func, xmin, xmax, npts, optimize_pwm, optimize_spat, symmetrize_spat));
    return rcpp_result_gen;
END_RCPP
}
// screen_kmers_cpp
Rcpp::DataFrame screen_kmers_cpp(const Rcpp::StringVector& sequences, const Rcpp::DataFrame& response, const Rcpp::LogicalVector& is_train_logical, const int& L, const int& from_range, const int& to_range, const float& min_cor, const int& min_gap, const int& max_gap, const int& n_in_train, const int& seed, const bool& verbose);
RcppExport SEXP _prego_screen_kmers_cpp(SEXP sequencesSEXP, SEXP responseSEXP, SEXP is_train_logicalSEXP, SEXP LSEXP, SEXP from_rangeSEXP, SEXP to_rangeSEXP, SEXP min_corSEXP, SEXP min_gapSEXP, SEXP max_gapSEXP, SEXP n_in_trainSEXP, SEXP seedSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type is_train_logical(is_train_logicalSEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int& >::type from_range(from_rangeSEXP);
    Rcpp::traits::input_parameter< const int& >::type to_range(to_rangeSEXP);
    Rcpp::traits::input_parameter< const float& >::type min_cor(min_corSEXP);
    Rcpp::traits::input_parameter< const int& >::type min_gap(min_gapSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_gap(max_gapSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_in_train(n_in_trainSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(screen_kmers_cpp(sequences, response, is_train_logical, L, from_range, to_range, min_cor, min_gap, max_gap, n_in_train, seed, verbose));
    return rcpp_result_gen;
END_RCPP
}
// interpolateFunction
Rcpp::NumericVector interpolateFunction(Rcpp::Function func, float xmin, float xmax, int npts, Rcpp::NumericVector x);
RcppExport SEXP _prego_interpolateFunction(SEXP funcSEXP, SEXP xminSEXP, SEXP xmaxSEXP, SEXP nptsSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type func(funcSEXP);
    Rcpp::traits::input_parameter< float >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< float >::type xmax(xmaxSEXP);
    Rcpp::traits::input_parameter< int >::type npts(nptsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(interpolateFunction(func, xmin, xmax, npts, x));
    return rcpp_result_gen;
END_RCPP
}
// n_nuc_distribution
Rcpp::DataFrame n_nuc_distribution(Rcpp::StringVector sequences, int n, int size);
RcppExport SEXP _prego_n_nuc_distribution(SEXP sequencesSEXP, SEXP nSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type sequences(sequencesSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(n_nuc_distribution(sequences, n, size));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_prego_kmer_matrix_cpp", (DL_FUNC) &_prego_kmer_matrix_cpp, 7},
    {"_prego_get_consensus_cpp", (DL_FUNC) &_prego_get_consensus_cpp, 3},
    {"_prego_compute_pwm_cpp", (DL_FUNC) &_prego_compute_pwm_cpp, 8},
    {"_prego_compute_local_pwm_cpp", (DL_FUNC) &_prego_compute_local_pwm_cpp, 7},
    {"_prego_mask_sequences_cpp", (DL_FUNC) &_prego_mask_sequences_cpp, 9},
    {"_prego_regress_pwm_cpp", (DL_FUNC) &_prego_regress_pwm_cpp, 28},
    {"_prego_screen_kmers_cpp", (DL_FUNC) &_prego_screen_kmers_cpp, 12},
    {"_prego_interpolateFunction", (DL_FUNC) &_prego_interpolateFunction, 5},
    {"_prego_n_nuc_distribution", (DL_FUNC) &_prego_n_nuc_distribution, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_prego(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
