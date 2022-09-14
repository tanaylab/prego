// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_pwm_cpp
Rcpp::NumericVector compute_pwm_cpp(const Rcpp::StringVector& sequences, const Rcpp::NumericMatrix& pssm_mat, const bool& is_bidirect, const int& spat_min, const int& spat_max, const Rcpp::NumericVector& spat_factor, const int& bin_size);
RcppExport SEXP _prego_compute_pwm_cpp(SEXP sequencesSEXP, SEXP pssm_matSEXP, SEXP is_bidirectSEXP, SEXP spat_minSEXP, SEXP spat_maxSEXP, SEXP spat_factorSEXP, SEXP bin_sizeSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(compute_pwm_cpp(sequences, pssm_mat, is_bidirect, spat_min, spat_max, spat_factor, bin_size));
    return rcpp_result_gen;
END_RCPP
}
// regress_pwm_cpp
Rcpp::List regress_pwm_cpp(const Rcpp::StringVector& sequences, const Rcpp::DataFrame& response, const Rcpp::LogicalVector& is_train_logical, const std::string& motif, const int& spat_min, const int& spat_max, const float& min_nuc_prob, const int& spat_bin, const float& improve_epsilon, const bool& is_bidirect, const float& unif_prior, const std::string& score_metric, const int& verbose, const int& seed, const Rcpp::NumericMatrix& pssm_mat);
RcppExport SEXP _prego_regress_pwm_cpp(SEXP sequencesSEXP, SEXP responseSEXP, SEXP is_train_logicalSEXP, SEXP motifSEXP, SEXP spat_minSEXP, SEXP spat_maxSEXP, SEXP min_nuc_probSEXP, SEXP spat_binSEXP, SEXP improve_epsilonSEXP, SEXP is_bidirectSEXP, SEXP unif_priorSEXP, SEXP score_metricSEXP, SEXP verboseSEXP, SEXP seedSEXP, SEXP pssm_matSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(regress_pwm_cpp(sequences, response, is_train_logical, motif, spat_min, spat_max, min_nuc_prob, spat_bin, improve_epsilon, is_bidirect, unif_prior, score_metric, verbose, seed, pssm_mat));
    return rcpp_result_gen;
END_RCPP
}
// screen_kmers_cpp
Rcpp::DataFrame screen_kmers_cpp(const Rcpp::StringVector& sequences, const Rcpp::DataFrame& response, const Rcpp::LogicalVector& is_train_logical, const int& L, const int& from_range, const int& to_range, const float& min_cor, const int& min_n, const int& min_gap, const int& max_gap, const int& n_in_train, const int& seed, const bool& verbose);
RcppExport SEXP _prego_screen_kmers_cpp(SEXP sequencesSEXP, SEXP responseSEXP, SEXP is_train_logicalSEXP, SEXP LSEXP, SEXP from_rangeSEXP, SEXP to_rangeSEXP, SEXP min_corSEXP, SEXP min_nSEXP, SEXP min_gapSEXP, SEXP max_gapSEXP, SEXP n_in_trainSEXP, SEXP seedSEXP, SEXP verboseSEXP) {
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
    Rcpp::traits::input_parameter< const int& >::type min_n(min_nSEXP);
    Rcpp::traits::input_parameter< const int& >::type min_gap(min_gapSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_gap(max_gapSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_in_train(n_in_trainSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(screen_kmers_cpp(sequences, response, is_train_logical, L, from_range, to_range, min_cor, min_n, min_gap, max_gap, n_in_train, seed, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_prego_compute_pwm_cpp", (DL_FUNC) &_prego_compute_pwm_cpp, 7},
    {"_prego_regress_pwm_cpp", (DL_FUNC) &_prego_regress_pwm_cpp, 15},
    {"_prego_screen_kmers_cpp", (DL_FUNC) &_prego_screen_kmers_cpp, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_prego(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
