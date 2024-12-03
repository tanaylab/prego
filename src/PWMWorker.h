#ifndef PWMWORKER_H
#define PWMWORKER_H

#include <Rcpp.h>
#include <RcppParallel.h>
#include "logSumExp.h"
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

class PWMWorker : public Worker {
private:
    const RMatrix<double> sequences;
    const RMatrix<double> combined_pwm;
    const int D_min;
    const Rcpp::IntegerVector motif_lengths;
    const bool bidirect;
    RMatrix<double> output;
    
    const int seq_length;
    const int motif_length;
    const int num_motifs;
    const int max_windows;
    
    // Spatial binning parameters
    const bool use_spatial;
    const RMatrix<double> spat_factors; 
    const int spat_bin_size;

    // Helper functions
    std::pair<int, int> calculate_window_counts(int seq_length, 
                                              int motif_length, 
                                              int D_min) const;

    void fill_full_windows(std::vector<double>& windows_mat,
                          size_t sequence_idx,
                          int n_full_windows) const;

    void fill_partial_windows(std::vector<double>& windows_mat,
                            size_t sequence_idx,
                            int n_full_windows,
                            int n_partial_windows,
                            int window_size) const;

    void compute_matrix_scores(std::vector<double>& window_scores,
                             const std::vector<double>& windows_mat,
                             int total_windows) const;

    void process_sequence(size_t sequence_idx,
                         std::vector<double>& windows_mat,
                         std::vector<double>& window_scores);

    void process_window_scores(RMatrix<double>& output,
                             const std::vector<double>& window_scores,
                             size_t sequence_idx,
                             int total_windows) const;

public:
    PWMWorker(const Rcpp::NumericMatrix sequences,
              const Rcpp::NumericMatrix combined_pwm,
              const int D_min,
              const Rcpp::IntegerVector motif_lengths,
              const bool bidirect,
              Rcpp::NumericMatrix output,
              const Rcpp::NumericMatrix spat_factors = Rcpp::NumericMatrix(0),
              const int spat_bin_size = 1);

    void operator()(std::size_t begin, std::size_t end);
};


#endif // PWMWORKER_H