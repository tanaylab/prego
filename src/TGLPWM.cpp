#include <Rcpp.h>
#include <RcppParallel.h>
#include "logSumExp.h"
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// Forward declarations
class PWMWorker;

// BLAS routine declaration
extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA, const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA, double* C, const int* LDC);
}

// Window calculation functions
std::pair<int, int> calculate_window_counts(int seq_length, 
                                          int motif_length, 
                                          int D_min) {
    int n_full_windows = (seq_length >= motif_length) ? 
        seq_length - motif_length + 1 : 0;
    int n_partial_windows = (seq_length > D_min) ? 
        std::min(motif_length - D_min, 
                seq_length - std::max(motif_length, D_min) + 1) : 0;
    
    return {n_full_windows, n_partial_windows};
}

void fill_full_windows(std::vector<double>& windows_mat,
                      const RMatrix<double>& sequences,
                      size_t sequence_idx,
                      int n_full_windows,
                      int motif_length) {
    for (int i = 0; i < n_full_windows; i++) {
        for (int d = 0; d < motif_length; d++) {
            for (int j = 0; j < 4; j++) {
                size_t win_idx = d * 4 + j;
                size_t seq_idx = (i + d) * 4 + j;
                if (win_idx < windows_mat.size() && seq_idx < sequences.ncol()) {
                    windows_mat[win_idx + i * (4 * motif_length)] = 
                        sequences(sequence_idx, seq_idx);
                }
            }
        }
    }
}

void fill_partial_windows(std::vector<double>& windows_mat,
                         const RMatrix<double>& sequences,
                         size_t sequence_idx,
                         int n_full_windows,
                         int n_partial_windows,
                         int seq_length,
                         int motif_length,
                         int D_min) {
    if (n_partial_windows <= 0) return;

    int col_idx = n_full_windows;
    for (int window_size = motif_length - 1; 
         window_size >= motif_length - n_partial_windows; 
         window_size--) {
        if (window_size >= D_min) {
            int pos = seq_length - window_size;
            for (int d = 0; d < window_size; d++) {
                for (int j = 0; j < 4; j++) {
                    size_t win_idx = d * 4 + j;
                    size_t seq_idx = (pos + d) * 4 + j;
                    if (win_idx < windows_mat.size() && seq_idx < sequences.ncol()) {
                        windows_mat[win_idx + col_idx * (4 * motif_length)] = 
                            sequences(sequence_idx, seq_idx);
                    }
                }
            }
            col_idx++;
        }
    }
}

// Score processing functions
void compute_matrix_scores(std::vector<double>& window_scores,
                         const std::vector<double>& windows_mat,
                         const RMatrix<double>& combined_pwm,
                         int total_windows,
                         int motif_length) {
    char trans1 = 'T';
    char trans2 = 'N';
    double alpha = 1.0;
    double beta = 0.0;
    int m = total_windows;
    int n = combined_pwm.ncol();
    int k = 4 * motif_length;

    if (m > 0 && n > 0 && k > 0) {
        dgemm_(&trans1, &trans2,
               &m, &n, &k,
               &alpha,
               windows_mat.data(), &k,
               combined_pwm.begin(), &k,
               &beta,
               window_scores.data(), &m);
    }
}


class PWMWorker : public Worker {
private:
    const RMatrix<double> sequences;
    const RMatrix<double> combined_pwm;
    const int D_min;
    const IntegerVector motif_lengths;
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

    void process_sequence(size_t sequence_idx,
                        std::vector<double>& windows_mat,
                        std::vector<double>& window_scores) {
        // Handle short sequences
        if (seq_length < D_min) {
            std::fill(output.row(sequence_idx).begin(), 
                     output.row(sequence_idx).end(), 
                     R_NegInf);
            return;
        }

        // Calculate window counts
        auto [n_full_windows, n_partial_windows] = 
            calculate_window_counts(seq_length, motif_length, D_min);
        int total_windows = n_full_windows + n_partial_windows;

        if (total_windows <= 0) {
            std::fill(output.row(sequence_idx).begin(), 
                     output.row(sequence_idx).end(), 
                     R_NegInf);
            return;
        }

        // Clear working memory
        std::fill(windows_mat.begin(), windows_mat.end(), 0.0);
        std::fill(window_scores.begin(), window_scores.end(), 0.0);

        // Fill windows matrix
        fill_full_windows(windows_mat, sequences, sequence_idx,
                         n_full_windows, motif_length);
        fill_partial_windows(windows_mat, sequences, sequence_idx,
                           n_full_windows, n_partial_windows,
                           seq_length, motif_length, D_min);

        // Compute and process scores
        compute_matrix_scores(window_scores, windows_mat, combined_pwm,
                            total_windows, motif_length);
        process_window_scores(output, window_scores, sequence_idx,
                            num_motifs, bidirect, total_windows,
                            seq_length, motif_lengths);
    }

public:
    PWMWorker(const NumericMatrix sequences,
              const NumericMatrix combined_pwm,
              const int D_min,
              const IntegerVector motif_lengths,
              const bool bidirect,
              NumericMatrix output,
              const NumericMatrix spat_factors = NumericMatrix(0),
              const int spat_bin_size = 1) 
        : sequences(sequences), combined_pwm(combined_pwm), D_min(D_min), 
          motif_lengths(motif_lengths), bidirect(bidirect),
          output(output),
          seq_length(sequences.ncol() / 4),
          motif_length(combined_pwm.nrow() / 4),
          num_motifs(bidirect ? combined_pwm.ncol() / 2 : combined_pwm.ncol()),
          max_windows(std::max(1, seq_length - std::min(motif_length, D_min) + 1)),
          use_spatial(spat_factors.ncol() > 1 || spat_bin_size > 1),
          spat_factors(spat_factors),
          spat_bin_size(spat_bin_size)
    {
        // if (use_spatial) {
        //     Rcpp::Rcout << "Using spatial factors: motifs = " << spat_factors.nrow() 
        //                << ", bins = " << spat_factors.ncol() 
        //                << ", bin size = " << spat_bin_size << "\n";
        // } else {
        //     Rcpp::Rcout << "Not using spatial factors\n";
        // }
    }

    void operator()(std::size_t begin, std::size_t end) {
        std::vector<double> windows_mat(4 * motif_length * max_windows, 0.0);
        std::vector<double> window_scores(max_windows * (bidirect ? 2*num_motifs : num_motifs), 0.0);
        
        for (std::size_t i = begin; i < end; i++) {
            process_sequence(i, windows_mat, window_scores);
        }
    }

private:
    void process_window_scores(RMatrix<double>& output,
                             const std::vector<double>& window_scores,
                             size_t sequence_idx,
                             int num_motifs,
                             bool bidirect,
                             int total_windows,
                             int seq_length,
                             const IntegerVector& motif_lengths) {
        for (int j = 0; j < num_motifs; j++) {
            int motif_length = motif_lengths[j];
            int relevant_windows = seq_length - motif_length + 1;

            if (!use_spatial) {                
                if (bidirect) {
                    const double* fwd_scores = &window_scores[j * total_windows];
                    const double* rev_scores = &window_scores[(j + num_motifs) * total_windows];
                    double fwd = log_sum_exp(fwd_scores, relevant_windows);
                    double rev = log_sum_exp(rev_scores, relevant_windows);
                    log_sum_log(fwd, rev);
                    output(sequence_idx, j) = fwd;
                } else {
                    const double* scores = &window_scores[j * total_windows];
                    output(sequence_idx, j) = log_sum_exp(scores, relevant_windows);
                }
                continue;
            }
            
            if (bidirect) {
                // Forward direction scores
                double fwd = -R_PosInf;
                for (int w = 0; w < relevant_windows; w++) {
                    double score = window_scores[j * total_windows + w];
                    // Calculate position relative to sequence start
                    size_t absolute_pos = w;  // w is already relative to min_range
                    size_t bin = absolute_pos / spat_bin_size;
                    
                    if (bin < static_cast<size_t>(spat_factors.ncol())) {
                        score += log(spat_factors(j, bin));
                        log_sum_log(fwd, score);
                    }
                }

                // Reverse direction scores  
                double rev = -R_PosInf;
                for (int w = 0; w < relevant_windows; w++) {
                    double score = window_scores[(j + num_motifs) * total_windows + w];
                    // Use same position calculation for reverse complement
                    size_t absolute_pos = w;
                    size_t bin = absolute_pos / spat_bin_size;
                    
                    if (bin < static_cast<size_t>(spat_factors.ncol())) {
                        score += log(spat_factors(j, bin));
                        log_sum_log(rev, score);
                    }
                }

                log_sum_log(fwd, rev);
                output(sequence_idx, j) = fwd;
            } else {
                double total = -R_PosInf;
                for (int w = 0; w < relevant_windows; w++) {
                    double score = window_scores[j * total_windows + w];
                    size_t absolute_pos = w;
                    size_t bin = absolute_pos / spat_bin_size;
                    
                    if (bin < static_cast<size_t>(spat_factors.ncol())) {
                        score += log(spat_factors(j, bin));
                        log_sum_log(total, score);
                    }
                }
                output(sequence_idx, j) = total;
            }
        }
    }
};


void validate_inputs(const NumericMatrix& sequences,
                    const NumericMatrix& pwm,
                    const IntegerVector& motif_lengths,
                    const int D_min) {
    if (sequences.ncol() % 4 != 0) {
        stop("Number of columns in sequences must be divisible by 4");
    }
    if (pwm.nrow() % 4 != 0) {
        stop("Number of rows in PWM must be divisible by 4");
    }
    if (motif_lengths.length() != pwm.ncol()) {
        stop("Length of motif_lengths must match number of PWM columns");
    }
    if (D_min < 1) {
        stop("D_min must be positive");
    }
    if (sequences.nrow() < 1) {
        stop("No sequences provided");
    }
}

// [[Rcpp::export]]
NumericMatrix calc_seq_pwm_parallel_cpp(const NumericMatrix& sequences,
                                      const NumericMatrix& pwm,
                                      const NumericMatrix& pwm_rc,  
                                      const IntegerVector& motif_lengths,
                                      const int D_min = 1,
                                      const bool bidirect = false,
                                      const NumericMatrix& spat_factors = NumericMatrix(0),  // Changed to matrix
                                      const int spat_bin_size = 1) {
    // Validate inputs
    validate_inputs(sequences, pwm, motif_lengths, D_min);
    
    // Validate spatial factors if provided
    if (spat_factors.nrow() > 0 && spat_factors.nrow() != pwm.ncol()) {
        stop("Number of rows in spatial factors matrix must match number of motifs");
    }
    
    // Prepare combined PWM matrix
    NumericMatrix combined_pwm;
    if (bidirect) {
        combined_pwm = NumericMatrix(pwm.nrow(), 2 * pwm.ncol());
        std::copy(pwm.begin(), pwm.end(), combined_pwm.begin());
        std::copy(pwm_rc.begin(), pwm_rc.end(), 
                 combined_pwm.begin() + pwm.nrow() * pwm.ncol());
    } else {
        combined_pwm = pwm;
    }
    
    // Create output matrix
    NumericMatrix output(sequences.nrow(), pwm.ncol());
    
    try {
        PWMWorker worker(sequences, combined_pwm, D_min, motif_lengths, 
                        bidirect, output, spat_factors, spat_bin_size);
        parallelFor(0, sequences.nrow(), worker);
    } catch (std::exception& e) {
        stop("Error in parallel processing: %s", e.what());
    }
  
    return output;
}