#include <Rcpp.h>
#include <RcppParallel.h>
#include "PWMWorker.h"
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;


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