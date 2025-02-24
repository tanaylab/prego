#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// Structure to hold a PSSM
struct PSSM {
    std::vector<double> data; // Flattened matrix: [A1,C1,G1,T1,A2,C2,G2,T2,...]
    int nrows;
    int ncols; // Always 4 for DNA

    PSSM(const NumericMatrix &mat) : data(mat.nrow() * 4), nrows(mat.nrow()), ncols(4) {
        // Copy and flatten the matrix
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                data[i * ncols + j] = mat(i, j);
            }
        }
    }
};

// Function to compute correlation between two vectors
double compute_correlation(const std::vector<double> &x, const std::vector<double> &y,
                           const std::string &method) {
    int n = x.size();
    if (method == "pearson") {
        double sum_x = 0, sum_y = 0, sum_xy = 0;
        double sum_x2 = 0, sum_y2 = 0;

        for (int i = 0; i < n; i++) {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += x[i] * y[i];
            sum_x2 += x[i] * x[i];
            sum_y2 += y[i] * y[i];
        }

        double numerator = n * sum_xy - sum_x * sum_y;
        double denominator = std::sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y));

        return (denominator > 0) ? numerator / denominator : 0.0;
    } else { // Spearman
        // Create vectors of pairs (value, index)
        std::vector<std::pair<double, int>> x_pairs, y_pairs;
        for (int i = 0; i < n; i++) {
            x_pairs.push_back({x[i], i});
            y_pairs.push_back({y[i], i});
        }

        // Sort pairs by value
        std::sort(x_pairs.begin(), x_pairs.end());
        std::sort(y_pairs.begin(), y_pairs.end());

        // Assign ranks, handling ties properly
        std::vector<double> x_rank(n), y_rank(n);

        // Process x ranks
        for (int i = 0; i < n;) {
            int start = i;
            // Find range of tied values
            while (i < n - 1 && x_pairs[i].first == x_pairs[i + 1].first) {
                i++;
            }
            // Calculate average rank for tied values
            double avg_rank = (start + i) / 2.0 + 1; // +1 because ranks start at 1
            // Assign average rank to all tied values
            for (int j = start; j <= i; j++) {
                x_rank[x_pairs[j].second] = avg_rank;
            }
            i++;
        }

        // Process y ranks
        for (int i = 0; i < n;) {
            int start = i;
            // Find range of tied values
            while (i < n - 1 && y_pairs[i].first == y_pairs[i + 1].first) {
                i++;
            }
            // Calculate average rank for tied values
            double avg_rank = (start + i) / 2.0 + 1;
            // Assign average rank to all tied values
            for (int j = start; j <= i; j++) {
                y_rank[y_pairs[j].second] = avg_rank;
            }
            i++;
        }

        // Compute Pearson correlation on the ranks
        return compute_correlation(x_rank, y_rank, "pearson");
    }
}

// Function to compute KL divergence between two vectors
// The vectors should contain probability distributions
double compute_kl_divergence(const std::vector<double> &p, const std::vector<double> &q) {
    double kl_div = 0.0;
    for (size_t i = 0; i < p.size(); i++) {
        // Avoid log(0) by ensuring both p and q have positive values
        if (p[i] > 0 && q[i] > 0) {
            kl_div += p[i] * std::log(p[i] / q[i]);
        }
    }
    return kl_div;
}

// Function to compute maximum correlation between two PSSMs with sliding window
double compute_max_pssm_correlation(const PSSM &pssm1, const PSSM &pssm2, const std::string &method,
                                    double prior) {
    int window_size = std::min(pssm1.nrows, pssm2.nrows);
    int max_pos = std::max(pssm1.nrows, pssm2.nrows);

    const PSSM &pssm_s = (pssm1.nrows <= pssm2.nrows) ? pssm1 : pssm2;
    const PSSM &pssm_l = (pssm1.nrows <= pssm2.nrows) ? pssm2 : pssm1;

    double best_score;
    // For Pearson and Spearman correlations, higher is better
    // For KL divergence, lower is better
    if (method == "kl") {
        best_score = std::numeric_limits<double>::max(); // Initialize to a high value
    } else {
        best_score = -1.0; // Initialize to a low value for correlations
    }

    std::vector<double> vec_s(window_size * 4), vec_l(window_size * 4);

    // Apply prior and normalize
    auto normalize_with_prior = [prior](std::vector<double> &vec) {
        for (size_t i = 0; i < vec.size(); i += 4) {
            // Add prior to each value
            for (int j = 0; j < 4; j++) {
                vec[i + j] += prior;
            }

            // Calculate sum in a single pass
            double sum = vec[i] + vec[i + 1] + vec[i + 2] + vec[i + 3];

            // Normalize
            for (int j = 0; j < 4; j++) {
                vec[i + j] /= sum;
            }
        }
    };

    // Copy shorter PSSM
    std::copy(pssm_s.data.begin(), pssm_s.data.end(), vec_s.begin());
    normalize_with_prior(vec_s);

    // Slide window and compute score
    for (int start = 0; start <= max_pos - window_size; start++) {
        // Copy window from longer PSSM
        for (int i = 0; i < window_size; i++) {
            for (int j = 0; j < 4; j++) {
                vec_l[i * 4 + j] = pssm_l.data[(start + i) * 4 + j];
            }
        }
        normalize_with_prior(vec_l);

        double score;
        if (method == "kl") {
            // For KL divergence, compute both directions and take average
            double kl_s_l = compute_kl_divergence(vec_s, vec_l);
            double kl_l_s = compute_kl_divergence(vec_l, vec_s);
            score = (kl_s_l + kl_l_s) / 2.0; // Symmetric KL divergence
        } else {
            score = compute_correlation(vec_s, vec_l, method);
        }

        // Update best score (min for KL, max for correlations)
        if (method == "kl") {
            best_score = std::min(best_score, score);
        } else {
            best_score = std::max(best_score, score);
        }
    }

    return best_score;
}

// Parallel worker for computing correlation matrix
struct PSSMCorrelationWorker : public Worker {
    // Input data
    const std::vector<PSSM> &pssms;
    const std::string method;
    const double prior;
    const int n_motifs;

    // Output correlation matrix
    std::vector<double> &correlations;

    PSSMCorrelationWorker(const std::vector<PSSM> &pssms, std::vector<double> &correlations,
                          const std::string &method, double prior)
        : pssms(pssms), method(method), prior(prior), n_motifs(pssms.size()),
          correlations(correlations) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t k = begin; k < end; k++) {
            int i = k / n_motifs;
            int j = k % n_motifs;
            if (i < j) { // Only compute upper triangle
                double corr = compute_max_pssm_correlation(pssms[i], pssms[j], method, prior);
                correlations[i * n_motifs + j] = corr;
                correlations[j * n_motifs + i] = corr; // Fill lower triangle
            }
        }
    }
};

// [[Rcpp::export]]
NumericMatrix pssm_dataset_cor_parallel(const List &pssm_list,
                                        const std::string &method = "spearman",
                                        double prior = 0.01) {
    int n_motifs = pssm_list.length();

    // Convert list of PSSMs to vector
    std::vector<PSSM> pssms;
    pssms.reserve(n_motifs);
    for (int i = 0; i < n_motifs; i++) {
        NumericMatrix mat = as<NumericMatrix>(pssm_list[i]);
        pssms.emplace_back(mat);
    }

    // Initialize correlation matrix
    std::vector<double> correlations(n_motifs * n_motifs, 0.0);

    // Set diagonal to 1
    for (int i = 0; i < n_motifs; i++) {
        correlations[i * n_motifs + i] = 1.0;
    }

    // Create parallel worker
    PSSMCorrelationWorker worker(pssms, correlations, method, prior);

    // Run parallel computation
    parallelFor(0, n_motifs * n_motifs, worker);

    // Convert to R matrix
    NumericMatrix result(n_motifs, n_motifs);
    for (int i = 0; i < n_motifs; i++) {
        for (int j = 0; j < n_motifs; j++) {
            result(i, j) = correlations[i * n_motifs + j];
        }
    }

    return result;
}