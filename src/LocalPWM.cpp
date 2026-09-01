#include <Rcpp.h>
#include <RcppParallel.h>
#include "blas_pin.h"
#include "logSumExp.h"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// Expected local PWM score of a motif against a per-position base frequency
// matrix, for every motif and every start position.
//
// For a motif of length L placed at start j, with motif column p_l and
// frequency column q_{j+l}, there are two ways to combine a pair of
// distributions, and they differ only in where the log sits:
//
//   multiply = false   score(j) = sum_l  q_{j+l} . log p_l    = E[log P]
//   multiply = true    score(j) = sum_l  log(q_{j+l} . p_l)   = log E[P]
//
// The first is linear in the frequencies, so it is the exact mean score over
// the ensemble whatever the joint distribution of positions is. The second
// assumes the positions are independent, but has a motif-independent floor on
// a flat ensemble (L*log(0.25)), which makes motifs comparable to each other.
// Both reduce to compute_local_pwm() when the frequency matrix is one-hot.
//
// Either way the shape is: one dgemm, an optional elementwise log, then a
// banded sum. The log between the two linear steps is why this cannot collapse
// into a single matrix multiplication.

// Cap on the m x (D*block) intermediate, in doubles (~20MB). Blocking over
// motifs keeps a full motif database off the >1GB an unblocked pass would need.
static const size_t BLOCK_TARGET_ELEMS = 2500000;
static const int MAX_MOTIF_BLOCK = 64;
// Floor for the shrink below. A very skinny dgemm, and a fold whose innermost
// axis is this short, cost more than the extra threads win back.
static const int MIN_MOTIF_BLOCK = 8;

struct FreqTask {
    int mat;
    int motif_from;  // index into the length-sorted motif order
    int n_block;
    int block_len;   // longest motif in this block; the dgemm runs at this width
    bool direct;     // block is contiguous at full width: feed the PWM to dgemm as-is
};

class FreqLocalPWMWorker : public Worker {
  private:
    const std::vector<const double *> &q;   // each m x 4, column major
    const std::vector<int> &q_len;          // m per matrix
    const std::vector<double *> &out;       // each n_motifs x m, column major
    const double *pwm;                      // 4D x n_motifs, log scale
    const double *pwm_rc;
    const int *motif_lengths;
    const int *order;                       // motifs sorted by length
    const int n_motifs;
    const int D;
    const bool multiply;
    const bool bidirect;
    const std::vector<FreqTask> &tasks;
    const int max_m;
    const int max_block;

    // One strand: dgemm the frequencies against a block of motif columns, then
    // fold the D columns of each motif back onto its start positions.
    void score_block(const double *Q, int m, const double *P, const FreqTask &task,
                     std::vector<double> &S, std::vector<double> &pbuf,
                     std::vector<double> &acc) const {
        const int B = task.n_block;
        const int Db = task.block_len;
        const double *A;

        if (task.direct) {
            // Uniform-length database: the sort is a no-op, the block is already
            // contiguous and full width, so the PWM is its own dgemm operand.
            A = P + (size_t)task.motif_from * 4 * D;
        } else {
            // Gather the block's motifs. A motif's first Db positions are contiguous
            // in the PWM, so each is one copy; the length-sorted order means Db is
            // the longest motif in this block rather than the longest in the whole
            // database, which is what keeps the dgemm off the padded columns.
            for (int b = 0; b < B; b++) {
                const double *src = P + (size_t)order[task.motif_from + b] * 4 * D;
                double *dst = pbuf.data() + (size_t)b * 4 * Db;
                if (multiply) {
                    for (int k = 0; k < 4 * Db; k++) {
                        dst[k] = std::exp(src[k]);
                    }
                } else {
                    std::copy(src, src + 4 * Db, dst);
                }
            }
            A = pbuf.data();
        }

        char no_trans = 'N';
        double alpha = 1.0, beta = 0.0;
        int mm = m, nn = Db * B, kk = 4;
        dgemm_(&no_trans, &no_trans, &mm, &nn, &kk, &alpha, Q, &mm, A, &kk, &beta,
               S.data(), &mm);

        std::fill(acc.begin(), acc.begin() + (size_t)m * B, 0.0);
        for (int b = 0; b < B; b++) {
            const int len = motif_lengths[order[task.motif_from + b]];
            const int last = m - len;
            double *a = acc.data() + (size_t)b * m;
            for (int l = 0; l < len; l++) {
                const double *col = S.data() + (size_t)(b * Db + l) * m;
                // Every (motif, offset) column is read exactly once, so the log
                // for multiply mode is fused in here rather than run over all of S.
                if (multiply) {
                    for (int j = 0; j <= last; j++) {
                        a[j] += std::log(col[j + l]);
                    }
                } else {
                    for (int j = 0; j <= last; j++) {
                        a[j] += col[j + l];
                    }
                }
            }
        }
    }

    void run_task(const FreqTask &task, std::vector<double> &S, std::vector<double> &pbuf,
                  std::vector<double> &acc, std::vector<double> &acc_rc) const {
        const int m = q_len[task.mat];
        const double *Q = q[task.mat];
        double *O = out[task.mat];

        score_block(Q, m, pwm, task, S, pbuf, acc);
        if (bidirect) {
            score_block(Q, m, pwm_rc, task, S, pbuf, acc_rc);
        }

        for (int b = 0; b < task.n_block; b++) {
            const int i = order[task.motif_from + b];
            const int last = m - motif_lengths[i];
            for (int j = 0; j <= last; j++) {
                double v = acc[(size_t)b * m + j];
                if (bidirect) {
                    log_sum_log(v, acc_rc[(size_t)b * m + j]);
                }
                O[(size_t)j * n_motifs + i] = v;
            }
            for (int j = last + 1; j < m; j++) {
                O[(size_t)j * n_motifs + i] = NA_REAL;
            }
        }
    }

  public:
    FreqLocalPWMWorker(const std::vector<const double *> &q, const std::vector<int> &q_len,
                       const std::vector<double *> &out, const double *pwm,
                       const double *pwm_rc, const int *motif_lengths, const int *order,
                       int n_motifs, int D,
                       bool multiply, bool bidirect, const std::vector<FreqTask> &tasks,
                       int max_m, int max_block)
        : q(q), q_len(q_len), out(out), pwm(pwm), pwm_rc(pwm_rc),
          motif_lengths(motif_lengths), order(order), n_motifs(n_motifs), D(D), multiply(multiply),
          bidirect(bidirect), tasks(tasks), max_m(max_m), max_block(max_block) {}

    void operator()(std::size_t begin, std::size_t end) {
        pin_worker_blas_to_one_thread();

        std::vector<double> S((size_t)max_m * D * max_block);
        std::vector<double> acc((size_t)max_m * max_block);
        std::vector<double> acc_rc(bidirect ? (size_t)max_m * max_block : 0);
        std::vector<double> pbuf((size_t)4 * D * max_block);

        for (std::size_t t = begin; t < end; t++) {
            run_task(tasks[t], S, pbuf, acc, acc_rc);
        }
    }
};

// [[Rcpp::export]]
Rcpp::List calc_freq_local_pwm_cpp(const Rcpp::List &freqs, const Rcpp::NumericMatrix &pwm,
                                   const Rcpp::NumericMatrix &pwm_rc,
                                   const Rcpp::IntegerVector &motif_lengths,
                                   const bool multiply = true, const bool bidirect = true,
                                   const int n_threads = 1) {
    if (pwm.nrow() < 4 || pwm.nrow() % 4 != 0) {
        stop("Number of rows in PWM must be a positive multiple of 4");
    }
    if (pwm.ncol() != motif_lengths.length()) {
        stop("Length of motif_lengths must match number of PWM columns");
    }
    if (pwm_rc.nrow() != pwm.nrow() || pwm_rc.ncol() != pwm.ncol()) {
        stop("Reverse complement PWM must have the same dimensions as the PWM");
    }

    const int n_motifs = pwm.ncol();
    const int D = pwm.nrow() / 4;
    const int n_mats = freqs.size();

    int max_len = 0;
    for (int i = 0; i < n_motifs; i++) {
        if (motif_lengths[i] < 1 || motif_lengths[i] > D) {
            stop("Motif lengths must be between 1 and the number of PWM positions");
        }
        max_len = std::max(max_len, motif_lengths[i]);
    }

    // Hold both the (possibly coerced) inputs and the freshly allocated outputs
    // in lists, and take the pointers back out of those lists - a NumericMatrix
    // built by coercing a list element is owned only by the local variable, so
    // its buffer would be freed before the workers ran.
    Rcpp::List q_list(n_mats);
    Rcpp::List out_list(n_mats);
    std::vector<const double *> q(n_mats);
    std::vector<double *> out(n_mats);
    std::vector<int> q_len(n_mats);
    int max_m = 0;

    for (int k = 0; k < n_mats; k++) {
        q_list[k] = Rcpp::NumericMatrix(freqs[k]);
        Rcpp::NumericMatrix qk = q_list[k];
        if (qk.ncol() != 4) {
            stop("Each frequency matrix must have 4 columns");
        }
        if (qk.nrow() < max_len) {
            stop("Frequency matrix %d is shorter than the longest motif", k + 1);
        }
        out_list[k] = Rcpp::NumericMatrix(n_motifs, qk.nrow());
        Rcpp::NumericMatrix ok = out_list[k];
        q[k] = qk.begin();
        out[k] = ok.begin();
        q_len[k] = qk.nrow();
        max_m = std::max(max_m, (int)qk.nrow());
    }

    if (n_mats == 0) {
        return out_list;
    }

    // Size the motif block so the intermediate stays near BLOCK_TARGET_ELEMS
    // regardless of how long the frequency matrices are.
    int max_block = (int)(BLOCK_TARGET_ELEMS / ((size_t)max_m * D));
    max_block = std::max(1, std::min(max_block, MAX_MOTIF_BLOCK));

    // Then shrink it so that one frequency matrix on its own already splits into
    // enough tasks to fill the workers. Without this, a single matrix scored
    // against a database that fits in one block is a single task and runs on one
    // thread - measured at 38ms against 61 motifs, where the same matrix inside a
    // batch costs 1ms.
    //
    // Deliberately derived from the database and the worker count only, never
    // from the batch size: the block width changes the dgemm's inner dimension,
    // and a BLAS that switches kernel on it returns results differing in the last
    // bit or two. Keeping the block independent of how many matrices were passed
    // means a region scores identically whether it is handed over on its own or
    // inside a batch of a thousand.
    if (n_threads > 1) {
        max_block = std::max(MIN_MOTIF_BLOCK, std::min(max_block, n_motifs / n_threads));
    }

    // Group motifs of similar length together, so a block's dgemm runs at the
    // longest motif in that block rather than the longest in the database. On a
    // database whose motifs vary in length this is most of the padded work, and
    // it costs nothing numerically - the fold already stops at each motif's own
    // length, so the columns dropped were never read.
    std::vector<int> order(n_motifs);
    for (int i = 0; i < n_motifs; i++) {
        order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return motif_lengths[a] < motif_lengths[b]; });

    bool identity = true;
    for (int i = 0; i < n_motifs; i++) {
        if (order[i] != i) {
            identity = false;
            break;
        }
    }

    std::vector<FreqTask> tasks;
    for (int k = 0; k < n_mats; k++) {
        for (int i = 0; i < n_motifs; i += max_block) {
            const int nb = std::min(max_block, n_motifs - i);
            const int bl = motif_lengths[order[i + nb - 1]];  // sorted: last is longest
            // The gather is only needed when the sort actually moved something or
            // the block is narrower than the database; "sum" would otherwise pay a
            // copy it never needed. "multiply" always needs its own buffer for exp().
            tasks.push_back({k, i, nb, bl, identity && bl == D && !multiply});
        }
    }

    FreqLocalPWMWorker worker(q, q_len, out, pwm.begin(), pwm_rc.begin(), motif_lengths.begin(),
                              order.data(), n_motifs, D, multiply, bidirect, tasks, max_m,
                              max_block);
    parallelFor(0, tasks.size(), worker);

    return out_list;
}
