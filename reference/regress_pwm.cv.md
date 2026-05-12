# Cross-validate a PWM regression model

Perform cross-validation on a PWM regression model. You can either
provide explicit folds, or use the `nfolds` argument to set the number
of folds. If the response is binary (0 and 1) or a `categories` vector
is given, the folds would be stratified by the response/categories.

## Usage

``` r
regress_pwm.cv(
  sequences,
  response,
  nfolds = NULL,
  metric = NULL,
  folds = NULL,
  categories = NULL,
  use_sample = FALSE,
  seed = 60427,
  parallel = getOption("prego.parallel", FALSE),
  fold_parallel = !parallel && getOption("prego.parallel", FALSE),
  add_full_model = TRUE,
  alternative = "less",
  ...
)
```

## Arguments

- sequences:

  A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through
  `toupper`). Sequences must be long enough to cover `spat_num_bins` \*
  `spat_bin_size` bp, and should be centered around the motif/signal.

- response:

  A matrix (or vector) of response variables. The number of rows must
  equal the number of sequences. A single binary vector (0/1) is
  required when `score_metric = "ks"`.

- nfolds:

  number of folds for cross-validation. Can be NULL if `folds` are
  provided.

- metric:

  metric to use for cross-validation. One of 'ks' or 'r2'. If NULL -
  'ks' would be set for binary response and 'r2' for continuous
  response.

- folds:

  vector of fold numbers for each sequence (optional)

- categories:

  vector of categories for each sequence (optional)

- use_sample:

  whether to use sampled optimization or not.

- seed:

  Random seed.

- parallel:

  whether to run the cross-validation in parallel.

- fold_parallel:

  whether to run the optimization in each fold in parallel. It is
  recommended to set this to FALSE if `parallel` is TRUE.

- add_full_model:

  whether to add the full model (without cross-validation) to the
  results.

- alternative:

  Alternative hypothesis for `ks.test`. One of `"two.sided"`, `"less"`,
  `"greater"`.

- ...:

  Arguments passed on to
  [`regress_pwm`](https://tanaylab.github.io/prego/reference/regress_pwm.md),
  [`regress_pwm.sample`](https://tanaylab.github.io/prego/reference/regress_pwm.sample.md)

  `motif`

  :   Initial seed for the regression. Either a kmer string ("\*"
      denotes a wildcard) or a data frame with a pre-computed PSSM (same
      shape as the `pssm` slot of the return value). If `NULL`
      (default), a kmer screen finds a seed automatically. Ignored when
      `init_from_dataset = TRUE`.

  `init_from_dataset`

  :   If TRUE, initialize from the best-matching PSSM in `motif_dataset`
      (chosen by `final_metric`; see
      [`screen_pwm`](https://tanaylab.github.io/prego/reference/screen_pwm.md)).
      Overrides `motif`.

  `motif_length`

  :   Length of the seed motif. Shorter seeds are padded with wildcards;
      longer seeds are *not* truncated.

  `score_metric`

  :   Metric optimized during PWM fitting. One of `"r2"` or `"ks"`.
      `"ks"` requires a single binary response.

  `bidirect`

  :   Whether the motif is bi-directional. If TRUE, the
      reverse-complement is also scored.

  `spat_bin_size`

  :   Size of each spatial bin (bp).

  `spat_num_bins`

  :   Number of spatial bins. Sequence length must be at least
      `spat_bin_size * spat_num_bins`; positions outside this window are
      ignored. If `bidirect = TRUE`, this should be odd (prego
      symmetrizes the spatial model around the center bin).

  `spat_model`

  :   A previously computed spatial model (the `spat` slot of a prior
      result) to reuse.

  `improve_epsilon`

  :   Minimum improvement in the objective for the optimizer to
      continue.

  `min_nuc_prob`

  :   Floor on per-position nucleotide probabilities at every iteration.

  `unif_prior`

  :   Uniform prior added to nucleotide probabilities.

  `include_response`

  :   If TRUE (default), the response matrix is stored in the result.

  `verbose`

  :   Show verbose messages from the optimizer.

  `consensus_single_thresh,consensus_double_thresh`

  :   Probability thresholds for calling a single-nucleotide or
      two-nucleotide IUPAC consensus, respectively.

  `match_with_db`

  :   If TRUE, match the resulting PSSM to the closest motif in
      `motif_dataset` using `pssm_match`. The closest match is always
      returned, even if similarity is low.

  `screen_db`

  :   If TRUE, screen `motif_dataset` via
      [`screen_pwm`](https://tanaylab.github.io/prego/reference/screen_pwm.md)
      and add the best database motif to the result (`db_motif`,
      `db_motif_pred`, `db_motif_pssm`, `db_motif_score`).

  `motif_dataset`

  :   A data frame with PSSMs (columns `A`, `C`, `G`, `T`, `motif`),
      e.g. `HOMER_motifs`, `JASPAR_motifs`, or
      [`all_motif_datasets()`](https://tanaylab.github.io/prego/reference/all_motif_datasets.md)
      (default).

  `multi_kmers`

  :   If TRUE (default) and `motif = NULL`, generate multiple kmer
      candidates and pick the best by `final_metric`. If `motif` is
      provided, this is ignored. See also `return_all`.

  `return_all`

  :   If TRUE, return *all* candidate-kmer regressions instead of only
      the best one. Requires `multi_kmers = TRUE`, `motif = NULL`, and
      `motif_num = 1`. Use this when you want N independent kmer-seeded
      motifs without the residual-rounds approach of `motif_num > 1`.
      Default: FALSE. See *Value* for the returned structure.

  `final_metric`

  :   Metric used to pick the best motif across candidates (and to score
      the final fit). One of `"ks"` or `"r2"`. Distinct from
      `score_metric`, which is used *inside* the optimizer. If NULL:
      `"ks"` for binary response, `"r2"` otherwise.

  `kmer_length`

  :   Vector of kmer lengths to try in the seed screen.

  `max_cands`

  :   Maximum number of candidate kmers to regress.

  `motif_num`

  :   Number of motifs to infer. If \> 1, dispatches to
      [`regress_multiple_motifs`](https://tanaylab.github.io/prego/reference/regress_pwm.md)
      (residual rounds; see Details).

  `smooth_k`

  :   Window size (in observations) for the running-mean residual
      computed between rounds when `motif_num > 1`.

  `min_kmer_cor`

  :   Minimum kmer-response correlation for a kmer to qualify as a seed
      candidate.

  `internal_num_folds`

  :   Number of internal cross-validation folds during optimization.

  `sample_for_kmers`

  :   If TRUE, run the kmer screen on a random sample of sequences
      (final regression is still on the full data). Useful for large
      datasets. Only relevant when `multi_kmers = TRUE`.

  `sample_frac`

  :   Fraction of sequences sampled for the kmer screen. Default: 0.1.

  `sample_idxs`

  :   Explicit sequence indices to use for the kmer screen. If NULL, a
      random sample is drawn.

  `sample_ratio`

  :   For binary response, ratio of '1' to '0' in the sampled subset.
      Relevant only when `sample_frac` is NULL.

  `val_frac`

  :   Fraction of sequences held out for internal validation in the
      multi-kmer mode. Default: 0.1.

  `log_energy`

  :   If TRUE, the energy is log-transformed at every iteration.

  `energy_func`

  :   Optional function applied to the energy at every iteration; takes
      and returns a numeric vector (e.g. `log` or `function(x) x^2`).
      Input energies are in the range `[0, 1]`; if your function was fit
      on log energies, use `log_energy = TRUE`.

  `xmin,xmax,npts`

  :   Range and resolution for interpolating `energy_func`.

  `energy_func_generator`

  :   A function that, given the previous-iteration result and the
      original response, returns an `energy_func`. Used in
      `motif_num > 1`: each round, the generator builds an energy
      function (e.g. a GAM-based one for non-monotonic relationships)
      and a second pass is run with it. Example:
      ` function(prev_reg, resp) { df <- data.frame(x = prev_reg$pred, y = resp) model <- mgcv::gam(y ~ s(x, k = 3, bs = "cr"), family = binomial("logit"), data = df, method = "REML") function(z) mgcv::predict.gam(model, newdata = data.frame(x = z)) }`

  `optimize_pwm`

  :   If FALSE, the PWM is held fixed and only the spatial model is
      optimized.

  `optimize_spat`

  :   If FALSE, the spatial model is held fixed and only the PWM is
      optimized.

  `kmer_sequence_length`

  :   Length of the central window used for the kmer screen. If NULL,
      the full sequence is used.

  `symmetrize_spat`

  :   If TRUE (default), the spatial model is symmetrized around the
      center bin.

  `min_gap,max_gap`

  :   the length of a gap to be considered in the pattern. Only one gap,
      of length min_gap:max_gap, is being used, and is located anywhere
      in the motif. Note that this greatly expand the search space (and
      increase multiple testing severely).

## Value

a list with the following elements:

- cv_models: :

  a list of models, one for each fold.

- cv_pred: :

  a vector of predictions for each sequence.

- score: :

  score of the model on the cross-validated predictions.

- cv_scores: :

  a vector of scores for each fold.

- folds: :

  a vector with the fold number for each sequence.

- full_model: :

  The full model (without cross-validation), if `add_full_model` is
  TRUE.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm.cv(
    cluster_sequences_example, cluster_mat_example[, 1],
    nfolds = 5, use_sample = TRUE, sample_frac = c(0.1, 1)
)
res$score
res$cv_scores

plot(
    res$cv_pred,
    res$full_model$pred,
    xlab = "CV predictions", ylab = "Full model predictions", cex = 0.1
)
plot_regression_prediction_binary(res$cv_pred, cluster_mat_example[, 1])
plot_regression_prediction_binary(res$full_model$pred, cluster_mat_example[, 1])

# without sampling
res <- regress_pwm.cv(
    cluster_sequences_example, cluster_mat_example[, 1],
    nfolds = 5, use_sample = FALSE
)
res$score
res$cv_scores
plot(res$cv_pred,
    res$full_model$pred,
    xlab = "CV predictions", ylab = "Full model predictions", cex = 0.1
)
} # }
```
