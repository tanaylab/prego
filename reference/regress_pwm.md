# Perform a PWM regression

Fits a PWM (and a spatial model) that predicts `response` from
`sequences`. The function operates in one of three modes depending on
its arguments:

## Usage

``` r
regress_multiple_motifs(
  sequences,
  response,
  motif = NULL,
  motif_length = 15,
  score_metric = "r2",
  bidirect = TRUE,
  spat_bin_size = NULL,
  spat_num_bins = NULL,
  spat_model = NULL,
  improve_epsilon = 0.0001,
  min_nuc_prob = 0.001,
  unif_prior = 0.05,
  include_response = TRUE,
  seed = 60427,
  verbose = FALSE,
  kmer_length = 8,
  multi_kmers = FALSE,
  final_metric = NULL,
  max_cands = 10,
  min_gap = 0,
  max_gap = 1,
  min_kmer_cor = 0.08,
  motif_num = 2,
  smooth_k = 100,
  consensus_single_thresh = 0.6,
  consensus_double_thresh = 0.85,
  internal_num_folds = 1,
  match_with_db = TRUE,
  motif_dataset = all_motif_datasets(),
  parallel = getOption("prego.parallel", FALSE),
  alternative = "less",
  sample_for_kmers = FALSE,
  sample_frac = NULL,
  sample_idxs = NULL,
  sample_ratio = 1,
  log_energy = FALSE,
  energy_func_generator = NULL,
  energy_func = NULL,
  optimize_pwm = TRUE,
  optimize_spat = TRUE,
  ...
)

regress_pwm(
  sequences,
  response,
  motif = NULL,
  motif_length = 15,
  init_from_dataset = FALSE,
  score_metric = "r2",
  bidirect = TRUE,
  spat_bin_size = NULL,
  spat_num_bins = NULL,
  spat_model = NULL,
  improve_epsilon = 0.0001,
  min_nuc_prob = 0.001,
  unif_prior = 0.05,
  include_response = TRUE,
  seed = 60427,
  verbose = FALSE,
  kmer_length = 8,
  multi_kmers = TRUE,
  final_metric = NULL,
  max_cands = 10,
  min_gap = 0,
  max_gap = 1,
  min_kmer_cor = 0.08,
  motif_num = 1,
  smooth_k = 100,
  consensus_single_thresh = 0.6,
  consensus_double_thresh = 0.85,
  internal_num_folds = 1,
  match_with_db = TRUE,
  screen_db = FALSE,
  motif_dataset = all_motif_datasets(),
  parallel = getOption("prego.parallel", FALSE),
  alternative = "less",
  sample_for_kmers = FALSE,
  sample_frac = NULL,
  sample_idxs = NULL,
  sample_ratio = 1,
  val_frac = 0.1,
  log_energy = FALSE,
  energy_func = NULL,
  xmin = -100,
  xmax = 100,
  npts = 10000,
  energy_func_generator = NULL,
  optimize_pwm = TRUE,
  optimize_spat = TRUE,
  kmer_sequence_length = NULL,
  symmetrize_spat = TRUE,
  return_all = FALSE,
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

- motif:

  Initial seed for the regression. Either a kmer string ("\*" denotes a
  wildcard) or a data frame with a pre-computed PSSM (same shape as the
  `pssm` slot of the return value). If `NULL` (default), a kmer screen
  finds a seed automatically. Ignored when `init_from_dataset = TRUE`.

- motif_length:

  Length of the seed motif. Shorter seeds are padded with wildcards;
  longer seeds are *not* truncated.

- score_metric:

  Metric optimized during PWM fitting. One of `"r2"` or `"ks"`. `"ks"`
  requires a single binary response.

- bidirect:

  Whether the motif is bi-directional. If TRUE, the reverse-complement
  is also scored.

- spat_bin_size:

  Size of each spatial bin (bp).

- spat_num_bins:

  Number of spatial bins. Sequence length must be at least
  `spat_bin_size * spat_num_bins`; positions outside this window are
  ignored. If `bidirect = TRUE`, this should be odd (prego symmetrizes
  the spatial model around the center bin).

- spat_model:

  A previously computed spatial model (the `spat` slot of a prior
  result) to reuse.

- improve_epsilon:

  Minimum improvement in the objective for the optimizer to continue.

- min_nuc_prob:

  Floor on per-position nucleotide probabilities at every iteration.

- unif_prior:

  Uniform prior added to nucleotide probabilities.

- include_response:

  If TRUE (default), the response matrix is stored in the result.

- seed:

  Random seed.

- verbose:

  Show verbose messages from the optimizer.

- kmer_length:

  Vector of kmer lengths to try in the seed screen.

- multi_kmers:

  If TRUE (default) and `motif = NULL`, generate multiple kmer
  candidates and pick the best by `final_metric`. If `motif` is
  provided, this is ignored. See also `return_all`.

- final_metric:

  Metric used to pick the best motif across candidates (and to score the
  final fit). One of `"ks"` or `"r2"`. Distinct from `score_metric`,
  which is used *inside* the optimizer. If NULL: `"ks"` for binary
  response, `"r2"` otherwise.

- max_cands:

  Maximum number of candidate kmers to regress.

- min_gap, max_gap:

  the length of a gap to be considered in the pattern. Only one gap, of
  length min_gap:max_gap, is being used, and is located anywhere in the
  motif. Note that this greatly expand the search space (and increase
  multiple testing severely).

- min_kmer_cor:

  Minimum kmer-response correlation for a kmer to qualify as a seed
  candidate.

- motif_num:

  Number of motifs to infer. If \> 1, dispatches to
  `regress_multiple_motifs` (residual rounds; see Details).

- smooth_k:

  Window size (in observations) for the running-mean residual computed
  between rounds when `motif_num > 1`.

- consensus_single_thresh, consensus_double_thresh:

  Probability thresholds for calling a single-nucleotide or
  two-nucleotide IUPAC consensus, respectively.

- internal_num_folds:

  Number of internal cross-validation folds during optimization.

- match_with_db:

  If TRUE, match the resulting PSSM to the closest motif in
  `motif_dataset` using `pssm_match`. The closest match is always
  returned, even if similarity is low.

- motif_dataset:

  A data frame with PSSMs (columns `A`, `C`, `G`, `T`, `motif`), e.g.
  `HOMER_motifs`, `JASPAR_motifs`, or
  [`all_motif_datasets()`](https://tanaylab.github.io/prego/reference/all_motif_datasets.md)
  (default).

- parallel:

  Whether to parallelize. Use `set_parallel` to set the number of cores.

- alternative:

  Alternative hypothesis for `ks.test`. One of `"two.sided"`, `"less"`,
  `"greater"`.

- sample_for_kmers:

  If TRUE, run the kmer screen on a random sample of sequences (final
  regression is still on the full data). Useful for large datasets. Only
  relevant when `multi_kmers = TRUE`.

- sample_frac:

  Fraction of sequences sampled for the kmer screen. Default: 0.1.

- sample_idxs:

  Explicit sequence indices to use for the kmer screen. If NULL, a
  random sample is drawn.

- sample_ratio:

  For binary response, ratio of '1' to '0' in the sampled subset.
  Relevant only when `sample_frac` is NULL.

- log_energy:

  If TRUE, the energy is log-transformed at every iteration.

- energy_func_generator:

  A function that, given the previous-iteration result and the original
  response, returns an `energy_func`. Used in `motif_num > 1`: each
  round, the generator builds an energy function (e.g. a GAM-based one
  for non-monotonic relationships) and a second pass is run with it.
  Example:
  ` function(prev_reg, resp) { df <- data.frame(x = prev_reg$pred, y = resp) model <- mgcv::gam(y ~ s(x, k = 3, bs = "cr"), family = binomial("logit"), data = df, method = "REML") function(z) mgcv::predict.gam(model, newdata = data.frame(x = z)) }`

- energy_func:

  Optional function applied to the energy at every iteration; takes and
  returns a numeric vector (e.g. `log` or `function(x) x^2`). Input
  energies are in the range `[0, 1]`; if your function was fit on log
  energies, use `log_energy = TRUE`.

- optimize_pwm:

  If FALSE, the PWM is held fixed and only the spatial model is
  optimized.

- optimize_spat:

  If FALSE, the spatial model is held fixed and only the PWM is
  optimized.

- ...:

  Arguments passed on to
  [`screen_kmers`](https://tanaylab.github.io/prego/reference/screen_kmers.md)

  `min_cor`

  :   Only patterns for which the maximum correlation to one of the
      response variable is larger than min_cor will be reported

  `is_train`

  :   a boolean vector that determine which subset of sequences to use
      when screening

  `min_gap,max_gap`

  :   the length of a gap to be considered in the pattern. Only one gap,
      of length min_gap:max_gap, is being used, and is located anywhere
      in the motif. Note that this greatly expand the search space (and
      increase multiple testing severely).

  `from_range`

  :   Sequences will be considered only from position from_range
      (default 0)

  `to_range`

  :   Sequences will be considered only up to position to_range (default
      NULL - using the length of the sequences)

  `return_mat`

  :   Return a matrix of patterns and their correlation to the response
      variables instead of a data frame. (default: FALSE)

- init_from_dataset:

  If TRUE, initialize from the best-matching PSSM in `motif_dataset`
  (chosen by `final_metric`; see
  [`screen_pwm`](https://tanaylab.github.io/prego/reference/screen_pwm.md)).
  Overrides `motif`.

- screen_db:

  If TRUE, screen `motif_dataset` via
  [`screen_pwm`](https://tanaylab.github.io/prego/reference/screen_pwm.md)
  and add the best database motif to the result (`db_motif`,
  `db_motif_pred`, `db_motif_pssm`, `db_motif_score`).

- val_frac:

  Fraction of sequences held out for internal validation in the
  multi-kmer mode. Default: 0.1.

- xmin, xmax, npts:

  Range and resolution for interpolating `energy_func`.

- kmer_sequence_length:

  Length of the central window used for the kmer screen. If NULL, the
  full sequence is used.

- symmetrize_spat:

  If TRUE (default), the spatial model is symmetrized around the center
  bin.

- return_all:

  If TRUE, return *all* candidate-kmer regressions instead of only the
  best one. Requires `multi_kmers = TRUE`, `motif = NULL`, and
  `motif_num = 1`. Use this when you want N independent kmer-seeded
  motifs without the residual-rounds approach of `motif_num > 1`.
  Default: FALSE. See *Value* for the returned structure.

## Value

The returned structure depends on the mode:

**Single-motif mode** (default, or with `motif` provided) - a list with:

- pssm:

  Data frame of the inferred PSSM: rows are positions, columns are
  nucleotides.

- spat:

  Data frame of the inferred spatial model (one row per bin).

- pred:

  Numeric vector of predictions, one per input sequence.

- consensus:

  IUPAC consensus derived from `pssm`.

- response:

  The response matrix (omitted when `include_response = FALSE`).

- r2:

  Per-response-column \\r^2\\ of `pred` vs the response.

- ks:

  Kolmogorov-Smirnov test of `pred` for response=1 vs response=0 (binary
  response only).

- score:

  The KS statistic (binary) or \\r^2\\ (continuous).

- seed_motif:

  The kmer/string used as the seed.

- kmers:

  The screened candidate kmers (only when the seed was found
  automatically).

- predict:

  A closure: `predict(new_sequences)` returns predictions.

When `match_with_db = TRUE`, also: `motif_db`, `db_match_cor`,
`db_match_pssm`, `db_match_pred`, `db_match_r2`, and (binary)
`db_match_ks`.

When `screen_db = TRUE`, also: `db_motif`, `db_motif_pred`,
`db_motif_pssm`, `db_motif_score`.

**Multi-kmer mode with `return_all = TRUE`** - a *named list* of
single-motif results (named by `seed_motif`), one per candidate kmer,
sorted by validation score (descending). Each element has the structure
above, plus a `val_score` field. Refitting behavior:

- `sample_for_kmers = FALSE` (default): each element's PSSM was fit on
  the train portion (`1 - val_frac` of the full data); `pred` therefore
  covers the train subset only. Call `res[[i]]$predict(sequences)` to
  score the full data, or refit explicitly via
  `regress_pwm(sequences, response, motif = res[[i]]$seed_motif, multi_kmers = FALSE)`.

- `sample_for_kmers = TRUE`: each element is refit on the full data
  (analogous to the single-best path).

**Multi-motif mode (`motif_num > 1`)** - a list with:

- models:

  List of single-motif results, one per round.

- multi_stats:

  Data frame with columns `model`, `score` (per-motif KS or \\r^2\\),
  `comb_score` (linear-combination score using models 1:i), `diff`,
  `consensus`, `seed_motif`, plus db-match columns when
  `match_with_db = TRUE`.

- pred:

  Predictions from the combined linear model.

- model:

  The fitted combined `lm` object.

- predict:

  Closure for predictions from the combined model.

- predict_multi:

  Closure returning a data frame of per-motif predictions (columns `e1`,
  `e2`, ...).

## Details

- Single motif from a fixed seed (`motif` given)::

  One PWM is optimized starting from the supplied kmer or PSSM. Returns
  a single regression result.

- Single motif from a kmer screen (`motif = NULL`, `multi_kmers = TRUE`,
  default)::

  Candidate kmer seeds are generated by `screen_kmers`, each is
  regressed on a train/validation split, the best by validation score is
  refit on the full data. Returns a single regression result. Set
  `return_all = TRUE` to get *all* candidate-kmer regressions instead of
  just the best - useful when you want the top N motifs without the
  residual-rounds approach used by `motif_num > 1`.

- Multiple motifs via residual rounds (`motif_num > 1`)::

  Delegates to `regress_multiple_motifs`. Each successive motif is fit
  against the residual of a linear model of the previous motifs. Returns
  a multi-motif result (see *Value*).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(sequences_example, response_mat_example)
res$pssm
res$spat
head(res$pred)

plot_regression_qc(res)

# intialize with a pre-computed PSSM
res1 <- regress_pwm(sequences_example, response_mat_example, motif = res$pssm)

# intialize with a pre-computed PSSM and spatial model
res2 <- regress_pwm(
    sequences_example,
    response_mat_example,
    motif = res$pssm,
    spat_model = res$spat
)

# binary response
res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
plot_regression_qc(res_binary)

# match with db
res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], match_with_db = TRUE)
plot_regression_qc(res_binary)

# Screen for best db motif
res_binary <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1], screen_db = TRUE)
plot_regression_qc(res_binary)

# initialize with a motif from the database
res_binary <- regress_pwm(
    cluster_sequences_example,
    cluster_mat_example[, 1],
    init_from_dataset = TRUE
)

# use multiple kmer seeds (returns the single best motif)
res_multi <- regress_pwm(
    cluster_sequences_example,
    cluster_mat_example[, 1],
    multi_kmers = TRUE,
    kmer_length = 6:8,
    final_metric = "ks"
)
plot_regression_qc(res_multi)

# return all candidate-kmer motifs (sorted by validation score)
res_all <- regress_pwm(
    cluster_sequences_example,
    cluster_mat_example[, 1],
    multi_kmers = TRUE,
    return_all = TRUE,
    max_cands = 5
)
length(res_all)
sapply(res_all, function(x) x$seed_motif)
sapply(res_all, function(x) x$val_score)

# Screen for multiple motifs via residual rounds
res_multi <- regress_pwm(
    cluster_sequences_example,
    cluster_mat_example[, 1],
    motif_num = 3,
    match_with_db = TRUE
)
res_multi$multi_stats
plot_regression_qc_multi(res_multi)
} # }

# regress_multiple_motifs() can also be called directly (equivalent to
# regress_pwm() with motif_num > 1, which internally dispatches to it)
if (FALSE) { # \dontrun{
res_multi2 <- regress_multiple_motifs(
    cluster_sequences_example,
    cluster_mat_example[, 1],
    motif_num = 5,
    match_with_db = TRUE
)
} # }
```
