# Run PWM regression on a sample of the data

The optimization would be performed with a sampled dataset of size
`sample_frac`, or explicit sampled indices `sample_idxs`.

## Usage

``` r
regress_pwm.sample(
  sequences,
  response,
  spat_bin_size = NULL,
  spat_num_bins = NULL,
  bidirect = TRUE,
  include_response = TRUE,
  motif_num = 1,
  multi_kmers = TRUE,
  sample_frac = NULL,
  sample_idxs = NULL,
  sample_ratio = 1,
  parallel = getOption("prego.parallel", TRUE),
  match_with_db = TRUE,
  screen_db = FALSE,
  motif_dataset = all_motif_datasets(),
  seed = 60427,
  final_metric = NULL,
  unif_prior = 0.05,
  alternative = "two.sided",
  energy_func = NULL,
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

- spat_bin_size:

  Size of each spatial bin (bp).

- spat_num_bins:

  Number of spatial bins. Sequence length must be at least
  `spat_bin_size * spat_num_bins`; positions outside this window are
  ignored. If `bidirect = TRUE`, this should be odd (prego symmetrizes
  the spatial model around the center bin).

- bidirect:

  Whether the motif is bi-directional. If TRUE, the reverse-complement
  is also scored.

- include_response:

  If TRUE (default), the response matrix is stored in the result.

- motif_num:

  Number of motifs to infer. If \> 1, dispatches to
  [`regress_multiple_motifs`](https://tanaylab.github.io/prego/reference/regress_pwm.md)
  (residual rounds; see Details).

- multi_kmers:

  If TRUE (default) and `motif = NULL`, generate multiple kmer
  candidates and pick the best by `final_metric`. If `motif` is
  provided, this is ignored. See also `return_all`.

- sample_frac:

  fraction of the dataset to sample. When `response` is categorical (0
  and 1), the sampling would be stratified by the category, i.e.
  `sample_frac` can be a vector of length 2 with the fraction of 0 and 1
  responses to sample respectively. If NULL - the default would be 0.1
  for continuous variables, and for binary variables - the number of 0
  responses would be equal to `sample_ratio` times the number of 1
  responses.

- sample_idxs:

  indices of the sequences to use. If NULL, the indices would be sampled
  using `sample_frac`.

- sample_ratio:

  ratio between the '1' category and the '0' category in the sampled
  dataset. Relevant only when `sample_frac` is NULL.

- parallel:

  Whether to parallelize. Use `set_parallel` to set the number of cores.

- match_with_db:

  If TRUE, match the resulting PSSM to the closest motif in
  `motif_dataset` using `pssm_match`. The closest match is always
  returned, even if similarity is low.

- screen_db:

  If TRUE, screen `motif_dataset` via
  [`screen_pwm`](https://tanaylab.github.io/prego/reference/screen_pwm.md)
  and add the best database motif to the result (`db_motif`,
  `db_motif_pred`, `db_motif_pssm`, `db_motif_score`).

- motif_dataset:

  A data frame with PSSMs (columns `A`, `C`, `G`, `T`, `motif`), e.g.
  `HOMER_motifs`, `JASPAR_motifs`, or
  [`all_motif_datasets()`](https://tanaylab.github.io/prego/reference/all_motif_datasets.md)
  (default).

- seed:

  Random seed.

- final_metric:

  Metric used to pick the best motif across candidates (and to score the
  final fit). One of `"ks"` or `"r2"`. Distinct from `score_metric`,
  which is used *inside* the optimizer. If NULL: `"ks"` for binary
  response, `"r2"` otherwise.

- unif_prior:

  Uniform prior added to nucleotide probabilities.

- alternative:

  Alternative hypothesis for `ks.test`. One of `"two.sided"`, `"less"`,
  `"greater"`.

- energy_func:

  Optional function applied to the energy at every iteration; takes and
  returns a numeric vector (e.g. `log` or `function(x) x^2`). Input
  energies are in the range `[0, 1]`; if your function was fit on log
  energies, use `log_energy = TRUE`.

- ...:

  Arguments passed on to
  [`regress_pwm`](https://tanaylab.github.io/prego/reference/regress_pwm.md),
  [`screen_kmers`](https://tanaylab.github.io/prego/reference/screen_kmers.md)

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

  `spat_model`

  :   A previously computed spatial model (the `spat` slot of a prior
      result) to reuse.

  `improve_epsilon`

  :   Minimum improvement in the objective for the optimizer to
      continue.

  `min_nuc_prob`

  :   Floor on per-position nucleotide probabilities at every iteration.

  `verbose`

  :   Show verbose messages from the optimizer.

  `consensus_single_thresh,consensus_double_thresh`

  :   Probability thresholds for calling a single-nucleotide or
      two-nucleotide IUPAC consensus, respectively.

  `return_all`

  :   If TRUE, return *all* candidate-kmer regressions instead of only
      the best one. Requires `multi_kmers = TRUE`, `motif = NULL`, and
      `motif_num = 1`. Use this when you want N independent kmer-seeded
      motifs without the residual-rounds approach of `motif_num > 1`.
      Default: FALSE. See *Value* for the returned structure.

  `kmer_length`

  :   Vector of kmer lengths to try in the seed screen.

  `max_cands`

  :   Maximum number of candidate kmers to regress.

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

  `val_frac`

  :   Fraction of sequences held out for internal validation in the
      multi-kmer mode. Default: 0.1.

  `log_energy`

  :   If TRUE, the energy is log-transformed at every iteration.

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

  `min_cor`

  :   Only patterns for which the maximum correlation to one of the
      response variable is larger than min_cor will be reported

  `is_train`

  :   a boolean vector that determine which subset of sequences to use
      when screening

  `from_range`

  :   Sequences will be considered only from position from_range
      (default 0)

  `to_range`

  :   Sequences will be considered only up to position to_range (default
      NULL - using the length of the sequences)

  `return_mat`

  :   Return a matrix of patterns and their correlation to the response
      variables instead of a data frame. (default: FALSE)

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

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm.sample(
    cluster_sequences_example,
    cluster_mat_example[, 1],
    final_metric = "ks",
    screen_db = TRUE
)

res$pssm
res$spat
head(res$pred)

plot_regression_qc(res)
} # }
```
