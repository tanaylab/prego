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
  `toupper`). Please make sure that the sequences are long enough to
  cover `spat_num_bins` \* `spat_bin_size` bp, and that they are
  centered around the motif/signal.

- response:

  A matrix of response variables - number of rows should equal the
  number of sequences

- spat_bin_size:

  size of the spatial bin (in bp).

- spat_num_bins:

  number of spatial bins. Please make sure that the sequences are long
  enough to cover this number of bins. bp outside of spat_bin_size \*
  spat_num_bins would be ignored. If `bidirect` is TRUE, the number of
  bins should be odd as 'prego' symmetrizes the motif around the center
  bin.

- bidirect:

  is the motif bi-directional. If TRUE, the reverse-complement of the
  motif will be used as well.

- include_response:

  include the response in the resulting list (default: TRUE)

- motif_num:

  Number of motifs to infer. When `motif_num` \> 1, the function would
  run `motif_num` times, each time on the residuals of a linear model of
  all the previous runs (see `smooth_k` parameter). The best motif is
  then returned, while all the others are stored at 'models' in the
  return value.

- multi_kmers:

  if TRUE, different candidates of kmers would be regressed in order to
  find the best seed according to `final_metric`.

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

  whether to run optimization in parallel. use `set_parallel` to set the
  number of cores to use.

- match_with_db:

  match the resulting PWMs with motif databases using `pssm_match`. Note
  that the closest match is returned, even if it is not similar enough
  in absolute terms.

- screen_db:

  Screen `motif_dataset` using `screen_pwm` and use the best motif as
  the initial motif. If TRUE, the following fields would be added to the
  return value: "db_motif", "db_motif_pred", "db_motif_pssm" and
  "db_motif_score".

- motif_dataset:

  a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an
  additional column 'motif' containing the motif name, for example
  `HOMER_motifs`, `JASPAR_motifs` or all_motif_datasets(). By default
  all_motif_datasets() would be used.

- seed:

  random seed

- final_metric:

  metric to use in order to choose the best motif. One of 'ks' or 'r2'.
  Note that unlike `score_metric` which is used in the regression
  itself, this metric is used only for choosing the best motif out of
  all the runs on the sampled dataset. If NULL - 'ks' would be used for
  binary response and 'r2' for continuous response.

- unif_prior:

  uniform prior for nucleotide probabilities

- alternative:

  alternative hypothesis for the p-value calculation when using
  `ks.test`. One of "two.sided", "less" or "greater".

- energy_func:

  a function to transform the energy at each iteration. Should accept a
  numeric vector and return a numeric vector. e.g. `log` or
  `function(x) x^2`. Note that the range of the input energies is
  between 0 and 1 (the probability of the motif in the sequence), so if
  you inferred the function using the the returned energies (which are
  in log scale) you should make sure that the function first log
  transforms using `log_energy=TRUE`.

- ...:

  Arguments passed on to
  [`regress_pwm`](https://tanaylab.github.io/prego/reference/regress_pwm.md),
  [`screen_kmers`](https://tanaylab.github.io/prego/reference/screen_kmers.md)

  `motif`

  :   Initial motif to start the regression from. Can be either a string
      with a kmer where the character "\*" indicates a wildcard or a
      data frame with a pre-computed PSSM (see the slot `pssm` in the
      return value of this function). If NULL - a K-mer screen would be
      performed in order to find the best kmer for initialization. If
      `init_from_dataset` is TRUE, the regression would be initialized
      from the PSSM of the best motif in the dataset.

  `init_from_dataset`

  :   initialize the regression from the PSSM of the best motif in
      `motif_dataset`, using `final_metric` as the metric. If TRUE, the
      `motif` parameter would be ignored. See
      [`screen_pwm`](https://tanaylab.github.io/prego/reference/screen_pwm.md)
      for more details.

  `motif_length`

  :   Length of the seed motif. If the motif is shorter than this, it
      will be extended by wildcards (stars). Note that If the motif is
      longer than this, it will *not* be truncated.

  `score_metric`

  :   metric to use for optimizing the PWM. One of "r2" or "ks". When
      using "ks" the response variable should be a single vector of 0
      and 1.

  `spat_model`

  :   a previously computed spatial model (see `spat`) in the return
      value of this function.

  `improve_epsilon`

  :   minimum improve in the objective function to continue the
      optimization

  `min_nuc_prob`

  :   minimum nucleotide probability in every iteration

  `verbose`

  :   show verbose messages.

  `consensus_single_thresh,consensus_double_thresh`

  :   thresholds for the consensus sequence calculation (single and
      double nucleotides)

  `kmer_length`

  :   a vector of kmer lengths to screen in order to find the best seed
      motif.

  `max_cands`

  :   maximum number of kmer candidates to try.

  `smooth_k`

  :   k for smoothing the predictions of each model in order to compute
      the residuals when `motif_num` \> 1. The residuals are computed as
      `response` - running mean of size 'k' of the current model.

  `min_kmer_cor`

  :   minimal correlation between the kmer and the response in order to
      use it as a seed.

  `internal_num_folds`

  :   number of folds to use in the internal cross-validation.

  `sample_for_kmers`

  :   Use a random sample of the dataset in order to find the best kmer.
      This is useful when the dataset is very large and the kmer screen
      would take a long time. Note that the final regression would be
      performed on the entire dataset. Only relevant when `multi_kmers`
      is TRUE.

  `val_frac`

  :   fraction of the dataset to use for the internal validation. when
      using multiple kmers. Default: 0.1.

  `log_energy`

  :   transform the energy to log scale on each iteration.

  `xmin,xmax,npts`

  :   range for the energy function and the number of points to use for
      its interpolation.

  `energy_func_generator`

  :   a function to generate the energy function when regressing
      multiple motifs. Should accept the result of the previous
      iteration + the original response and return a function similar to
      `energy_func`. e.g.
      ` function(prev_reg, resp) { df <- data.frame(x = prev_reg$pred, y = resp) fn_gam <- as.formula("y ~ s(x, k=3, bs='cr')") model <- mgcv::gam(fn_gam, family = binomial(link = "logit"), data = df, method="REML") function(z){ mgcv::predict.gam(object = model, newdata = data.frame(x = z)) }}`.
      When this parameter is not NULL, energy_func_generator would
      create an energy function and then run another step of regression
      initialized with the previous motif with `energy_func` as the
      energy function. This is useful when the energy function is not
      monotonic, for example - one might want to use a gam model to fit
      the energy function like in the example above.

  `optimize_pwm`

  :   optimize the PWM model (Default: TRUE). If FALSE, the PWM model
      would be used as the initial model for the spatial model.

  `optimize_spat`

  :   optimize the spatial model (Default: TRUE). If FALSE, the spatial
      model would be used as the initial model for the PWM model.

  `kmer_sequence_length`

  :   the length of the sequence to use for the kmer screen. If NULL,
      the entire sequence would be used.

  `symmetrize_spat`

  :   if TRUE, the spatial model would be symmetrized around the center
      bin. Default: TRUE.

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

a list with the following elements:

- pssm: :

  data frame with the pssm matrix with the inferred motif, where rows
  are positions and columns are nucleotides.

- spat: :

  a data frame with the inferred spatial model, with the spatial factor
  for each bin.

- pred: :

  a vector with the predicted pwm for each sequence.

- consensus: :

  Consensus sequence based on the PSSM.

- response: :

  The response matrix. If `include_response` is FALSE, the response
  matrix is not included in the list.

- r2: :

  \\r^2\\ of the prediction with respect to the each response variable.

- ks: :

  If response is binary, Kolmogorov-Smirnov test results of the
  predictions where the response was 1 vs the predictions where the
  response was 0.

- seed_motif: :

  The seed motif that started the regression.

- kmers: :

  The k-mers that were screened in order to find the best seed motif (if
  motif was NULL).

- sample_idxs: :

  The indices of the sequences that were used for the regression (only
  for `regress_pwm.sample`).

- predict: :

  a function that can be used to predict the PWM for a new sequence.

When `match_with_db` is TRUE, the following additional elements are
returned:

- motif_db: :

  The motif database that the most similar to the resulting PSSM.

- db_match_cor: :

  The correlation between the resulting PSSM and the closest match in
  the motif database.

- db_match_pssm: :

  The PSSM of the closest match in the motif database.

- db_match_pred: :

  The predicted PWM of the closest match in the motif database.

- db_match_r2: :

  The \\r^2\\ of the predicted PWM of the closest match in the motif
  database and the response

- db_match_ks: :

  If response is binary, the Kolmogorov-Smirnov test results of the
  predicted PWM of the closest match in the motif database where the
  response was 1 vs the predictions where the response was 0.

When `screen_db` is TRUE, the following additional elements are
returned:

- db_motif: :

  The best motif from the motif database.

- db_motif_pred: :

  The predicted PWM of the best motif from the motif database.

- db_motif_pssm: :

  The PSSM of the best motif from the motif database.

- db_motif_score: :

  The score of the best motif from the motif database.

When `n_motifs` is greater than 1, a list with the following elements is
returned:

- models: :

  A list (as above) of each inferred model

- multi_stats: :

  A data frame with the following columns: `model`, `score` (KS for
  binary, r^2 otherwise), `comb_score` (score for the combined linear
  model for models 1:i) and additional statistics per model

- pred: :

  a vector with the predicted pwm for using a linear model of the
  combined scores.

- comb_modle: :

  a linear model of the combined scores.

- predict: :

  a function that can be used to predict the PWM for a new sequence.

- predict_multi: :

  a function that can be used to predict the PWM for the different
  models for a new sequence

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
