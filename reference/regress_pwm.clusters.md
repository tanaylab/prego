# Run PWM regression on clusters

Run PWM regression on clusters

## Usage

``` r
regress_pwm.clusters(
  sequences,
  clusters,
  use_sample = TRUE,
  match_with_db = TRUE,
  screen_db = FALSE,
  sample_frac = NULL,
  sample_ratio = 1,
  final_metric = "ks",
  parallel = getOption("prego.parallel", TRUE),
  use_sge = FALSE,
  dataset = all_motif_datasets(),
  motifs = NULL,
  min_D = 0,
  prior = 0.01,
  alternative = "two.sided",
  ...
)
```

## Arguments

- sequences:

  A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through
  `toupper`). Sequences must be long enough to cover `spat_num_bins` \*
  `spat_bin_size` bp, and should be centered around the motif/signal.

- clusters:

  a vector with the cluster assignments for each sequence

- use_sample:

  whether to use sampled optimization or not (default: FALSE).

- match_with_db:

  match the resulting PWMs with motif databases using `pssm_match`. This
  would add a column named 'db_match' to the stats data frame, together
  with 'pred_mat_db' with the database motif predictions, and
  'db_dataset' which is similar to 'motif_dataset' for the database
  motifs. Note that the closest match is returned, even if it is not
  similar enough in absolute terms. Also, the match is done between the
  resulting regression *pssm* and the pssms in the database - in order
  to find the best motif in the database set `screen_db=TRUE`.

- screen_db:

  screen for the best motif in the database which explains the clusters.
  See `screen_pwm.clusters`.

- sample_frac:

  a vector of two numbers, specifying the fraction of sequences to use
  in when sampling for the sequences which are not in the cluster (first
  number) and in the cluster (second number). If NULL -

- sample_ratio:

  When `sample_frac` is NULL, the number of sequences not in the cluster
  would be equal to `sample_ratio` times the number of sequences in the
  cluster.

- final_metric:

  Metric used to pick the best motif across candidates (and to score the
  final fit). One of `"ks"` or `"r2"`. Distinct from `score_metric`,
  which is used *inside* the optimizer. If NULL: `"ks"` for binary
  response, `"r2"` otherwise.

- parallel:

  Whether to parallelize. Use `set_parallel` to set the number of cores.

- use_sge:

  use the function `gcluster.run2` from the misha.ext package to run the
  optimization on a SGE cluster. Only relevant if the `misha.ext`
  package is installed. Note that `gcluster.run2` writes the current
  environment before starting the parallelization, so it is better to
  run this function in a clean environment. Also, Note that 'prego'
  needs to be installed in order for this to work, i.e. you cannot use
  `devtools::load_all()` or
  [`pkgload::load_all()`](https://pkgload.r-lib.org/reference/load_all.html)
  to load the package.

- dataset:

  a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an
  additional column 'motif' containing the motif name, for example
  `HOMER_motifs` or `JASPAR_motifs`, or
  [`all_motif_datasets()`](https://tanaylab.github.io/prego/reference/all_motif_datasets.md),
  or a MotifDB object.

- motifs:

  names of specific motifs to extract from the dataset

- min_D:

  minimum distance to consider a match

- prior:

  a prior probability for each nucleotide.

- alternative:

  alternative hypothesis for the KS test. Can be "two.sided", "less" or
  "greater"

- ...:

  Arguments passed on to
  [`regress_pwm`](https://tanaylab.github.io/prego/reference/regress_pwm.md),
  [`regress_pwm.sample`](https://tanaylab.github.io/prego/reference/regress_pwm.sample.md)

  `response`

  :   A matrix (or vector) of response variables. The number of rows
      must equal the number of sequences. A single binary vector (0/1)
      is required when `score_metric = "ks"`.

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

  `seed`

  :   Random seed.

  `consensus_single_thresh,consensus_double_thresh`

  :   Probability thresholds for calling a single-nucleotide or
      two-nucleotide IUPAC consensus, respectively.

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

  `sample_idxs`

  :   Explicit sequence indices to use for the kmer screen. If NULL, a
      random sample is drawn.

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

- models: :

  a list with the models for each cluster

- cluster_mat: :

  an indicator matrix with the cluster assignments

- pred_mat: :

  a matrix of the energies of the predicted motifs per cluster (columns)
  in each sequence (rows)

- motif_dataset: :

  a data frame with the PSSMs for each cluster

- spat_dataset: :

  a data frame with the spatial model for each cluster

- stats: :

  a data frame with statistics for each cluster

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm.clusters(cluster_sequences_example, clusters_example)
head(res$pred_mat)
res$stats

plot_regression_qc(res$models[[1]], title = names(res$models)[1])

# multiple motifs per cluster
res_multi <- regress_pwm.clusters(cluster_sequences_example, clusters_example, motif_num = 3)
res_multi$multi_stats
plot_regression_qc_multi(res_multi$models[[1]], title = names(res_multi$models)[1])
} # }

# screen also for the best motif in the database
if (FALSE) { # \dontrun{
res_screen <- regress_pwm.clusters(cluster_sequences_example, clusters_example, screen_db = TRUE)
res_screen$stats

plot_regression_qc(res_screen$models[[1]], title = names(res_screen$models)[1])
} # }
```
