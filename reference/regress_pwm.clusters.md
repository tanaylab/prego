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
  `toupper`). Please make sure that the sequences are long enough to
  cover `spat_num_bins` \* `spat_bin_size` bp, and that they are
  centered around the motif/signal.

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

  metric to use in order to choose the best motif. One of 'ks' or 'r2'.
  Note that unlike `score_metric` which is used in the regression
  itself, this metric is used only for choosing the best motif out of
  all the runs on the sampled dataset. If NULL - 'ks' would be used for
  binary response and 'r2' for continuous response.

- parallel:

  whether to run optimization in parallel. use `set_parallel` to set the
  number of cores to use.

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

  :   A matrix of response variables - number of rows should equal the
      number of sequences

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

  `bidirect`

  :   is the motif bi-directional. If TRUE, the reverse-complement of
      the motif will be used as well.

  `spat_bin_size`

  :   size of the spatial bin (in bp).

  `spat_num_bins`

  :   number of spatial bins. Please make sure that the sequences are
      long enough to cover this number of bins. bp outside of
      spat_bin_size \* spat_num_bins would be ignored. If `bidirect` is
      TRUE, the number of bins should be odd as 'prego' symmetrizes the
      motif around the center bin.

  `spat_model`

  :   a previously computed spatial model (see `spat`) in the return
      value of this function.

  `improve_epsilon`

  :   minimum improve in the objective function to continue the
      optimization

  `min_nuc_prob`

  :   minimum nucleotide probability in every iteration

  `unif_prior`

  :   uniform prior for nucleotide probabilities

  `include_response`

  :   include the response in the resulting list (default: TRUE)

  `verbose`

  :   show verbose messages.

  `seed`

  :   random seed

  `consensus_single_thresh,consensus_double_thresh`

  :   thresholds for the consensus sequence calculation (single and
      double nucleotides)

  `motif_dataset`

  :   a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an
      additional column 'motif' containing the motif name, for example
      `HOMER_motifs`, `JASPAR_motifs` or all_motif_datasets(). By
      default all_motif_datasets() would be used.

  `multi_kmers`

  :   if TRUE, different candidates of kmers would be regressed in order
      to find the best seed according to `final_metric`.

  `kmer_length`

  :   a vector of kmer lengths to screen in order to find the best seed
      motif.

  `max_cands`

  :   maximum number of kmer candidates to try.

  `motif_num`

  :   Number of motifs to infer. When `motif_num` \> 1, the function
      would run `motif_num` times, each time on the residuals of a
      linear model of all the previous runs (see `smooth_k` parameter).
      The best motif is then returned, while all the others are stored
      at 'models' in the return value.

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

  `sample_idxs`

  :   indices of the sequences to use for the kmer screen. If NULL, a
      random sample would be used.

  `val_frac`

  :   fraction of the dataset to use for the internal validation. when
      using multiple kmers. Default: 0.1.

  `log_energy`

  :   transform the energy to log scale on each iteration.

  `energy_func`

  :   a function to transform the energy at each iteration. Should
      accept a numeric vector and return a numeric vector. e.g. `log` or
      `function(x) x^2`. Note that the range of the input energies is
      between 0 and 1 (the probability of the motif in the sequence), so
      if you inferred the function using the the returned energies
      (which are in log scale) you should make sure that the function
      first log transforms using `log_energy=TRUE`.

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
