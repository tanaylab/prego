# Load a multiple motif regression model from a file

Load a multiple motif regression model from a file

## Usage

``` r
load_multi_regression(
  fn,
  response = NULL,
  sequences = NULL,
  motif_dataset = all_motif_datasets(),
  parallel = getOption("prego.parallel", FALSE),
  alternative = "two.sided"
)
```

## Arguments

- fn:

  file name or a list with the model

- response:

  A matrix (or vector) of response variables. The number of rows must
  equal the number of sequences. A single binary vector (0/1) is
  required when `score_metric = "ks"`.

- sequences:

  A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through
  `toupper`). Sequences must be long enough to cover `spat_num_bins` \*
  `spat_bin_size` bp, and should be centered around the motif/signal.

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

## Value

a list with the following elements:

- models: :

  a list of models.

- model: :

  the combined model.

- spat_min: :

  the minimum spatial position.

- spat_max: :

  the maximum spatial position.

- bidirect: :

  whether the model is bidirectional.

- spat_bin_size: :

  the spatial bin size.

- seq_length: :

  the sequence length.

- motif_num: :

  the number of motifs.

- predict: :

  a function to predict the response.

- predict_multi: :

  a function to predict the response for each motif.

## Examples

``` r
if (FALSE) { # \dontrun{
res_multi <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1],
    final_metric = "ks", spat_bin_size = 40,
    spat_num_bins = 7,
    motif_num = 2
)
tmp <- tempfile()
res_multi$export(tmp)
r <- load_multi_regression(tmp)
} # }
```
