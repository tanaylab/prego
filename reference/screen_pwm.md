# Screen for motifs in a database given a response variable

Screen for motifs in a database given a response variable

## Usage

``` r
screen_pwm(
  sequences,
  response,
  metric = NULL,
  dataset = all_motif_datasets(),
  motifs = NULL,
  parallel = getOption("prego.parallel", TRUE),
  only_best = FALSE,
  prior = 0.01,
  alternative = "two.sided",
  ...
)
```

## Arguments

- sequences:

  a vector with the sequences

- response:

  a vector of response variable for each sequence. If the response is a
  matrix, the average will be used.

- metric:

  metric to use in order to choose the best motif. One of 'ks' or 'r2'.
  If NULL - the default would be 'ks' for binary variables, and 'r2' for
  continuous variables.

- dataset:

  a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an
  additional column 'motif' containing the motif name, for example
  `HOMER_motifs` or `JASPAR_motifs`, or
  [`all_motif_datasets()`](https://tanaylab.github.io/prego/reference/all_motif_datasets.md),
  or a MotifDB object.

- motifs:

  names of specific motifs to extract from the dataset

- parallel:

  logical, whether to use parallel processing

- only_best:

  return only the best motif (the one with the highest score). If FALSE,
  all the motifs will be returned.

- prior:

  a prior probability for each nucleotide.

- alternative:

  alternative hypothesis for the KS test. One of 'two.sided', 'less' or
  'greater'.

- ...:

  Arguments passed on to
  [`compute_pwm`](https://tanaylab.github.io/prego/reference/compute_pwm.md)

  `pssm`

  :   a PSSM matrix or data frame. The columns of the matrix or data
      frame should be named with the nucleotides ('A', 'C', 'G' and
      'T').

  `spat`

  :   a data frame with the spatial model (as returned from the `$spat`
      slot from the regression). Should contain a column called 'bin'
      and a column called 'spat_factor'.

  `spat_min`

  :   the minimum position to use from the sequences. The default is 1.

  `spat_max`

  :   the maximum position to use from the sequences. The default is the
      length of the sequences.

  `bidirect`

  :   is the motif bi-directional. If TRUE, the reverse-complement of
      the motif will be used as well.

  `func`

  :   the function to use to combine the PWMs for each sequence. Either
      'logSumExp' or 'max'. The default is 'logSumExp'.

## Value

a data frame with the following columns:

- motif: :

  the motif name.

- score: :

  the score of the motif (depending on `metric`).

if `only_best` is TRUE, only the best motif would be returned (a data
framw with a single row).

## Examples

``` r
if (FALSE) { # \dontrun{
res_screen <- screen_pwm(cluster_sequences_example, cluster_mat_example[, 1])
head(res_screen)

# only best match
screen_pwm(cluster_sequences_example, cluster_mat_example[, 1])

# with r^2 metric
res_screen <- screen_pwm(sequences_example, response_mat_example[, 1], metric = "r2")
head(res_screen)
} # }
```
