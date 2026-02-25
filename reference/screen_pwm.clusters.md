# Screen for motifs in a database for every cluster

Screen for motifs in a database for every cluster

## Usage

``` r
screen_pwm.clusters(
  sequences,
  clusters,
  dataset = all_motif_datasets(),
  motifs = NULL,
  parallel = getOption("prego.parallel", TRUE),
  min_D = 0.4,
  only_best = FALSE,
  prior = 0.01,
  alternative = "two.sided",
  ...
)
```

## Arguments

- sequences:

  a vector with the sequences

- clusters:

  a vector with the cluster assignments

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

- min_D:

  minimum distance to consider a match

- only_best:

  if TRUE, only return the best match for each cluster

- prior:

  a prior probability for each nucleotide.

- alternative:

  alternative hypothesis for the KS test. Can be "two.sided", "less" or
  "greater"

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

a matrix with the KS D statistics for each cluster (columns) and every
motif (rows) that had at least one cluster with D \>= min_D. If
`only_best` is TRUE, a named vector with the name of best motif match
for each cluster is returned (regardless of `min_D`).

## Examples

``` r
if (FALSE) { # \dontrun{
D_mat <- screen_pwm.clusters(cluster_sequences_example, clusters_example)
dim(D_mat)
D_mat[1:5, 1:5]

# return only the best match
screen_pwm.clusters(cluster_sequences_example, clusters_example, only_best = TRUE)
} # }
```
