# Compute PWMs for a set of sequences given a PSSM matrix

Compute PWMs for a set of sequences given a PSSM matrix

## Usage

``` r
compute_pwm(
  sequences,
  pssm,
  spat = NULL,
  spat_min = 1,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  func = "logSumExp"
)
```

## Arguments

- sequences:

  a vector of sequences

- pssm:

  a PSSM matrix or data frame. The columns of the matrix or data frame
  should be named with the nucleotides ('A', 'C', 'G' and 'T').

- spat:

  a data frame with the spatial model (as returned from the `$spat` slot
  from the regression). Should contain a column called 'bin' and a
  column called 'spat_factor'.

- spat_min:

  the minimum position to use from the sequences. The default is 1.

- spat_max:

  the maximum position to use from the sequences. The default is the
  length of the sequences.

- bidirect:

  is the motif bi-directional. If TRUE, the reverse-complement of the
  motif will be used as well.

- prior:

  a prior probability for each nucleotide.

- func:

  the function to use to combine the PWMs for each sequence. Either
  'logSumExp' or 'max'. The default is 'logSumExp'.

## Value

a vector with the predicted pwm for each sequence.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])

pwm <- compute_pwm(cluster_sequences_example, res$pssm, res$spat)
head(pwm)

# this is similar to the prediction in the regression
head(res$pred)
} # }
```
