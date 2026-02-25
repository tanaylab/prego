# Compute local PWMs for a set of sequences given a PSSM matrix

compute the local PWM for each position in every sequence. The edges of
each sequences would become NA.

## Usage

``` r
compute_local_pwm(
  sequences,
  pssm,
  spat = NULL,
  spat_min = 0,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  return_list = FALSE
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

- return_list:

  Logical. If TRUE, returns a list with one vector per sequence (useful
  for sequences of different lengths). If FALSE, returns a matrix
  (requires all sequences to have the same length). Default is FALSE.

## Value

If return_list is FALSE: a matrix with `length(sequences)` rows and
`max(nchar(sequences))` columns with the local PWM for each sequence in
each position. If return_list is TRUE: a list with `length(sequences)`
elements, where each element is a numeric vector of local PWM scores for
the corresponding sequence.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])

# Return matrix (sequences must have same length)
pwm <- compute_local_pwm(cluster_sequences_example, res$pssm, res$spat)
head(pwm)

# Return list (allows sequences of different lengths)
pwm_list <- compute_local_pwm(cluster_sequences_example, res$pssm, res$spat, return_list = TRUE)
pwm_list[[1]]

# Using a motif from MOTIF_DB
hnf1a_pwm <- compute_local_pwm(
    cluster_sequences_example,
    as.matrix(MOTIF_DB["JASPAR.HNF1A"]),
    return_list = TRUE
)
hnf1a_pwm[[1]]
} # }
```
