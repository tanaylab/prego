# Screen sequences for positions with PWM scores meeting a threshold condition

This function screens sequences to find positions where the local PWM
score meets a specified threshold condition using various operators (\>,
\<, \>=, \<=, ==).

## Usage

``` r
screen_local_pwm(
  sequences,
  pssm,
  operator,
  threshold,
  spat = NULL,
  spat_min = 0,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01
)
```

## Arguments

- sequences:

  a vector of sequences

- pssm:

  a PSSM matrix or data frame. The columns of the matrix or data frame
  should be named with the nucleotides ('A', 'C', 'G' and 'T').

- operator:

  A character string specifying the comparison operator. One of: "\>",
  "\<", "\>=", "\<=", "=="

- threshold:

  A numeric value specifying the threshold for comparison

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

## Value

A list with one element per sequence, where each element is a numeric
vector containing the 1-indexed positions that meet the threshold
condition.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find positions where HNF1A motif scores above -19
hnf1a_positions <- screen_local_pwm(
    cluster_sequences_example,
    as.matrix(MOTIF_DB["JASPAR.HNF1A"]),
    operator = ">",
    threshold = -37
)

which(purrr::map_dbl(hnf1a_positions, length) > 0) # which sequences have positions?
hnf1a_positions[[4]] # positions for the 4th sequence

hnf1a_energy <- compute_local_pwm(cluster_sequences_example, as.matrix(MOTIF_DB["JASPAR.HNF1A"]))
hnf1a_energy[4, hnf1a_positions[[4]]]
} # }
```
