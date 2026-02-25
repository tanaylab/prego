# Calculate the frequency of a position weight matrix (PWM) in a given set of intervals

Calculate the frequency of a position weight matrix (PWM) in a given set
of intervals

## Usage

``` r
gextract.local_pwm_freq(
  intervals,
  pssm,
  q_threshold,
  bg_intervals = NULL,
  spat = NULL,
  spat_min = 0,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  n_sequences = 10000,
  dist_from_edge = 3000000,
  chromosomes = NULL
)
```

## Arguments

- intervals:

  The intervals to extract

- pssm:

  a PSSM matrix or data frame. The columns of the matrix or data frame
  should be named with the nucleotides ('A', 'C', 'G' and 'T').

- q_threshold:

  The quantile threshold of the PWM (e.g. 0.99 for the top percentile)

- bg_intervals:

  (optional) an intervals set for the background. If not provided,
  random intervals will be used

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

- n_sequences:

  number of sequences to sample in order to compute the quantiles. The
  default is 1e4.

- dist_from_edge:

  The minimum distance from the edge of the chromosome for a region to
  start or end(default: 3e6)

- chromosomes:

  The chromosomes to sample from (default: all chromosomes)

## Value

a matrix with `nrow(intervals)` rows and `ncol(pssm)` columns with the
TRUE if the PWM is above the threshold for each sequence in each
position.
