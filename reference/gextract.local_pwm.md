# Extracts local position weight matrix (PWM) scores for given intervals and a PWM.

Extracts local position weight matrix (PWM) scores for given intervals
and a PWM.

## Usage

``` r
gextract.local_pwm(
  intervals,
  pssm,
  spat = NULL,
  spat_min = 0,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01
)
```

## Arguments

- intervals:

  The intervals to extract

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

## Value

A matrix with `nrow(intervals)` rows and `ncol(pssm)` columns with the
local PWM for each sequence in each position.
