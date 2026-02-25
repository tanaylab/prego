# Center intervals by PSSM

This function takes a set of intervals and a position-specific scoring
matrix (PSSM) and centers the intervals based on the maximum score
position in the PSSM. The intervals are shifted so that the maximum
score position becomes the center of each interval.

## Usage

``` r
gintervals.center_by_pssm(
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

A data frame containing the centered intervals. The intervals will have
the same columns as the input intervals, but the start and end positions
will be adjusted to center the intervals based on the maximum score
position in the PSSM.
