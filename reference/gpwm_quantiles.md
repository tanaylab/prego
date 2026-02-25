# Compute quantile of pwm for a given interval size

Computes the quantile of the pwm for a given interval size by sampling
random intervals from the genome, or using given intervals. The number
of sequences to sample can be specified with `n_sequences`.

## Usage

``` r
gpwm_quantiles(
  size,
  quantiles,
  pssm,
  bg_intervals = NULL,
  spat = NULL,
  spat_min = 1,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  n_sequences = 10000,
  dist_from_edge = 3000000,
  chromosomes = NULL,
  func = "logSumExp"
)
```

## Arguments

- size:

  size of the intervals to sample

- quantiles:

  quantiles to compute. See `quantile` for more details.

- pssm:

  PSSM matrix or data frame

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

- func:

  the function to use to combine the PWMs for each sequence. Either
  'logSumExp' or 'max'. The default is 'logSumExp'.

## Value

a named vector with the quantiles of the pwm for the given interval
size.

## Examples

``` r
if (FALSE) { # \dontrun{
library(misha)
library(dplyr)
gdb.init_examples()
pssm <- JASPAR_motifs %>%
    filter(motif == "JASPAR.CDX1") %>%
    select(-motif)
gpwm_quantiles(1000, seq(0, 1, 0.1), pssm, dist_from_edge = 100)
} # }
```
