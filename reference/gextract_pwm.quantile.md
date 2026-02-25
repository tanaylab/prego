# Extract quantiles of pwm of intervals from a motif database

Extract for each interval its quantile in the genome for each motif
given its length. Note that the quantiles are computed for each motif
separately, and therefore this might be slow for intervals with
un-normalized lengths.

## Usage

``` r
gextract_pwm.quantile(
  intervals,
  motifs = NULL,
  dataset = MOTIF_DB,
  percision = 0.01,
  spat = NULL,
  spat_min = 1,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  func = "logSumExp",
  n_sequences = 10000,
  dist_from_edge = 3000000,
  chromosomes = NULL,
  parallel = getOption("prego.parallel", TRUE)
)
```

## Arguments

- intervals:

  misha intervals set

- motifs:

  names of specific motifs to extract from the dataset

- dataset:

  a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an
  additional column 'motif' containing the motif name, for example
  `HOMER_motifs` or `JASPAR_motifs`, or
  [`all_motif_datasets()`](https://tanaylab.github.io/prego/reference/all_motif_datasets.md),
  or a MotifDB object.

- percision:

  the percision of the quantiles. Default is 0.01, which means that the
  quantiles will be computed for every 1% of the pwm.

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

- n_sequences:

  number of sequences to sample in order to compute the quantiles. The
  default is 1e4.

- dist_from_edge:

  The minimum distance from the edge of the chromosome for a region to
  start or end(default: 3e6)

- chromosomes:

  The chromosomes to sample from (default: all chromosomes)

- parallel:

  logical, whether to use parallel processing

## Value

a data frame with the quantiles of the pwm for each interval and motif.
The quantiles columns would be of the form {motif}.q

## Examples

``` r
if (FALSE) { # \dontrun{
library(misha)
gdb.init_examples()
gextract_pwm.quantile("annotations", motifs = c("JASPAR.CDX1", "JASPAR.CDX2"), dist_from_edge = 100)
} # }
```
