# Extract pwm of intervals from a motif database

Extract the pwm of each interval for each motif from a motif database.
`gextract_pwm_old` is an older version of this function, which is
slower, and returns slightly different results due to float percision
instead of double.

## Usage

``` r
gextract_pwm_old(
  intervals,
  motifs = NULL,
  dataset = all_motif_datasets(),
  spat = NULL,
  spat_min = 1,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  func = "logSumExp",
  parallel = getOption("prego.parallel", TRUE)
)

gextract_pwm(
  intervals,
  motifs = NULL,
  dataset = MOTIF_DB,
  spat = NULL,
  spat_min = 1,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  func = "logSumExp",
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

- parallel:

  logical, whether to use parallel processing

## Value

The intervals set with additional columns per motif, containing the pwm
of each interval for each motif

## Examples

``` r
if (FALSE) { # \dontrun{
library(misha)
gdb.init_examples()
pwms <- gextract_pwm(gintervals.load("annotations"))
pwms[, 1:20]
} # }
```
