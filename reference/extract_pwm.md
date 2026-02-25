# Extract pwm of sequences from a motif database

Extracts the pwm of a motif from a motif database. `extract_pwm_old` is
a deprecated version of this function, which is slower, and returns
slightly different results due to float percision instead of double. If
the sequences are not of the same length, the old version will be used.

## Usage

``` r
extract_pwm_old(
  sequences,
  motifs = NULL,
  dataset = all_motif_datasets(),
  spat = NULL,
  spat_min = 0,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  func = "logSumExp",
  parallel = getOption("prego.parallel", TRUE)
)

extract_pwm(
  sequences,
  motifs = NULL,
  dataset = MOTIF_DB,
  spat = NULL,
  spat_min = NULL,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01,
  func = "logSumExp",
  parallel = getOption("prego.parallel", TRUE)
)
```

## Arguments

- sequences:

  a vector of sequences

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

a matrix size of \# of sequences x \# of motifs with the pwm of each
sequence for each motif

## Details

Unlike the old extraction functions (`extract_pwm_old`, `compute_pwm`,
and `compute_local_pwm`), `extract_pwm` returns NA values when the
sequence contains 'N' or '\*' characters.

## Examples

``` r
if (FALSE) { # \dontrun{
pwms <- extract_pwm(
    cluster_sequences_example,
    motifs = c("JASPAR.CDX1", "HOMER.Hnf1", "HOMER.GATA3_2")
)
head(pwms)

# all motifs
all_pwms <- extract_pwm(cluster_sequences_example, prior = 0.01)
dim(all_pwms)
all_pwms[1:5, 1:5]

# for a specific dataset
pwms_jaspar <- extract_pwm(cluster_sequences_example, dataset = JASPAR_motifs, prior = 0.01)
head(pwms_jaspar)

# for specific motifs
pwms_jaspar <- extract_pwm(
    cluster_sequences_example,
    motifs = c("JASPAR.CDX1", "JASPAR.CDX2"),
    prior = 0.01
)
} # }
```
