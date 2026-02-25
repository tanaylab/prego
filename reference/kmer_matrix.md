# Generate a kmer Matrix

This function calculates the frequency of each kmer for each DNA
sequence.

## Usage

``` r
kmer_matrix(
  sequences,
  kmer_length,
  max_gap = 0,
  mask = NULL,
  add_mask = FALSE,
  from_range = 1,
  to_range = NULL,
  set_rownames = FALSE
)
```

## Arguments

- sequences:

  A vector of strings with DNA sequences ('T', 'C', 'G', 'A' or 'N').

- kmer_length:

  The length of the kmers to be considered.

- max_gap:

  The maximum length of a gap to be considered in the pattern. Default:
  0

- mask:

  a string the length of `kmer_length` where 'N' indicates a wildcard
  position (default NULL - no mask).

- add_mask:

  if TRUE, the result of the mask will be added to the non-masked kmers.
  Otherwise - only the masked kmers would be returned.

- from_range:

  Sequences will be considered only from position from_range.

- to_range:

  Sequences will be considered only up to position to_range (default
  NULL - using the length of the sequences).

- set_rownames:

  If TRUE, the rownames of the matrix will be set to the sequences
  (default FALSE).

## Value

A matrix where rows are the number of sequences, columns are the kmers
and the values are the number of occurrences of each kmer.

## Examples

``` r
kmer_matrix(c("ATCG", "TCGA", "ATAT"), 2)
#>      CG TC AT GA TA
#> [1,]  1  1  1  0  0
#> [2,]  1  1  0  1  0
#> [3,]  0  0  2  0  1
kmer_matrix(c("ATCG", "TCGA", "ATAT"), 3)
#>      TCG ATC CGA TAT ATA
#> [1,]   1   1   0   0   0
#> [2,]   1   0   1   0   0
#> [3,]   0   0   0   1   1
kmer_matrix(c("ATCG", "TCGA", "ATAT"), 3, mask = "ATN")
#>      TCN ATN CGN TAN
#> [1,]   1   1   0   0
#> [2,]   1   0   1   0
#> [3,]   0   1   0   1
```
