# Transform PSSM (Position-Specific Scoring Matrix) to a KMER

This function transforms a PSSM into a k-mer of a given length.

## Usage

``` r
pssm_to_kmer(pssm, kmer_length = NULL, pos_bits_thresh = NULL, prior = 0.01)
```

## Arguments

- pssm:

  PSSM matrix or data frame. The PSSM must have at least kmer_length
  rows.

- kmer_length:

  The length of the k-mer to return. If NULL - the length of the k-mer
  is equal to the number of rows in the PSSM.

- pos_bits_thresh:

  A numeric value indicating the minimum number of bits per position to
  include the nucleotide in the k-mer. If the nucleotide does not meet
  this threshold, it is replaced with 'N'. Default is NULL.

- prior:

  A numeric value indicating the prior probability for each nucleotide.
  Default is 0.01.

## Value

A character vector of length 1 containing the k-mer.

## Examples

``` r
pssm_to_kmer(get_motif_pssm("HOMER.AP_1"))
#> [1] "ATGACTCATC"
plot_pssm_logo_dataset("HOMER.AP_1")

```
