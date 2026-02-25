# Transform k-mers to PSSM (Position-Specific Scoring Matrix)

This function transforms a vector of k-mers into a position-specific
scoring matrix (PSSM). A PSSM represents the frequency of each
nucleotide at each position in the k-mers. If a nucleotide is 'N', it is
treated as equal probabilities for 'A', 'C', 'G', and 'T'. The result is
returned as a data frame with columns for the k-mer, position, and
nucleotide frequencies.

## Usage

``` r
kmers_to_pssm(kmers, prior = 0.01)
```

## Arguments

- kmers:

  A character vector of k-mers.

- prior:

  A numeric value indicating the prior probability for each nucleotide.
  Default is 0.01.

## Value

A data frame with columns for the k-mer, position, and nucleotide
frequencies, 'kmer', 'pos', 'A', 'C', 'G', 'T'.

## Examples

``` r
kmers_to_pssm(c("ACGTN", "TGCAN"), prior = 0.01)
#>     kmer pos           A           C           G           T
#> 1  ACGTN   1 0.970873786 0.009708738 0.009708738 0.009708738
#> 2  ACGTN   2 0.009708738 0.970873786 0.009708738 0.009708738
#> 3  ACGTN   3 0.009708738 0.009708738 0.970873786 0.009708738
#> 4  ACGTN   4 0.009708738 0.009708738 0.009708738 0.970873786
#> 5  ACGTN   5 0.250000000 0.250000000 0.250000000 0.250000000
#> 6  TGCAN   1 0.009708738 0.009708738 0.009708738 0.970873786
#> 7  TGCAN   2 0.009708738 0.009708738 0.970873786 0.009708738
#> 8  TGCAN   3 0.009708738 0.970873786 0.009708738 0.009708738
#> 9  TGCAN   4 0.970873786 0.009708738 0.009708738 0.009708738
#> 10 TGCAN   5 0.250000000 0.250000000 0.250000000 0.250000000
```
