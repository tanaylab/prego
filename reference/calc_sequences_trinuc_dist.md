# Calculate Trinucleotide Distribution in Sequences

Calculate Trinucleotide Distribution in Sequences

## Usage

``` r
calc_sequences_trinuc_dist(sequences, size = NULL)
```

## Arguments

- sequences:

  a character vector containing the sequences to analyze. Each element
  of the vector should be a single sequence.

- size:

  an integer specifying the size to consider for the analysis. If NULL
  (default), the maximum length of the sequences in the `sequences`
  vector is used.

## Value

a data frame with columns 'pos' and 64 columns representing each
possible trinucleotide. Each row represents a position in the sequences
(from 1 to `size`), and contains the fraction of each trinucleotide at
that position across all sequences.

## Examples

``` r
# Generate some random sequences for testing
set.seed(60427)
sequences <- sapply(1:100, function(x) {
    paste0(sample(c("A", "C", "G", "T"), 1000, replace = TRUE), collapse = "")
})
sequences <- as.character(sequences)

# Calculate the trinucleotide distribution
result <- calc_sequences_trinuc_dist(sequences)

head(result)
#>   pos  AAA  AAC  AAG  AAT  ACA  ACC  ACG  ACT  AGA  AGC  AGG  AGT  ATA  ATC
#> 1   1 0.02 0.02 0.04 0.02 0.01 0.04 0.01 0.00 0.02 0.02 0.00 0.00 0.03 0.04
#> 2   2 0.01 0.01 0.02 0.03 0.03 0.01 0.01 0.01 0.01 0.01 0.04 0.04 0.01 0.00
#> 3   3 0.00 0.01 0.01 0.03 0.02 0.04 0.01 0.02 0.02 0.01 0.01 0.01 0.03 0.01
#> 4   4 0.02 0.00 0.00 0.01 0.00 0.03 0.01 0.02 0.01 0.00 0.04 0.01 0.01 0.02
#> 5   5 0.02 0.03 0.04 0.02 0.00 0.01 0.02 0.01 0.01 0.01 0.01 0.00 0.02 0.00
#> 6   6 0.00 0.02 0.01 0.03 0.00 0.05 0.02 0.01 0.02 0.01 0.02 0.01 0.01 0.02
#>    ATG  ATT  CAA  CAC  CAG  CAT  CCA  CCC  CCG  CCT  CGA  CGC  CGG  CGT  CTA
#> 1 0.02 0.01 0.02 0.01 0.03 0.01 0.00 0.02 0.03 0.01 0.01 0.04 0.00 0.03 0.01
#> 2 0.02 0.01 0.02 0.00 0.00 0.00 0.02 0.02 0.03 0.01 0.01 0.02 0.02 0.02 0.00
#> 3 0.03 0.02 0.01 0.03 0.04 0.01 0.01 0.03 0.03 0.02 0.01 0.01 0.01 0.01 0.01
#> 4 0.02 0.01 0.04 0.01 0.01 0.01 0.01 0.02 0.05 0.02 0.00 0.00 0.02 0.03 0.01
#> 5 0.01 0.00 0.01 0.01 0.02 0.00 0.01 0.02 0.02 0.04 0.00 0.03 0.03 0.02 0.03
#> 6 0.04 0.00 0.00 0.00 0.01 0.02 0.01 0.02 0.02 0.02 0.01 0.02 0.01 0.01 0.03
#>    CTC  CTG  CTT  GAA  GAC  GAG  GAT  GCA  GCC  GCG  GCT  GGA  GGC  GGG  GGT
#> 1 0.01 0.03 0.02 0.00 0.02 0.01 0.00 0.00 0.00 0.02 0.01 0.02 0.01 0.03 0.00
#> 2 0.02 0.00 0.00 0.02 0.06 0.00 0.01 0.02 0.03 0.00 0.03 0.02 0.00 0.01 0.01
#> 3 0.01 0.02 0.02 0.02 0.01 0.00 0.02 0.03 0.02 0.00 0.01 0.02 0.02 0.03 0.03
#> 4 0.03 0.01 0.01 0.01 0.03 0.02 0.00 0.02 0.02 0.00 0.01 0.03 0.00 0.01 0.01
#> 5 0.00 0.06 0.00 0.02 0.01 0.00 0.03 0.00 0.01 0.00 0.00 0.02 0.03 0.03 0.02
#> 6 0.02 0.00 0.03 0.02 0.01 0.00 0.01 0.01 0.03 0.02 0.02 0.04 0.00 0.03 0.02
#>    GTA  GTC  GTG  GTT  TAA  TAC  TAG  TAT  TCA  TCC  TCG  TCT  TGA  TGC  TGG
#> 1 0.04 0.00 0.03 0.00 0.03 0.01 0.02 0.01 0.01 0.02 0.01 0.00 0.04 0.01 0.01
#> 2 0.01 0.01 0.01 0.02 0.00 0.02 0.03 0.05 0.02 0.03 0.00 0.01 0.01 0.03 0.03
#> 3 0.01 0.05 0.01 0.02 0.00 0.01 0.01 0.00 0.01 0.01 0.01 0.01 0.01 0.01 0.00
#> 4 0.02 0.00 0.01 0.04 0.04 0.00 0.00 0.01 0.01 0.02 0.02 0.04 0.02 0.01 0.03
#> 5 0.03 0.01 0.01 0.00 0.01 0.03 0.00 0.02 0.02 0.03 0.01 0.03 0.01 0.01 0.02
#> 6 0.00 0.02 0.02 0.01 0.00 0.02 0.04 0.02 0.00 0.01 0.00 0.00 0.02 0.04 0.05
#>    TGT  TTA  TTC  TTG  TTT
#> 1 0.02 0.02 0.01 0.01 0.00
#> 2 0.02 0.00 0.01 0.01 0.01
#> 3 0.02 0.00 0.02 0.00 0.02
#> 4 0.00 0.02 0.04 0.01 0.01
#> 5 0.01 0.00 0.00 0.03 0.04
#> 6 0.00 0.01 0.02 0.00 0.01
```
