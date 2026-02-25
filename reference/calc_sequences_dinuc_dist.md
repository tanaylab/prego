# Calculate Dinucleotide Distribution in Sequences

Calculate Dinucleotide Distribution in Sequences

## Usage

``` r
calc_sequences_dinuc_dist(sequences, size = NULL)
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

a data frame with columns 'pos' and 16 columns representing each
possible dinucleotide. Each row represents a position in the sequences
(from 1 to `size`), and contains the fraction of each dinucleotide at
that position across all sequences.

## Examples

``` r
# Generate some random sequences for testing
set.seed(60427)
sequences <- sapply(1:100, function(x) {
    paste0(sample(c("A", "C", "G", "T"), 1000, replace = TRUE), collapse = "")
})
sequences <- as.character(sequences)

# Calculate the dinucleotide distribution
result <- calc_sequences_dinuc_dist(sequences)

head(result)
#>   pos   AA   AC   AG   AT   CA   CC   CG   CT   GA   GC   GG   GT   TA   TC
#> 1   1 0.10 0.06 0.04 0.10 0.07 0.06 0.08 0.07 0.03 0.03 0.06 0.07 0.07 0.04
#> 2   2 0.07 0.06 0.10 0.04 0.02 0.08 0.07 0.02 0.09 0.08 0.04 0.05 0.10 0.06
#> 3   3 0.05 0.09 0.05 0.09 0.09 0.09 0.04 0.06 0.05 0.06 0.10 0.09 0.02 0.04
#> 4   4 0.03 0.06 0.06 0.06 0.07 0.10 0.05 0.06 0.06 0.05 0.05 0.07 0.05 0.09
#> 5   5 0.11 0.04 0.03 0.03 0.04 0.09 0.08 0.09 0.06 0.01 0.10 0.05 0.06 0.09
#> 6   6 0.06 0.08 0.06 0.07 0.03 0.07 0.05 0.08 0.04 0.08 0.09 0.05 0.08 0.01
#>     TG   TT
#> 1 0.08 0.04
#> 2 0.09 0.03
#> 3 0.04 0.04
#> 4 0.06 0.08
#> 5 0.05 0.07
#> 6 0.11 0.04
```
