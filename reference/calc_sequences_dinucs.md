# Calculate Dinucleotide Counts for Sequences

This function calculates the total count of each dinucleotide for each
sequence in a vector of DNA sequences.

## Usage

``` r
calc_sequences_dinucs(sequences)
```

## Arguments

- sequences:

  A character vector of DNA sequences. Each element should be a string
  representing a DNA sequence composed of A, T, C, and G.

## Value

A numeric matrix where:

- Each row corresponds to a sequence in the input vector.

- Each column represents a specific dinucleotide (AA, AC, AG, AT, CA,
  CC, etc.).

- The values in the matrix are the counts of each dinucleotide in each
  sequence.

- Column names are set to the corresponding dinucleotides.

## Examples

``` r
sequences <- c("ATCG", "GCTA", "AATT")
result <- calc_sequences_dinucs(sequences)
print(result)
#>      AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
#> [1,]  0  0  0  1  0  0  1  0  0  0  0  0  0  1  0  0
#> [2,]  0  0  0  0  0  0  0  1  0  1  0  0  1  0  0  0
#> [3,]  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1
```
