# Calculate the number of bits per position in a Position-Specific Scoring Matrix (PSSM).

This function takes a PSSM as input and calculates the number of bits
per position. The PSSM should be a data frame or matrix with columns
representing the nucleotides A, C, G, and T. The function first
normalizes the PSSM by dividing each element by the sum of its row.
Then, it calculates the entropy for each position using the formula:
bits = log2(4) + sum(p \* log2(p)), where p is the probability of each
nucleotide at the position. Finally, it sets any negative values to zero
and returns the resulting bits per position.

## Usage

``` r
bits_per_pos(pssm, prior = 0.01)
```

## Arguments

- pssm:

  A data frame or matrix representing the Position-Specific Scoring
  Matrix (PSSM).

- prior:

  A numeric value indicating the prior probability for each nucleotide.
  Default is 0.01.

## Value

A numeric vector representing the number of bits per position in the
PSSM.

## Examples

``` r
pssm <- data.frame(
    A = c(0.2, 0.3, 0.1, 0.4),
    C = c(0.1, 0.2, 0.3, 0.4),
    G = c(0.4, 0.3, 0.2, 0.1),
    T = c(0.3, 0.2, 0.4, 0.1)
)
bits_per_pos(pssm)
#> [1] 0.14121539 0.02684396 0.14121539 0.25558682
```
