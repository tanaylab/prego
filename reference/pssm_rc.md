# Reverse complement a PSSM

Reverse complement a PSSM

## Usage

``` r
pssm_rc(pssm)
```

## Arguments

- pssm:

  A PSSM. Data frame with columns 'A', 'C', 'G', 'T' and 'pos' or a
  matrix with columns 'A', 'C', 'G', 'T'

## Value

A PSSM with the same format, but reverse complemented.

## Examples

``` r
# Create simulated PSSM data frame
pssm <- data.frame(
    pos = 1:4,
    A = c(0.1, 0.2, 0.3, 0.1),
    C = c(0.1, 0.3, 0.2, 0.1),
    G = c(0.1, 0.3, 0.3, 0.7),
    T = c(0.7, 0.2, 0.2, 0.1)
)

# Reverse complement the PSSM
rc_pssm <- pssm_rc(pssm)
```
