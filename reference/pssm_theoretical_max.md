# Calculate Theoretical Scores for a Position-Specific Scoring Matrix (PSSM)

These functions compute theoretical scores (maximum, minimum, and
quantiles) that can be achieved by sequences matching a given
Position-Specific Scoring Matrix (PSSM).

## Usage

``` r
pssm_theoretical_max(pssm, prior = 0.01, regularization = 0.01)

pssm_theoretical_min(pssm, prior = 0.01, regularization = 0.01)

pssm_quantile(pssm, q, prior = 0.01, regularization = 0.01)
```

## Arguments

- pssm:

  A matrix or data frame containing position-specific probabilities for
  nucleotides A, C, G, T

- prior:

  A numeric value (default: 0.01) added to each probability to avoid
  zero probabilities

- regularization:

  A numeric value (default: 0.01) added inside the log function to
  prevent -Inf values

- q:

  A numeric value between 0 and 1 specifying the quantile to calculate
  (only for pssm_quantile)

## Value

- pssm_theoretical_max: Maximum possible score for the given PSSM

- pssm_theoretical_min: Minimum possible score for the given PSSM

- pssm_quantile: Score at the specified quantile between min and max
  scores

## Details

The functions normalize the probabilities and calculate scores based on
logarithms with regularization. The quantile function interpolates
linearly between min and max scores.

## Examples

``` r
pssm <- as.matrix(MOTIF_DB["HOMER.CTCF"])
pssm_theoretical_max(pssm)
#> [1] -17.47226
pssm_theoretical_min(pssm)
#> [1] -123.9643
pssm_quantile(pssm, 0.85)
#> [1] -33.44606
```
