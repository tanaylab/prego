# Convert PSSM to consensus sequence

Convert PSSM to consensus sequence

## Usage

``` r
consensus_from_pssm(pssm, single_thresh = 0.4, double_thresh = 0.6)
```

## Arguments

- pssm:

  A PSSM matrix

- single_thresh, double_thresh:

  thresholds for the consensus sequence calculation (single and double
  nucleotides)

## Value

A consensus sequence for the PSSM. If no consensus sequence can be
found, the function returns NA.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
consensus_from_pssm(res$pssm)
} # }
```
