# Compute the KL divergence between two given PSSMs

The KL divergence is computed by shifting the shorter PSSM along the
longer one and computing the divergence at each position. The minimum
divergence is returned.

## Usage

``` r
pssm_diff(pssm1, pssm2, prior = 0.01)
```

## Arguments

- pssm1:

  first PSSM matrix or data frame

- pssm2:

  second PSSM matrix or data frame

- prior:

  a prior probability for each nucleotide.

## Value

KL divergence between the two PSSMs (lower values indicate more similar
PSSMs)

## Examples

``` r
if (FALSE) { # \dontrun{
res1 <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
pssm_diff(res1$pssm, JASPAR_motifs[JASPAR_motifs$motif == "HNF1A", ])
} # }

# Compare KL divergence between two motifs from MOTIF_DB
pssm_diff(as.matrix(MOTIF_DB["HOMER.GATA3_2"]), as.matrix(MOTIF_DB["JASPAR.CDX1"]))
#> [1] 8.327341

# Lower values indicate more similar motifs
pssm_diff(as.matrix(MOTIF_DB["HOMER.GATA3_2"]), as.matrix(MOTIF_DB["JASPAR.GATA3"]))
#> [1] 9.37645
```
