# Compute the correlation between two given PSSMs

The correlation is computed by shifting the shorter PSSM along the
longer one and computing the correlation at each position. The maximum
correlation is returned.

## Usage

``` r
pssm_cor(pssm1, pssm2, method = c("spearman", "pearson"), prior = 0.01)
```

## Arguments

- pssm1:

  first PSSM matrix or data frame

- pssm2:

  second PSSM matrix or data frame

- method:

  method to use for computing the correlation. See
  [`cor`](https://rdrr.io/r/stats/cor.html) for details.

- prior:

  a prior probability for each nucleotide.

## Value

Correlation between the two PSSMs

## Examples

``` r
if (FALSE) { # \dontrun{
res1 <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
pssm_cor(res1$pssm, JASPAR_motifs[JASPAR_motifs$motif == "HNF1A", ])
} # }

# Compare two motifs from MOTIF_DB
pssm_cor(as.matrix(MOTIF_DB["HOMER.GATA3_2"]), as.matrix(MOTIF_DB["JASPAR.CDX1"]))
#> [1] 0.3846587

# Compare using different correlation methods
pssm_cor(
    as.matrix(MOTIF_DB["HOMER.GATA3_2"]),
    as.matrix(MOTIF_DB["JASPAR.GATA3"]),
    method = "pearson"
)
#> [1] 0.2425705
pssm_cor(
    as.matrix(MOTIF_DB["HOMER.GATA3_2"]),
    as.matrix(MOTIF_DB["JASPAR.GATA3"]),
    method = "spearman"
)
#> [1] 0.3971397
```
