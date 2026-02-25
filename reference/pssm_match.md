# Match PSSM to a directory of motifs

Match a PSSM to a directory of motifs. The PSSM is matched to each motif
in the directory by computing the specified similarity measure.

## Usage

``` r
pssm_match(
  pssm,
  motifs,
  best = FALSE,
  method = c("spearman", "pearson", "kl"),
  prior = 0.01
)
```

## Arguments

- pssm:

  PSSM matrix or data frame with columns 'A', 'C', 'G', 'T'

- motifs:

  A data frame with PSSMs ('A', 'C', 'G', 'T' columns), with an
  additional column 'motif' containing the motif name

- best:

  Whether to return only the best match (default: FALSE)

- method:

  Method to use for matching: "spearman", "pearson", or "kl" for KL
  divergence. Note that the KL divergence is a measure of dissimilarity,
  so lower values indicate better matches.

- prior:

  A prior probability added to each nucleotide frequency

## Value

If best is TRUE, returns a string with the best matching motif name. If
best is FALSE, returns a data frame with columns 'motif' and either
'cor' or 'kl', sorted by decreasing similarity (for correlations) or
increasing divergence (for KL).

## Examples

``` r
if (FALSE) { # \dontrun{
res1 <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])

# Find the best match using Spearman correlation
pssm_match(res1$pssm, JASPAR_motifs, best = TRUE)

# Find all matches using Spearman correlation
matches <- pssm_match(res1$pssm, all_motif_datasets(), method = "spearman")
head(matches)

# Find best match using KL divergence
best_match <- pssm_match(res1$pssm, all_motif_datasets(), best = TRUE, method = "kl")
} # }

# Match a motif against a small motif subset
motif_subset <- as.data.frame(MOTIF_DB[c(
    "HOMER.GATA3_2",
    "JASPAR.GATA3",
    "JASPAR.CDX1",
    "JOLMA.GATA3_mono_DBD"
)])
matches <- pssm_match(as.matrix(MOTIF_DB["HOMER.GATA3_2"]), motif_subset)
head(matches)
#>                  motif       cor
#> 1        HOMER.GATA3_2 1.0000000
#> 2 JOLMA.GATA3_mono_DBD 0.9607420
#> 3         JASPAR.Gata3 0.3971397
#> 4          JASPAR.CDX1 0.3846587

if (FALSE) { # \dontrun{
# Find the best match for a GATA3 motif
best_match <- pssm_match(as.matrix(MOTIF_DB["HOMER.GATA3_2"]), as.data.frame(MOTIF_DB), best = TRUE)

# Use KL divergence for matching
kl_matches <- pssm_match(as.matrix(MOTIF_DB["JASPAR.CDX1"]), as.data.frame(MOTIF_DB), method = "kl")
head(kl_matches)
} # }
```
