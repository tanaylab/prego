# Compute correlation matrix for PSSM datasets

Compute correlation matrix for PSSM datasets

## Usage

``` r
pssm_dataset_cor(
  dataset1,
  dataset2 = NULL,
  method = c("spearman", "pearson"),
  prior = 0.01
)
```

## Arguments

- dataset1:

  First PSSM dataset. A data frame with columns 'motif', 'A', 'C', 'G',
  'T'

- dataset2:

  Optional second PSSM dataset with the same structure as dataset1. If
  provided, computes correlations between motifs in dataset1 and
  dataset2. If NULL (default), computes correlations between all motifs
  in dataset1.

- method:

  Method to use for correlation calculation ("spearman" or "pearson")

- prior:

  A prior probability added to each nucleotide frequency (default: 0.01)

## Value

If dataset2 is NULL, returns a symmetric square matrix of correlations
between all motifs in dataset1. If dataset2 is provided, returns a
matrix with rows corresponding to motifs in dataset1 and columns to
motifs in dataset2.

## Examples

``` r
if (FALSE) { # \dontrun{
# Correlations within a single dataset
cm <- pssm_dataset_cor(JASPAR_motifs)

# Correlations between two datasets
cm2 <- pssm_dataset_cor(JASPAR_motifs, HOMER_motifs)
} # }

# Correlations within MOTIF_DB (subset)
gata_motifs <- as.data.frame(MOTIF_DB["GATA", pattern = TRUE])
cm_gata <- pssm_dataset_cor(gata_motifs)

# Compare GATA motifs with CDX motifs
cdx_motifs <- as.data.frame(MOTIF_DB["CDX", pattern = TRUE])
cm_cross <- pssm_dataset_cor(gata_motifs, cdx_motifs)
```
