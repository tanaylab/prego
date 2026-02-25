# Compute KL divergence matrix for PSSM datasets

Compute KL divergence matrix for PSSM datasets

## Usage

``` r
pssm_dataset_diff(dataset1, dataset2 = NULL, prior = 0.01)
```

## Arguments

- dataset1:

  First PSSM dataset. A data frame with columns 'motif', 'A', 'C', 'G',
  'T'

- dataset2:

  Optional second PSSM dataset with the same structure as dataset1. If
  provided, computes divergences between motifs in dataset1 and
  dataset2. If NULL (default), computes divergences between all motifs
  in dataset1.

- prior:

  A prior probability added to each nucleotide frequency (default: 0.01)

## Value

If dataset2 is NULL, returns a symmetric square matrix of KL divergences
between all motifs in dataset1. If dataset2 is provided, returns a
matrix with rows corresponding to motifs in dataset1 and columns to
motifs in dataset2.

## Examples

``` r
if (FALSE) { # \dontrun{
# KL divergences within a single dataset
dm <- pssm_dataset_diff(JASPAR_motifs)

# KL divergences between two datasets
dm2 <- pssm_dataset_diff(JASPAR_motifs, HOMER_motifs)
} # }

# KL divergences within MOTIF_DB (subset)
gata_motifs <- as.data.frame(MOTIF_DB["GATA", pattern = TRUE])
dm_gata <- pssm_dataset_diff(gata_motifs)

# Compare GATA motifs with CDX motifs using KL divergence
cdx_motifs <- as.data.frame(MOTIF_DB["CDX", pattern = TRUE])
dm_cross <- pssm_dataset_diff(gata_motifs, cdx_motifs)
```
