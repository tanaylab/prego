# Plot LOGO of pssm from dataset (e.g. "HOMER" or "JASPAR")

Plot LOGO of pssm from dataset (e.g. "HOMER" or "JASPAR")

## Usage

``` r
plot_pssm_logo_dataset(
  motif,
  dataset = all_motif_datasets(),
  title = motif,
  subtitle = ggplot2::waiver(),
  pos_bits_thresh = NULL,
  revcomp = FALSE,
  method = "bits"
)
```

## Arguments

- motif:

  the motif name (e.g. "GATA4")

- dataset:

  a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an
  additional column 'motif' containing the motif name, for example
  `HOMER_motifs`, `JASPAR_motifs` or all_motif_datasets()

- title:

  title of the plot

- subtitle:

  subtitle of the plot

- pos_bits_thresh:

  Positions with bits above this threshold would be highlighted in red.
  If `NULL`, no positions would be highlighted.

- revcomp:

  whether to plot the reverse complement of the PSSM

- method:

  Height method, can be one of "bits" or "probability" (default:"bits")

## Value

a ggplot object

## Examples

``` r
plot_pssm_logo_dataset("JASPAR.Brachyury")


plot_pssm_logo_dataset("GATA5", JASPAR_motifs)

```
