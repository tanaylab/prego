# Plot a motif from a MotifDB object

Plot a motif from a MotifDB object

## Usage

``` r
# S4 method for class 'MotifDB'
plot(
  x,
  title = "",
  subtitle = ggplot2::waiver(),
  revcomp = FALSE,
  method = "bits",
  force = FALSE,
  ...
)
```

## Arguments

- x:

  MotifDB object

- title:

  title of the plot

- subtitle:

  subtitle of the plot

- revcomp:

  whether to plot the reverse complement of the PSSM

- method:

  Height method, can be one of "bits" or "probability" (default:"bits")

- force:

  force plotting more than 30 motifs

- ...:

  additional arguments passed to
  [`ggseqlogo::ggseqlogo()`](https://rdrr.io/pkg/ggseqlogo/man/ggseqlogo.html)

## Value

a ggplot object

## Examples

``` r
plot(MOTIF_DB["HOMER.GATA3_2"])
#> Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
#> of ggplot2 3.3.4.
#> ℹ The deprecated feature was likely used in the ggseqlogo package.
#>   Please report the issue at <https://github.com/omarwagih/ggseqlogo/issues>.
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the ggseqlogo package.
#>   Please report the issue at <https://github.com/omarwagih/ggseqlogo/issues>.

plot(MOTIF_DB["HNF1", pattern = TRUE])

plot(MOTIF_DB[c("HOMER.GATA3_2", "JASPAR.CDX1")])

```
