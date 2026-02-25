# Sample rows respecting quantiles of a reference distribution

This function randomly samples rows from a data frame in such a way that
the quantiles of the selected data match as closely as possible those of
the full data.

## Usage

``` r
sample_quantile_matched_rows(
  data_frame,
  reference,
  sample_fraction,
  num_quantiles = 10,
  seed = 60427,
  verbose = TRUE
)
```

## Arguments

- data_frame:

  A data frame from which to sample rows.

- reference:

  A numeric vector of the same length as the number of rows in the data
  frame.

- sample_fraction:

  A fraction specifying the proportion of rows to sample from the data
  frame.

- num_quantiles:

  An integer specifying the number of quantiles, default is 10.

- seed:

  An integer specifying the random seed to use.

- verbose:

  A logical specifying whether to print messages.

## Value

A data frame of sampled rows.

## Examples

``` r
sampled <- sample_quantile_matched_rows(mtcars, mtcars$mpg, sample_fraction = 0.1)
#> ✔ Sampled 3 rows from the data frame.
plot(quantile(mtcars$mpg), quantile(sampled$mpg))
abline(0, 1)

```
