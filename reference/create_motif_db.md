# Create a MotifDB object from a tidy data frame

Create a MotifDB object from a tidy data frame

## Usage

``` r
create_motif_db(
  motif_db,
  prior = 0.01,
  spat_factors = NULL,
  spat_bin_size = 1,
  spat_min = NA_real_,
  spat_max = NA_real_
)
```

## Arguments

- motif_db:

  A tidy data frame containing motif information

- prior:

  Pseudocount prior to add to probabilities (default: 0.01)

- spat_factors:

  Matrix of spatial factors (rows=motifs, cols=bins) or NULL

- spat_bin_size:

  Size of spatial bins (default: 1)

- spat_min:

  Starting position of sequence or NA (default: NA)

- spat_max:

  Ending position of sequence or NA (default: NA)

## Value

A MotifDB object

## Examples

``` r
create_motif_db(all_motif_datasets())
#> MotifDB object with 3867 motifs and prior 0.01
#> Slots include: @mat, @rc_mat, @motif_lengths, @prior, @spat_factors,
#> @spat_bin_size, @spat_min, @spat_max
```
