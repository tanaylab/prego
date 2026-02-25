# Get specific motifs from the MotifDB

Get specific motifs from the MotifDB

## Usage

``` r
# S4 method for class 'MotifDB'
x[i, j, ..., pattern = TRUE, drop = TRUE]
```

## Arguments

- x:

  MotifDB object

- i:

  Character vector of motif names, numeric indices, or regex pattern(s)

- j:

  Not used

- ...:

  Not used

- pattern:

  Logical indicating whether to treat character input as regex pattern
  (default: TRUE)

- drop:

  Not used

## Value

MotifDB object containing the specified motifs

## Examples

``` r
MOTIF_DB["HOMER.GATA3_2"]
#> MotifDB object with 1 motifs and prior 0.01
#> Slots include: @mat, @rc_mat, @motif_lengths, @prior, @spat_factors,
#> @spat_bin_size, @spat_min, @spat_max
MOTIF_DB[c("HOMER.GATA3_2", "JASPAR.CDX1")]
#> MotifDB object with 2 motifs and prior 0.01
#> Slots include: @mat, @rc_mat, @motif_lengths, @prior, @spat_factors,
#> @spat_bin_size, @spat_min, @spat_max
MOTIF_DB["GATA", pattern = TRUE]
#> MotifDB object with 44 motifs and prior 0.01
#> Slots include: @mat, @rc_mat, @motif_lengths, @prior, @spat_factors,
#> @spat_bin_size, @spat_min, @spat_max
```
