# Trim PSSM

This function trims a Position-Specific Scoring Matrix (PSSM) by
removing positions with low information content at the beginning and end
of the motif.

## Usage

``` r
trim_pssm(pssm, bits_thresh = 0.1)

pssm_trim(pssm, bits_thresh = 0.1)
```

## Arguments

- pssm:

  A data frame representing the PSSM, with columns for position (pos)
  and bits per position (bits).

- bits_thresh:

  The threshold value for bits per position. Positions with bits above
  this threshold will be kept, while positions with bits below this
  threshold at the beginning and the end of the motif will be removed.
  The default value is 0.1.

## Value

A trimmed PSSM data frame, with positions filtered based on the bits
threshold.
