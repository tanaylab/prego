# Initialize method for MotifDB objects

Initialize method for MotifDB objects

## Usage

``` r
# S4 method for class 'MotifDB'
initialize(
  .Object,
  mat = matrix(0, 4, 1),
  rc_mat = matrix(0, 4, 1),
  motif_lengths = setNames(as.integer(1), colnames(mat)[1]),
  prior = 0.01,
  spat_factors = matrix(1, nrow = 1, ncol = 1),
  spat_bin_size = 1,
  spat_min = NA_real_,
  spat_max = NA_real_
)
```

## Arguments

- .Object:

  MotifDB object

- mat:

  Position weight matrix

- rc_mat:

  Reverse complement matrix

- motif_lengths:

  Named vector of motif lengths

- prior:

  Pseudocount prior value

- spat_factors:

  Matrix of spatial factors

- spat_bin_size:

  Size of spatial bins

- spat_min:

  Starting position of sequence or NULL

- spat_max:

  Ending position of sequence or NULL
