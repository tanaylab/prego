# Expected local PWM scores over a base frequency matrix

Score every motif in a database against a per-position base frequency
matrix at every start position. Where
[`compute_local_pwm`](https://tanaylab.github.io/prego/reference/compute_local_pwm.md)
scores one concrete sequence, this scores an ensemble of sequences
summarised by its per-position nucleotide distribution.

## Usage

``` r
calc_freq_local_pwm(
  freqs,
  motifs,
  combine = c("multiply", "sum"),
  bidirect = TRUE
)
```

## Arguments

- freqs:

  Base frequencies, in any of:

  - a single matrix, either `m x 4` or `4 x m`, in A,C,G,T order.
    Returns a single matrix rather than a list.

  - a list of such matrices, which may differ in length.

  - a 3-D array, `n_matrices x m x 4` or `n_matrices x 4 x m`.

  - a stacked matrix, `n_matrices x 4m`, the layout
    [`calc_seq_pwm`](https://tanaylab.github.io/prego/reference/calc_seq_pwm.md)
    takes. To stack exactly 4 matrices use a list instead - a 4-row
    matrix is read as a single `4 x m` one.

  Each position must be a distribution, i.e. sum to 1.

- motifs:

  A `MotifDB` object, or a tidy motif data frame with columns `motif`,
  `pos`, `A`, `C`, `G`, `T`.

- combine:

  How to combine the per-column products, `"multiply"` (the default) or
  `"sum"`. See Details.

- bidirect:

  Score both strands, combining them at each position as
  `log(exp(fwd) + exp(rev))` - the same convention as
  [`compute_local_pwm`](https://tanaylab.github.io/prego/reference/compute_local_pwm.md).
  `FALSE` scores the forward strand only.

## Value

A numeric matrix with one row per motif and one column per position of
the input, holding the score of the motif *starting* at that position.
The last `L - 1` entries of each row are `NA`, where `L` is that motif's
length. A list of such matrices if `freqs` held more than one.

## Details

For a motif of length \\L\\ placed at start \\j\\, with motif column
\\p_l\\ and frequency column \\q\_{j+l}\\, there are two ways to combine
each pair of distributions. They differ only in where the log sits:

- `combine = "multiply"`:

  \\\sum_l \log(q\_{j+l} \cdot p_l)\\. Multiplies the per-column
  probabilities, i.e. the log of the expected likelihood. On a flat
  ensemble every motif gets the same floor, \\L \log(0.25)\\, so scores
  are comparable across motifs. This is exact only if the positions of
  the ensemble are independent - the usual PWM assumption, but one that
  an ensemble with covarying positions violates.

- `combine = "sum"`:

  \\\sum_l q\_{j+l} \cdot \log p_l\\. Sums the per-column
  log-probabilities, i.e. the expected log-likelihood - the mean score
  you would get by drawing sequences from the ensemble and running
  [`compute_local_pwm`](https://tanaylab.github.io/prego/reference/compute_local_pwm.md)
  on each. Being linear in the frequencies, it is exact whatever the
  joint distribution of positions is. Its downside is that a flat
  ensemble gives each motif a different floor, strongly anti-correlated
  with the motif's information content, so scores are not comparable
  across motifs without normalisation.

Both reduce to
[`compute_local_pwm`](https://tanaylab.github.io/prego/reference/compute_local_pwm.md)
when the frequency matrix is one-hot, i.e. when the ensemble is a single
sequence.

## See also

[`compute_local_pwm`](https://tanaylab.github.io/prego/reference/compute_local_pwm.md),
[`calc_seq_pwm`](https://tanaylab.github.io/prego/reference/calc_seq_pwm.md)

## Examples

``` r
db <- all_motif_datasets()
mdb <- create_motif_db(db[db$motif %in% head(unique(db$motif), 20), ])

# An ensemble that is certain at every position is just a sequence, and
# scores the same as compute_local_pwm() does on that sequence.
seq <- "ACGTACGTTGCAAGGTCCATACGT"
onehot <- diag(4)[match(strsplit(seq, "")[[1]], c("A", "C", "G", "T")), ]
scores <- calc_freq_local_pwm(onehot, mdb)
dim(scores)
#> [1] 20 24
scores[1:3, 1:5]
#>                      [,1]      [,2]      [,3]      [,4]      [,5]
#> HOMER.AP_1      -25.29433 -22.01231 -22.92823 -21.10492 -21.10696
#> HOMER.AP_2gamma -31.01041 -27.36552 -24.02633 -23.31629 -24.95362
#> HOMER.AP_2alpha -30.23389 -27.48152 -24.45515 -23.79706 -25.29995

# Under "multiply", a flat ensemble scores L * log(0.25) for every motif -
# the floor depends on the motif's length, not on the motif.
flat <- matrix(0.25, nrow = 24, ncol = 4)
calc_freq_local_pwm(flat, mdb, bidirect = FALSE)[1:3, 1]
#>      HOMER.AP_1 HOMER.AP_2gamma HOMER.AP_2alpha 
#>       -13.86294       -16.63553       -16.63553 
mdb@motif_lengths[1:3] * log(0.25)
#>      HOMER.AP_1 HOMER.AP_2gamma HOMER.AP_2alpha 
#>       -13.86294       -16.63553       -16.63553 

# Several ensembles in one call
res <- calc_freq_local_pwm(list(certain = onehot, flat = flat), mdb)
sapply(res, dim)
#>      certain flat
#> [1,]      20   20
#> [2,]      24   24
```
