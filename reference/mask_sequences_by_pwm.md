# Mask sequences by thresholding the PWM

Mask sequences by thresholding the PWM. Sequences with a PWM above the
threshold will be masked by 'N'. Sequences at the edges of the sequences
will also be masked by 'N'.

## Usage

``` r
mask_sequences_by_pwm(
  sequences,
  pssm,
  mask_thresh,
  pos_bits_thresh = 0.2,
  spat = NULL,
  spat_min = 0,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01
)
```

## Arguments

- sequences:

  a vector of sequences

- pssm:

  a PSSM matrix or data frame. The columns of the matrix or data frame
  should be named with the nucleotides ('A', 'C', 'G' and 'T').

- mask_thresh:

  Threshold for masking. Sequences with a PWM above this threshold will
  be masked by 'N'.

- pos_bits_thresh:

  Mask only positions with amount of information contributed (Shannon
  entropy, measured in bits) above this threshold. The scale is the same
  as the y axis in the pssm logo plots.

- spat:

  a data frame with the spatial model (as returned from the `$spat` slot
  from the regression). Should contain a column called 'bin' and a
  column called 'spat_factor'.

- spat_min:

  the minimum position to use from the sequences. The default is 1.

- spat_max:

  the maximum position to use from the sequences. The default is the
  length of the sequences.

- bidirect:

  is the motif bi-directional. If TRUE, the reverse-complement of the
  motif will be used as well.

- prior:

  a prior probability for each nucleotide.

## Value

A vector with the masked sequences.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
new_sequences <- mask_sequences_by_pwm(
    cluster_sequences_example,
    res$pssm,
    quantile(res$pred, 0.95),
    spat = res$spat
)

head(new_sequences)
} # }
```
