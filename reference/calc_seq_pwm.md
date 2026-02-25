# Calculate Position Weight Matrix (PWM) Scores for DNA Sequences

Calculate Position Weight Matrix (PWM) Scores for DNA Sequences

## Usage

``` r
calc_seq_pwm(sequences, mdb, bidirect = TRUE)
```

## Arguments

- sequences:

  Character vector of DNA sequences.

- mdb:

  MotifDB object containing PWMs.

- bidirect:

  is the motif bi-directional. If TRUE, the reverse-complement of the
  motif will be used as well.

## Value

A numeric matrix with sequences as rows and motifs as columns,
containing PWM scores. Row names are preserved from input sequences if
they exist. Column names are preserved from the PWM matrix if they
exist.

## Examples

``` r
sequences <- c("ACGTACGT", "TGCATGCA")
scores <- calc_seq_pwm(sequences, MOTIF_DB)
```
