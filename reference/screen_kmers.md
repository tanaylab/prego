# Screen for kmers

Screen for kmers

## Usage

``` r
screen_kmers(
  sequences,
  response,
  kmer_length = 6,
  min_cor = 0.08,
  is_train = NULL,
  min_gap = 0,
  max_gap = 0,
  from_range = 0,
  to_range = NULL,
  return_mat = FALSE,
  seed = 60427,
  verbose = FALSE
)
```

## Arguments

- sequences:

  A vector of DNA sequences ('A', 'T', 'C' or 'G'. Will go through
  `toupper`)

- response:

  A matrix of response variables - number of rows should equal the
  number of sequences

- kmer_length:

  The number of non-gap characters in motifs that will be screened

- min_cor:

  Only patterns for which the maximum correlation to one of the response
  variable is larger than min_cor will be reported

- is_train:

  a boolean vector that determine which subset of sequences to use when
  screening

- min_gap, max_gap:

  the length of a gap to be considered in the pattern. Only one gap, of
  length min_gap:max_gap, is being used, and is located anywhere in the
  motif. Note that this greatly expand the search space (and increase
  multiple testing severely).

- from_range:

  Sequences will be considered only from position from_range (default 0)

- to_range:

  Sequences will be considered only up to position to_range (default
  NULL - using the length of the sequences)

- return_mat:

  Return a matrix of patterns and their correlation to the response
  variables instead of a data frame. (default: FALSE)

- seed:

  random seed

- verbose:

  show verbose messages

## Value

A data frame with the following columns, together with a column for each
response variable with the correlation of the kmers to the response
variable:

- kmer: :

  the kmer pattern, where "\*" indicates a wildcard

- max_r2: :

  the maximum R^2 to one of the response variables

- avg_n: :

  the average number of times the kmer appears in the sequences

- avg_var: :

  the variance of the number of times the kmer appears in the sequences

if `return_mat` is TRUE, a matrix with correlations to the response
variables (where rows are the kmers) is returned instead of a data
frame. If no kmer is found, an empty data frame is returned.

## Examples

``` r
kmers <- screen_kmers(sequences_example, response_mat_example)
#> ℹ Number of response variables: 5
#> ℹ Screening kmers of length 6, from position 0 to position 300
#> ℹ minimal correlation: 0.08
#> ✔ Found 575 kmers in 1000 sequences.
head(kmers)
#>     kmer     max_r2 avg_n    avg_var          c1          c2         c3
#> 1 AGATAA 0.04998882 0.140 0.15240000 -0.01025383 -0.03046725 -0.2235818
#> 2 TTATCT 0.02975253 0.128 0.13761599  0.03136108  0.04927359 -0.1724892
#> 3 CTTATC 0.02680831 0.080 0.07960001  0.02861430  0.01967959 -0.1637324
#> 4 GGGGAG 0.02300140 0.243 0.28195098  0.04611732  0.02960520  0.1516621
#> 5 GGGCGG 0.02280306 0.062 0.06015600  0.06576437  0.02662074  0.1510068
#> 6 TAACTG 0.02215597 0.071 0.06995900 -0.01145993 -0.02249238 -0.1488488
#>             c4          c5
#> 1  0.004382457 -0.03498838
#> 2  0.035843253  0.02917432
#> 3  0.059437271  0.03772539
#> 4  0.044013925  0.03089626
#> 5  0.104473867  0.12732512
#> 6 -0.037136175 -0.01858661

kmers <- screen_kmers(sequences_example, response_mat_example, return_mat = TRUE)
#> ℹ Number of response variables: 5
#> ℹ Screening kmers of length 6, from position 0 to position 300
#> ℹ minimal correlation: 0.08
#> ✔ Found 575 kmers in 1000 sequences.
head(kmers)
#>                 c1          c2         c3           c4          c5
#> AGATAA -0.01025383 -0.03046725 -0.2235818  0.004382457 -0.03498838
#> TTATCT  0.03136108  0.04927359 -0.1724892  0.035843253  0.02917432
#> CTTATC  0.02861430  0.01967959 -0.1637324  0.059437271  0.03772539
#> GGGGAG  0.04611732  0.02960520  0.1516621  0.044013925  0.03089626
#> GGGCGG  0.06576437  0.02662074  0.1510068  0.104473867  0.12732512
#> TAACTG -0.01145993 -0.02249238 -0.1488488 -0.037136175 -0.01858661

kmers <- screen_kmers(sequences_example, response_mat_example, max_gap = 3)
#> ℹ Number of response variables: 5
#> ℹ Screening kmers of length 6, from position 0 to position 300
#> ℹ Gaps of length 0:3 are allowed
#> ℹ minimal correlation: 0.08
#> ✔ Found 4884 kmers in 1000 sequences.
head(kmers)
#>       kmer     max_r2 avg_n    avg_var           c1           c2         c3
#> 1   AGATAA 0.05012437 0.139 0.15167901 -0.011159671 -0.033565160 -0.2238847
#> 2  AGA*AAG 0.03404951 0.220 0.23760001  0.062884174  0.063842162 -0.1845251
#> 3   TTATCT 0.02915060 0.127 0.13487099  0.036021627  0.051995803 -0.1707355
#> 4 ACAT**CT 0.02848161 0.090 0.09190000 -0.002383306 -0.023468828 -0.1687650
#> 5  AG*TAAG 0.02829135 0.118 0.12007600  0.034630746 -0.005794165 -0.1682003
#> 6   CTTATC 0.02680831 0.080 0.07960001  0.028614303  0.019679585 -0.1637324
#>            c4          c5
#> 1 0.001539792 -0.03616555
#> 2 0.025056362  0.02373981
#> 3 0.032191906  0.02568413
#> 4 0.019714575 -0.02059141
#> 5 0.046150818  0.02593128
#> 6 0.059437271  0.03772539
```
