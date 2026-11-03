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
#> 1 AGATAA 0.04998422 0.140 0.15240000 -0.01025364 -0.03046969 -0.2235715
#> 2 TTATCT 0.02974978 0.128 0.13761599  0.03136203  0.04927634 -0.1724813
#> 3 CTTATC 0.02680594 0.080 0.07960001  0.02861515  0.01968119 -0.1637252
#> 4 GGGGAG 0.02299936 0.243 0.28195098  0.04611842  0.02960637  0.1516554
#> 5 GGGCGG 0.02280105 0.062 0.06015600  0.06576587  0.02662213  0.1510002
#> 6 TAACTG 0.02215385 0.071 0.06995900 -0.01146036 -0.02249384 -0.1488417
#>            c4          c5
#> 1  0.00438225 -0.03498891
#> 2  0.03583921  0.02917477
#> 3  0.05942981  0.03772579
#> 4  0.04400812  0.03089627
#> 5  0.10446008  0.12732546
#> 6 -0.03713173 -0.01858625

kmers <- screen_kmers(sequences_example, response_mat_example, return_mat = TRUE)
#> ℹ Number of response variables: 5
#> ℹ Screening kmers of length 6, from position 0 to position 300
#> ℹ minimal correlation: 0.08
#> ✔ Found 575 kmers in 1000 sequences.
head(kmers)
#>                 c1          c2         c3          c4          c5
#> AGATAA -0.01025364 -0.03046969 -0.2235715  0.00438225 -0.03498891
#> TTATCT  0.03136203  0.04927634 -0.1724813  0.03583921  0.02917477
#> CTTATC  0.02861515  0.01968119 -0.1637252  0.05942981  0.03772579
#> GGGGAG  0.04611842  0.02960637  0.1516554  0.04400812  0.03089627
#> GGGCGG  0.06576587  0.02662213  0.1510002  0.10446008  0.12732546
#> TAACTG -0.01146036 -0.02249384 -0.1488417 -0.03713173 -0.01858625

kmers <- screen_kmers(sequences_example, response_mat_example, max_gap = 3)
#> ℹ Number of response variables: 5
#> ℹ Screening kmers of length 6, from position 0 to position 300
#> ℹ Gaps of length 0:3 are allowed
#> ℹ minimal correlation: 0.08
#> ✔ Found 4884 kmers in 1000 sequences.
head(kmers)
#>       kmer     max_r2 avg_n    avg_var           c1          c2         c3
#> 1   AGATAA 0.05012024 0.139 0.15167901 -0.011159840 -0.03356719 -0.2238755
#> 2  AGA*AAG 0.03404636 0.220 0.23760001  0.062885828  0.06384566 -0.1845165
#> 3   TTATCT 0.02914828 0.127 0.13487099  0.036022298  0.05199933 -0.1707287
#> 4 ACAT**CT 0.02847921 0.090 0.09190000 -0.002383035 -0.02347000 -0.1687578
#> 5  AG*TAAG 0.02828897 0.118 0.12007600  0.034631759 -0.00579447 -0.1681933
#> 6   CTTATC 0.02680594 0.080 0.07960001  0.028615145  0.01968119 -0.1637252
#>           c4          c5
#> 1 0.00153976 -0.03616545
#> 2 0.02505355  0.02373969
#> 3 0.03218812  0.02568405
#> 4 0.01971230 -0.02059132
#> 5 0.04614470  0.02593124
#> 6 0.05942981  0.03772579
```
