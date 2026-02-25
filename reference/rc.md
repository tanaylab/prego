# Reverse Complement DNA Sequences

This function takes a character vector of DNA sequences and returns
their reverse complements. It uses an efficient C++ implementation via
Rcpp for improved performance.

## Usage

``` r
rc(dna)
```

## Arguments

- dna:

  A character vector of DNA sequences. Can be a single sequence or
  multiple sequences. The sequences can be in upper or lower case.

## Value

A character vector of the same length as the input, where each element
is the reverse complement of the corresponding input sequence.

## Details

The function performs the following operations on each sequence: 1.
Converts the sequence to uppercase. 2. Reverses the sequence. 3.
Complements each base (A\<-\>T, C\<-\>G). Non-standard characters (not
A, T, C, or G) are preserved in their reversed positions.

## Examples

``` r
rc("ATCG") # Returns "CGAT"
#> [1] "CGAT"
rc(c("ATCG", "GGCC", "TATA")) # Returns c("CGAT", "GGCC", "TATA")
#> [1] "CGAT" "GGCC" "TATA"
```
