# Convert intervals to sequences

This function takes a set of intervals and converts them into sequences.
It requires the 'misha' package to be installed. If the package is not
installed, it will display an error message with instructions on how to
install it.

## Usage

``` r
intervals_to_seq(intervals, size = NULL)
```

## Arguments

- intervals:

  The intervals set as a data frame with 'chrom', 'start', and 'end'
  columns. Can be a string with misha intervals name.

- size:

  The size to normalize the intervals to. If NULL, the intervals will
  not be normalized.

## Value

A character vector of sequences.

## Examples

``` r
if (FALSE) { # \dontrun{
library(misha)
gdb.init_examples()
intervals_to_seq("annotations")
intervals_to_seq("annotations", 20)
} # }
```
