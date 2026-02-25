# Set parallel threads

Set the number of parallel threads to use. prego uses
[`doMC::registerDoMC`](https://rdrr.io/pkg/doMC/man/registerDoMC.html)
(when available) to register the parallelization backend for plyr. By
default, prego uses 80% of the number of available cores. The options
are saved under 'prego.parallel' (should we use parallelization,
logical) and 'prego.parallel.nc' (number of cores to use, integer).

## Usage

``` r
set_parallel(thread_num = max(1, round(parallel::detectCores() * 0.8)))
```

## Arguments

- thread_num:

  number of threads. use '1' for non parallel behavior

## Value

None

## Examples

``` r
# \donttest{
set_parallel(8)
# }
```
