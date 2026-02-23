## R CMD check results

`R CMD check --as-cran --no-manual prego_0.0.9.tar.gz`

0 errors | 0 warnings | 2 notes

### Notes

1. `checking CRAN incoming feasibility ... NOTE`
   - New submission (expected for first CRAN submission).
2. `checking installed package size ... NOTE`
   - Installed size is 19.4Mb (`data` 3.9Mb, `libs` 13.2Mb).
   - The package includes bundled motif datasets and compiled C++ code used by core functionality.

### Additional checks

1. All examples complete under `--as-cran` and `--run-donttest`.
2. Test suite passes during `R CMD check`.
