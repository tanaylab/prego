## R CMD check results

`R CMD check --as-cran --no-manual prego_0.0.9.tar.gz`

0 errors | 0 warnings | 3 notes

### Notes

1. `checking CRAN incoming feasibility ... NOTE`
   - New submission (expected for first CRAN submission).
2. `checking installed package size ... NOTE`
   - Installed size is 19.4Mb (`data` 3.9Mb, `libs` 13.2Mb).
   - The package includes bundled motif datasets and compiled C++ code used by core functionality.
3. `checking for future file timestamps ... NOTE`
   - Environment-specific on our shared filesystem (`unable to verify current time`).

### Additional checks

1. All examples complete under `--as-cran` and `--run-donttest`.
2. Test suite passes during `R CMD check`.
3. CRAN incoming policy-related documentation notes were addressed:
   - Updated moved source URLs.
   - Replaced DOI URLs with `\\doi{...}` where required.
   - Updated DESCRIPTION text opening sentence.
