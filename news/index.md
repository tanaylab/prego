# Changelog

## prego 0.0.10

- Added `return_all` parameter to `regress_pwm` (multi-kmer path). When
  TRUE, returns every candidate-kmer regression (sorted by validation
  score) instead of just the best one - useful for getting N independent
  motifs without the residual-rounds approach used by `motif_num > 1`.
  When `sample_for_kmers = TRUE`, each candidate is refit on the full
  data.
- Improved docs for `regress_pwm` (clarified the three operating modes,
  fixed `n_motifs`/`comb_modle` typos in the return-value section).

## prego 0.0.9

- Fix: Floating point discrepancies between `predict` and `regress_pwm`
  output.

## prego 0.0.8

- Fix: `compute_local_pwm` now returns NA when sequence length \< motif
  length (issue \#35).

## prego 0.0.7

- Fix: `regress_multiple_motifs` now works with multiple response
  variables.

## prego 0.0.6

- Added `screen_local_pwm` to find positions in sequences that match a
  PSSM.
- Added `return_list` parameter to `compute_local_pwm`.

## prego 0.0.5

- Faster pssm correlation computation using `RcppParallel`.
- `pssm_match` when `best=FALSE` now returns a `score` field instead of
  `cor`.
- Implemented computation of KL divergence between two PSSMs. Note that
  spearman correlation is still the best way to match PSSMs.

## prego 0.0.4

- Added `size` paramter to `intervals_to_seq`.

## prego 0.0.3

- Added MotifDB object to store motif information.
- Implmented a faster energy computation method, which is now used by
  default.

## prego 0.0.0.9000

- Added a `NEWS.md` file to track changes to the package.
