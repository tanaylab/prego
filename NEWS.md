# prego 0.0.10

* Fix: `calc_seq_pwm()` / `extract_pwm()` now error clearly when given sequences
  of unequal length instead of silently recycling the shorter ones (via `rbind`)
  and returning wrong scores. These functions build a single rectangular one-hot
  matrix and require equal-length sequences; for variable lengths use
  `compute_pwm()`.
* Fix: PWM scoring no longer opens thousands of threads / fails on core-limited
  machines. `compute_pwm()`, `compute_local_pwm()` and `calc_seq_pwm()` (hence
  `extract_pwm()`) run a small per-sequence BLAS `dgemm` inside an
  RcppParallel/TBB `parallelFor`. With a multi-threaded BLAS (MKL, OpenBLAS) each
  TBB worker spawned its own BLAS thread team, so one call opened
  `n_threads * n_blas_threads` OS threads (thousands on a many-core node) - which
  oversubscribed the machine and could fail outright on a cluster job with a
  thread/process (cgroup pids) limit. The inner BLAS is now pinned to a single
  thread for the duration of these calls (the TBB loop over sequences is the
  parallelism layer), and the previous BLAS thread count is restored afterwards.
  Results are unchanged; throughput is the same or better. Adds a dependency on
  the `RhpcBLASctl` package.
* Fix: `compute_pwm()` gave a sequence a different score depending on the other
  sequences in the batch. The motif-scan window was capped to the length of the
  *first* sequence for the whole batch, so a longer sequence sitting behind a
  shorter one was only scanned over its first `nchar(sequences[1])` positions and
  missed motif hits further along. Each sequence is now scored on its own full
  length, independently of its batch companions.
* Added `return_all` parameter to `regress_pwm` (multi-kmer path). When TRUE, returns every candidate-kmer regression (sorted by validation score) instead of just the best one - useful for getting N independent motifs without the residual-rounds approach used by `motif_num > 1`. When `sample_for_kmers = TRUE`, each candidate is refit on the full data.
* Improved docs for `regress_pwm` (clarified the three operating modes, fixed `n_motifs`/`comb_modle` typos in the return-value section).

# prego 0.0.9

* Fix: Floating point discrepancies between `predict` and `regress_pwm` output.

# prego 0.0.8 

* Fix: `compute_local_pwm` now returns NA when sequence length < motif length (issue #35).

# prego 0.0.7

* Fix: `regress_multiple_motifs` now works with multiple response variables.

# prego 0.0.6

* Added `screen_local_pwm` to find positions in sequences that match a PSSM.
* Added `return_list` parameter to `compute_local_pwm`.

# prego 0.0.5 

* Faster pssm correlation computation using `RcppParallel`.
* `pssm_match` when `best=FALSE` now returns a `score` field instead of `cor`.
* Implemented computation of KL divergence between two PSSMs. Note that spearman correlation is still the best way to match PSSMs.

# prego 0.0.4

* Added `size` paramter to `intervals_to_seq`.

# prego 0.0.3

* Added MotifDB object to store motif information.
* Implmented a faster energy computation method, which is now used by default.

# prego 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
