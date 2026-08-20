# prego 0.0.11

* Fix: the thread-explosion guard added in 0.0.10 did not actually work on an
  OpenMP-linked OpenBLAS - the `libopenblas-*-openmp` build conda ships, and
  therefore what most lab environments run. `RhpcBLASctl::blas_set_num_threads(1)`
  reports success there while changing nothing: `openblas_set_num_threads()` is a
  no-op for that build and `omp_get_max_threads()` is what sizes the thread team.
  `local_serial_blas()` now pins the OpenMP thread count as well. Measured on a
  128-core node, `extract_pwm()` over 3,000 x 500bp sequences and 20 motifs at
  `set_parallel(16)`: unguarded and blas-pin-only both exceeded 600s at 2,049 OS
  threads and ~11,000% CPU (~80% of it system time); with the OpenMP pin, 4.8s.

  This does not close the gap completely. `omp_set_num_threads()` writes the
  *calling* thread's libgomp ICV, while the parallel workers that actually run
  the BLAS call are separate threads which still inherit the untouched global
  one. The same benchmark runs in 0.47s with 17 threads when `OMP_NUM_THREADS=1`
  is set in the environment *before* R starts, so on a many-core machine that is
  still worth setting: `OMP_NUM_THREADS=1 Rscript your_script.R`. Note that
  `Sys.setenv(OMP_NUM_THREADS = 1)` inside the script cannot work - libgomp reads
  the variable in an ELF constructor at process load, long before R runs line 1.
  Pinning the count from inside the worker threads themselves would remove the
  need for the variable; that is worth doing but needs validating on a machine
  where the affected BLAS can be rebuilt against.

# prego 0.0.10

* Fix: `compute_pwm()` scored `N` (and the `*` wildcard) inconsistently between
  strands - the forward strand used the column's average log-probability while
  the reverse strand used a flat `log(0.25)`. With `bidirect = TRUE` this made a
  sequence and its reverse-complement score differently whenever an `N` fell on
  an informative position. Both strands now use the column average
  (`get_avg_log_prob()`), matching the other likelihood routines, so scoring is
  strand-symmetric again. Only affects sequences containing `N`/`*`.
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
