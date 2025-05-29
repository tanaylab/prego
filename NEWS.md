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
