# MotifDB Class

S4 class to store position weight matrices and their properties

## Slots

- `mat`:

  A numeric matrix containing position weight matrices in log scale

- `rc_mat`:

  A numeric matrix containing reverse complement of position weight
  matrices

- `motif_lengths`:

  A named integer vector containing the length of each motif

- `prior`:

  The pssm prior probability

- `spat_factors`:

  A numeric matrix containing spatial factors

- `spat_bin_size`:

  The size of spatial bins

- `spat_min`:

  The starting position of the sequence or NA

- `spat_max`:

  The ending position of the sequence or NA
