calculate_bins <- function(max_seq_len, spat_num_bins = NULL, spat_bin_size = NULL, default_bin_size = 40) {
    if (!is.null(spat_num_bins) && !is.null(spat_bin_size)) {
        return(list(spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size))
    }

    if (!is.null(spat_num_bins)) {
        spat_bin_size <- floor(max_seq_len / spat_num_bins)
        return(list(spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size))
    }

    if (!is.null(spat_bin_size)) {
        spat_num_bins <- floor(max_seq_len / spat_bin_size)
        if (spat_num_bins %% 2 == 0) {
            spat_num_bins <- spat_num_bins - 1
        }
        return(list(spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size))
    }

    spat_bin_size <- default_bin_size
    spat_num_bins <- floor(max_seq_len / spat_bin_size)

    if (spat_num_bins %% 2 == 0) {
        spat_num_bins <- spat_num_bins - 1
    }

    if (spat_num_bins < 3) {
        cli_abort("Calculated spat_num_bins is less than 3, which is not permissible.")
    }

    return(list(spat_num_bins = spat_num_bins, spat_bin_size = spat_bin_size))
}

calc_spat_min_max <- function(spat_num_bins, max_seq_len, spat_bin_size) {
    if (spat_num_bins %% 2 != 1) {
        cli_abort("The {.field spat_num_bins} must be an odd number")
    }
    if (spat_bin_size * spat_num_bins > max_seq_len) {
        cli_abort("The {.field spat_bin_size} ({.val {spat_bin_size}}) times the {.field spat_num_bins} ({.val {spat_num_bins}}) must be smaller than the maximum sequence length ({.val {max_seq_len}})")
    }

    center <- round(max_seq_len / 2)

    if (spat_num_bins == 1) {
        spat_min <- center - spat_bin_size / 2
        spat_max <- center + spat_bin_size / 2
    } else {
        # position one bin at the center, and then add bins to the left and to the right
        spat_min <- center - ((spat_num_bins - 1) / 2) * spat_bin_size - spat_bin_size / 2
        spat_max <- center + ((spat_num_bins - 1) / 2) * spat_bin_size + spat_bin_size / 2
    }


    return(list(spat_min = round(spat_min), spat_max = round(spat_max)))
}
