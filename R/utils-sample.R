#' Sample rows respecting quantiles of a reference distribution
#'
#' This function randomly samples rows from a data frame in such a way that the
#' quantiles of the selected data match as closely as possible those of the full data.
#'
#' @param data_frame A data frame from which to sample rows.
#' @param reference A numeric vector of the same length as the number of rows in the data frame.
#' @param sample_fraction A fraction specifying the proportion of rows to sample from the data frame.
#' @param num_quantiles An integer specifying the number of quantiles, default is 10.
#' @param seed An integer specifying the random seed to use.
#'
#' @return A data frame of sampled rows.
#'
#' @examples
#' sampled <- sample_quantile_matched_rows(mtcars, mtcars$mpg, sample_fraction = 0.1)
#' plot(quantile(mtcars$mpg), quantile(sampled$mpg))
#' abline(0, 1)
#'
#' @export
sample_quantile_matched_rows <- function(data_frame, reference, sample_fraction, num_quantiles = 10, seed = 60427) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    columns <- colnames(data_frame)

    sample_size <- floor(nrow(data_frame) * sample_fraction)

    if (sample_size < 1) {
        cli::cli_abort("Sample fraction {.val {sample_fraction}} is too small for the given data frame size {.val {nrow(data_frame)}}.")
    }

    # Adjust num_quantiles if it's larger than sample_size
    num_quantiles <- min(sample_size, num_quantiles)

    # Calculate quantiles of the reference
    ref_quantiles <- quantile(reference, probs = seq(0, 1, length.out = num_quantiles + 1))

    # Create a data frame with reference and its corresponding quantile
    df_with_quantiles <- data_frame %>%
        mutate(quantile = cut(reference, breaks = unique(ref_quantiles), labels = FALSE, include.lowest = TRUE)) %>%
        group_by(quantile)

    # Sample rows such that the quantiles of the sampled data matches as closely as possible those of the full data
    quantile_matched_sampled_df <- df_with_quantiles %>%
        dplyr::sample_n(size = floor(sample_size / num_quantiles), replace = FALSE) %>%
        ungroup()

    quantile_matched_sampled_df <- quantile_matched_sampled_df %>%
        select(any_of(columns))

    cli::cli_alert_success("Sampled {.val {nrow(quantile_matched_sampled_df)}} rows from the data frame.")
    return(quantile_matched_sampled_df)
}

sample_response <- function(response, sample_frac = NULL, sample_ratio = 1, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (is.null(sample_frac)) {
        if (is_binary_response(response)) {
            sample_frac <- c(pmin(1, sample_ratio * sum(response[, 1] == 1) / sum(response[, 1] == 0)), 1)
        } else {
            cli_alert_info("Using {.code sample_frac = 0.1}")
            sample_frac <- 0.1
        }
    }
    cli_alert_info("Sampling {.val {round(sample_frac, digits = 2)}} of the dataset")
    categorical <- ncol(response) == 1 && all(response[, 1] == 0 | response[, 1] == 1)
    if (categorical) {
        cli_alert_info("Stratified sampling")
        if (length(sample_frac) == 1) {
            sample_frac <- c(sample_frac, sample_frac)
        }
        samp_idx_0 <- sample(which(response[, 1] == 0), size = round(sample_frac[1] * sum(response[, 1] == 0)))
        samp_idx_1 <- sample(which(response[, 1] == 1), size = round(sample_frac[2] * sum(response[, 1] == 1)))
        sample_idxs <- c(samp_idx_0, samp_idx_1)
    } else {
        sample_idxs <- sample_quantile_matched_rows(data.frame(i = 1:nrow(response)), response[, 1], sample_frac, seed = seed) %>%
            pull(i) %>%
            sort()        
    }

    if (is_binary_response(response)) {
        cli_alert_info("sampled {.val {sum(response[sample_idxs, 1] == 0)}} zeros and {.val {sum(response[sample_idxs, 1] == 1)}} ones")
    }
    return(sample_idxs)
}
