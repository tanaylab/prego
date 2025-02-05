#' @title MotifDB Class
#' @slot mat A numeric matrix containing position weight matrices in log scale
#' @slot rc_mat A numeric matrix containing reverse complement of position weight matrices
#' @slot motif_lengths A named integer vector containing the length of each motif
#' @slot prior The pssm prior probability
#' @slot spat_factors A numeric matrix containing spatial factors
#' @slot spat_bin_size The size of spatial bins
#' @slot spat_min The starting position of the sequence or NA
#' @slot spat_max The ending position of the sequence or NA
#'
#' @description S4 class to store position weight matrices and their properties
setClass("MotifDB",
    slots = list(
        mat = "matrix",
        rc_mat = "matrix",
        motif_lengths = "integer",
        prior = "numeric",
        spat_factors = "matrix", # Matrix of spatial factors (motifs x bins)
        spat_bin_size = "numeric", # Size of each spatial bin
        spat_min = "numeric", # Starting position of sequence
        spat_max = "numeric" # Ending position of sequence
    ),
    validity = function(object) {
        errors <- character()

        # Check if mat exists and has the right structure
        if (nrow(object@mat) %% 4 != 0) {
            errors <- c(errors, "Matrix rows must be multiple of 4 (A,C,G,T)")
        }

        # Check if rc_mat has same dimensions as mat
        if (!identical(dim(object@mat), dim(object@rc_mat))) {
            errors <- c(errors, "Reverse complement matrix must have same dimensions as main matrix")
        }

        # Check if motif_lengths matches mat columns
        if (length(object@motif_lengths) != ncol(object@mat)) {
            errors <- c(errors, "Length of motif_lengths must match number of matrix columns")
        }

        # Check if motif_lengths are positive
        if (any(object@motif_lengths <= 0)) {
            errors <- c(errors, "All motif lengths must be positive")
        }

        # Check if motif_lengths has names matching mat columns
        if (!identical(names(object@motif_lengths), colnames(object@mat))) {
            errors <- c(errors, "Names of motif_lengths must match column names of matrix")
        }

        # Check if prior is valid
        if (object@prior <= 0 || object@prior >= 1) {
            errors <- c(errors, "Prior must be between 0 and 1")
        }

        # Check that spatial factors matrix has correct dimensions and values
        if (nrow(object@spat_factors) > 0) { # Only check if not empty
            # Check number of rows matches number of motifs
            if (nrow(object@spat_factors) != ncol(object@mat)) {
                errors <- c(errors, "Number of rows in spatial factors matrix must match number of motifs")
            }

            # Check for non-negative values
            if (any(object@spat_factors < 0, na.rm = TRUE)) {
                errors <- c(errors, "Spatial factors must be non-negative")
            }

            # Check row names match motif names
            if (!identical(rownames(object@spat_factors), colnames(object@mat))) {
                errors <- c(errors, "Row names of spatial factors matrix must match motif names")
            }
        }

        # Check that spat_bin_size is positive
        if (object@spat_bin_size <= 0) {
            errors <- c(errors, "Spatial bin size must be positive")
        }

        # Check spatial min/max values if they exist
        if (!is.na(object@spat_min) && !is.na(object@spat_max)) {
            if (object@spat_min > object@spat_max) {
                errors <- c(errors, "Spatial min must be less than or equal to spatial max")
            }
            if (object@spat_min < 0) {
                errors <- c(errors, "Spatial min must be non-negative")
            }
        }

        if (length(errors) == 0) TRUE else errors
    }
)

#' Create a MotifDB object from a tidy data frame
#'
#' @param motif_db A tidy data frame containing motif information
#' @param prior Pseudocount prior to add to probabilities (default: 0.01)
#' @param spat_factors Matrix of spatial factors (rows=motifs, cols=bins) or NULL
#' @param spat_bin_size Size of spatial bins (default: 1)
#' @param spat_min Starting position of sequence or NA (default: NA)
#' @param spat_max Ending position of sequence or NA (default: NA)
#' @return A MotifDB object
#' @examples
#' create_motif_db(all_motif_datasets())
#' @export
create_motif_db <- function(motif_db, prior = 0.01, spat_factors = NULL,
                            spat_bin_size = 1, spat_min = NA_real_, spat_max = NA_real_) {
    # Calculate matrices using modified function
    matrices <- motif_db_to_mat(motif_db, prior)

    # Calculate motif lengths
    motif_lengths <- motif_db %>%
        dplyr::count(motif) %>%
        select(motif, n) %>%
        tibble::deframe()

    # Transform matrices to log scale
    matrices$mat <- log(matrices$mat)
    matrices$rc_mat <- log(matrices$rc_mat)

    motif_names <- unique(motif_db$motif)
    matrices$mat <- matrices$mat[, motif_names, drop = FALSE]
    matrices$rc_mat <- matrices$rc_mat[, motif_names, drop = FALSE]
    motif_lengths <- motif_lengths[motif_names]

    # Zero out positions after motif length for both matrices
    for (i in 1:ncol(matrices$mat)) {
        if (motif_lengths[i] * 4 < nrow(matrices$mat)) {
            zero_rows <- (motif_lengths[i] * 4 + 1):nrow(matrices$mat)
            matrices$mat[zero_rows, i] <- 0
            matrices$rc_mat[zero_rows, i] <- 0
        }
    }

    # Create default spatial factors if none provided
    if (is.null(spat_factors)) {
        spat_factors <- matrix(1,
            nrow = length(motif_lengths), ncol = 1,
            dimnames = list(names(motif_lengths), NULL)
        )
    } else {
        # Ensure row names match motif names
        if (any(!rownames(spat_factors) %in% names(motif_lengths))) {
            cli::cli_abort("Some motifs not found in spatial factors matrix")
        }
        spat_factors <- spat_factors[names(motif_lengths), , drop = FALSE]
    }

    # Create and return MotifDB object
    new("MotifDB",
        mat = matrices$mat,
        rc_mat = matrices$rc_mat,
        motif_lengths = motif_lengths,
        prior = prior,
        spat_factors = spat_factors,
        spat_bin_size = as.numeric(spat_bin_size),
        spat_min = spat_min,
        spat_max = spat_max
    )
}

#' Show method for MotifDB objects
#'
#' @param object MotifDB object
setMethod(
    "show", "MotifDB",
    function(object) {
        cli::cli_text("{.cl MotifDB} object with {.val {ncol(object@mat)}} motifs and prior {.val {object@prior}}")
        if (ncol(object@spat_factors) > 1 || object@spat_bin_size > 1) {
            cli::cli_text("Spatial factors: bin size {.val {object@spat_bin_size}} with {.val {ncol(object@spat_factors)}} bins per motif")
        }
        if (!is.na(object@spat_min) && !is.na(object@spat_max)) {
            cli::cli_text("Spatial range: {.val {object@spat_min}} to {.val {object@spat_max}}")
        }
        cli::cli_text("Slots include: {.field @mat}, {.field @rc_mat}, {.field @motif_lengths}, {.field @prior}, {.field @spat_factors}, {.field @spat_bin_size}, {.field @spat_min}, {.field @spat_max}")
    }
)

#' Initialize method for MotifDB objects
#'
#' @param .Object MotifDB object
#' @param mat Position weight matrix
#' @param rc_mat Reverse complement matrix
#' @param motif_lengths Named vector of motif lengths
#' @param prior Pseudocount prior value
#' @param spat_factors Matrix of spatial factors
#' @param spat_bin_size Size of spatial bins
#' @param spat_min Starting position of sequence or NULL
#' @param spat_max Ending position of sequence or NULL
setMethod(
    "initialize", "MotifDB",
    function(.Object,
             mat = matrix(0, 4, 1),
             rc_mat = matrix(0, 4, 1),
             motif_lengths = setNames(as.integer(1), colnames(mat)[1]),
             prior = 0.01,
             spat_factors = matrix(1, nrow = 1, ncol = 1),
             spat_bin_size = 1,
             spat_min = NA_real_,
             spat_max = NA_real_) {
        .Object@mat <- mat
        .Object@rc_mat <- rc_mat
        .Object@motif_lengths <- motif_lengths
        .Object@prior <- prior
        .Object@spat_factors <- spat_factors
        .Object@spat_bin_size <- as.numeric(spat_bin_size)
        .Object@spat_min <- as.numeric(spat_min)
        .Object@spat_max <- as.numeric(spat_max)
        validObject(.Object)
        return(.Object)
    }
)

#' Get specific motifs from the MotifDB
#'
#' @param x MotifDB object
#' @param i Character vector of motif names, numeric indices, or regex pattern(s)
#' @param j Not used
#' @param drop Not used
#' @param ... Not used
#' @param pattern Logical indicating whether to treat character input as regex pattern (default: FALSE)
#' @return MotifDB object containing the specified motifs
#' @examples
#' MOTIF_DB["HOMER.GATA3_2"]
#' MOTIF_DB[c("HOMER.GATA3_2", "JASPAR.CDX1")]
#' MOTIF_DB["GATA", pattern = T]
#' @export
setMethod(
    "[", "MotifDB",
    function(x, i, j, ..., pattern = TRUE, drop = TRUE) {
        if (is.character(i)) {
            if (pattern) {
                # pattern matching
                matching_motifs <- unique(unlist(lapply(i, function(pat) {
                    matches <- grep(pat, colnames(x@mat), ignore.case = TRUE, value = TRUE)
                    if (length(matches) == 0) {
                        cli::cli_warn("Pattern {.val {pat}} matched no motifs")
                    }
                    return(matches)
                })))

                if (length(matching_motifs) == 0) {
                    cli::cli_abort("No motifs matched any of the provided patterns", call = parent.frame())
                }

                idx <- match(matching_motifs, colnames(x@mat))
            } else {
                # exact matching
                missing_motifs <- i[!i %in% colnames(x@mat)]
                if (length(missing_motifs) > 0) {
                    cli::cli_abort(c(
                        "Some motifs not found in MotifDB:",
                        "x" = "{.val {missing_motifs}}"
                    ), call = parent.frame())
                }
                idx <- match(i, colnames(x@mat))
            }
        } else {
            # numeric indices
            if (any(i > ncol(x@mat) | i < 1)) {
                cli::cli_abort("Index out of bounds")
            }
            idx <- i
        }

        new("MotifDB",
            mat = x@mat[, idx, drop = FALSE],
            rc_mat = x@rc_mat[, idx, drop = FALSE],
            motif_lengths = x@motif_lengths[idx],
            prior = x@prior,
            spat_factors = x@spat_factors[idx, , drop = FALSE],
            spat_bin_size = x@spat_bin_size,
            spat_min = x@spat_min,
            spat_max = x@spat_max
        )
    }
)


motif_db_to_mat <- function(motif_db, prior = 0.01) {
    # Add position column
    motif_db <- motif_db %>%
        group_by(motif) %>%
        mutate(pos = dplyr::row_number())

    D <- max(motif_db$pos)

    # Process for forward matrix
    forward_mat <- motif_db %>%
        select(motif, pos, A, C, G, T) %>%
        tidyr::pivot_longer(
            cols = c(A, C, G, T),
            names_to = "nucleotide",
            values_to = "value"
        ) %>%
        group_by(motif, pos) %>%
        mutate(value = value / sum(value)) %>%
        mutate(value = value + prior) %>%
        mutate(value = value / sum(value))

    # Create row IDs and arrange for forward matrix
    forward_df <- forward_mat %>%
        mutate(
            row_id = paste0(nucleotide, "_", pos)
        ) %>%
        arrange(motif, pos, match(nucleotide, c("A", "C", "G", "T"))) %>%
        tidyr::pivot_wider(
            id_cols = row_id,
            names_from = motif,
            values_from = value,
            values_fill = 0
        ) %>%
        tibble::column_to_rownames("row_id")

    # Create reverse complement matrix
    rc_df <- forward_mat %>%
        group_by(motif) %>%
        mutate(
            # Reverse the position order
            rc_pos = max(pos) - pos + 1,
            # Complement the nucleotides
            rc_nucleotide = dplyr::case_when(
                nucleotide == "A" ~ "T",
                nucleotide == "T" ~ "A",
                nucleotide == "C" ~ "G",
                nucleotide == "G" ~ "C"
            )
        ) %>%
        mutate(
            row_id = paste0(rc_nucleotide, "_", rc_pos)
        ) %>%
        arrange(motif, rc_pos, match(rc_nucleotide, c("A", "C", "G", "T"))) %>%
        tidyr::pivot_wider(
            id_cols = row_id,
            names_from = motif,
            values_from = value,
            values_fill = 0
        ) %>%
        tibble::column_to_rownames("row_id")

    # Convert to matrices and return as list
    list(
        mat = as.matrix(forward_df),
        rc_mat = as.matrix(rc_df)
    )
}

#' Get or set the prior for a MotifDB object
#'
#' @param object MotifDB object
#' @return The prior value
setGeneric("prior", function(object) standardGeneric("prior"))

#' Set the prior for a MotifDB object
#'
#' @param object MotifDB object
#' @param value New prior value
#' @return Updated MotifDB object
setGeneric("prior<-", function(object, value) standardGeneric("prior<-"))

#' Get the prior value from a MotifDB object
#'
#' @param object MotifDB object#'
#' @return The prior value
#' @examples
#' prior(MOTIF_DB)
#' @export
setMethod(
    "prior", "MotifDB",
    function(object) {
        return(object@prior)
    }
)

#' Convert a MotifDB object back to a tidy data frame
#'
#' @param motif_db A MotifDB object
#' @return A tidy data frame with columns for motif, position, and nucleotide probabilities
#' @examples
#' head(motif_db_to_dataframe(MOTIF_DB))
#' @export
motif_db_to_dataframe <- function(motif_db) {
    # Convert log values back to probabilities
    prob_mat <- exp(motif_db@mat)

    result <- prob_mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column("p") %>%
        separate(p, c("nuc", "pos"), sep = "_") %>%
        gather("motif", "prob", -nuc, -pos) %>%
        group_by(motif, pos) %>%
        mutate(prob = sum(prob + motif_db@prior) * prob - motif_db@prior) %>%
        ungroup() %>%
        spread("nuc", "prob") %>%
        select(motif, pos, A, C, G, T) %>%
        mutate(pos = as.numeric(pos)) %>%
        arrange(motif, pos)

    motif_len <- tibble::enframe(motif_db@motif_lengths, name = "motif", value = "max_len")

    result <- result %>%
        left_join(motif_len, by = "motif") %>%
        filter(pos <= max_len) %>%
        select(-max_len)

    return(result)
}

#' Convert a MotifDB object to a data frame
#'
#' @param x A MotifDB object
#' @param row.names NULL or a character vector giving the row names for the data frame
#' @param optional logical. If TRUE, setting row names and converting column names (to syntactic names: see make.names) is optional
#' @param ... additional arguments to be passed to or from methods
#' @return A data frame containing the motif probabilities
#' @examples
#' dataset <- as.data.frame(MOTIF_DB)
#' head(dataset)
#' nrow(dataset)
#' length(unique(dataset$motif))
#' @export
setMethod(
    "as.data.frame", "MotifDB",
    function(x, row.names = NULL, optional = FALSE, ...) {
        df <- motif_db_to_dataframe(x)
        return(as.data.frame(df, row.names = row.names, optional = optional, ...))
    }
)

#' Convert a MotifDB object to a matrix
#' @param x MotifDB object
#' @param ... ignored arguments
#' @return A matrix containing the motif probabilities, rownames are motif_position, colnames are nucleotides
#' @examples
#' as.matrix(MOTIF_DB["HOMER.GATA3_2"])
#' @export
setMethod(
    "as.matrix", "MotifDB",
    function(x, ...) {
        return(as.data.frame(x) %>% tidyr::unite("pos", motif, pos) %>% tibble::column_to_rownames("pos") %>% as.matrix())
    }
)

#' Plot a motif from a MotifDB object
#' @param x MotifDB object
#' @param title title of the plot
#' @param subtitle subtitle of the plot
#' @param revcomp whether to plot the reverse complement of the PSSM
#' @param method Height method, can be one of "bits" or "probability" (default:"bits")
#' @param force force plotting more than 30 motifs
#' @inheritDotParams ggseqlogo::ggseqlogo
#' @return a ggplot object
#' @examples
#' plot(MOTIF_DB["HOMER.GATA3_2"])
#' plot(MOTIF_DB["HNF1", pattern = T])
#' plot(MOTIF_DB[c("HOMER.GATA3_2", "JASPAR.CDX1")])
#'
#' @export
setMethod(
    "plot", "MotifDB",
    function(x, title = "", subtitle = ggplot2::waiver(), revcomp = FALSE, method = "bits", force = FALSE, ...) {
        if (length(x) > 30) {
            if (!force) {
                cli::cli_abort("Plotting more than 30 motifs is not recommended. Please subset the MotifDB object first.", call = parent.frame())
            } else {
                cli::cli_warn("Plotting more than 30 motifs is not recommended.")
            }
        }
        x <- as.data.frame(x)
        x <- split(x, x$motif)
        pfm <- purrr::map(x, ~ {
            pssm <- .x %>% select(-motif)
            if (revcomp) {
                pssm <- pssm_rc(pssm)
            }
            t(pssm_to_mat(pssm))
        })
        ggseqlogo::ggseqlogo(pfm, method = method) +
            ggtitle(title, subtitle = subtitle)
    }
)

#' Set a new prior for a MotifDB object
#'
#' @param object MotifDB object
#' @param value New prior value between 0 and 1
#' @return Updated MotifDB object with new prior
#' @examples
#' prior(MOTIF_DB)
#' prior(MOTIF_DB) <- 0.2
#' prior(MOTIF_DB)
#' @export
setMethod(
    "prior<-", "MotifDB",
    function(object, value) {
        # Convert to data frame
        df <- motif_db_to_dataframe(object)

        # Create new MotifDB with new prior
        new_db <- create_motif_db(
            df,
            prior = value,
            spat_factors = object@spat_factors,
            spat_bin_size = object@spat_bin_size,
            spat_min = object@spat_min,
            spat_max = object@spat_max
        )

        # Return new object
        return(new_db)
    }
)

#' Get the length of a MotifDB object
#' @param x MotifDB object
#' @return The number of motifs in the object
#' @examples
#' length(MOTIF_DB)
#' length(MOTIF_DB[c("HOMER.GATA3_2", "JASPAR.CDX1")])
#' @export
setMethod(
    "length", "MotifDB",
    function(x) {
        return(ncol(x@mat))
    }
)

#' Get the names of motifs in a MotifDB object
#' @param x MotifDB object
#' @return The names of motifs in the object
#' @examples
#' names(MOTIF_DB)
#' names(MOTIF_DB[c("HOMER.GATA3_2", "JASPAR.CDX1")])
#' @export
setMethod(
    "names", "MotifDB",
    function(x) {
        return(colnames(x@mat))
    }
)
