#' @title MotifDB Class
#' @slot mat A numeric matrix containing position weight matrices in log scale
#' @slot rc_mat A numeric matrix containing reverse complement of position weight matrices
#' @slot motif_lengths A named integer vector containing the length of each motif
#' @slot prior The pssm prior probability
#'
#' @description S4 class to store position weight matrices and their properties
setClass("MotifDB",
    slots = list(
        mat = "matrix",
        rc_mat = "matrix",
        motif_lengths = "integer",
        prior = "numeric"
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

        if (length(errors) == 0) TRUE else errors
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

#' Create a MotifDB object from a tidy data frame
#'
#' @param motif_db A tidy data frame containing motif information
#' @param prior Pseudocount prior to add to probabilities (default: 0.01)
#' @return A MotifDB object
#' @export
create_motif_db <- function(motif_db, prior = 0.01) {
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

    # Zero out positions after motif length for both matrices
    for (i in 1:ncol(matrices$mat)) {
        if (motif_lengths[i] * 4 < nrow(matrices$mat)) {
            zero_rows <- (motif_lengths[i] * 4 + 1):nrow(matrices$mat)
            matrices$mat[zero_rows, i] <- 0
            matrices$rc_mat[zero_rows, i] <- 0
        }
    }

    # Create and return MotifDB object
    new("MotifDB",
        mat = matrices$mat,
        rc_mat = matrices$rc_mat,
        motif_lengths = motif_lengths,
        prior = prior
    )
}

#' Show method for MotifDB objects
#'
#' @param object MotifDB object
setMethod(
    "show", "MotifDB",
    function(object) {
        cli::cli_text("{.cl MotifDB} object with {.val {ncol(object@mat)}} motifs and prior {.val {object@prior}}")
        cli::cli_text("Slots include: {.field @mat}, {.field @rc_mat} {.field @motif_lengths}, {.field @prior}")
    }
)

#' Initialize method for MotifDB objects
#'
#' @param .Object MotifDB object
#' @param mat Position weight matrix
#' @param rc_mat Reverse complement matrix
#' @param motif_lengths Named vector of motif lengths
#' @param prior Pseudocount prior value
setMethod(
    "initialize", "MotifDB",
    function(.Object,
             mat = matrix(0, 4, 1),
             rc_mat = matrix(0, 4, 1),
             motif_lengths = setNames(as.integer(1), colnames(mat)[1]),
             prior = 0.01) {
        .Object@mat <- mat
        .Object@rc_mat <- rc_mat
        .Object@motif_lengths <- motif_lengths
        .Object@prior <- prior
        validObject(.Object)
        return(.Object)
    }
)

#' Get specific motifs from the MotifDB
#'
#' @param object MotifDB object
#' @param motif_names Character vector of motif names or numeric indices
#' @return MotifDB object containing the specified motifs
setMethod(
    "[", "MotifDB",
    function(x, i) {
        if (is.character(i)) {
            missing_motifs <- i[!i %in% colnames(x@mat)]
            if (length(missing_motifs) > 0) {
                cli::cli_abort(c(
                    "Some motifs not found in MotifDB:",
                    "x" = "{.val {missing_motifs}}"
                ))
            }
            idx <- match(i, colnames(x@mat))
        } else {
            idx <- i
        }

        new("MotifDB",
            mat = x@mat[, idx, drop = FALSE],
            rc_mat = x@rc_mat[, idx, drop = FALSE],
            motif_lengths = x@motif_lengths[idx],
            prior = x@prior
        )
    }
)
