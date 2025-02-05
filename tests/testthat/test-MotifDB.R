create_test_motif_db <- function(prior = 0.01, spat_min = NA_real_, spat_max = NA_real_) {
    # Create a simple motif database with two motifs
    motif_data <- tibble(
        motif = rep(c("motif1", "motif2"), each = 4),
        pos = rep(1, 8),
        A = c(0.7, 0.1, 0.2, 0.3, 0.1, 0.6, 0.2, 0.1),
        C = c(0.1, 0.7, 0.3, 0.2, 0.2, 0.1, 0.6, 0.2),
        G = c(0.1, 0.1, 0.3, 0.3, 0.6, 0.2, 0.1, 0.1),
        T = c(0.1, 0.1, 0.2, 0.2, 0.1, 0.1, 0.1, 0.6)
    )

    create_motif_db(motif_data, prior = prior, spat_min = spat_min, spat_max = spat_max)
}

# Test MotifDB class creation and validation
test_that("MotifDB object creation works with valid input", {
    motif_db <- create_test_motif_db()
    expect_s4_class(motif_db, "MotifDB")
    expect_true(validObject(motif_db))

    # Test with spatial boundaries
    motif_db_with_bounds <- create_test_motif_db(spat_min = 0, spat_max = 100)
    expect_s4_class(motif_db_with_bounds, "MotifDB")
    expect_true(validObject(motif_db_with_bounds))
})

test_that("MotifDB validates spatial boundaries correctly", {
    # Test valid spatial boundaries
    expect_error(create_test_motif_db(spat_min = 0, spat_max = 100), NA)
    expect_error(create_test_motif_db(spat_min = 0, spat_max = 0), NA)

    # Test invalid spatial boundaries
    expect_error(create_test_motif_db(spat_min = 100, spat_max = 0))
    expect_error(create_test_motif_db(spat_min = -1, spat_max = 100))

    # Test NA values (should be valid)
    expect_error(create_test_motif_db(spat_min = NA_real_, spat_max = NA_real_), NA)
    expect_error(create_test_motif_db(spat_min = 0, spat_max = NA_real_), NA)
    expect_error(create_test_motif_db(spat_min = NA_real_, spat_max = 100), NA)
})

test_that("MotifDB validates matrix dimensions", {
    motif_data <- tibble(
        motif = c("motif1"),
        pos = c(1),
        A = c(0.7),
        C = c(0.1),
        G = c(0.1),
        T = c(0.1)
    )
    expect_error(create_motif_db(motif_data), NA)
})

test_that("MotifDB validates prior constraints", {
    expect_error(create_test_motif_db(prior = 0))
    expect_error(create_test_motif_db(prior = 1))
    expect_error(create_test_motif_db(prior = -0.1))
})

test_that("MotifDB subsetting preserves spatial boundaries", {
    motif_db <- create_test_motif_db(spat_min = 0, spat_max = 100)

    # Test character subsetting
    subset_db <- motif_db["motif1"]
    expect_equal(subset_db@spat_min, 0)
    expect_equal(subset_db@spat_max, 100)

    # Test numeric subsetting
    subset_db2 <- motif_db[1]
    expect_equal(subset_db2@spat_min, 0)
    expect_equal(subset_db2@spat_max, 100)

    # Test multiple motif selection
    multi_subset <- motif_db[c("motif1", "motif2")]
    expect_equal(multi_subset@spat_min, 0)
    expect_equal(multi_subset@spat_max, 100)
})

test_that("motif_db_to_dataframe conversion works with spatial boundaries", {
    motif_db <- create_test_motif_db(spat_min = 0, spat_max = 100)
    df <- motif_db_to_dataframe(motif_db)

    # Conversion should work the same regardless of spatial boundaries
    motif_db_no_bounds <- create_test_motif_db()
    df_no_bounds <- motif_db_to_dataframe(motif_db_no_bounds)

    expect_equal(df, df_no_bounds)
})

test_that("spatial factors validation works", {
    motif_db <- create_test_motif_db()

    # Test with valid spatial factors
    spat_factors <- matrix(1,
        nrow = 2,
        ncol = 3,
        dimnames = list(c("motif1", "motif2"), NULL)
    )
    expect_error(
        create_motif_db(
            motif_db_to_dataframe(motif_db),
            spat_factors = spat_factors,
            spat_min = 0,
            spat_max = 100
        ),
        NA
    )

    # Test with invalid dimensions
    invalid_spat_factors <- matrix(1, nrow = 3, ncol = 3)
    expect_error(
        create_motif_db(
            motif_db_to_dataframe(motif_db),
            spat_factors = invalid_spat_factors
        )
    )

    # Test with negative values
    invalid_values_spat_factors <- spat_factors
    invalid_values_spat_factors[1, 1] <- -1
    expect_error(
        create_motif_db(
            motif_db_to_dataframe(motif_db),
            spat_factors = invalid_values_spat_factors
        )
    )
})

test_that("reverse complement matrix is correctly computed with spatial boundaries", {
    motif_db <- create_test_motif_db(spat_min = 0, spat_max = 100)

    # For a single position, check if reverse complement is correct
    check_rc <- function(pos, nuc_idx) {
        forward <- motif_db@mat[pos * 4 - (4 - nuc_idx), ]
        rc_pos <- motif_db@motif_lengths[1] - pos + 1
        rc_nuc_idx <- c(4, 3, 2, 1)[nuc_idx] # A->T, C->G, G->C, T->A
        reverse <- motif_db@rc_mat[rc_pos * 4 - (4 - rc_nuc_idx), ]
        expect_equal(forward, reverse)
    }

    # Check all positions and nucleotides
    for (pos in 1:4) {
        for (nuc_idx in 1:4) {
            check_rc(pos, nuc_idx)
        }
    }
})

test_that("prior modification preserves spatial boundaries", {
    # Create original object with spatial boundaries
    motif_db <- create_test_motif_db(prior = 0.01, spat_min = 0, spat_max = 100)

    # Modify prior
    new_prior <- 0.02
    prior(motif_db) <- new_prior

    # Check spatial boundaries are preserved
    expect_equal(motif_db@spat_min, 0)
    expect_equal(motif_db@spat_max, 100)
})

test_that("motif_lengths validation works", {
    motif_db <- create_test_motif_db()

    # Test that motif lengths match matrix dimensions
    expect_equal(
        length(motif_db@motif_lengths),
        ncol(motif_db@mat)
    )

    # Test that names match
    expect_equal(
        names(motif_db@motif_lengths),
        colnames(motif_db@mat)
    )

    # Test that lengths are positive
    expect_true(all(motif_db@motif_lengths > 0))
})

test_that("as.data.frame works correctly with spatial boundaries", {
    motif_db <- create_test_motif_db(spat_min = 0, spat_max = 100)

    # Test direct conversion
    df1 <- as.data.frame(motif_db)
    df2 <- as.data.frame(motif_db_to_dataframe(motif_db))
    # Results should be identical
    expect_equal(df1, df2)

    # Compare with non-bounded version
    motif_db_no_bounds <- create_test_motif_db()
    df_no_bounds <- as.data.frame(motif_db_no_bounds)
    expect_equal(df1, df_no_bounds)
})
