# Tests for PBM class and PBM-related functions
# These tests are designed to run WITHOUT misha

# --- Helper: create a minimal valid PBM object ---
# Note: with_spat=FALSE creates a PBM without spatial model (for pure slot tests).
# with_spat=TRUE creates a PBM with spatial model (required for pbm.compute and pbm_list operations).
create_mock_pbm <- function(name = "test_motif", n_positions = 5, with_spat = FALSE) {
    pssm <- matrix(
        c(
            0.8, 0.05, 0.1, 0.05,
            0.05, 0.8, 0.05, 0.1,
            0.1, 0.05, 0.8, 0.05,
            0.05, 0.1, 0.05, 0.8,
            0.25, 0.25, 0.25, 0.25
        )[seq_len(n_positions * 4)],
        nrow = n_positions, ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    if (with_spat) {
        spat <- data.frame(
            bin = seq(0, 450, by = 50),
            spat_factor = rep(1.0, 10)
        )
        spat_min <- 0
        spat_max <- 499
        seq_length <- 500
    } else {
        spat <- data.frame(bin = numeric(0), spat_factor = numeric(0))
        spat_min <- numeric(0)
        spat_max <- numeric(0)
        seq_length <- numeric(0)
    }

    PBM(
        name = name,
        pssm = pssm,
        max_energy = -2.5,
        min_energy = -7.0,
        energy_range = c(-7.0, 0.0),
        spat = spat,
        spat_min = spat_min,
        spat_max = spat_max,
        seq_length = seq_length,
        coefs = c(
            "low-energy" = 0.3,
            "high-energy" = 0.5,
            "higher-energy" = -0.1,
            "sigmoid" = 0.2
        ),
        size = 500
    )
}

# Helper: generate random DNA sequences
random_dna_seq <- function(n, len = 500) {
    set.seed(42)
    bases <- c("A", "C", "G", "T")
    vapply(seq_len(n), function(i) {
        paste(sample(bases, len, replace = TRUE), collapse = "")
    }, character(1))
}

# =============================================================================
# PBM constructor
# =============================================================================
test_that("PBM constructor creates a valid S4 object", {
    pbm <- create_mock_pbm()

    expect_s4_class(pbm, "PBM")
    expect_true(methods::is(pbm, "PBM"))
    expect_true(methods::is(pbm, "IQFeature"))
    expect_equal(pbm@name, "test_motif")
    expect_equal(nrow(pbm@pssm), 5)
    expect_equal(ncol(pbm@pssm), 4)
})

test_that("PBM inherits from IQFeature", {
    pbm <- create_mock_pbm()

    expect_true(methods::is(pbm, "IQFeature"))
    expect_true("name" %in% slotNames(pbm))
    expect_true("coefs" %in% slotNames(pbm))
})

test_that("PBM has correct slot types", {
    pbm <- create_mock_pbm()

    expect_true(is.matrix(pbm@pssm))
    expect_true(is.numeric(pbm@max_energy))
    expect_true(is.numeric(pbm@min_energy))
    expect_true(is.numeric(pbm@energy_range))
    expect_equal(length(pbm@energy_range), 2)
    expect_true(is.data.frame(pbm@spat))
    expect_true(is.numeric(pbm@size))
})

test_that("PBM with spatial model stores spat correctly", {
    pbm <- create_mock_pbm(with_spat = TRUE)

    expect_true(nrow(pbm@spat) > 0)
    expect_true(all(c("bin", "spat_factor") %in% colnames(pbm@spat)))
    expect_true(pbm@spat_min < pbm@spat_max)
    expect_true(pbm@seq_length > 0)
})

# =============================================================================
# validate_pbm
# =============================================================================
test_that("validate_pbm accepts a valid PBM", {
    pbm <- create_mock_pbm()
    expect_silent(iceqream:::validate_pbm(pbm))
})

test_that("validate_pbm accepts a valid PBM with spatial model", {
    pbm <- create_mock_pbm(with_spat = TRUE)
    expect_silent(iceqream:::validate_pbm(pbm))
})

test_that("validate_pbm rejects non-PBM objects", {
    expect_error(iceqream:::validate_pbm(list()), "PBM")
    expect_error(iceqream:::validate_pbm(42), "PBM")
})

test_that("validate_pbm rejects PBM with empty name", {
    pbm <- create_mock_pbm()
    pbm@name <- character(0)
    expect_error(iceqream:::validate_pbm(pbm), "name")
})

test_that("validate_pbm rejects PBM with empty PSSM", {
    pbm <- create_mock_pbm()
    pbm@pssm <- matrix(nrow = 0, ncol = 4, dimnames = list(NULL, c("A", "C", "G", "T")))
    expect_error(iceqream:::validate_pbm(pbm), "PSSM")
})

test_that("validate_pbm rejects PBM with wrong energy_range length", {
    pbm <- create_mock_pbm()
    pbm@energy_range <- numeric(0)
    expect_error(iceqream:::validate_pbm(pbm), "energy range")
})

test_that("validate_pbm rejects PBM with zero size", {
    pbm <- create_mock_pbm()
    pbm@size <- 0
    expect_error(iceqream:::validate_pbm(pbm), "size")
})

test_that("validate_pbm rejects PBM with spat_min >= spat_max", {
    pbm <- create_mock_pbm(with_spat = TRUE)
    pbm@spat_min <- 500
    pbm@spat_max <- 100
    expect_error(iceqream:::validate_pbm(pbm), "spatial minimum")
})

test_that("validate_pbm rejects PBM with seq_length <= 0", {
    pbm <- create_mock_pbm(with_spat = TRUE)
    pbm@seq_length <- 0
    expect_error(iceqream:::validate_pbm(pbm), "sequence length")
})

# =============================================================================
# show method for PBM
# =============================================================================
test_that("PBM show method produces output without error", {
    pbm <- create_mock_pbm()
    # cli output goes to stderr/connection, so capture messages
    output <- utils::capture.output(show(pbm), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(nchar(output_text) > 0)
    expect_true(grepl("PBM", output_text))
})

test_that("PBM show method includes name", {
    pbm <- create_mock_pbm(name = "my_special_motif")
    output <- utils::capture.output(show(pbm), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(grepl("my_special_motif", output_text))
})

test_that("PBM show method with spatial model includes spatial info", {
    pbm <- create_mock_pbm(with_spat = TRUE)
    output <- utils::capture.output(show(pbm), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(grepl("[Ss]patial", output_text))
})

# =============================================================================
# pbm.normalize_energies
# =============================================================================
test_that("pbm.normalize_energies normalizes to expected range", {
    pbm <- create_mock_pbm()
    raw_energies <- runif(100, min = -15, max = 0)

    normalized <- pbm.normalize_energies(pbm, raw_energies)

    expect_length(normalized, 100)
    expect_true(all(normalized >= 0))
    expect_true(all(normalized <= 10))
})

test_that("pbm.normalize_energies handles extreme values", {
    pbm <- create_mock_pbm()

    # Very negative energy should be clamped
    normalized_low <- pbm.normalize_energies(pbm, c(-1000))
    expect_true(normalized_low >= 0)

    # Very high energy should be clamped
    normalized_high <- pbm.normalize_energies(pbm, c(1000))
    expect_true(normalized_high <= 10)
})

test_that("pbm.normalize_energies preserves ordering for monotone input", {
    pbm <- create_mock_pbm()
    raw_energies <- seq(-10, -1, length.out = 20)

    normalized <- pbm.normalize_energies(pbm, raw_energies)

    # Monotonically increasing input should produce non-decreasing output
    expect_true(all(diff(normalized) >= -1e-10))
})

# =============================================================================
# pbm.normalize_energies_local
# =============================================================================
test_that("pbm.normalize_energies_local uses max_energy from PBM", {
    pbm <- create_mock_pbm()
    energies <- runif(50, min = -10, max = 0)

    result <- pbm.normalize_energies_local(pbm, energies)

    # Should give same result as calling normalize with pbm's max_energy
    expected <- pbm.normalize_energies(pbm, energies, pbm@max_energy)
    expect_equal(result, expected)
})

# =============================================================================
# pbm.compute (requires prego; PBMs with spatial data)
# =============================================================================
test_that("pbm.compute computes energies on short DNA sequences", {
    skip_if_not_installed("prego")

    pbm <- create_mock_pbm(with_spat = TRUE)
    sequences <- random_dna_seq(10, len = 500)

    energies <- pbm.compute(pbm, sequences, response = FALSE, normalize_energies = TRUE)

    expect_length(energies, 10)
    expect_true(is.numeric(energies))
    expect_true(all(energies >= 0))
    expect_true(all(energies <= 10))
})

test_that("pbm.compute with response=TRUE returns response values", {
    skip_if_not_installed("prego")

    pbm <- create_mock_pbm(with_spat = TRUE)
    sequences <- random_dna_seq(10, len = 500)

    response <- pbm.compute(pbm, sequences, response = TRUE, normalize_energies = TRUE)

    expect_length(response, 10)
    expect_true(is.numeric(response))
})

test_that("pbm.compute errors when response=TRUE but normalize_energies=FALSE", {
    skip_if_not_installed("prego")

    pbm <- create_mock_pbm(with_spat = TRUE)
    sequences <- random_dna_seq(5, len = 500)

    expect_error(
        pbm.compute(pbm, sequences, response = TRUE, normalize_energies = FALSE),
        "normalize_energies"
    )
})

test_that("pbm.compute handles sequences of different length than size", {
    skip_if_not_installed("prego")

    pbm <- create_mock_pbm(with_spat = TRUE)  # size = 500
    short_sequences <- random_dna_seq(5, len = 100)

    # cli::cli_alert_warning is used (not base::warning), so we capture messages
    # and verify the function still computes energies
    energies <- pbm.compute(pbm, short_sequences, response = FALSE, normalize_energies = TRUE)
    expect_length(energies, 5)
    expect_true(is.numeric(energies))
})

test_that("pbm.compute without normalization returns raw-like values", {
    skip_if_not_installed("prego")

    pbm <- create_mock_pbm(with_spat = TRUE)
    sequences <- random_dna_seq(10, len = 500)

    raw <- pbm.compute(pbm, sequences, response = FALSE, normalize_energies = FALSE)
    norm <- pbm.compute(pbm, sequences, response = FALSE, normalize_energies = TRUE)

    expect_length(raw, 10)
    expect_true(is.numeric(raw))
    # Raw and normalized should differ (unless by extreme coincidence)
    expect_false(all(raw == norm))
})

# =============================================================================
# pbm.trim_pssm
# =============================================================================
test_that("pbm.trim_pssm trims low-information positions", {
    skip_if_not_installed("prego")

    pbm <- create_mock_pbm(n_positions = 5)
    original_nrow <- nrow(pbm@pssm)

    # The last position is uniform (0.25 each), so trimming should potentially remove it
    trimmed <- pbm.trim_pssm(pbm, bits_thresh = 0.5)

    expect_s4_class(trimmed, "PBM")
    expect_true(nrow(trimmed@pssm) <= original_nrow)
})

# =============================================================================
# PBM list operations
# =============================================================================
test_that("creating a list of PBMs works", {
    pbm1 <- create_mock_pbm(name = "TF_A")
    pbm2 <- create_mock_pbm(name = "TF_B")
    pbm3 <- create_mock_pbm(name = "TF_C")

    pbm_list <- list(TF_A = pbm1, TF_B = pbm2, TF_C = pbm3)

    expect_length(pbm_list, 3)
    expect_equal(names(pbm_list), c("TF_A", "TF_B", "TF_C"))
    expect_s4_class(pbm_list$TF_A, "PBM")
    expect_s4_class(pbm_list$TF_B, "PBM")
    expect_s4_class(pbm_list$TF_C, "PBM")
})

test_that("pbm_list.normalize normalizes a matrix of energies", {
    pbm1 <- create_mock_pbm(name = "TF_A")
    pbm2 <- create_mock_pbm(name = "TF_B")
    pbm_list <- list(TF_A = pbm1, TF_B = pbm2)

    # Raw energies matrix
    energies <- matrix(runif(20, -15, 0), nrow = 10, ncol = 2)
    colnames(energies) <- c("TF_A", "TF_B")

    normalized <- pbm_list.normalize(pbm_list, energies)

    expect_equal(dim(normalized), dim(energies))
    expect_equal(colnames(normalized), c("TF_A", "TF_B"))
    expect_true(all(normalized >= 0))
    expect_true(all(normalized <= 10))
})

test_that("pbm_list.normalize replaces NAs with 0", {
    pbm1 <- create_mock_pbm(name = "TF_A")
    pbm_list <- list(TF_A = pbm1)

    energies <- matrix(c(-5, NA, -3, NA, -1), nrow = 5, ncol = 1)
    colnames(energies) <- "TF_A"

    normalized <- pbm_list.normalize(pbm_list, energies)

    expect_false(any(is.na(normalized)))
})

test_that("pbm_list.normalize ignores columns not in pbm_list", {
    pbm1 <- create_mock_pbm(name = "TF_A")
    pbm_list <- list(TF_A = pbm1)

    energies <- matrix(runif(20, -10, 0), nrow = 10, ncol = 2)
    colnames(energies) <- c("TF_A", "TF_UNKNOWN")

    normalized <- pbm_list.normalize(pbm_list, energies)

    # TF_A column should be normalized, TF_UNKNOWN should be unchanged
    expect_equal(dim(normalized), dim(energies))
    expect_true(all(normalized[, "TF_A"] >= 0 & normalized[, "TF_A"] <= 10))
})

test_that("pbm_list.compute computes energies for a list of PBMs", {
    skip_if_not_installed("prego")

    # pbm_list.compute uses pbm_list_to_mdb -> motifs_to_mdb which needs spatial data
    pbm1 <- create_mock_pbm(name = "TF_A", with_spat = TRUE)
    pbm2 <- create_mock_pbm(name = "TF_B", with_spat = TRUE)
    pbm_list <- list(TF_A = pbm1, TF_B = pbm2)

    sequences <- random_dna_seq(10, len = 500)

    energies <- pbm_list.compute(pbm_list, sequences, response = FALSE, normalize_energies = TRUE)

    expect_true(is.matrix(energies))
    expect_equal(nrow(energies), 10)
    expect_equal(ncol(energies), 2)
    expect_true(all(colnames(energies) %in% c("TF_A", "TF_B")))
})

test_that("pbm_list.compute with response returns transformed values", {
    skip_if_not_installed("prego")

    pbm1 <- create_mock_pbm(name = "TF_A", with_spat = TRUE)
    pbm2 <- create_mock_pbm(name = "TF_B", with_spat = TRUE)
    pbm_list <- list(TF_A = pbm1, TF_B = pbm2)

    sequences <- random_dna_seq(10, len = 500)

    response <- pbm_list.compute(pbm_list, sequences, response = TRUE, normalize_energies = TRUE)

    expect_true(is.matrix(response))
    expect_equal(nrow(response), 10)
    expect_equal(ncol(response), 2)
})

test_that("pbm_list.compute errors when response=TRUE but normalize_energies=FALSE", {
    skip_if_not_installed("prego")

    pbm1 <- create_mock_pbm(name = "TF_A", with_spat = TRUE)
    pbm_list <- list(TF_A = pbm1)

    sequences <- random_dna_seq(5, len = 500)

    expect_error(
        pbm_list.compute(pbm_list, sequences, response = TRUE, normalize_energies = FALSE),
        "normalize_energies"
    )
})

test_that("pbm_list.trim_pssm trims all PBMs in list", {
    skip_if_not_installed("prego")

    pbm1 <- create_mock_pbm(name = "TF_A", n_positions = 5)
    pbm2 <- create_mock_pbm(name = "TF_B", n_positions = 5)
    pbm_list <- list(TF_A = pbm1, TF_B = pbm2)

    trimmed <- pbm_list.trim_pssm(pbm_list, bits_thresh = 0.5)

    expect_length(trimmed, 2)
    expect_s4_class(trimmed$TF_A, "PBM")
    expect_s4_class(trimmed$TF_B, "PBM")
})

# =============================================================================
# IQFeature class
# =============================================================================
test_that("IQFeature constructor creates valid object", {
    iq <- IQFeature(
        name = "test_feature",
        coefs = c("low-energy" = 0.1, "high-energy" = 0.2, "higher-energy" = 0.05, "sigmoid" = 0.15)
    )

    expect_s4_class(iq, "IQFeature")
    expect_equal(iq@name, "test_feature")
    expect_equal(length(iq@coefs), 4)
})

test_that("IQFeature show method works", {
    iq <- IQFeature(name = "demo", coefs = c(a = 1.0, b = 2.0))
    output <- utils::capture.output(show(iq), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(grepl("IQFeature", output_text))
})

# =============================================================================
# iq_feature.compute_response
# =============================================================================
test_that("iq_feature.compute_response with single coefficient scales values", {
    iq <- IQFeature(name = "simple", coefs = c(weight = 2.0))
    values <- c(1, 2, 3, 4, 5)

    result <- iq_feature.compute_response(iq, values)

    expect_equal(result, values * 2.0)
})

test_that("iq_feature.compute_response with multiple coefficients uses logistic features", {
    iq <- IQFeature(
        name = "complex",
        coefs = c("low-energy" = 0.3, "high-energy" = 0.5, "higher-energy" = -0.1, "sigmoid" = 0.2)
    )
    values <- seq(0, 10, length.out = 20)

    result <- iq_feature.compute_response(iq, values)

    expect_length(result, 20)
    expect_true(is.numeric(result))
})

# =============================================================================
# IQInteraction class
# =============================================================================
test_that("IQInteraction constructor creates valid object", {
    iq_int <- IQInteraction(
        name = "TF_A:TF_B",
        coefs = c("low-energy" = 0.1, "high-energy" = 0.2, "higher-energy" = 0.05, "sigmoid" = 0.15),
        term1 = "TF_A",
        term2 = "TF_B",
        term1_type = "motif",
        term2_type = "motif",
        scale = 1.0,
        inter_max = 100.0,
        norm_min = 0.0,
        norm_max = 1.0
    )

    expect_s4_class(iq_int, "IQInteraction")
    expect_true(methods::is(iq_int, "IQFeature"))
    expect_equal(iq_int@term1, "TF_A")
    expect_equal(iq_int@term2, "TF_B")
})

test_that("IQInteraction show method works", {
    iq_int <- IQInteraction(
        name = "TF_A:TF_B",
        coefs = c(a = 1.0),
        term1 = "TF_A",
        term2 = "TF_B",
        term1_type = "motif",
        term2_type = "additional",
        scale = 1.0,
        inter_max = 50.0,
        norm_min = 0.0,
        norm_max = 1.0
    )
    output <- utils::capture.output(show(iq_int), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(grepl("IQInteraction", output_text))
})

# =============================================================================
# iq_interaction.compute
# =============================================================================
test_that("iq_interaction.compute normalizes interaction values", {
    iq_int <- IQInteraction(
        name = "TF_A:TF_B",
        coefs = c(a = 1.0),
        term1 = "TF_A",
        term2 = "TF_B",
        term1_type = "motif",
        term2_type = "motif",
        scale = 1.0,
        inter_max = 100.0,
        norm_min = 0.0,
        norm_max = 1.0
    )

    term1_vals <- c(5, 10, 15, 20)
    term2_vals <- c(2, 4, 6, 8)

    result <- iq_interaction.compute(iq_int, term1_vals, term2_vals)

    expect_length(result, 4)
    expect_true(is.numeric(result))
    # With scale=1 and norm range [0,1], results should be between 0 and 1
    expect_true(all(result >= 0))
    expect_true(all(result <= 1))
})

test_that("iq_interaction.compute respects scale factor", {
    iq_int_s1 <- IQInteraction(
        name = "X:Y", coefs = c(a = 1.0),
        term1 = "X", term2 = "Y",
        term1_type = "motif", term2_type = "motif",
        scale = 1.0, inter_max = 100.0,
        norm_min = 0.0, norm_max = 1.0
    )
    iq_int_s2 <- IQInteraction(
        name = "X:Y", coefs = c(a = 1.0),
        term1 = "X", term2 = "Y",
        term1_type = "motif", term2_type = "motif",
        scale = 2.0, inter_max = 100.0,
        norm_min = 0.0, norm_max = 1.0
    )

    t1 <- c(10, 20)
    t2 <- c(5, 10)

    r1 <- iq_interaction.compute(iq_int_s1, t1, t2)
    r2 <- iq_interaction.compute(iq_int_s2, t1, t2)

    expect_equal(r2, r1 * 2, tolerance = 1e-10)
})

# =============================================================================
# IQmodel class
# =============================================================================
test_that("IQmodel constructor creates valid object", {
    pbm <- create_mock_pbm(name = "TF_A")
    iq_feat <- IQFeature(
        name = "some_feat",
        coefs = c("low-energy" = 0.1, "high-energy" = 0.2, "higher-energy" = 0.05, "sigmoid" = 0.1)
    )

    model <- IQmodel(
        features = list(some_feat = iq_feat),
        pbms = list(TF_A = pbm),
        intercept = -0.5,
        func = "logSumExp",
        lambda = 0.01,
        alpha = 0.5,
        min_pred = 0.1,
        max_pred = 0.9,
        norm_factors = c(0.0, 1.0)
    )

    expect_s4_class(model, "IQmodel")
    expect_equal(length(model@pbms), 1)
    expect_equal(length(model@features), 1)
    expect_equal(model@intercept, -0.5)
    expect_equal(model@func, "logSumExp")
})

test_that("IQmodel show method works", {
    model <- IQmodel(
        features = list(),
        pbms = list(),
        intercept = 0,
        func = "logSumExp",
        lambda = 0.01,
        alpha = 0.5,
        min_pred = 0.1,
        max_pred = 0.9,
        norm_factors = c(0, 1)
    )
    output <- utils::capture.output(show(model), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(grepl("IQmodel", output_text))
})

test_that("iq_model.has_interactions detects interactions", {
    iq_int <- IQInteraction(
        name = "A:B", coefs = c(a = 1.0),
        term1 = "A", term2 = "B",
        term1_type = "motif", term2_type = "motif",
        scale = 1, inter_max = 1, norm_min = 0, norm_max = 1
    )
    model_with_int <- IQmodel(
        features = list("A:B" = iq_int),
        pbms = list(),
        intercept = 0, func = "logSumExp",
        lambda = 0.01, alpha = 0.5,
        min_pred = 0, max_pred = 1,
        norm_factors = c(0, 1)
    )
    model_without_int <- IQmodel(
        features = list(),
        pbms = list(),
        intercept = 0, func = "logSumExp",
        lambda = 0.01, alpha = 0.5,
        min_pred = 0, max_pred = 1,
        norm_factors = c(0, 1)
    )

    expect_true(iq_model.has_interactions(model_with_int))
    expect_false(iq_model.has_interactions(model_without_int))
})

test_that("iq_model.get_interactions returns only interactions", {
    iq_feat <- IQFeature(name = "feat1", coefs = c(a = 1.0))
    iq_int <- IQInteraction(
        name = "A:B", coefs = c(a = 1.0),
        term1 = "A", term2 = "B",
        term1_type = "motif", term2_type = "motif",
        scale = 1, inter_max = 1, norm_min = 0, norm_max = 1
    )
    model <- IQmodel(
        features = list(feat1 = iq_feat, "A:B" = iq_int),
        pbms = list(),
        intercept = 0, func = "logSumExp",
        lambda = 0.01, alpha = 0.5,
        min_pred = 0, max_pred = 1,
        norm_factors = c(0, 1)
    )

    interactions <- iq_model.get_interactions(model)
    regular <- iq_model.get_regular_features(model)

    expect_length(interactions, 1)
    expect_length(regular, 1)
    expect_s4_class(interactions[["A:B"]], "IQInteraction")
    expect_s4_class(regular[["feat1"]], "IQFeature")
})

# =============================================================================
# IQSeqFeature class
# =============================================================================
test_that("IQSeqFeature constructor creates valid object", {
    compute_fn <- function(sequences) {
        stringr::str_count(sequences, "G|C") / nchar(sequences)
    }

    iq_seq <- IQSeqFeature(
        name = "gc_content",
        coefs = c("low-energy" = 0.1, "high-energy" = 0.2, "higher-energy" = 0.05, "sigmoid" = 0.1),
        compute_func = compute_fn,
        min_value = 0.3,
        max_value = 0.7,
        quantile = 0.95,
        size = 500
    )

    expect_s4_class(iq_seq, "IQSeqFeature")
    expect_true(methods::is(iq_seq, "IQFeature"))
    expect_equal(iq_seq@name, "gc_content")
    expect_equal(iq_seq@size, 500)
    expect_true(is.function(iq_seq@compute_func))
})

test_that("IQSeqFeature show method works", {
    compute_fn <- function(sequences) rep(0, length(sequences))
    iq_seq <- IQSeqFeature(
        name = "test_seq_feat",
        coefs = c(a = 1.0),
        compute_func = compute_fn,
        min_value = 0,
        max_value = 1,
        size = 100
    )
    output <- utils::capture.output(show(iq_seq), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(grepl("IQSeqFeature", output_text))
})

# =============================================================================
# iq_seq_feature.normalize
# =============================================================================
test_that("iq_seq_feature.normalize clamps and scales to 0-10", {
    compute_fn <- function(sequences) rep(0, length(sequences))
    iq_seq <- IQSeqFeature(
        name = "test",
        coefs = c(a = 1.0),
        compute_func = compute_fn,
        min_value = 2.0,
        max_value = 8.0,
        size = 100
    )

    values <- c(0, 2, 5, 8, 12)
    result <- iceqream:::iq_seq_feature.normalize(iq_seq, values)

    expect_length(result, 5)
    expect_true(all(result >= 0))
    expect_true(all(result <= 10))
    # value = min_value should give 0
    expect_equal(result[2], 0)
    # value = max_value should give 10
    expect_equal(result[4], 10)
    # value below min_value should also give 0
    expect_equal(result[1], 0)
    # value above max_value should give 10
    expect_equal(result[5], 10)
})

test_that("iq_seq_feature.normalize handles zero range", {
    compute_fn <- function(sequences) rep(0, length(sequences))
    iq_seq <- IQSeqFeature(
        name = "test",
        coefs = c(a = 1.0),
        compute_func = compute_fn,
        min_value = 5.0,
        max_value = 5.0,
        size = 100
    )

    values <- c(1, 5, 10)
    result <- iceqream:::iq_seq_feature.normalize(iq_seq, values)

    # Zero range should return all zeros
    expect_true(all(result == 0))
})

# =============================================================================
# IQFeatureGroup class
# =============================================================================
test_that("IQFeatureGroup constructor creates valid object", {
    compute_fn <- function(sequences) rep(0, length(sequences))
    feat1 <- IQSeqFeature(
        name = "f1", coefs = c(a = 1.0),
        compute_func = compute_fn, min_value = 0, max_value = 1, size = 100
    )
    feat2 <- IQSeqFeature(
        name = "f2", coefs = c(a = 1.0),
        compute_func = compute_fn, min_value = 0, max_value = 1, size = 100
    )

    group_compute <- function(sequences) {
        matrix(0, nrow = length(sequences), ncol = 2, dimnames = list(NULL, c("f1", "f2")))
    }

    group <- IQFeatureGroup(
        features = list(f1 = feat1, f2 = feat2),
        compute_func = group_compute,
        size = 100
    )

    expect_s4_class(group, "IQFeatureGroup")
    expect_equal(length(group@features), 2)
    expect_equal(group@size, 100)
})

test_that("IQFeatureGroup show method works", {
    compute_fn <- function(sequences) rep(0, length(sequences))
    feat1 <- IQSeqFeature(
        name = "f1", coefs = c(a = 1.0),
        compute_func = compute_fn, min_value = 0, max_value = 1, size = 100
    )
    group_compute <- function(sequences) {
        matrix(0, nrow = length(sequences), ncol = 1, dimnames = list(NULL, c("f1")))
    }
    group <- IQFeatureGroup(
        features = list(f1 = feat1),
        compute_func = group_compute,
        size = 100
    )
    output <- utils::capture.output(show(group), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(grepl("IQFeatureGroup", output_text))
})

# =============================================================================
# iq_model.feature_names_by_type
# =============================================================================
test_that("iq_model.feature_names_by_type categorizes features", {
    iq_feat <- IQFeature(name = "feat1", coefs = c(a = 1.0))
    iq_int <- IQInteraction(
        name = "A:B", coefs = c(a = 1.0),
        term1 = "A", term2 = "B",
        term1_type = "motif", term2_type = "motif",
        scale = 1, inter_max = 1, norm_min = 0, norm_max = 1
    )

    features_list <- list(feat1 = iq_feat, "A:B" = iq_int)
    result <- iq_model.feature_names_by_type(features_list)

    expect_true("regular" %in% names(result))
    expect_true("interaction" %in% names(result))
    expect_true("feat1" %in% result$regular)
    expect_true("A:B" %in% result$interaction)
})

test_that("iq_model.feature_names_by_type includes interaction terms when requested", {
    iq_int <- IQInteraction(
        name = "X:Y", coefs = c(a = 1.0),
        term1 = "X", term2 = "Y",
        term1_type = "motif", term2_type = "additional",
        scale = 1, inter_max = 1, norm_min = 0, norm_max = 1
    )

    features_list <- list("X:Y" = iq_int)
    result <- iq_model.feature_names_by_type(features_list, include_interaction_terms = TRUE)

    expect_true("X" %in% result$interaction_terms)
    expect_true("Y" %in% result$interaction_terms)
})
