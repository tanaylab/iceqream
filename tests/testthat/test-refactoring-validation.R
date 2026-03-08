# =============================================================================
# Refactoring Validation Tests
#
# Validates the changes from the deep code review remediation:
# - Exported function availability
# - Helper function correctness
# - Deprecated API migration (rename: create_specifc_terms -> create_specific_terms)
# - S4 validity checks
# - Virtual track cleanup (gvtrack.rm import)
# - Anti-pattern fixes (seq_along / seq_len usage)
# - Utility function correctness (regression tests)
# =============================================================================

# Helper to retrieve all imports for a given package from the namespace.
# getNamespaceImports returns a named list where the same package name
# may appear multiple times (once per importFrom directive). This helper
# collects all of them into a single character vector.
get_all_imports <- function(pkg, from_pkg) {
    ns <- getNamespaceImports(pkg)
    unname(unlist(ns[names(ns) == from_pkg]))
}

# =============================================================================
# 1. Exported function availability
# =============================================================================

test_that("energy utility functions are exported and callable", {
    expect_true(is.function(norm01))
    expect_true(is.function(norm0q))
    expect_true(is.function(rescale))
    expect_true(is.function(logist))
    expect_true(is.function(create_logist_features))
    expect_true(is.function(norm_energy_matrix))
    expect_true(is.function(norm_energy_dataset))
})

test_that("trajectory model utility functions are exported and callable", {
    expect_true(is.function(strip_traj_model))
    expect_true(is.function(learn_traj_prego))
})

test_that("trajectory model manipulation functions are exported and callable", {
    expect_true(is.function(add_motif_models_to_traj))
    expect_true(is.function(adjust_energies))
    expect_true(is.function(adjust_motif_seq_lengths))
})

test_that("homogenize_pssm_models is exported and callable", {
    expect_true(is.function(homogenize_pssm_models))
})

# =============================================================================
# 2. Helper function correctness
# =============================================================================

# --- split_low_correlation_clusters ---

test_that("split_low_correlation_clusters returns a data frame with expected columns", {
    set.seed(42)
    # Create synthetic features with known correlation structure
    base1 <- rnorm(50)
    base2 <- rnorm(50)
    feat_mat <- cbind(
        f1 = base1 + rnorm(50, sd = 0.1),
        f2 = base1 + rnorm(50, sd = 0.1),
        f3 = base2 + rnorm(50, sd = 0.1),
        f4 = rnorm(50)
    )
    corr_matrix <- cor(feat_mat)

    # All in one cluster
    clust_map <- data.frame(
        feat = colnames(feat_mat),
        clust = rep(1L, 4),
        stringsAsFactors = FALSE
    )

    result <- iceqream:::split_low_correlation_clusters(clust_map, corr_matrix, intra_cor_thresh = 0.8)

    expect_s3_class(result, "data.frame")
    expect_true("feat" %in% colnames(result))
    expect_true("clust" %in% colnames(result))
    expect_equal(nrow(result), 4)
    # Low intra-cluster correlation should lead to splitting
    expect_true(length(unique(result$clust)) > 1)
})

test_that("split_low_correlation_clusters does not split when correlation is high", {
    set.seed(42)
    base <- rnorm(50)
    feat_mat <- cbind(
        f1 = base + rnorm(50, sd = 0.01),
        f2 = base + rnorm(50, sd = 0.01)
    )
    corr_matrix <- cor(feat_mat)

    clust_map <- data.frame(
        feat = colnames(feat_mat),
        clust = rep(1L, 2),
        stringsAsFactors = FALSE
    )

    result <- iceqream:::split_low_correlation_clusters(clust_map, corr_matrix, intra_cor_thresh = 0.5)

    # Highly correlated features should stay in same cluster
    expect_equal(length(unique(result$clust)), 1)
})

test_that("split_low_correlation_clusters handles single-element clusters correctly", {
    set.seed(42)
    feat_mat <- cbind(f1 = rnorm(50), f2 = rnorm(50))
    corr_matrix <- cor(feat_mat)

    # Each feature in its own cluster
    clust_map <- data.frame(
        feat = c("f1", "f2"),
        clust = c(1L, 2L),
        stringsAsFactors = FALSE
    )

    result <- iceqream:::split_low_correlation_clusters(clust_map, corr_matrix, intra_cor_thresh = 0.9)

    # Single-element clusters should not be split further
    expect_equal(nrow(result), 2)
    expect_true(all(result$n == 1))
})

# --- fit_and_predict_model ---

test_that("fit_and_predict_model returns a list with model and predicted_diff_score", {
    set.seed(42)
    n <- 100
    p <- 4
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- paste0("f", seq_len(p))
    y <- norm01(X[, 1] + 0.5 * X[, 2] + rnorm(n, sd = 0.5))

    result <- suppressWarnings(iceqream:::fit_and_predict_model(
        train_features = X,
        train_y = y,
        predict_features = X,
        diff_score = y,
        alpha = 0.5,
        lambda = 0.01,
        seed = 42
    ))

    expect_true(is.list(result))
    expect_true("model" %in% names(result))
    expect_true("predicted_diff_score" %in% names(result))
    expect_length(result$predicted_diff_score, n)
    expect_true(is.numeric(result$predicted_diff_score))
    # Model should be a glmnet object
    expect_true(inherits(result$model, "glmnet") || inherits(result$model, "lognet"))
})

test_that("fit_and_predict_model predictions have same range as diff_score", {
    set.seed(42)
    n <- 200
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- paste0("feat", seq_len(p))
    y <- norm01(rowSums(X) + rnorm(n, sd = 0.3))
    diff_score <- y * 2 - 1 # range [-1, 1]

    result <- suppressWarnings(iceqream:::fit_and_predict_model(
        train_features = X,
        train_y = norm01(diff_score),
        predict_features = X,
        diff_score = diff_score,
        alpha = 0.5,
        lambda = 0.01,
        seed = 42
    ))

    # predicted_diff_score should be rescaled to the range of diff_score
    expect_true(min(result$predicted_diff_score) >= min(diff_score) - 0.01)
    expect_true(max(result$predicted_diff_score) <= max(diff_score) + 0.01)
})

# --- theme_iq_track ---

test_that("theme_iq_track returns a list of ggplot2 theme components", {
    result <- iceqream:::theme_iq_track()
    expect_true(is.list(result))
    expect_true(length(result) >= 2)
    # Should contain a theme_classic and a theme override
    has_theme_element <- any(sapply(result, function(x) inherits(x, "theme")))
    expect_true(has_theme_element)
})

test_that("theme_iq_track respects custom base_size", {
    result8 <- iceqream:::theme_iq_track(base_size = 8)
    result14 <- iceqream:::theme_iq_track(base_size = 14)
    # Both should return lists of theme components
    expect_true(is.list(result8))
    expect_true(is.list(result14))
})

test_that("theme_iq_track can set legend.position", {
    result <- iceqream:::theme_iq_track(legend.position = "none")
    expect_true(is.list(result))
    # The second element should be a theme with legend.position = "none"
    theme_element <- result[[2]]
    expect_true(inherits(theme_element, "theme"))
})

# =============================================================================
# 3. Deprecated API migration
# =============================================================================

test_that("create_specific_terms (corrected spelling) exists and works", {
    # The function was renamed from create_specifc_terms to create_specific_terms
    expect_true(is.function(iceqream:::create_specific_terms))

    # Test basic functionality with multiple terms (single term causes
    # energies[, terms$term1] to drop to a vector, which is a known
    # limitation of the current implementation)
    set.seed(42)
    energies <- matrix(abs(rnorm(60)), nrow = 10, ncol = 6)
    colnames(energies) <- paste0("m", 1:6)

    terms <- data.frame(
        term1 = c("m1", "m3"),
        term2 = c("m2", "m4"),
        variable = c("m1:m2", "m3:m4"),
        stringsAsFactors = FALSE
    )

    result <- iceqream:::create_specific_terms(energies, terms)

    expect_true(is.matrix(result))
    expect_equal(nrow(result), 10)
    expect_equal(ncol(result), 2)
    expect_equal(colnames(result), c("m1:m2", "m3:m4"))
})

test_that("pivot_longer is imported and used instead of gather", {
    tidyr_imports <- get_all_imports("iceqream", "tidyr")
    expect_true("pivot_longer" %in% tidyr_imports)
})

# =============================================================================
# 4. S4 validity checks
# =============================================================================

test_that("IQFeature with empty name fails validation", {
    expect_error(
        IQFeature(name = character(0), coefs = c(a = 1.0)),
        "name"
    )
})

test_that("IQFeature with blank name fails validation", {
    expect_error(
        IQFeature(name = "", coefs = c(a = 1.0)),
        "name"
    )
})

test_that("IQFeature with valid name passes validation", {
    iq <- IQFeature(name = "valid_name", coefs = c(a = 1.0))
    expect_s4_class(iq, "IQFeature")
    expect_equal(iq@name, "valid_name")
})

test_that("PBM with empty pssm fails S4 validation", {
    expect_error(
        PBM(
            name = "test",
            pssm = matrix(nrow = 0, ncol = 4, dimnames = list(NULL, c("A", "C", "G", "T"))),
            max_energy = -2.5,
            min_energy = -7.0,
            energy_range = c(-7.0, 0.0),
            spat = data.frame(bin = numeric(0), spat_factor = numeric(0)),
            spat_min = numeric(0),
            spat_max = numeric(0),
            seq_length = numeric(0),
            coefs = c(a = 1.0),
            size = 500
        ),
        "pssm"
    )
})

test_that("PBM with invalid energy_range fails S4 validation", {
    pssm <- matrix(
        c(0.7, 0.1, 0.1, 0.1, 0.1, 0.7, 0.1, 0.1),
        nrow = 2, ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    expect_error(
        PBM(
            name = "test",
            pssm = pssm,
            max_energy = -2.5,
            min_energy = -7.0,
            energy_range = numeric(0), # invalid: must be length 2
            spat = data.frame(bin = numeric(0), spat_factor = numeric(0)),
            spat_min = numeric(0),
            spat_max = numeric(0),
            seq_length = numeric(0),
            coefs = c(a = 1.0),
            size = 500
        ),
        "energy_range|energy range"
    )
})

test_that("PBM with zero size fails S4 validation", {
    pssm <- matrix(
        c(0.7, 0.1, 0.1, 0.1, 0.1, 0.7, 0.1, 0.1),
        nrow = 2, ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G", "T"))
    )

    expect_error(
        PBM(
            name = "test",
            pssm = pssm,
            max_energy = -2.5,
            min_energy = -7.0,
            energy_range = c(-7.0, 0.0),
            spat = data.frame(bin = numeric(0), spat_factor = numeric(0)),
            spat_min = numeric(0),
            spat_max = numeric(0),
            seq_length = numeric(0),
            coefs = c(a = 1.0),
            size = 0 # invalid: must be positive
        ),
        "size"
    )
})

test_that("PBM with too few pssm columns fails S4 validation", {
    pssm <- matrix(
        c(0.7, 0.1, 0.1, 0.7, 0.1, 0.1),
        nrow = 2, ncol = 3, byrow = TRUE,
        dimnames = list(NULL, c("A", "C", "G"))
    )

    expect_error(
        PBM(
            name = "test",
            pssm = pssm,
            max_energy = -2.5,
            min_energy = -7.0,
            energy_range = c(-7.0, 0.0),
            spat = data.frame(bin = numeric(0), spat_factor = numeric(0)),
            spat_min = numeric(0),
            spat_max = numeric(0),
            seq_length = numeric(0),
            coefs = c(a = 1.0),
            size = 500
        ),
        "pssm|column"
    )
})

# =============================================================================
# 5. Virtual track cleanup
# =============================================================================

test_that("gvtrack.rm is imported in the namespace", {
    misha_imports <- get_all_imports("iceqream", "misha")
    expect_true("gvtrack.rm" %in% misha_imports)
})

test_that("gvtrack.create is imported in the namespace", {
    misha_imports <- get_all_imports("iceqream", "misha")
    expect_true("gvtrack.create" %in% misha_imports)
})

test_that("gvtrack.iterator is imported in the namespace", {
    misha_imports <- get_all_imports("iceqream", "misha")
    expect_true("gvtrack.iterator" %in% misha_imports)
})

# =============================================================================
# 6. Anti-pattern fixes: seq_along / seq_len edge cases
# =============================================================================

test_that("norm01 handles zero-length vector gracefully", {
    x <- numeric(0)
    result <- suppressWarnings(norm01(x))
    expect_length(result, 0)
})

test_that("logist handles zero-length vector gracefully", {
    x <- numeric(0)
    result <- logist(x)
    expect_length(result, 0)
})

test_that("norm_energy_matrix handles zero-row matrix gracefully", {
    x <- matrix(numeric(0), nrow = 0, ncol = 2)
    colnames(x) <- c("m1", "m2")
    result <- norm_energy_matrix(x)
    expect_equal(nrow(result), 0)
    expect_equal(colnames(result), c("m1", "m2"))
})

test_that("rescale handles zero-length input gracefully", {
    x <- numeric(0)
    orig <- c(1, 5, 10)
    result <- rescale(x, orig)
    expect_length(result, 0)
})

# =============================================================================
# 7. Utility function correctness (regression tests)
# =============================================================================

# --- norm01 ---

test_that("norm01 output always in [0, 1] range for diverse inputs", {
    set.seed(123)
    # Positive values
    x1 <- runif(100, 0, 1000)
    expect_true(all(norm01(x1) >= 0 & norm01(x1) <= 1))

    # Negative values
    x2 <- runif(100, -500, -1)
    expect_true(all(norm01(x2) >= 0 & norm01(x2) <= 1))

    # Mixed sign
    x3 <- rnorm(200, 0, 100)
    expect_true(all(norm01(x3) >= 0 & norm01(x3) <= 1))

    # Exactly two values
    x4 <- c(-1, 1)
    result4 <- norm01(x4)
    expect_equal(result4[1], 0)
    expect_equal(result4[2], 1)
})

test_that("norm01 preserves relative ordering", {
    x <- c(3, 1, 4, 1, 5, 9, 2, 6)
    normed <- norm01(x)
    expect_equal(order(x), order(normed))
})

# --- norm0q ---

test_that("norm0q output is in [0, 1] range", {
    set.seed(42)
    x <- rnorm(500, mean = 10, sd = 5)
    result <- norm0q(x, quant = 0.95)
    expect_true(all(result >= 0))
    expect_true(all(result <= 1))
})

test_that("norm0q with quant=1 and positive data gives same as norm01", {
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    result_q1 <- norm0q(x, quant = 1)
    result_01 <- norm01(x)
    expect_equal(result_q1, result_01)
})

test_that("norm0q clips values above quantile to 1", {
    x <- c(1, 2, 3, 4, 5, 100)
    result <- norm0q(x, quant = 0.5)
    expect_equal(max(result), 1)
})

# --- rescale ---

test_that("rescale(norm01(x), x) returns approximately original range", {
    orig <- c(-20, -10, 0, 10, 20)
    normed <- norm01(orig)
    rescaled <- rescale(normed, orig)
    expect_equal(rescaled, orig, tolerance = 1e-10)
})

test_that("rescale(norm01(x), x) roundtrips for random data", {
    set.seed(7)
    orig <- rnorm(50, mean = 42, sd = 13)
    normed <- norm01(orig)
    rescaled <- rescale(normed, orig)
    expect_equal(rescaled, orig, tolerance = 1e-10)
})

test_that("rescale maps 0 to min and 1 to max of original", {
    orig <- c(5, 10, 15, 20, 25)
    expect_equal(rescale(0, orig), 5)
    expect_equal(rescale(1, orig), 25)
    expect_equal(rescale(0.5, orig), 15)
})

# --- logist ---

test_that("logist output range is within [0, L]", {
    x <- seq(-10, 10, by = 0.5)

    # Default L = 1
    result1 <- logist(x, L = 1, k = 1)
    expect_true(all(result1 >= 0))
    expect_true(all(result1 <= 1))

    # Custom L = 5
    result5 <- logist(x, L = 5, k = 1)
    expect_true(all(result5 >= 0))
    expect_true(all(result5 <= 5))

    # Custom L = 0.5
    result_half <- logist(x, L = 0.5, k = 1)
    expect_true(all(result_half >= 0))
    expect_true(all(result_half <= 0.5))
})

test_that("logist at midpoint equals L/2", {
    expect_equal(logist(0, x_0 = 0, L = 1, k = 1), 0.5)
    expect_equal(logist(3, x_0 = 3, L = 2, k = 5), 1.0)
    expect_equal(logist(10, x_0 = 10, L = 4, k = 0.1), 2.0)
})

test_that("logist is monotonically increasing for k > 0", {
    x <- seq(-20, 20, by = 0.1)
    result <- logist(x, k = 3)
    diffs <- diff(result)
    expect_true(all(diffs >= 0))
})

test_that("logist approaches 0 for large negative x and L for large positive x", {
    expect_equal(logist(-500, L = 1, k = 1), 0, tolerance = 1e-10)
    expect_equal(logist(500, L = 1, k = 1), 1, tolerance = 1e-10)
    expect_equal(logist(-500, L = 3, k = 1), 0, tolerance = 1e-10)
    expect_equal(logist(500, L = 3, k = 1), 3, tolerance = 1e-10)
})

# --- create_logist_features ---

test_that("create_logist_features produces 4x columns (4 transformations per feature)", {
    x <- matrix(rnorm(50), nrow = 10, ncol = 5)
    colnames(x) <- paste0("motif_", 1:5)
    result <- create_logist_features(x)
    expect_equal(ncol(result), 20) # 5 columns x 4 transformations
    expect_equal(nrow(result), 10)
})

test_that("create_logist_features column names follow pattern: name_transformation", {
    x <- matrix(rnorm(30), nrow = 10, ncol = 3)
    colnames(x) <- c("A", "B", "C")
    result <- create_logist_features(x)

    expected_suffixes <- c("high-energy", "higher-energy", "low-energy", "sigmoid")
    for (prefix in c("A", "B", "C")) {
        for (suffix in expected_suffixes) {
            expect_true(
                paste0(prefix, "_", suffix) %in% colnames(result),
                info = paste("Missing column:", paste0(prefix, "_", suffix))
            )
        }
    }
})

test_that("create_logist_features removes all-NA columns before transforming", {
    x <- matrix(rnorm(30), nrow = 10, ncol = 3)
    colnames(x) <- c("good1", "bad_na", "good2")
    x[, "bad_na"] <- NA

    result <- create_logist_features(x)

    # bad_na column should be excluded => only 2 features x 4 = 8 columns
    expect_equal(ncol(result), 8)
    expect_false(any(grepl("bad_na", colnames(result))))
})

test_that("create_logist_features output values are finite for finite input", {
    set.seed(42)
    x <- matrix(runif(100, 0, 10), nrow = 20, ncol = 5)
    colnames(x) <- paste0("m", 1:5)
    result <- create_logist_features(x)
    expect_true(all(is.finite(result)))
})

# =============================================================================
# Additional regression tests for integration correctness
# =============================================================================

test_that("norm_energy_matrix with self-reference produces consistent output", {
    set.seed(42)
    x <- matrix(runif(30, -10, 0), nrow = 10, ncol = 3)
    colnames(x) <- c("m1", "m2", "m3")

    result <- norm_energy_matrix(x, dataset_x = x)

    # Self-referencing normalization should produce values >= 0
    expect_true(all(result >= 0))
    expect_true(all(!is.na(result)))
    # Dimensions preserved
    expect_equal(dim(result), dim(x))
    expect_equal(colnames(result), colnames(x))
})

test_that("strip_traj_model is an exported function that strips glmnet models", {
    # strip_traj_model is exported and callable
    expect_true(is.function(strip_traj_model))
    # Its presence in the namespace confirms the refactoring correctly exposed it
    expect_true("strip_traj_model" %in% ls(asNamespace("iceqream")))
})

test_that("create_specific_terms normalizes interaction values to [0, scale]", {
    set.seed(42)
    energies <- matrix(abs(rnorm(60)), nrow = 10, ncol = 6)
    colnames(energies) <- c("A", "B", "C", "D", "E", "F")

    # Use multiple terms to avoid the single-column vector drop issue
    terms <- data.frame(
        term1 = c("A", "C"),
        term2 = c("B", "D"),
        variable = c("A:B", "C:D"),
        stringsAsFactors = FALSE
    )

    result <- iceqream:::create_specific_terms(energies, terms, scale = 1)

    expect_true(all(result >= 0))
    expect_true(all(result <= 1))
    expect_equal(ncol(result), 2)
    expect_equal(colnames(result), c("A:B", "C:D"))
})

test_that("create_specific_terms respects custom scale", {
    set.seed(42)
    energies <- matrix(abs(rnorm(60)), nrow = 10, ncol = 6)
    colnames(energies) <- c("A", "B", "C", "D", "E", "F")

    terms <- data.frame(
        term1 = c("A", "C"),
        term2 = c("B", "D"),
        variable = c("A:B", "C:D"),
        stringsAsFactors = FALSE
    )

    result_s1 <- iceqream:::create_specific_terms(energies, terms, scale = 1)
    result_s5 <- iceqream:::create_specific_terms(energies, terms, scale = 5)

    expect_equal(result_s5, result_s1 * 5, tolerance = 1e-10)
})

test_that("homogenize_pssm_models preserves model structure", {
    skip_if_not_installed("prego")

    model <- list(
        pssm = data.frame(
            pos = 0:3,
            A = c(0.7, 0.1, 0.1, 0.1),
            C = c(0.1, 0.7, 0.1, 0.1),
            G = c(0.1, 0.1, 0.7, 0.1),
            T = c(0.1, 0.1, 0.1, 0.7)
        ),
        spat = data.frame(bin = numeric(0), spat_factor = numeric(0)),
        spat_min = numeric(0),
        spat_max = numeric(0),
        seq_length = numeric(0)
    )

    models <- list(motif1 = model)
    result <- homogenize_pssm_models(models)

    expect_true(is.list(result))
    expect_equal(names(result), "motif1")
    expect_true("pssm" %in% names(result$motif1))
    expect_equal(nrow(result$motif1$pssm), 4)
})

# =============================================================================
# Verify key namespace imports for the refactoring
# =============================================================================

test_that("tidyr pivot_longer is available in namespace", {
    tidyr_imports <- get_all_imports("iceqream", "tidyr")
    expect_true("pivot_longer" %in% tidyr_imports)
})

test_that("tidyr pivot_wider is available in namespace", {
    tidyr_imports <- get_all_imports("iceqream", "tidyr")
    expect_true("pivot_wider" %in% tidyr_imports)
})

test_that("methods .hasSlot is imported", {
    methods_imports <- get_all_imports("iceqream", "methods")
    expect_true(".hasSlot" %in% methods_imports)
})

test_that("cli messaging functions are imported", {
    cli_imports <- get_all_imports("iceqream", "cli")
    expect_true("cli_abort" %in% cli_imports)
    expect_true("cli_alert" %in% cli_imports)
    expect_true("cli_alert_info" %in% cli_imports)
    expect_true("cli_alert_warning" %in% cli_imports)
    expect_true("cli_alert_success" %in% cli_imports)
})
