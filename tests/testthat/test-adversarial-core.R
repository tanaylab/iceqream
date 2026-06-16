# Adversarial tests for core numeric helpers and the partial-response path,
# focused on degenerate inputs and the interaction path (which the B5 report
# optimization, compute_partial_response(vars=), must handle correctly).

# ---------------------------------------------------------------------------
# Core numeric utilities: norm01 / norm0q / rescale / logist
# ---------------------------------------------------------------------------

test_that("norm01 maps a constant vector to all zeros (no divide-by-zero)", {
    expect_equal(norm01(c(5, 5, 5)), c(0, 0, 0))
    expect_equal(norm01(rep(-3.2, 4)), rep(0, 4))
})

test_that("norm01 propagates NA but normalizes the rest", {
    expect_equal(norm01(c(NA, 1, 3)), c(NA, 0, 1))
    expect_equal(norm01(c(-2, -1, 0)), c(0, 0.5, 1))
})

test_that("norm01 of a single element is 0", {
    expect_equal(norm01(5), 0)
})

test_that("norm0q output stays within [0, 1] and clips above the quantile", {
    out <- norm0q(c(0, 1, 2, 100), quant = 0.5)
    expect_true(all(out >= 0 & out <= 1))
    expect_equal(max(out), 1)
})

test_that("rescale handles a constant orig_x without NaN/Inf", {
    expect_equal(rescale(c(0, 0.5, 1), c(7, 7, 7)), c(7, 7, 7))
})

test_that("rescale maps [0,1] onto [min(orig_x), max(orig_x)]", {
    expect_equal(rescale(c(0, 0.5, 1), c(2, 4)), c(2, 3, 4))
})

test_that("logist does not overflow to NaN on extreme inputs", {
    out <- logist(c(-1000, 0, 1000))
    expect_true(all(is.finite(out)))
    expect_equal(out, c(0, 0.5, 1))
})

# ---------------------------------------------------------------------------
# create_logist_features
# ---------------------------------------------------------------------------

test_that("create_logist_features produces 4 named transforms per input column", {
    m <- matrix(c(0, 5, 10), ncol = 1, dimnames = list(NULL, "m1"))
    out <- create_logist_features(m)
    expect_equal(ncol(out), 4)
    expect_setequal(
        colnames(out),
        c("m1_low-energy", "m1_high-energy", "m1_higher-energy", "m1_sigmoid")
    )
    expect_true(all(is.finite(out)))
})

test_that("create_logist_features coerces a bare vector to a single named column", {
    out <- create_logist_features(c(1, 2, 3))
    expect_equal(ncol(out), 4)
    expect_true(all(grepl("^V1_", colnames(out))))
    expect_equal(nrow(out), 3)
})

test_that("create_logist_features drops all-NA columns but keeps partial-NA ones", {
    m <- matrix(c(NA, NA, NA, 1, 2, 3), ncol = 2, dimnames = list(NULL, c("a", "b")))
    out <- create_logist_features(m)
    expect_equal(ncol(out), 4) # only "b" survives -> 4 transforms
    expect_true(all(grepl("^b_", colnames(out))))
})

test_that("create_logist_features propagates NA entries without erroring", {
    m <- matrix(c(1, NA, 3), ncol = 1, dimnames = list(NULL, "m1"))
    out <- create_logist_features(m)
    expect_true(any(is.na(out)))
    expect_equal(nrow(out), 3)
})

# ---------------------------------------------------------------------------
# norm_energy_matrix degenerate inputs (no NaN / Inf escapes)
# ---------------------------------------------------------------------------

test_that("norm_energy_matrix on a single-column input stays finite in [0, 10]", {
    set.seed(101)
    E <- matrix(runif(20, -10, 0), ncol = 1, dimnames = list(NULL, "m"))
    out <- norm_energy_matrix(E, E)
    expect_true(all(is.finite(out)))
    expect_true(all(out >= 0 & out <= 10))
})

test_that("norm_energy_matrix with a degenerate q=0 reference yields finite output (no NaN)", {
    set.seed(102)
    E <- matrix(runif(30, -10, 0), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
    out <- norm_energy_matrix(E, E, q = 0)
    expect_false(any(is.na(out)))
    expect_true(all(is.finite(out)))
})

test_that("norm_energy_matrix on a single row does not produce NaN", {
    E <- matrix(c(-1, -3), nrow = 1, dimnames = list(NULL, c("a", "b")))
    out <- norm_energy_matrix(E, E)
    expect_false(any(is.na(out)))
})

# ---------------------------------------------------------------------------
# compute_partial_response on a model WITH interactions (the complex path the
# B5 vars= optimization has to handle: variable names contain ":").
# ---------------------------------------------------------------------------

build_interaction_model <- function(logist_interactions = FALSE) {
    tm <- create_interaction_traj_model(n_peaks = 200, n_motifs = 8, seed = 5)
    suppressWarnings(suppressMessages(
        add_interactions(
            tm,
            interaction_threshold = 0.001,
            logist_interactions = logist_interactions,
            seed = 5
        )
    ))
}

test_that("compute_partial_response(vars=) is value-identical on interaction variables (logist_interactions = FALSE)", {
    tm <- build_interaction_model(logist_interactions = FALSE)
    skip_if(ncol(tm@interactions) < 2, "fixture produced too few interactions")

    full <- iceqream:::compute_partial_response(tm)
    inter_vars <- colnames(tm@interactions)[1:2]
    motif_var <- colnames(tm@normalized_energies)[1]
    target <- c(motif_var, inter_vars)

    sub <- iceqream:::compute_partial_response(tm, vars = target)
    expect_true(all(target %in% colnames(sub)))
    expect_equal(sub[, colnames(sub), drop = FALSE], full[, colnames(sub), drop = FALSE])
})

test_that("compute_partial_response(vars=) is value-identical on interaction variables (logist_interactions = TRUE)", {
    tm <- build_interaction_model(logist_interactions = TRUE)
    skip_if(ncol(tm@interactions) < 2, "fixture produced too few interactions")

    full <- iceqream:::compute_partial_response(tm)
    inter_vars <- colnames(tm@interactions)[1:2]

    sub <- iceqream:::compute_partial_response(tm, vars = inter_vars)
    expect_true(all(inter_vars %in% colnames(sub)))
    expect_equal(sub[, colnames(sub), drop = FALSE], full[, colnames(sub), drop = FALSE])
})

test_that("compute_partial_response(reverse=TRUE) with vars= matches the full reverse computation", {
    tm <- build_interaction_model(logist_interactions = FALSE)
    skip_if(ncol(tm@interactions) < 1, "fixture produced no interactions")

    full_rev <- iceqream:::compute_partial_response(tm, reverse = TRUE)
    target <- c(colnames(tm@normalized_energies)[1], colnames(tm@interactions)[1])

    sub_rev <- iceqream:::compute_partial_response(tm, vars = target, reverse = TRUE)
    expect_true(all(target %in% colnames(sub_rev)))
    expect_equal(sub_rev[, colnames(sub_rev), drop = FALSE], full_rev[, colnames(sub_rev), drop = FALSE])
})

test_that("compute_partial_response is order- and duplicate-insensitive in vars=", {
    tm <- create_mock_traj_model(n_peaks = 60, n_motifs = 4)
    a <- iceqream:::compute_partial_response(tm, vars = c("motif2", "motif1"))
    b <- iceqream:::compute_partial_response(tm, vars = c("motif1", "motif2", "motif1"))
    # same set of columns, same values regardless of request order / duplicates
    expect_setequal(colnames(a), c("motif1", "motif2"))
    expect_setequal(colnames(b), c("motif1", "motif2"))
    expect_equal(a[, c("motif1", "motif2")], b[, c("motif1", "motif2")])
})
