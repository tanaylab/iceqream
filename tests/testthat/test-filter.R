# Adversarial tests for filter_traj_model() and its helpers (all genome-free:
# they operate on the in-memory model, not on extracted sequences). These cover
# edge cases that previously crashed or silently mis-behaved.

test_that("filter_traj_model returns a valid model on a normal multi-motif model", {
    tm <- create_interaction_traj_model(n_peaks = 200, n_motifs = 10, seed = 1)
    out <- suppressWarnings(suppressMessages(
        filter_traj_model(tm, r2_threshold = 1e-6, bits_threshold = 0, sample_frac = NULL)
    ))
    expect_s4_class(out, "TrajectoryModel")
    expect_lte(length(out@motif_models), length(tm@motif_models))
    expect_gte(length(out@motif_models), 1)
})

test_that("filter_traj_model never removes every motif (keeps at least one)", {
    tm <- create_interaction_traj_model(n_peaks = 200, n_motifs = 10, seed = 1)
    # absurd thresholds select every motif for removal
    out <- suppressWarnings(suppressMessages(
        filter_traj_model(tm, r2_threshold = 1e9, bits_threshold = 1e9, sample_frac = NULL)
    ))
    expect_s4_class(out, "TrajectoryModel")
    expect_gte(length(out@motif_models), 1)
    # the kept model is still fittable (>= 2 feature columns)
    expect_gte(ncol(out@model_features), 2)
})

test_that("filter_traj_model is a no-op on a single-motif model (no crash)", {
    tm <- create_mock_traj_model(n_peaks = 80, n_motifs = 1)
    out <- suppressWarnings(suppressMessages(
        filter_traj_model(tm, r2_threshold = 1e-6, bits_threshold = 0, sample_frac = NULL)
    ))
    expect_s4_class(out, "TrajectoryModel")
    expect_equal(length(out@motif_models), 1)
    # even with absurd thresholds it cannot remove the only motif
    out2 <- suppressWarnings(suppressMessages(
        filter_traj_model(tm, r2_threshold = 1e9, bits_threshold = 1e9, sample_frac = NULL)
    ))
    expect_equal(length(out2@motif_models), 1)
})

test_that("filter_traj_model reduces a two-motif model to one without crashing", {
    tm <- create_mock_traj_model(n_peaks = 80, n_motifs = 2)
    out <- suppressWarnings(suppressMessages(
        filter_traj_model(tm, r2_threshold = 1e9, bits_threshold = 1e9, sample_frac = NULL)
    ))
    expect_gte(length(out@motif_models), 1)
})

test_that("filter_traj_model applies the bits filter even when no feature fails the R^2 threshold", {
    # Regression test: previously the bits-based removal was silently discarded
    # unless some feature also failed the R^2 threshold.
    tm <- create_mock_traj_model(n_peaks = 150, n_motifs = 3)
    # Make motif1 information-poor (uniform PSSM -> ~0 bits); motif2/3 stay informative.
    tm@motif_models[["motif1"]]$pssm[, c("A", "C", "G", "T")] <- 0.25
    out <- suppressWarnings(suppressMessages(
        filter_traj_model(tm, r2_threshold = 1e-9, bits_threshold = 1, sample_frac = NULL)
    ))
    expect_false("motif1" %in% names(out@motif_models))
    expect_true(all(c("motif2", "motif3") %in% names(out@motif_models)))
})

test_that("add_features_r2 handles few-motif models without crashing", {
    # 2 motifs: leave-one-out leaves 1 motif (4 logist features) -> fittable.
    tm2 <- create_mock_traj_model(n_peaks = 80, n_motifs = 2)
    r2 <- suppressWarnings(suppressMessages(add_features_r2(tm2, sample_frac = NULL)))
    expect_length(r2@features_r2, 2)
    expect_true(all(is.finite(r2@features_r2)))

    # 1 motif: leave-one-out would leave 0 features -> guarded to R^2-added of 0.
    tm1 <- create_mock_traj_model(n_peaks = 80, n_motifs = 1)
    r2_1 <- suppressWarnings(suppressMessages(add_features_r2(tm1, sample_frac = NULL)))
    expect_length(r2_1@features_r2, 1)
    expect_true(is.finite(r2_1@features_r2))
})

test_that("filter_traj_model_intervals preserves matrix slots for a single-motif model", {
    tm <- create_mock_traj_model(n_peaks = 80, n_motifs = 1)
    sub <- iceqream:::filter_traj_model_intervals(tm, c(rep(TRUE, 40), rep(FALSE, 40)))
    expect_true(is.matrix(sub@normalized_energies))
    expect_equal(ncol(sub@normalized_energies), 1L)
    expect_equal(nrow(sub@normalized_energies), 40L)
    expect_true(is.matrix(sub@model_features))
})
