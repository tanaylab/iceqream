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

# --- filter_traj_model_by_beta (the by-beta variant used by the manifold
#     filter): same robustness guards as filter_traj_model. ---

test_that("filter_traj_model_by_beta keeps a fittable model under all thresholds", {
    tm <- create_mock_traj_model(n_peaks = 120, n_motifs = 4)
    # low threshold -> nothing removed
    keep <- suppressWarnings(suppressMessages(iceqream:::filter_traj_model_by_beta(tm, threshold = 1e-9)))
    expect_equal(length(keep@motif_models), 4)
    # absurd threshold would remove everything -> kept >= 1 (no glmnet crash)
    rm_all <- suppressWarnings(suppressMessages(iceqream:::filter_traj_model_by_beta(tm, threshold = 1e9)))
    expect_gte(length(rm_all@motif_models), 1)
})

test_that("filter_traj_model_by_beta is a no-op on a single-motif model", {
    tm <- create_mock_traj_model(n_peaks = 80, n_motifs = 1)
    out <- suppressWarnings(suppressMessages(iceqream:::filter_traj_model_by_beta(tm, threshold = 1e9)))
    expect_equal(length(out@motif_models), 1)
})

test_that("create_logist_features returns an empty matrix on a 0-column input", {
    # Regression test: a model with no additional features (a 0-column
    # data.frame) used to crash create_logist_features at colnames<-. It should
    # now return an empty matrix that still carries the rows for downstream cbind.
    empty <- data.frame(row.names = paste0("peak", 1:5))
    out <- iceqream:::create_logist_features(empty)
    expect_true(is.matrix(out))
    expect_equal(nrow(out), 5L)
    expect_equal(ncol(out), 0L)

    # An all-NA single column collapses to 0 columns too -> same empty result.
    out2 <- iceqream:::create_logist_features(matrix(NA_real_, nrow = 5, ncol = 1))
    expect_equal(ncol(out2), 0L)
})

test_that("relearn_traj_model(new_logist=TRUE) works on a model with no additional features", {
    # The multi by-beta unify block relearns each trajectory from scratch; with
    # no additional features this previously crashed inside create_logist_features.
    tm <- create_mock_traj_model(n_peaks = 120, n_motifs = 4)
    expect_equal(ncol(tm@additional_features), 0L)
    out <- suppressWarnings(suppressMessages(
        iceqream:::relearn_traj_model(tm, new_logist = TRUE)
    ))
    expect_s4_class(out, "TrajectoryModel")
    expect_equal(length(out@motif_models), 4L)
    expect_gte(ncol(out@model_features), 2L)
})

test_that("filter_multi_traj_model_by_beta unifies without crashing on no-additional-feature models", {
    # Regression test for the unify block: each trajectory is relearned with
    # new_logist=TRUE over the unified motif set. Models with no additional
    # features (the common case) used to crash there.
    m1 <- iceqream:::add_traj_model_stats(create_mock_traj_model(n_peaks = 150, n_motifs = 5))
    m2 <- iceqream:::add_traj_model_stats(create_mock_traj_model(n_peaks = 150, n_motifs = 5))
    # Give the two trajectories distinct motif names so the unified set is wider
    # than either model (exercises the cross-model normalized_energies cbind).
    names(m1@motif_models) <- paste0("a", seq_along(m1@motif_models))
    colnames(m1@normalized_energies) <- names(m1@motif_models)
    names(m2@motif_models) <- paste0("b", seq_along(m2@motif_models))
    colnames(m2@normalized_energies) <- names(m2@motif_models)
    mm <- c(m1@motif_models, m2@motif_models)
    multi <- TrajectoryModelMulti(
        models = list(t1 = m1, t2 = m2), models_full = list(t1 = m1, t2 = m2),
        motif_models = mm, cluster_map = data.frame(), stats = data.frame(), params = list()
    )
    out <- suppressWarnings(suppressMessages(
        iceqream:::filter_multi_traj_model_by_beta(multi, beta_threshold = 1e-9)
    ))
    expect_s4_class(out, "TrajectoryModelMulti")
    expect_length(out@models_full, 2)
    # unify makes each full model carry the union of motifs (10 unique here)
    expect_true(all(purrr::map_lgl(out@models_full, ~ length(.x@motif_models) >= 1)))
})

test_that("filter_multi_traj_model maps the (fixed) single-trajectory filter over a multi model", {
    m1 <- iceqream:::add_traj_model_stats(create_interaction_traj_model(n_peaks = 150, n_motifs = 5, seed = 1))
    m2 <- iceqream:::add_traj_model_stats(create_interaction_traj_model(n_peaks = 150, n_motifs = 5, seed = 2))
    mm <- c(m1@motif_models, m2@motif_models)
    mm <- mm[unique(names(mm))]
    multi <- TrajectoryModelMulti(
        models = list(t1 = m1, t2 = m2), models_full = list(t1 = m1, t2 = m2),
        motif_models = mm, cluster_map = data.frame(), stats = data.frame(), params = list()
    )
    out <- suppressWarnings(suppressMessages(
        filter_multi_traj_model(multi, r2_threshold = 1e-6, bits_threshold = 0, sample_frac = NULL)
    ))
    expect_s4_class(out, "TrajectoryModelMulti")
    expect_length(out@models_full, 2)
    expect_true(all(purrr::map_lgl(out@models_full, ~ length(.x@motif_models) >= 1)))
})
