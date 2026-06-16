# Adversarial tests for the 2026-06 API/efficiency refactor.
# These deliberately probe edge cases and use differential testing (new code
# vs. a faithful reimplementation of the pre-refactor code) to prove the
# vectorizations / extractions did not change behavior.

# ---------------------------------------------------------------------------
# B4: create_features_terms vectorization (matrixStats + sweep) must be
# numerically identical to the old map_dfc -> apply(,2,norm01) implementation.
# ---------------------------------------------------------------------------

# Faithful copy of the pre-refactor implementation.
old_create_features_terms <- function(energies, features, data, scale = 1) {
    interactions <- purrr::map_dfc(features, function(.x) {
        inter <- energies[, setdiff(colnames(energies), .x), drop = FALSE] * data[, .x]
        max_vals <- apply(inter, 2, max, na.rm = TRUE)
        max_vals[max_vals == 0] <- 1
        inter <- t(t(inter) / max_vals)
        colnames(inter) <- paste0(.x, ":", colnames(inter))
        inter
    })
    interactions <- apply(interactions, 2, norm01) * scale
    interactions
}

make_energies <- function(n, p, seed) {
    set.seed(seed)
    E <- matrix(runif(n * p, 0, 10), nrow = n, ncol = p)
    colnames(E) <- paste0("m", seq_len(p))
    rownames(E) <- paste0("peak", seq_len(n))
    E
}

# Compare values + column names but not row names: the new (vectorized)
# create_features_terms preserves input row names while the old apply(,2,norm01)
# implementation dropped them. This is benign because the production caller
# create_interaction_terms() overwrites rownames(interactions) <-
# rownames(energies) immediately afterward (interactions.R:93-94), so the final
# rownames are identical either way. We assert the numbers and column layout match.
expect_same_values <- function(a, b, ...) {
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b, ...)
}

test_that("create_features_terms matches the old implementation on random energies", {
    E <- make_energies(150, 6, 11)
    feats <- colnames(E)[1:3]
    new <- suppressWarnings(iceqream:::create_features_terms(E, feats, E))
    old <- suppressWarnings(old_create_features_terms(E, feats, E))
    expect_same_values(new, old)
})

test_that("create_features_terms matches old impl with a non-default scale", {
    E <- make_energies(120, 5, 12)
    feats <- colnames(E)[c(1, 4)]
    new <- suppressWarnings(iceqream:::create_features_terms(E, feats, E, scale = 2.5))
    old <- suppressWarnings(old_create_features_terms(E, feats, E, scale = 2.5))
    expect_same_values(new, old)
})

test_that("create_features_terms matches old impl with NAs, an all-zero column, and a constant column", {
    E <- make_energies(100, 5, 13)
    E[c(3, 17, 42), 2] <- NA   # scattered NAs
    E[, 3] <- 0                # all-zero motif column (product -> 0 -> norm01 -> 0)
    E[, 4] <- 5                # constant motif column (zero-range product after scaling)
    feats <- colnames(E)[c(1, 2, 5)]
    new <- suppressWarnings(iceqream:::create_features_terms(E, feats, E))
    old <- suppressWarnings(old_create_features_terms(E, feats, E))
    expect_same_values(new, old)
})

test_that("create_features_terms matches old impl with negative and extreme-magnitude values", {
    set.seed(14)
    E <- matrix(rnorm(80 * 4, sd = 1e5), nrow = 80, ncol = 4)
    E[, 1] <- E[, 1] * 1e-6
    colnames(E) <- paste0("m", 1:4)
    rownames(E) <- paste0("peak", 1:80)
    feats <- colnames(E)[1:2]
    new <- suppressWarnings(iceqream:::create_features_terms(E, feats, E))
    old <- suppressWarnings(old_create_features_terms(E, feats, E))
    expect_same_values(new, old)
})

test_that("create_features_terms matches old impl when energies arrive as a data.frame", {
    E <- make_energies(90, 4, 15)
    df <- as.data.frame(E)
    feats <- colnames(E)[1:2]
    new <- suppressWarnings(iceqream:::create_features_terms(df, feats, df))
    old <- suppressWarnings(old_create_features_terms(df, feats, df))
    expect_same_values(new, old)
})

test_that("create_features_terms matches old impl for a single anchor feature", {
    E <- make_energies(70, 4, 16)
    feats <- colnames(E)[1]
    new <- suppressWarnings(iceqream:::create_features_terms(E, feats, E))
    old <- suppressWarnings(old_create_features_terms(E, feats, E))
    expect_same_values(new, old)
    # one anchor -> (p-1) interaction columns
    expect_equal(ncol(new), ncol(E) - 1)
})

test_that("create_features_terms output columns are bounded in [0, scale]", {
    E <- make_energies(100, 5, 17)
    feats <- colnames(E)[1:3]
    out <- suppressWarnings(iceqream:::create_features_terms(E, feats, E, scale = 3))
    finite <- out[is.finite(out)]
    expect_true(all(finite >= 0))
    expect_true(all(finite <= 3 + 1e-9))
})

# ---------------------------------------------------------------------------
# B4: min_signal_correlation filter behaves monotonically and never crashes.
# ---------------------------------------------------------------------------

test_that("add_interactions min_signal_correlation is monotone and never drops on the normal path", {
    tm <- create_interaction_traj_model(n_peaks = 200, n_motifs = 10, seed = 5)

    keep_count <- function(msc) {
        m <- suppressWarnings(suppressMessages(
            add_interactions(tm, interaction_threshold = 0.001, min_signal_correlation = msc, seed = 5)
        ))
        ncol(m@interactions)
    }

    n_none <- keep_count(NULL)   # no post-filter
    n_loose <- keep_count(1 / 8) # Akhiad-style loose filter
    n_tight <- keep_count(0.9)   # very aggressive

    expect_gt(n_none, 0)
    # tighter threshold keeps no more than a looser one
    expect_lte(n_loose, n_none)
    expect_lte(n_tight, n_loose)
})

# ---------------------------------------------------------------------------
# C3: fit_and_predict_model equivalence across alpha / lambda, and with a
# distinct predict set (train != predict).
# ---------------------------------------------------------------------------

manual_fit_predict <- function(Xtr, ytr, Xpr, diff, alpha, lambda, seed) {
    m <- iceqream:::strip_glmnet(suppressWarnings(glmnet::glmnet(
        Xtr, ytr, binomial(link = "logit"), alpha = alpha, lambda = lambda, seed = seed
    )))
    p <- iceqream:::logist(glmnet::predict.glmnet(m, newx = Xpr, type = "link", s = lambda))[, 1]
    list(model = m, pred = iceqream:::rescale(norm01(p), diff))
}

test_that("fit_and_predict_model equals the manual block across alpha in {0, 0.5, 1}", {
    set.seed(21)
    n <- 150
    X <- matrix(rnorm(n * 5), n, 5)
    colnames(X) <- paste0("f", 1:5)
    diff <- rnorm(n)
    y <- norm01(diff)
    for (a in c(0, 0.5, 1)) {
        res <- suppressWarnings(iceqream:::fit_and_predict_model(X, y, X, diff, alpha = a, lambda = 1e-4, seed = 21))
        man <- manual_fit_predict(X, y, X, diff, alpha = a, lambda = 1e-4, seed = 21)
        expect_equal(res$predicted_diff_score, man$pred, info = paste("alpha =", a))
        expect_equal(as.numeric(res$model$beta[, 1]), as.numeric(man$model$beta[, 1]), info = paste("alpha =", a))
    }
})

test_that("fit_and_predict_model handles train != predict features (split workflow)", {
    set.seed(22)
    n <- 200
    X <- matrix(rnorm(n * 4), n, 4)
    colnames(X) <- paste0("f", 1:4)
    diff <- rnorm(n)
    y <- norm01(diff)
    tr <- 1:140
    res <- suppressWarnings(iceqream:::fit_and_predict_model(
        X[tr, ], y[tr], X, diff, alpha = 1, lambda = 1e-4, seed = 22
    ))
    man <- manual_fit_predict(X[tr, ], y[tr], X, diff, alpha = 1, lambda = 1e-4, seed = 22)
    expect_length(res$predicted_diff_score, n)
    expect_equal(res$predicted_diff_score, man$pred)
})

# ---------------------------------------------------------------------------
# C3: rebuild_traj_model adversarial carryover.
# ---------------------------------------------------------------------------

test_that("rebuild_traj_model carries non-empty additional_features verbatim", {
    tm <- create_mock_traj_model(n_peaks = 40, n_motifs = 3)
    af <- data.frame(gc = runif(40), cpg = runif(40), row.names = rownames(tm@normalized_energies))
    tm@additional_features <- af

    ne <- tm@normalized_energies[, 1:2, drop = FALSE]
    mf <- create_logist_features(ne)
    fit <- suppressWarnings(iceqream:::fit_and_predict_model(
        mf, norm01(tm@diff_score), mf, tm@diff_score,
        alpha = tm@params$alpha, lambda = tm@params$lambda, seed = 1
    ))
    rebuilt <- iceqream:::rebuild_traj_model(
        tm, model = fit$model, motif_models = tm@motif_models[1:2],
        normalized_energies = ne, model_features = mf,
        predicted_diff_score = fit$predicted_diff_score
    )
    expect_identical(rebuilt@additional_features, af)
})

test_that("rebuild_traj_model coerces a data.frame energies input to a matrix and honors a params override", {
    tm <- create_mock_traj_model(n_peaks = 30, n_motifs = 3)
    ne_df <- as.data.frame(tm@normalized_energies[, 1:2, drop = FALSE])
    mf <- create_logist_features(as.matrix(ne_df))
    fit <- suppressWarnings(iceqream:::fit_and_predict_model(
        mf, norm01(tm@diff_score), mf, tm@diff_score,
        alpha = tm@params$alpha, lambda = tm@params$lambda, seed = 1
    ))
    new_params <- tm@params
    new_params$distilled_features <- c("a", "b")
    rebuilt <- iceqream:::rebuild_traj_model(
        tm, model = fit$model, motif_models = tm@motif_models[1:2],
        normalized_energies = ne_df, model_features = mf,
        predicted_diff_score = fit$predicted_diff_score, params = new_params
    )
    expect_true(is.matrix(rebuilt@normalized_energies))
    expect_identical(rebuilt@params$distilled_features, c("a", "b"))
})

# ---------------------------------------------------------------------------
# B5: compute_partial_response(vars=) is value-identical to subsetting the
# full computation (the report optimization must not change numbers).
# ---------------------------------------------------------------------------

test_that("compute_partial_response(vars=) equals the full computation restricted to those vars", {
    tm <- create_mock_traj_model(n_peaks = 60, n_motifs = 4)
    full <- iceqream:::compute_partial_response(tm)
    target <- colnames(tm@normalized_energies)[c(1, 3)]
    sub <- iceqream:::compute_partial_response(tm, vars = target)
    expect_setequal(colnames(sub), target)
    expect_equal(sub[, colnames(sub), drop = FALSE], full[, colnames(sub), drop = FALSE])
})

test_that("compute_partial_response with vars not in the model returns an empty frame, not an error", {
    tm <- create_mock_traj_model(n_peaks = 30, n_motifs = 2)
    res <- iceqream:::compute_partial_response(tm, vars = "does_not_exist")
    expect_s3_class(res, "data.frame")
    expect_equal(ncol(res), 0)
})

# ---------------------------------------------------------------------------
# C1: normalization edge cases.
# ---------------------------------------------------------------------------

test_that("norm_energy_matrix maps a zero-range (constant) column to 0 without NaN/Inf", {
    set.seed(31)
    E <- matrix(runif(50 * 3, -10, 0), 50, 3)
    E[, 2] <- -1 # constant column -> zero observed range
    colnames(E) <- c("a", "b", "c")
    out <- norm_energy_matrix(E, E)
    expect_true(all(is.finite(out)))
    expect_true(all(out[, "b"] == 0))
})

test_that("norm_energy_matrix tolerates NAs in the input", {
    set.seed(32)
    E <- matrix(runif(60 * 3, -10, 0), 60, 3)
    E[c(1, 5, 9), 1] <- NA
    colnames(E) <- c("a", "b", "c")
    out <- norm_energy_matrix(E, E)
    # NA fill at the end of norm_energy_matrix turns NA into 0
    expect_false(any(is.na(out)))
    expect_true(all(out >= 0 & out <= 10))
})

test_that("normalize_with_db_quantiles is invariant to extra motifs/quantiles in db_quantiles", {
    set.seed(33)
    me <- matrix(runif(20 * 2, -8, 0), 20, 2)
    colnames(me) <- c("m1", "m2")
    dbq <- matrix(c(-1, -2, -3, -1.5, -2.5, -3.5), nrow = 3, byrow = FALSE,
                  dimnames = list(c("m1", "m2", "m3_extra"), c("0.99", "1")))
    out_full <- iceqream:::normalize_with_db_quantiles(me, dbq, 1, -7, 10)
    # adding an unused motif row to db_quantiles must not change the result
    dbq2 <- rbind(dbq, m4_extra = c(-9, -9))
    out_extra <- iceqream:::normalize_with_db_quantiles(me, dbq2, 1, -7, 10)
    expect_equal(out_full, out_extra)
    expect_true(all(out_full >= 0 & out_full <= 10))
})
