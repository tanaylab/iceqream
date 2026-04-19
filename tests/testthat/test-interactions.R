# Tests for motif-motif interaction machinery in R/interactions.R
# Runs without misha.

# Helper: build a traj_model with pre-populated interaction columns so we can
# exercise remove_interactions / add_interactions branches without wiring the
# full glmnet-based selection path.
build_traj_model_with_logist_interactions <- function() {
    tm <- create_mock_traj_model(n_peaks = 50, n_motifs = 3)

    n <- nrow(tm@model_features)
    inter_mat <- matrix(
        runif(n * 2, 0, 1),
        nrow = n,
        ncol = 2,
        dimnames = list(rownames(tm@model_features), c("motif1:motif2", "motif1:motif3"))
    )
    tm@interactions <- inter_mat

    logist_inter <- create_logist_features(inter_mat)
    tm@model_features <- cbind(tm@model_features, logist_inter)
    tm@params$logist_interactions <- TRUE

    tm
}

test_that("remove_interactions strips logist-expanded columns when logist_interactions = TRUE", {
    tm <- build_traj_model_with_logist_interactions()

    cols_before <- colnames(tm@model_features)
    inter_cols <- colnames(tm@interactions)
    logist_inter_cols <- as.vector(outer(
        inter_cols,
        c("low-energy", "high-energy", "higher-energy", "sigmoid"),
        function(a, b) paste0(a, "_", b)
    ))

    expect_true(all(logist_inter_cols %in% cols_before))
    motif_only_cols <- setdiff(cols_before, logist_inter_cols)

    tm_clean <- iceqream:::remove_interactions(tm)

    expect_false(any(logist_inter_cols %in% colnames(tm_clean@model_features)))
    expect_setequal(colnames(tm_clean@model_features), motif_only_cols)
    expect_equal(ncol(tm_clean@interactions), 0)
    expect_false(isTRUE(tm_clean@params$logist_interactions))
})

test_that("remove_interactions strips raw product columns when logist_interactions = FALSE", {
    tm <- create_mock_traj_model(n_peaks = 50, n_motifs = 3)
    n <- nrow(tm@model_features)
    inter_mat <- matrix(
        runif(n * 2, 0, 1),
        nrow = n,
        ncol = 2,
        dimnames = list(rownames(tm@model_features), c("motif1:motif2", "motif1:motif3"))
    )
    tm@interactions <- inter_mat
    tm@model_features <- cbind(tm@model_features, inter_mat)
    tm@params$logist_interactions <- FALSE

    cols_before <- colnames(tm@model_features)
    inter_cols <- colnames(inter_mat)
    expect_true(all(inter_cols %in% cols_before))

    tm_clean <- iceqream:::remove_interactions(tm)

    expect_false(any(inter_cols %in% colnames(tm_clean@model_features)))
    expect_equal(ncol(tm_clean@interactions), 0)
})

test_that("remove_interactions is a no-op when there are no interactions", {
    tm <- create_mock_traj_model(n_peaks = 30, n_motifs = 2)
    cols_before <- colnames(tm@model_features)

    tm_clean <- iceqream:::remove_interactions(tm)

    expect_equal(colnames(tm_clean@model_features), cols_before)
})

test_that("relearn_traj_model use_cv path reuses cv.glmnet's full-path fit", {
    set.seed(1)
    n <- 80
    p <- 6
    X <- matrix(rnorm(n * p), n, p)
    y <- plogis(X[, 1] - X[, 2] + rnorm(n, sd = 0.3))

    cv_model <- glmnet::cv.glmnet(
        X, y,
        family = binomial(link = "logit"),
        alpha = 0.5, nfolds = 5, seed = 1
    )
    lambda <- cv_model$lambda.min

    old_path_model <- glmnet::glmnet(
        X, y,
        family = binomial(link = "logit"),
        alpha = 0.5, lambda = lambda, seed = 1
    )
    new_path_model <- cv_model$glmnet.fit

    pred_old <- glmnet::predict.glmnet(old_path_model, newx = X, type = "link", s = lambda)[, 1]
    pred_new <- glmnet::predict.glmnet(new_path_model, newx = X, type = "link", s = lambda)[, 1]

    expect_lt(max(abs(pred_old - pred_new)), 1e-3)
})
