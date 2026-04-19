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

# ------------------------------------------------------------------------------
# C0 baseline: single-pass add_interactions on the seeded synthetic fixture.
# These snapshot values are the parity target for the forthcoming progressive
# redesign. Any task C1/C2/C3/C4 commit that touches interaction selection
# MUST still pass add_interactions(..., interaction_threshold = 0.001,
# only_sig_motifs = FALSE, only_sig_add_motifs = TRUE) with the snapshot.
# ------------------------------------------------------------------------------

test_that("baseline: add_interactions(thr=0.001) on fixture produces snapshot", {
    skip_on_cran()
    tm <- create_interaction_traj_model(n_peaks = 200, n_motifs = 10, seed = 1)

    r2_train_before <- cor(
        tm@diff_score[tm@type == "train"],
        tm@predicted_diff_score[tm@type == "train"]
    )^2

    tm_with <- suppressMessages(suppressWarnings(
        add_interactions(
            tm,
            interaction_threshold = 0.001,
            only_sig_motifs = FALSE,
            only_sig_add_motifs = TRUE,
            seed = 1
        )
    ))

    r2_train_after <- cor(
        tm_with@diff_score[tm_with@type == "train"],
        tm_with@predicted_diff_score[tm_with@type == "train"]
    )^2

    # Snapshot values recorded on branch deep-code-review-remediation @ 0aa612c.
    # Do NOT change these without a conscious decision + NEWS entry.
    baseline <- list(
        n_interactions = ncol(tm_with@interactions),
        top_interactions = sort(colnames(tm_with@interactions))[seq_len(min(5, ncol(tm_with@interactions)))],
        r2_train_after = r2_train_after
    )

    # Soft checks: fixture must produce SOME interactions and R2 must improve.
    expect_gt(baseline$n_interactions, 0)
    expect_gte(baseline$r2_train_after, r2_train_before - 1e-6)

    # Hard snapshot: the seeded fixture must produce identical output on every
    # run of the current implementation.
    expect_snapshot(baseline)
})

test_that("progressive with length-1 thresholds matches single add_interactions", {
    tm <- create_interaction_traj_model(n_peaks = 100, n_motifs = 6, seed = 4)

    tm_single <- suppressMessages(suppressWarnings(
        add_interactions(tm, interaction_threshold = 0.001,
            only_sig_motifs = FALSE, only_sig_add_motifs = TRUE, seed = 4)
    ))
    tm_progressive <- suppressMessages(suppressWarnings(
        add_interactions_progressive(
            tm,
            thresholds = 0.001,
            only_sig_motifs = FALSE,
            only_sig_add_motifs = TRUE,
            seed = 4
        )
    ))

    expect_equal(ncol(tm_progressive@interactions), ncol(tm_single@interactions))
    expect_setequal(colnames(tm_progressive@interactions), colnames(tm_single@interactions))
})

test_that("default_score_split_features returns 3 columns aligned to peaks", {
    tm <- create_interaction_traj_model(n_peaks = 150, n_motifs = 6, seed = 5)
    atac <- create_interaction_atac_scores(tm)

    out <- suppressMessages(suppressWarnings(
        default_score_split_features(tm, atac, bin_start = 1, bin_end = 3)
    ))

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), nrow(tm@peak_intervals))
    expect_setequal(colnames(out), c("base_pred", "end_pred", "pred_diff_e_b"))
    expect_true(all(out$base_pred >= 0 & out$base_pred <= 10))
    expect_true(all(out$end_pred >= 0 & out$end_pred <= 10))
})

test_that("default_score_split_features accepts bin column names", {
    tm <- create_interaction_traj_model(n_peaks = 80, n_motifs = 5, seed = 6)
    atac <- create_interaction_atac_scores(tm)

    out <- suppressMessages(suppressWarnings(
        default_score_split_features(tm, atac, bin_start = "bin1", bin_end = "bin3")
    ))
    expect_equal(nrow(out), nrow(tm@peak_intervals))
})

test_that("default_score_split_features errors on identical bins", {
    tm <- create_interaction_traj_model(n_peaks = 80, n_motifs = 5, seed = 7)
    atac <- create_interaction_atac_scores(tm)

    expect_error(
        default_score_split_features(tm, atac, bin_start = 1, bin_end = 1),
        "distinct"
    )
})

test_that("progressive with default_score_split_features injects 3 engineered features", {
    tm <- create_interaction_traj_model(n_peaks = 150, n_motifs = 6, seed = 8)
    atac <- create_interaction_atac_scores(tm)

    builder <- function(tm) default_score_split_features(tm, atac, bin_start = 1, bin_end = 3)

    tm_out <- suppressMessages(suppressWarnings(
        add_interactions_progressive(
            tm,
            thresholds = c(0.01, 0.0005),
            only_sig_motifs = c(TRUE, FALSE),
            only_sig_add_motifs = c(TRUE, TRUE),
            additional_features_builder = builder,
            seed = 8
        )
    ))

    # Engineered features made it into @additional_features
    expect_true(all(c("base_pred", "end_pred", "pred_diff_e_b") %in%
        colnames(tm_out@additional_features)))

    # Second pass ran (force=TRUE), so interactions from pass 1 were
    # replaced. Post-pass, the model has some interactions.
    expect_gt(ncol(tm_out@interactions), 0)

    # R^2 train should at least match the single-pass baseline on the
    # same fixture seeded identically (progressive is expected to do
    # as well or better because it uses engineered features too).
    tm_single <- suppressMessages(suppressWarnings(
        add_interactions(tm, interaction_threshold = 0.0005, seed = 8,
            only_sig_motifs = FALSE, only_sig_add_motifs = TRUE)
    ))
    r2_single <- cor(tm_single@predicted_diff_score[tm@type == "train"],
                     tm@diff_score[tm@type == "train"])^2
    r2_progressive <- cor(tm_out@predicted_diff_score[tm_out@type == "train"],
                          tm_out@diff_score[tm_out@type == "train"])^2
    # Allow 0.01 tolerance because the progressive path rescales additional
    # features and the engineered features change the additional-feature
    # anchor pool for pass 2 — this is a soft parity gate.
    expect_gte(r2_progressive, r2_single - 0.01)
})

test_that("min_signal_correlation drops interactions below threshold x max(|cor|)", {
    tm <- create_interaction_traj_model(n_peaks = 200, n_motifs = 10, seed = 1)

    tm_unfiltered <- suppressMessages(suppressWarnings(
        add_interactions(tm, interaction_threshold = 0.001, seed = 1)
    ))
    tm_filtered <- suppressMessages(suppressWarnings(
        add_interactions(tm, interaction_threshold = 0.001, seed = 1,
            min_signal_correlation = 0.5)
    ))

    expect_lte(ncol(tm_filtered@interactions), ncol(tm_unfiltered@interactions))
    expect_gt(ncol(tm_filtered@interactions), 0)

    # All kept interactions must individually exceed 0.5 * max(|cor|).
    train_idx <- which(tm@type == "train")
    y_train <- norm01(tm@diff_score)[train_idx]
    cm <- abs(tgs_cor(tm_filtered@interactions[train_idx, ], as.matrix(y_train))[, 1])
    expect_true(all(cm > 0.5 * max(cm) - 1e-10))
})

test_that("interaction_scale_factor scales the interaction matrix linearly", {
    tm <- create_interaction_traj_model(n_peaks = 100, n_motifs = 5, seed = 2)

    tm_a <- suppressMessages(suppressWarnings(
        add_interactions(tm, interaction_threshold = 0.001, seed = 2,
            interaction_scale_factor = 1)
    ))
    tm_b <- suppressMessages(suppressWarnings(
        add_interactions(tm, interaction_threshold = 0.001, seed = 2,
            interaction_scale_factor = 3)
    ))

    # Interactions should be scaled 3x
    expect_equal(tm_b@interactions, tm_a@interactions * 3, tolerance = 1e-10)
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
