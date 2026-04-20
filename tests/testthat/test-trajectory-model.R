# Tests for TrajectoryModel S4 class and traj-model-utils functions
# These tests are designed to run WITHOUT misha
# create_mock_traj_model() is defined in helper-traj-model.R

# =============================================================================
# TrajectoryModel constructor
# =============================================================================
test_that("TrajectoryModel constructor creates a valid S4 object", {
    tm <- create_mock_traj_model()

    expect_s4_class(tm, "TrajectoryModel")
    expect_true(methods::is(tm, "TrajectoryModel"))
    expect_equal(length(tm@motif_models), 2)
    expect_equal(nrow(tm@peak_intervals), 50)
    expect_equal(length(tm@diff_score), 50)
    expect_equal(length(tm@predicted_diff_score), 50)
    expect_equal(ncol(tm@normalized_energies), 2)
    expect_equal(nrow(tm@normalized_energies), 50)
})

test_that("TrajectoryModel slots have correct types", {
    tm <- create_mock_traj_model()

    expect_true(is.character(tm@type))
    expect_true(is.numeric(tm@diff_score))
    expect_true(is.numeric(tm@predicted_diff_score))
    expect_true(is.matrix(tm@normalized_energies))
    expect_true(is.matrix(tm@model_features))
    expect_true(is.data.frame(tm@peak_intervals))
    expect_true(is.data.frame(tm@coefs))
    expect_true(is.list(tm@motif_models))
    expect_true(is.list(tm@params))
})

test_that("TrajectoryModel type vector has correct values", {
    tm <- create_mock_traj_model(n_peaks = 100, train_frac = 0.7)

    expect_equal(sum(tm@type == "train"), 70)
    expect_equal(sum(tm@type == "test"), 30)
    expect_true(all(tm@type %in% c("train", "test")))
})

# =============================================================================
# validate_traj_model
# =============================================================================
test_that("validate_traj_model accepts a valid TrajectoryModel", {
    tm <- create_mock_traj_model()
    expect_silent(iceqream:::validate_traj_model(tm))
})

test_that("validate_traj_model rejects non-TrajectoryModel objects", {
    expect_error(iceqream:::validate_traj_model(list()), "TrajectoryModel")
    expect_error(iceqream:::validate_traj_model(42), "TrajectoryModel")
    expect_error(iceqream:::validate_traj_model("not a model"), "TrajectoryModel")
})

# =============================================================================
# get_model_coefs
# =============================================================================
test_that("get_model_coefs returns a data frame with correct structure", {
    tm <- create_mock_traj_model()
    coefs <- iceqream:::get_model_coefs(tm@model, s = tm@params$lambda)

    expect_true(is.data.frame(coefs))
    expect_true("variable" %in% colnames(coefs))
    # The logistic feature types
    expected_types <- c("low-energy", "high-energy", "higher-energy", "sigmoid")
    expect_true(all(expected_types %in% colnames(coefs)))
})

test_that("get_model_coefs extracts variable names from feature columns", {
    tm <- create_mock_traj_model(n_motifs = 3)
    coefs <- iceqream:::get_model_coefs(tm@model, s = tm@params$lambda)

    expect_true(all(c("motif1", "motif2", "motif3") %in% coefs$variable))
})

test_that("get_model_coefs returns one row per variable", {
    tm <- create_mock_traj_model(n_motifs = 2)
    coefs <- iceqream:::get_model_coefs(tm@model, s = tm@params$lambda)

    # Should have 2 rows for 2 motifs
    expect_equal(nrow(coefs), 2)
})

test_that("get_model_coefs replaces NAs with zeros", {
    tm <- create_mock_traj_model()
    coefs <- iceqream:::get_model_coefs(tm@model, s = tm@params$lambda)

    # No NAs should be present
    expect_false(any(is.na(coefs)))
})

# =============================================================================
# show method
# =============================================================================
test_that("show method produces output without error", {
    tm <- create_mock_traj_model()
    # cli writes to stderr/connection, so capture both stdout and messages
    output <- utils::capture.output(show(tm), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(nchar(output_text) > 0)
    expect_true(grepl("TrajectoryModel", output_text))
})

test_that("show method reports correct number of motifs", {
    tm <- create_mock_traj_model(n_motifs = 3)
    output <- utils::capture.output(show(tm), type = "message")
    output_text <- paste(output, collapse = " ")
    expect_true(grepl("3", output_text))
})

# =============================================================================
# has_interactions / n_interactions
# =============================================================================
test_that("has_interactions returns FALSE for model without interactions", {
    tm <- create_mock_traj_model()
    expect_false(iceqream:::has_interactions(tm))
})

test_that("n_interactions returns 0 for model without interactions", {
    tm <- create_mock_traj_model()
    expect_equal(iceqream:::n_interactions(tm), 0)
})

# =============================================================================
# traj_model_has_test
# =============================================================================
test_that("traj_model_has_test returns TRUE when test samples exist", {
    tm <- create_mock_traj_model(train_frac = 0.8)
    expect_true(iceqream:::traj_model_has_test(tm))
})

test_that("traj_model_has_test returns FALSE when no test samples", {
    tm <- create_mock_traj_model(train_frac = 1.0)
    expect_false(iceqream:::traj_model_has_test(tm))
})

# =============================================================================
# has_additional_features
# =============================================================================
test_that("has_additional_features returns FALSE for model with empty additional_features", {
    tm <- create_mock_traj_model()
    expect_false(iceqream:::has_additional_features(tm))
})

# =============================================================================
# split_traj_model_to_train_test
# =============================================================================
test_that("split_traj_model_to_train_test splits correctly", {
    tm <- create_mock_traj_model(n_peaks = 100, train_frac = 0.7)

    result <- split_traj_model_to_train_test(tm)

    expect_true(is.list(result))
    expect_true("train" %in% names(result))
    expect_true("test" %in% names(result))

    train_model <- result$train
    test_model <- result$test

    expect_s4_class(train_model, "TrajectoryModel")
    expect_s4_class(test_model, "TrajectoryModel")

    # All types should be "train" in the train split
    expect_true(all(train_model@type == "train"))
    # All types should be "test" in the test split
    expect_true(all(test_model@type == "test"))

    # Sizes should add up
    expect_equal(
        nrow(train_model@peak_intervals) + nrow(test_model@peak_intervals),
        nrow(tm@peak_intervals)
    )
    expect_equal(
        length(train_model@diff_score) + length(test_model@diff_score),
        length(tm@diff_score)
    )
})

test_that("split_traj_model_to_train_test with no test returns NULL for test", {
    tm <- create_mock_traj_model(n_peaks = 50, train_frac = 1.0)

    result <- split_traj_model_to_train_test(tm)

    expect_null(result$test)
    expect_s4_class(result$train, "TrajectoryModel")
    # The train model should be the original
    expect_equal(nrow(result$train@peak_intervals), 50)
})

test_that("split preserves model_features dimensions", {
    tm <- create_mock_traj_model(n_peaks = 80, n_motifs = 3, train_frac = 0.75)

    result <- split_traj_model_to_train_test(tm)

    expect_equal(ncol(result$train@model_features), ncol(tm@model_features))
    expect_equal(ncol(result$test@model_features), ncol(tm@model_features))
    expect_equal(
        nrow(result$train@model_features) + nrow(result$test@model_features),
        nrow(tm@model_features)
    )
})

# =============================================================================
# filter_traj_model_intervals
# =============================================================================
test_that("filter_traj_model_intervals filters correctly", {
    tm <- create_mock_traj_model(n_peaks = 50)
    idxs <- c(1, 5, 10, 20)

    filtered <- filter_traj_model_intervals(tm, idxs)

    expect_equal(nrow(filtered@peak_intervals), 4)
    expect_equal(length(filtered@diff_score), 4)
    expect_equal(length(filtered@predicted_diff_score), 4)
    expect_equal(nrow(filtered@model_features), 4)
    expect_equal(nrow(filtered@normalized_energies), 4)
    expect_equal(length(filtered@type), 4)
})

test_that("filter_traj_model_intervals with logical index", {
    tm <- create_mock_traj_model(n_peaks = 20)
    keep <- tm@type == "train"

    filtered <- filter_traj_model_intervals(tm, keep)

    expect_equal(nrow(filtered@peak_intervals), sum(keep))
    expect_true(all(filtered@type == "train"))
})

test_that("filter_traj_model_intervals preserves column structure", {
    tm <- create_mock_traj_model(n_peaks = 30, n_motifs = 3)
    idxs <- 1:10

    filtered <- filter_traj_model_intervals(tm, idxs)

    expect_equal(ncol(filtered@model_features), ncol(tm@model_features))
    expect_equal(ncol(filtered@normalized_energies), ncol(tm@normalized_energies))
    expect_equal(colnames(filtered@normalized_energies), colnames(tm@normalized_energies))
})

# =============================================================================
# add_traj_model_stats
# =============================================================================
test_that("add_traj_model_stats computes R^2 stats", {
    tm <- create_mock_traj_model()

    tm_with_stats <- iceqream:::add_traj_model_stats(tm)

    expect_true("stats" %in% names(tm_with_stats@params))
    expect_true("r2_train" %in% names(tm_with_stats@params$stats))
    expect_true("r2_test" %in% names(tm_with_stats@params$stats))
    expect_true(is.numeric(tm_with_stats@params$stats$r2_train))
    expect_true(is.numeric(tm_with_stats@params$stats$r2_test))
    expect_true(tm_with_stats@params$stats$r2_train >= 0)
    expect_true(tm_with_stats@params$stats$r2_train <= 1)
})

test_that("add_traj_model_stats without test set omits r2_test", {
    tm <- create_mock_traj_model(train_frac = 1.0)

    tm_with_stats <- iceqream:::add_traj_model_stats(tm)

    expect_true("r2_train" %in% names(tm_with_stats@params$stats))
    expect_false("r2_test" %in% names(tm_with_stats@params$stats))
})

# =============================================================================
# feat_to_variable
# =============================================================================
test_that("feat_to_variable maps features to variables", {
    tm <- create_mock_traj_model(n_motifs = 2)

    ftv <- iceqream:::feat_to_variable(tm)

    expect_true(is.data.frame(ftv) || tibble::is_tibble(ftv))
    expect_true("feature" %in% colnames(ftv))
    expect_true("variable" %in% colnames(ftv))

    # Each feature should map to a motif variable
    expect_true(all(ftv$variable %in% c("motif1", "motif2")))

    # Number of features should match model_features columns
    expect_equal(nrow(ftv), ncol(tm@model_features))
})

test_that("feat_to_variable with add_types returns type column", {
    tm <- create_mock_traj_model(n_motifs = 2)

    ftv <- iceqream:::feat_to_variable(tm, add_types = TRUE)

    expect_true("type" %in% colnames(ftv))
    # All should be "motif" type since we have no additional features
    expect_true(all(ftv$type == "motif"))
})

# =============================================================================
# rename_motif_models (partial test -- full test requires relearning)
# =============================================================================
test_that("rename_motif_models validates names_map completeness", {
    tm <- create_mock_traj_model(n_motifs = 2)

    # Incomplete names_map should error
    incomplete_map <- c(motif1 = "TF_A")
    expect_error(
        rename_motif_models(tm, incomplete_map),
        "not found"
    )
})

test_that("rename_motif_models rejects duplicate new names", {
    tm <- create_mock_traj_model(n_motifs = 2)

    dup_map <- c(motif1 = "same_name", motif2 = "same_name")
    expect_error(
        rename_motif_models(tm, dup_map),
        "unique"
    )
})

# =============================================================================
# strip_glmnet
# =============================================================================
test_that("strip_glmnet retains essential components", {
    tm <- create_mock_traj_model()
    stripped <- iceqream:::strip_glmnet(tm@model)

    expect_true("beta" %in% names(stripped))
    expect_true("lambda" %in% names(stripped))
    expect_true("a0" %in% names(stripped))
    expect_true(inherits(stripped, "glmnet") || inherits(stripped, "lognet"))
})

test_that("strip_glmnet produces a smaller object", {
    tm <- create_mock_traj_model()
    original_size <- object.size(tm@model)
    stripped <- iceqream:::strip_glmnet(tm@model)
    stripped_size <- object.size(stripped)

    expect_true(stripped_size <= original_size)
})
