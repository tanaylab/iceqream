# Adversarial tests for rename_motif_models() (genome-free). Focus: renaming a
# model that has interactions, where the interaction feature names ("A:B") and
# the @interactions slot must be renamed too, plus the input guards.

test_that("rename_motif_models renames a plain model and stays self-consistent", {
    tm <- create_mock_traj_model(n_peaks = 80, n_motifs = 3)
    nm <- setNames(c("A", "B", "C"), names(tm@motif_models))
    out <- suppressWarnings(suppressMessages(rename_motif_models(tm, nm)))

    expect_setequal(names(out@motif_models), c("A", "B", "C"))
    expect_setequal(colnames(out@normalized_energies), c("A", "B", "C"))
    expect_false(any(is.na(colnames(out@model_features))))
    # the model_features carry the new names (4 logist features each)
    expect_true(all(grepl("^(A|B|C)_", colnames(out@model_features))))
    expect_equal(out@params$names_map, nm)
})

test_that("rename_motif_models renames interaction features on both @model_features and @interactions", {
    # Regression test: renaming an interaction model previously set the
    # interaction columns of @model_features to NA and left @interactions on the
    # old motif names, silently corrupting the model.
    tm <- create_interaction_traj_model(n_peaks = 200, n_motifs = 6, seed = 1)
    tm <- suppressWarnings(suppressMessages(add_interactions(tm, interaction_threshold = 0.001, seed = 1)))
    skip_if(ncol(tm@interactions) < 1, "fixture produced no interactions")

    nm <- setNames(paste0("TF", seq_along(names(tm@motif_models))), names(tm@motif_models))
    out <- suppressWarnings(suppressMessages(rename_motif_models(tm, nm)))

    expect_false(any(is.na(colnames(out@model_features))))
    # @interactions columns are renamed ("TFi:TFj") with no leftover old names
    expect_true(all(grepl("^TF[0-9]+:TF[0-9]+$", colnames(out@interactions))))
    expect_false(any(grepl("motif", colnames(out@interactions))))
    expect_false(any(grepl("motif", colnames(out@model_features))))
    expect_true(all(is.finite(out@predicted_diff_score)))

    # The renamed model is internally consistent: feature parsing and partial
    # responses work on the new names.
    expect_silent(ftv <- iceqream:::feat_to_variable(out, add_types = TRUE))
    pr <- iceqream:::compute_partial_response(out, vars = c("TF1", colnames(out@interactions)[1]))
    expect_setequal(colnames(pr), c("TF1", colnames(out@interactions)[1]))
})

test_that("rename_motif_models rejects duplicate target names and incomplete maps", {
    tm <- create_mock_traj_model(n_peaks = 60, n_motifs = 3)
    dup <- setNames(c("A", "A", "B"), names(tm@motif_models))
    expect_error(suppressMessages(rename_motif_models(tm, dup)), "unique")

    incomplete <- setNames(c("A", "B"), names(tm@motif_models)[1:2])
    expect_error(suppressMessages(rename_motif_models(tm, incomplete)), "not found")
})
