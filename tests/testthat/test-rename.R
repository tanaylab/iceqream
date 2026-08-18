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

test_that("rename_motif_models handles motif names containing '::' (JASPAR dimers)", {
    # Regression test: dimer motifs are named like "JASPAR.GATA1::TAL1". The
    # interaction branch used to claim any column with a ":" in it, so those
    # columns went through rename_interaction_names(), which splits on ":",
    # matches nothing and returns the name unchanged. @motif_models was renamed
    # while @model_features kept the old column, and inference then failed with
    # "Missing model feature columns" / "subscript out of bounds".
    tm <- create_mock_traj_model(n_peaks = 80, n_motifs = 3)

    # step 1: give one motif a dimer-style name (old names have no ":" here, so
    # this step is unaffected by the bug)
    to_dimer <- setNames(c("JASPAR.GATA1::TAL1", "B", "C"), names(tm@motif_models))
    tm <- suppressWarnings(suppressMessages(rename_motif_models(tm, to_dimer)))
    expect_true(any(grepl("^JASPAR.GATA1::TAL1_", colnames(tm@model_features))))

    # step 2: rename away from the dimer name - this is what used to break
    out <- suppressWarnings(suppressMessages(
        rename_motif_models(tm, setNames(c("TAL1", "B2", "C2"), names(tm@motif_models)))
    ))

    expect_setequal(names(out@motif_models), c("TAL1", "B2", "C2"))
    expect_false(any(is.na(colnames(out@model_features))))
    expect_false(any(grepl("GATA1", colnames(out@model_features))))
    # every motif has its 4 logist feature columns under the new name
    expect_true(all(paste0("TAL1_", c("low-energy", "high-energy", "higher-energy", "sigmoid")) %in%
        colnames(out@model_features)))
    # and the model stays usable: features present for every motif model
    expect_true(all(names(out@motif_models) %in% sub("_[^_]+$", "", colnames(out@model_features))))
})
