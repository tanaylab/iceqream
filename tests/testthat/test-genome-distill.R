# Genome-dependent adversarial tests for the distillation / merge /
# homogenization paths. These exercise functions that call prego on extracted
# sequences (merge_trajectory_motifs, distill_traj_model,
# distill_traj_model_multi) and cannot run on synthetic in-memory fixtures.
#
# Gated exactly like test-e2e.R via the shared genome_traj_models() fixture
# (see helper-genome.R): skipped on CRAN and when the misha genome or the cached
# vignette data is unavailable. ICEQREAM_GENOME_ROOT overrides the genome.

test_that("regress_trajectory_motifs accepts additional_features = NULL (normalized to a 0-column data.frame)", {
    m <- genome_traj_models()$base
    expect_s4_class(m, "TrajectoryModel")
    expect_equal(ncol(m@additional_features), 0L)
    expect_equal(nrow(m@additional_features), nrow(m@peak_intervals))
    expect_true(all(is.finite(m@predicted_diff_score)))
    expect_gt(length(m@motif_models), 0)
})

test_that("merge_trajectory_motifs replaces the merged motifs with the new one consistently across slots", {
    base <- genome_traj_models()$base
    mm <- names(base@motif_models)[1:2]
    merged <- suppressWarnings(suppressMessages(
        merge_trajectory_motifs(base, mm, "MERGED_TEST")
    ))

    expect_s4_class(merged, "TrajectoryModel")
    expect_true("MERGED_TEST" %in% names(merged@motif_models))
    expect_true("MERGED_TEST" %in% colnames(merged@normalized_energies))

    # The merged-away motifs are removed from every per-motif slot, not just
    # @model_features, so the model stays internally consistent.
    expect_false(any(mm %in% names(merged@motif_models)))
    expect_false(any(mm %in% colnames(merged@normalized_energies)))
    removed_feats <- as.vector(outer(
        mm, c("_low-energy", "_high-energy", "_higher-energy", "_sigmoid"), paste0
    ))
    expect_false(any(removed_feats %in% colnames(merged@model_features)))

    # Net effect: the two merged motifs are replaced by exactly one new motif.
    expect_equal(length(merged@motif_models), length(base@motif_models) - length(mm) + 1L)
    expect_equal(ncol(merged@normalized_energies), ncol(base@normalized_energies) - length(mm) + 1L)

    expect_true(all(is.finite(merged@predicted_diff_score)))
})

test_that("distill_traj_model returns a valid model no larger than its input", {
    base <- genome_traj_models()$base
    distilled <- suppressWarnings(suppressMessages(
        distill_traj_model(base, max_motif_num = 6)
    ))
    expect_s4_class(distilled, "TrajectoryModel")
    expect_lte(length(distilled@motif_models), length(base@motif_models))
    expect_equal(nrow(distilled@normalized_energies), nrow(base@peak_intervals))
    expect_true(all(is.finite(distilled@predicted_diff_score)))
})

test_that("distill_traj_model_multi homogenizes multiple trajectories into a shared motif set", {
    m <- genome_traj_models()
    multi <- suppressWarnings(suppressMessages(distill_traj_model_multi(
        list(t1 = m$base, t2 = m$base2),
        max_motif_num = 12, parallel = FALSE
    )))
    expect_s4_class(multi, "TrajectoryModelMulti")
    expect_length(multi@models, 2)
    expect_gt(length(multi@motif_models), 0)
    # every per-trajectory model uses the shared motif set and predicts finitely
    expect_true(all(purrr::map_lgl(
        multi@models, ~ all(is.finite(.x@predicted_diff_score))
    )))
})
