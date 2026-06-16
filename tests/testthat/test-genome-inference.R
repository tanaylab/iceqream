# Genome-dependent round-trip / consistency tests for the inference and export
# paths. These check that a trained model's stored predictions and per-motif
# energies are reproduced by the public inference / export entry points - the
# kind of cross-component agreement where subtle scale/size bugs hide.
#
# Built via the shared genome_traj_models() fixture (helper-genome.R); gated like
# test-e2e.R (skips on CRAN / without the genome or cached vignette data).

test_that("predict(TrajectoryModel) reproduces the training predictions exactly", {
    base <- genome_traj_models()$base
    p <- suppressWarnings(suppressMessages(predict(base, peak_intervals = base@peak_intervals)))
    expect_length(as.numeric(p), nrow(base@peak_intervals))
    expect_equal(as.numeric(p), base@predicted_diff_score, tolerance = 1e-8)
})

test_that("pbm_list.gextract reproduces the model's normalized energies (extracts at the PBM size)", {
    # Regression guard for the size bug: pbm_list.gextract must extract sequences
    # at the PBMs' trained size, not the input intervals' native width (the
    # gastrulation peaks are 300bp but the model is trained at peaks_size = 500).
    base <- genome_traj_models()$base
    pl <- traj_model_to_pbm_list(base)
    e <- suppressWarnings(suppressMessages(pbm_list.gextract(pl, base@peak_intervals)))
    mn <- names(base@motif_models)
    ec <- as.matrix(e[, mn, drop = FALSE])
    expect_equal(ec, base@normalized_energies[, mn], tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("pbm.gextract (single PBM) reproduces one motif's normalized energy", {
    base <- genome_traj_models()$base
    pl <- traj_model_to_pbm_list(base)
    m <- names(base@motif_models)[1]
    e <- suppressWarnings(suppressMessages(pbm.gextract(pl[[m]], base@peak_intervals)))
    # Looser tolerance than the pbm_list path: pbm.gextract uses prego::compute_pwm
    # (single motif) while training used the batched extract_pwm, which agree to
    # ~1e-5 rather than bit-for-bit.
    expect_equal(e[[m]], base@normalized_energies[, m], tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("create_iq_model + predict reproduces training predictions (no interactions)", {
    base <- genome_traj_models()$base
    iqm <- create_iq_model(base)
    p <- as.numeric(suppressWarnings(suppressMessages(
        predict(iqm, intervals = base@peak_intervals)
    )))
    expect_equal(p, base@predicted_diff_score, tolerance = 1e-6)
})

test_that("create_iq_model + predict reproduces an interaction model's predictions", {
    base <- genome_traj_models()$base
    bi <- suppressWarnings(suppressMessages(
        add_interactions(base, interaction_threshold = 0.001, seed = 60427)
    ))
    skip_if(ncol(bi@interactions) < 1, "fixture produced no interactions")
    iqm <- create_iq_model(bi)

    # add_interactions() rescales predicted_diff_score (rescale_pred = TRUE),
    # consistent with regress_trajectory_motifs(), so the IQmodel reproduces it
    # with the default predict() (rescale = TRUE) - exactly like the
    # no-interaction case above. (Before this was fixed, an interaction model's
    # predictions were left on the raw [0,1] logistic scale.)
    p <- as.numeric(suppressWarnings(suppressMessages(
        predict(iqm, intervals = base@peak_intervals)
    )))
    expect_equal(p, bi@predicted_diff_score, tolerance = 1e-6, ignore_attr = TRUE)
})
