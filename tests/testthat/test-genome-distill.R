# Genome-dependent adversarial tests for the distillation / merge /
# homogenization paths. These exercise functions that call prego on extracted
# sequences (merge_trajectory_motifs, distill_traj_model,
# distill_traj_model_multi) and cannot run on synthetic in-memory fixtures.
#
# Gated exactly like test-e2e.R: skipped on CRAN and when the misha genome or
# the cached vignette data is unavailable. Override the genome location with
# ICEQREAM_GENOME_ROOT (defaults to mm10). The base models are built once and
# memoized across the tests in this file.

.genome_root <- Sys.getenv("ICEQREAM_GENOME_ROOT", unset = "/home/aviezerl/mm10")
.cache_dir <- file.path(tools::R_user_dir("iceqream", "cache"), "gastrulation-example")

# Memoized builder: returns list(base = <traj 1->4>, base2 = <traj 1->2>),
# both regressed on a small subset of the cached gastrulation data with
# additional_features = NULL (which also exercises the NULL-handling fix).
.genome_models <- local({
    cache <- NULL
    function() {
        if (!is.null(cache)) {
            return(cache)
        }
        skip_on_cran()
        skip_if_not_installed("misha")
        skip_if_not_installed("prego")
        skip_if_not_installed("misha.ext")
        skip_if(!dir.exists(.genome_root), paste("misha genome not found at", .genome_root))
        needed <- file.path(.cache_dir, c("peak_intervals.rds", "atac_scores.rds", "motif_energies.rds"))
        skip_if(!all(file.exists(needed)), "cached vignette data not available")

        misha::gsetroot(.genome_root)
        pk <- readr::read_rds(needed[1])
        at <- as.matrix(readr::read_rds(needed[2]))
        me <- readr::read_rds(needed[3])

        set.seed(1)
        ri <- sort(sample(nrow(pk), 600))
        ci <- sort(sample(ncol(me), 80))
        pk2 <- pk[ri, c("chrom", "start", "end", "peak_name")]
        at2 <- as.matrix(at[ri, ])
        me2 <- me[ri, ci]
        rownames(me2) <- pk2$peak_name

        mk <- function(bin_start, bin_end) {
            suppressWarnings(suppressMessages(regress_trajectory_motifs(
                peak_intervals = pk2, atac_scores = at2, motif_energies = me2,
                additional_features = NULL, norm_intervals = pk2,
                min_tss_distance = NULL, max_motif_num = 10,
                bin_start = bin_start, bin_end = bin_end,
                n_prego_motifs = 0, peaks_size = 500, seed = 60427
            )))
        }

        cache <<- list(base = mk(1, 4), base2 = mk(1, 2))
        cache
    }
})

test_that("regress_trajectory_motifs accepts additional_features = NULL (normalized to a 0-column data.frame)", {
    m <- .genome_models()$base
    expect_s4_class(m, "TrajectoryModel")
    expect_equal(ncol(m@additional_features), 0L)
    expect_equal(nrow(m@additional_features), nrow(m@peak_intervals))
    expect_true(all(is.finite(m@predicted_diff_score)))
    expect_gt(length(m@motif_models), 0)
})

test_that("merge_trajectory_motifs replaces the merged motifs with the new one consistently across slots", {
    base <- .genome_models()$base
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
    base <- .genome_models()$base
    distilled <- suppressWarnings(suppressMessages(
        distill_traj_model(base, max_motif_num = 6)
    ))
    expect_s4_class(distilled, "TrajectoryModel")
    expect_lte(length(distilled@motif_models), length(base@motif_models))
    expect_equal(nrow(distilled@normalized_energies), nrow(base@peak_intervals))
    expect_true(all(is.finite(distilled@predicted_diff_score)))
})

test_that("distill_traj_model_multi homogenizes multiple trajectories into a shared motif set", {
    m <- .genome_models()
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
