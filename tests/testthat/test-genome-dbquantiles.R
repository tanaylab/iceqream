# Genome-backed regression test for the db_quantiles normalization path.
#
# Gated like the other genome tests (skip on CRAN / when the misha genome or the
# cached vignette data is missing; override the genome with ICEQREAM_GENOME_ROOT).
#
# Locks the conclusion of the db_quantiles correctness investigation: the
# fixed-range db_quantiles normalization in compute_motif_energies() feeds only
# the (scale-invariant) correlation-based motif selection; the trained model's
# energies are always recomputed observed-range during distillation, so a model
# built from db_quantiles energies stays normalization-self-consistent and its
# predictions round-trip through inference exactly - identical to the default
# path. (See dev review doc, "db_quantiles correctness investigation".)

test_that("a model trained on db_quantiles energies round-trips through inference exactly", {
    skip_on_cran()
    skip_if_not_installed("misha")
    skip_if_not_installed("prego")
    groot <- Sys.getenv("ICEQREAM_GENOME_ROOT", unset = "/home/aviezerl/mm10")
    skip_if(!dir.exists(groot), paste("misha genome not found at", groot))
    cdir <- file.path(tools::R_user_dir("iceqream", "cache"), "gastrulation-example")
    needed <- file.path(cdir, c("peak_intervals.rds", "atac_scores.rds"))
    skip_if(!all(file.exists(needed)), "cached vignette data not available")

    misha::gsetroot(groot)
    utils::data("motif_db", package = "iceqream", envir = environment())
    utils::data("mouse_db_quantiles", package = "iceqream", envir = environment())

    pk <- readr::read_rds(needed[1])
    at <- as.matrix(readr::read_rds(needed[2]))
    set.seed(1)
    ri <- sort(sample(nrow(pk), 400))
    pk2 <- pk[ri, c("chrom", "start", "end", "peak_name")]
    at2 <- at[ri, ]
    set.seed(2)
    sel <- sort(sample(unique(motif_db$motif), 50))
    db_sub <- motif_db[motif_db$motif %in% sel, ]

    e_dbq <- suppressWarnings(suppressMessages(compute_motif_energies(
        peak_intervals = pk2[, c("chrom", "start", "end")], db = db_sub,
        normalization_intervals = pk2[, c("chrom", "start", "end")],
        energy_norm_quantile = 1, db_quantiles = mouse_db_quantiles
    )))
    rownames(e_dbq) <- pk2$peak_name

    m <- suppressWarnings(suppressMessages(regress_trajectory_motifs(
        peak_intervals = pk2, atac_scores = at2, motif_energies = e_dbq,
        additional_features = NULL, norm_intervals = pk2, min_tss_distance = NULL,
        max_motif_num = 8, bin_start = 1, bin_end = ncol(at2), n_prego_motifs = 0,
        peaks_size = 500, seed = 60427
    )))

    # 1. The model's stored energies are observed-range (recomputed during
    #    distillation), so they exactly match what inference recomputes.
    ne_infer <- suppressWarnings(suppressMessages(
        iceqream:::calc_traj_model_energies(m, m@peak_intervals)
    ))
    common <- intersect(colnames(m@normalized_energies), colnames(ne_infer))
    expect_gt(length(common), 0)
    expect_equal(
        unname(m@normalized_energies[, common]),
        unname(ne_infer[, common]),
        tolerance = 1e-8
    )

    # 2. predict() reproduces the training predictions exactly despite the
    #    db_quantiles input (no train/infer mismatch).
    p <- suppressWarnings(suppressMessages(predict(m, peak_intervals = m@peak_intervals)))
    expect_equal(p, m@predicted_diff_score, tolerance = 1e-6)

    # 3. The exported IQ model reproduces them too.
    iqm <- suppressWarnings(suppressMessages(create_iq_model(m)))
    p_iqm <- suppressWarnings(suppressMessages(predict(iqm, intervals = m@peak_intervals)))
    expect_gt(cor(p_iqm, m@predicted_diff_score), 0.999)
})
