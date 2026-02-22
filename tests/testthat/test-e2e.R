# =============================================================================
# End-to-end test for iq_regression using vignette gastrulation data
# =============================================================================

# Helper: download and cache vignette data, returns path to extracted directory
get_vignette_data <- function() {
    cache_dir <- tools::R_user_dir("iceqream", "cache")
    data_dir <- file.path(cache_dir, "gastrulation-example")

    if (dir.exists(data_dir) &&
        file.exists(file.path(data_dir, "peak_intervals.rds")) &&
        file.exists(file.path(data_dir, "atac_scores.rds")) &&
        file.exists(file.path(data_dir, "additional_features.rds")) &&
        file.exists(file.path(data_dir, "gastrulation_intervals.tsv")) &&
        file.exists(file.path(data_dir, "motif_energies.rds"))) {
        return(data_dir)
    }

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    tarball <- file.path(cache_dir, "gastrulation-example.tar.gz")

    tryCatch(
        {
            download.file(
                "https://iceqream.s3.eu-west-1.amazonaws.com/gastrulation-example.tar.gz",
                tarball,
                mode = "wb",
                quiet = TRUE
            )
            untar(tarball, exdir = cache_dir)
        },
        error = function(e) {
            return(NULL)
        }
    )

    if (!dir.exists(data_dir)) {
        return(NULL)
    }
    return(data_dir)
}

# =============================================================================
# Full pipeline test
# =============================================================================
test_that("iq_regression produces valid TrajectoryModel on gastrulation data", {
    skip_on_cran()
    skip_if_not_installed("prego")
    skip_if_not_installed("misha")
    skip_if_not_installed("misha.ext")

    genome_root <- Sys.getenv("ICEQREAM_GENOME_ROOT", unset = "/home/aviezerl/mm10")
    skip_if(!dir.exists(genome_root), paste("misha genome not found at", genome_root))

    data_dir <- get_vignette_data()
    skip_if(is.null(data_dir), "Could not download vignette data")

    # Load data
    peak_intervals <- readr::read_rds(file.path(data_dir, "peak_intervals.rds"))
    atac_scores <- readr::read_rds(file.path(data_dir, "atac_scores.rds"))
    additional_features <- readr::read_rds(file.path(data_dir, "additional_features.rds"))
    normalization_intervals <- readr::read_tsv(
        file.path(data_dir, "gastrulation_intervals.tsv"),
        show_col_types = FALSE
    )
    motif_energies <- readr::read_rds(file.path(data_dir, "motif_energies.rds"))

    # Set genome
    misha::gsetroot(genome_root)

    # Run full pipeline
    traj_model <- iq_regression(
        peak_intervals = peak_intervals,
        atac_scores = atac_scores,
        motif_energies = motif_energies,
        normalize_energies = FALSE,
        additional_features = additional_features,
        norm_intervals = normalization_intervals,
        seed = 60427,
        n_prego_motifs = 0,
        frac_train = 0.8,
        max_motif_num = 30,
        plot_report = FALSE,
        rename_motifs = FALSE
    )

    # --- Validate S4 class ---
    expect_s4_class(traj_model, "TrajectoryModel")

    # --- Motif models ---
    expect_true(length(traj_model@motif_models) > 0, info = "Model should have at least one motif")

    # --- Predictions ---
    expect_true(length(traj_model@predicted_diff_score) > 0)
    expect_false(all(is.na(traj_model@predicted_diff_score)), info = "predicted_diff_score should not be all NA")
    expect_false(all(traj_model@predicted_diff_score == 0), info = "predicted_diff_score should not be all zero")

    # --- diff_score and predicted_diff_score have same length ---
    expect_equal(
        length(traj_model@diff_score),
        length(traj_model@predicted_diff_score),
        info = "diff_score and predicted_diff_score should have same length"
    )

    # --- Type vector ---
    expect_true(all(traj_model@type %in% c("train", "test")))
    expect_true(sum(traj_model@type == "train") > 0)
    expect_true(sum(traj_model@type == "test") > 0)

    # --- R² statistics (computed from model slots since stats may not include test R²) ---
    train_idx <- which(traj_model@type == "train")
    test_idx <- which(traj_model@type == "test")
    r2_train <- cor(
        traj_model@diff_score[train_idx],
        traj_model@predicted_diff_score[train_idx],
        use = "pairwise.complete.obs"
    )^2
    r2_test <- cor(
        traj_model@diff_score[test_idx],
        traj_model@predicted_diff_score[test_idx],
        use = "pairwise.complete.obs"
    )^2
    expect_gt(r2_train, 0.1, label = "Train R²")
    expect_gt(r2_test, 0, label = "Test R²")

    # --- Peak intervals consistency ---
    expect_true(nrow(traj_model@peak_intervals) > 0)
    expect_equal(
        nrow(traj_model@peak_intervals),
        length(traj_model@diff_score),
        info = "peak_intervals rows should match diff_score length"
    )

    # --- Normalized energies dimensions ---
    expect_equal(
        nrow(traj_model@normalized_energies),
        nrow(traj_model@peak_intervals),
        info = "normalized_energies rows should match number of peaks"
    )
    motif_names <- names(traj_model@motif_models)
    expect_true(
        all(motif_names %in% colnames(traj_model@normalized_energies)),
        info = "All motif model names should appear in normalized_energies columns"
    )

    # --- Coefficients ---
    expect_true(nrow(traj_model@coefs) > 0, info = "Should have model coefficients")

    # ==========================================================================
    # Extended tests: traj_model_to_pbm_list and pbm_list.compute
    # ==========================================================================

    pbm_list <- traj_model_to_pbm_list(traj_model)
    expect_true(is.list(pbm_list))
    expect_true(length(pbm_list) > 0)
    expect_true(all(sapply(pbm_list, function(x) methods::is(x, "PBM"))))
    expect_equal(
        sort(names(pbm_list)),
        sort(names(traj_model@motif_models)),
        info = "PBM list names should match motif model names"
    )

    # Compute on a few small sequences
    test_seqs <- c(
        paste0(rep("ACGT", 125), collapse = ""),
        paste0(rep("TGCA", 125), collapse = ""),
        paste0(rep("AATT", 125), collapse = "")
    )
    pbm_result <- pbm_list.compute(pbm_list, test_seqs)
    expect_true(is.matrix(pbm_result))
    expect_equal(nrow(pbm_result), 3)
    expect_equal(ncol(pbm_result), length(pbm_list))
})

# =============================================================================
# Test custom bin_start/bin_end: verifies test diff_score uses correct bins
# =============================================================================
test_that("iq_regression with custom bin_start/bin_end computes test diff correctly", {
    skip_on_cran()
    skip_if_not_installed("prego")
    skip_if_not_installed("misha")
    skip_if_not_installed("misha.ext")

    genome_root <- Sys.getenv("ICEQREAM_GENOME_ROOT", unset = "/home/aviezerl/mm10")
    skip_if(!dir.exists(genome_root), paste("misha genome not found at", genome_root))

    data_dir <- get_vignette_data()
    skip_if(is.null(data_dir), "Could not download vignette data")

    peak_intervals <- readr::read_rds(file.path(data_dir, "peak_intervals.rds"))
    atac_scores <- readr::read_rds(file.path(data_dir, "atac_scores.rds"))
    additional_features <- readr::read_rds(file.path(data_dir, "additional_features.rds"))
    normalization_intervals <- readr::read_tsv(
        file.path(data_dir, "gastrulation_intervals.tsv"),
        show_col_types = FALSE
    )
    motif_energies <- readr::read_rds(file.path(data_dir, "motif_energies.rds"))

    misha::gsetroot(genome_root)

    # Use bin1 -> bin2 (not default bin1 -> bin4) to exercise custom bin forwarding
    traj_model <- iq_regression(
        peak_intervals = peak_intervals,
        atac_scores = atac_scores,
        motif_energies = motif_energies,
        normalize_energies = FALSE,
        additional_features = additional_features,
        norm_intervals = normalization_intervals,
        seed = 60427,
        n_prego_motifs = 0,
        frac_train = 0.8,
        max_motif_num = 10,
        bin_start = 1,
        bin_end = 2,
        plot_report = FALSE,
        rename_motifs = FALSE
    )

    expect_s4_class(traj_model, "TrajectoryModel")

    test_idx <- which(traj_model@type == "test")
    train_idx <- which(traj_model@type == "train")

    # Test diff_score should be the bin2 - bin1 difference (after normalization),
    # not the default bin4 - bin1
    r2_test <- cor(
        traj_model@diff_score[test_idx],
        traj_model@predicted_diff_score[test_idx],
        use = "pairwise.complete.obs"
    )^2
    r2_train <- cor(
        traj_model@diff_score[train_idx],
        traj_model@predicted_diff_score[train_idx],
        use = "pairwise.complete.obs"
    )^2

    # Test R² should be in a reasonable range relative to train R²
    # (with wrong bins, test R² would be near 0)
    expect_gt(r2_test, r2_train * 0.3,
        label = "Test R² should be comparable to train R² when using correct bins"
    )
})
