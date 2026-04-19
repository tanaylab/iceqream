#' Perform IQ regression on peak intervals
#'
#' Perform IQ regression on peak intervals using the provided ATAC-seq scores, ATAC-seq score differences, normalized intervals, motif energies, and additional features, after dividing the intervals into training and testing sets.
#'
#' @param frac_train  A numeric value indicating the fraction of intervals to use for training (default is 0.8).
#' @param filter_model A logical value indicating whether to filter the model (default is TRUE).
#' @param filter_sample_frac The fraction of samples to use for computing the r2 without each model at the filtering step. When NULL, all samples are used.
#' @param n_cores The number of cores to use for parallel processing. When NULL, the number of threads is automatically determined as 80% of the available cores. See \code{prego::set_parallel()} for more details.
#' @param add_sequences_features Add CG content and dinuceotide content to the additional features.
#' @param train_idxs A vector of indices to use for training. If NULL, the training set is randomly selected.
#' @param test_idxs A vector of indices to use for testing. If NULL, the testing set is the complement of the training set.
#' @param output_dir A directory to save intermediate results. If not NULL, the train and test indices are saved to a CSV file, together with the models at each step (before filtering, after filtering, and after adding interactions). The models are saved as RDS files. If the directory exists, the files are overwritten.
#' @param plot_report A logical value indicating whether to plot the model report. Default is TRUE.
#' @param rename_motifs A logical value indicating whether to rename the motifs based on the HOMER database. Default is TRUE.
#' @param atac_scores Optional. A numeric matrix, representing mean ATAC score per bin per peak. When using \code{\link{preprocess_data}}, this should be the \code{atac_norm_prob} element of the returned list. Rows: peaks, columns: bins. By default iceqream would regress the last column minus the first column. If you want to regress something else, please either change bin_start or bin_end, or provide \code{atac_diff} instead. If \code{normalize_bins} is TRUE, the scores will be normalized to be between 0 and 1.
#' @param prego_min_diff minimal ATAC difference for a peak to participate in prego motif inference. Default: same as \code{min_diff}.
#' @param prego_energy_norm_quantile quantile of the energy used for normalization in prego motif inference. Default: 1
#' @param prego_spat_bin_size size of each spatial bin for prego motif inference. Default: NULL (uses prego default).
#' @param prego_spat_num_bins number of spatial bins for prego motif inference. Default: NULL (uses prego default).
#' @param max_n_interactions maximum number of interactions to consider. Default: NULL (all interactions).
#' @param strategy Interaction-selection strategy when `include_interactions = TRUE`. Either `"progressive"` (default) or `"single"`. Progressive runs [add_interactions_progressive()] — a two-pass selection that seeds with a tight threshold (Akhiad's pattern) and, when multi-bin `atac_scores` are available, injects `base_pred`/`end_pred`/`pred_diff_e_b` engineered additional features between passes via [default_score_split_features()] before a looser second pass. Single reproduces the 0.0.6 single-pass behavior at `interaction_threshold` — use this for numeric parity with the published paper or with the pycqream Phase 2 validation.
#' @param interaction_thresholds When `strategy = "progressive"`, the per-pass `interaction_threshold` values. Default `c(0.01, 0.0005)`.
#' @param interaction_only_sig_motifs When `strategy = "progressive"`, the per-pass `only_sig_motifs`. Default `c(TRUE, FALSE)`.
#' @param interaction_only_sig_add_motifs When `strategy = "progressive"`, the per-pass `only_sig_add_motifs`. Default `c(TRUE, TRUE)`.
#' @param interaction_scale_factor Forwarded to each [add_interactions()] call.
#' @param interaction_min_signal_correlation Forwarded to each [add_interactions()] call; if non-NULL, drops interactions whose training-set absolute correlation with `diff_score` is below this fraction of the best interaction's. `1/8` mirrors Akhiad's manual post-filter.
#'
#' @return An instance of \code{TrajectoryModel} with the final model.
#'
#'
#' @inheritParams regress_trajectory_motifs
#' @inheritParams filter_traj_model
#' @inheritDotParams regress_trajectory_motifs
#' @inherit regress_trajectory_motifs return
#' @export
iq_regression <- function(
  peak_intervals = NULL,
  peaks = NULL,
  atac_scores = NULL,
  atac_diff = NULL,
  normalize_bins = TRUE,
  norm_intervals = NULL,
  motif_energies = NULL,
  additional_features = NULL,
  min_tss_distance = 5000,
  add_sequences_features = TRUE,
  max_motif_num = 30,
  traj_prego = NULL,
  peaks_size = 500,
  bin_start = 1,
  bin_end = NULL,
  seed = 60427,
  frac_train = 0.8,
  train_idxs = NULL,
  test_idxs = NULL,
  filter_model = TRUE,
  min_diff = 0.1,
  prego_min_diff = min_diff,
  prego_sample_for_kmers = TRUE,
  prego_sample_fraction = 0.1,
  prego_energy_norm_quantile = 1,
  prego_spat_bin_size = NULL,
  prego_spat_num_bins = NULL,
  r2_threshold = 0.0005,
  bits_threshold = 1.75,
  filter_sample_frac = 0.1,
  include_interactions = FALSE,
  interaction_threshold = 0.001,
  max_motif_interaction_n = NULL,
  max_add_interaction_n = NULL,
  max_n_interactions = NULL,
  strategy = c("progressive", "single"),
  interaction_thresholds = c(0.01, 0.0005),
  interaction_only_sig_motifs = c(TRUE, FALSE),
  interaction_only_sig_add_motifs = c(TRUE, TRUE),
  interaction_scale_factor = 1,
  interaction_min_signal_correlation = NULL,
  n_prego_motifs = 0,
  n_cores = NULL,
  output_dir = NULL,
  plot_report = TRUE,
  rename_motifs = TRUE,
  ...
) {
    peak_intervals <- peak_intervals %||% peaks
    strategy <- match.arg(strategy)
    supplied_atac_diff <- !is.null(atac_diff)

    if (!is.null(n_cores)) {
        cli::cli_alert_info("Setting the number of cores to {.val {n_cores}}")
        if (!is.null(getOption("prego.parallel.nc"))) {
            cli::cli_alert("(Previous number of cores was {.val {getOption('prego.parallel.nc')}})")
            withr::defer(prego::set_parallel(getOption("prego.parallel.nc")))
        }
        prego::set_parallel(n_cores)
    }

    if (add_sequences_features) {
        cli::cli_alert("Computing sequence features")
        seq_feats <- create_sequence_features(peak_intervals, peaks_size)
        if (is.null(additional_features)) {
            additional_features <- seq_feats
        } else {
            seq_feats <- seq_feats[, setdiff(colnames(seq_feats), colnames(additional_features)), drop = FALSE]
            if (ncol(seq_feats) > 0) {
                cli::cli_alert("Added the following sequence features: {.val {colnames(seq_feats)}}")
                additional_features <- cbind(additional_features, seq_feats)
            }
        }
    }

    if (!is.null(atac_diff)) {
        # create fake atac scores with bin1 = 0
        atac_scores <- matrix(0, nrow = length(atac_diff), ncol = 2)
        atac_scores[, 2] <- atac_diff
        normalize_bins <- FALSE
    }

    if (!is.null(min_tss_distance)) {
        validate_misha()
    }
    tss_filter <- get_tss_distance_filter(peak_intervals, min_tss_distance)
    kept_idxs <- which(tss_filter)
    if (sum(!tss_filter) > 0) {
        cli::cli_alert_info("Filtering out {.val {sum(!tss_filter)}} peaks that are within {.val {min_tss_distance}}bp of a TSS")
    }
    peak_intervals <- peak_intervals[tss_filter, ]
    if (!is.null(atac_scores)) {
        atac_scores <- atac_scores[tss_filter, , drop = FALSE]
    }
    if (!is.null(atac_diff)) {
        atac_diff <- atac_diff[tss_filter]
    }
    if (!is.null(motif_energies)) {
        motif_energies <- motif_energies[tss_filter, , drop = FALSE]
    }
    if (!is.null(additional_features)) {
        additional_features <- additional_features[tss_filter, , drop = FALSE]
    }
    if (!is.null(train_idxs)) {
        train_idxs <- match(intersect(train_idxs, kept_idxs), kept_idxs)
        train_idxs <- train_idxs[!is.na(train_idxs)]
    }
    if (!is.null(test_idxs)) {
        test_idxs <- match(intersect(test_idxs, kept_idxs), kept_idxs)
        test_idxs <- test_idxs[!is.na(test_idxs)]
    }

    cli::cli_alert_info("Seed: {.val {seed}}")
    set.seed(seed)
    n_intervals <- nrow(peak_intervals)
    train_idxs <- train_idxs %||% sample(1:n_intervals, frac_train * n_intervals)
    test_idxs <- test_idxs %||% setdiff(1:n_intervals, train_idxs)

    # make sure that train and test indices are disjoint
    if (length(intersect(train_idxs, test_idxs)) > 0) {
        cli::cli_abort("Train and test indices must be disjoint")
    }

    frac_train <- length(train_idxs) / n_intervals

    cli::cli_alert_info("Training on {.val {length(train_idxs)}} intervals ({scales::percent(frac_train)}) and testing on {.val {length(test_idxs)}} intervals ({scales::percent(1 - frac_train)})")

    if (!is.null(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        cli::cli_alert_info("Saving the train and test indices to {.val {output_dir}}")
        train_test <- ifelse(1:n_intervals %in% train_idxs, "train", "test")
        readr::write_csv(peak_intervals %>% mutate(type = train_test), file.path(output_dir, "train_test_indices.csv"))
    }

    if (is.null(traj_prego) && n_prego_motifs > 0) {
        cli::cli_alert_info("Learning prego motifs (de-novo motifs)")
        if (is.null(atac_diff)) {
            atac_diff <- atac_scores[, 2] - atac_scores[, 1]
        }

        symmetrize_spat <- list(...)$symmetrize_spat %||% TRUE
        traj_prego <- learn_traj_prego(peak_intervals[train_idxs, ], atac_diff[train_idxs],
            n_motifs = n_prego_motifs, min_diff = prego_min_diff,
            sample_for_kmers = prego_sample_for_kmers,
            sample_fraction = prego_sample_fraction, energy_norm_quantile = prego_energy_norm_quantile, norm_intervals = norm_intervals %||% peak_intervals, seed = seed, spat_bin_size = prego_spat_bin_size, spat_num_bins = prego_spat_num_bins, peaks_size = peaks_size, symmetrize_spat = symmetrize_spat
        )
        if (!is.null(output_dir)) {
            out_file <- file.path(output_dir, "prego_model.rds")
            cli::cli_alert("Saving the prego model to {.val {out_file}}")
            readr::write_rds(traj_prego, out_file)
        }
    }


    infer_and_save <- function(traj_model, file = NULL) {
        # Compute test diff_score with correct bin selection and normalization
        test_atac <- as.matrix(atac_scores[test_idxs, , drop = FALSE])
        test_bin_start <- bin_start
        test_bin_end <- bin_end %||% ncol(test_atac)
        if (normalize_bins) {
            test_atac[, test_bin_start] <- norm01(test_atac[, test_bin_start])
            test_atac[, test_bin_end] <- norm01(test_atac[, test_bin_end])
        }
        test_diff <- test_atac[, test_bin_end] - test_atac[, test_bin_start]

        traj_model_all <- infer_trajectory_motifs(traj_model, peak_intervals[test_idxs, ],
            additional_features = additional_features[test_idxs, ],
            diff_score = test_diff
        )
        cli::cli_text("\n")
        cli::cli_text("Number of motifs: {.val {length(traj_model_all@motif_models)}}")
        cli::cli_text("R^2 train: {.val {round(cor(traj_model_all@diff_score[traj_model_all@type == 'train'], traj_model_all@predicted_diff_score[traj_model_all@type == 'train'], use = 'pairwise.complete.obs')^2, digits = 3)}}")
        cli::cli_text("R^2 test: {.val {round(cor(traj_model_all@diff_score[traj_model_all@type == 'test'], traj_model_all@predicted_diff_score[traj_model_all@type == 'test'], use = 'pairwise.complete.obs')^2, digits = 3)}}")
        cli::cli_text("\n")
        if (!is.null(file) && !is.null(output_dir)) {
            out_file <- file.path(output_dir, file)
            cli::cli_alert("Saving the model to {.val {out_file}}")
            traj_model_all <- strip_traj_model(traj_model_all)
            readr::write_rds(traj_model_all, out_file)

            if (plot_report) {
                plot_traj_model_report(traj_model_all, filename = file.path(output_dir, gsub("\\.rds$", ".pdf", file)))
            }
        }
        traj_model_all
    }

    cli::cli_alert("Regressing on train set")

    traj_model <- regress_trajectory_motifs(
        peak_intervals = peak_intervals[train_idxs, ],
        atac_scores = atac_scores[train_idxs, ],
        normalize_bins = normalize_bins,
        norm_intervals = norm_intervals,
        motif_energies = motif_energies[train_idxs, ],
        additional_features = additional_features[train_idxs, ],
        traj_prego = traj_prego,
        peaks_size = peaks_size,
        bin_start = bin_start,
        bin_end = bin_end,
        max_motif_num = max_motif_num,
        seed = seed,
        min_diff = min_diff,
        min_tss_distance = NULL,
        ...
    )

    if (rename_motifs) {
        nm <- match_traj_model_motif_names(traj_model)
        traj_model <- rename_motif_models(traj_model, nm)
    }

    final_model <- infer_and_save(traj_model, "iq_regression_model.rds")

    if (filter_model) {
        cli::cli_alert("Filtering the model")
        traj_model <- filter_traj_model(traj_model, r2_threshold = r2_threshold, bits_threshold = bits_threshold, sample_frac = filter_sample_frac, seed = seed)
        final_model <- infer_and_save(traj_model, "iq_regression_filtered_model.rds")
    }

    if (include_interactions) {
        if (strategy == "single") {
            cli::cli_alert("Using single-pass interaction strategy (paper / 0.0.6-compatible).")
            traj_model <- add_interactions(
                traj_model,
                interaction_threshold = interaction_threshold,
                max_motif_n = max_motif_interaction_n,
                max_add_n = max_add_interaction_n,
                max_n = max_n_interactions,
                interaction_scale_factor = interaction_scale_factor,
                min_signal_correlation = interaction_min_signal_correlation,
                seed = seed
            )
        } else {
            # Progressive. Use default_score_split_features as the between-pass
            # builder when multi-bin atac_scores are available; degrade to a
            # tight single-pass when they aren't (e.g. user supplied atac_diff).
            train_atac <- if (!supplied_atac_diff) atac_scores[train_idxs, , drop = FALSE] else NULL
            resolved_bin_end <- bin_end %||% ifelse(is.null(train_atac), 2L, ncol(train_atac))
            has_bins <- !is.null(train_atac) &&
                ncol(train_atac) >= 2 &&
                !identical(bin_start, resolved_bin_end)
            builder <- if (has_bins) {
                function(tm) {
                    default_score_split_features(
                        tm,
                        atac_scores = train_atac,
                        bin_start = bin_start,
                        bin_end = resolved_bin_end
                    )
                }
            } else {
                cli::cli_alert_info(
                    "Progressive interaction strategy: no multi-bin {.field atac_scores} available, running single tight pass at threshold {.val {interaction_thresholds[1]}}."
                )
                NULL
            }
            if (is.null(builder) && length(interaction_thresholds) > 1) {
                # Without a builder, the tail passes are redundant. Collapse to a single tight pass.
                thresholds <- interaction_thresholds[1]
                only_sig_motifs <- interaction_only_sig_motifs[1]
                only_sig_add_motifs <- interaction_only_sig_add_motifs[1]
            } else {
                thresholds <- interaction_thresholds
                only_sig_motifs <- interaction_only_sig_motifs
                only_sig_add_motifs <- interaction_only_sig_add_motifs
            }
            traj_model <- add_interactions_progressive(
                traj_model,
                thresholds = thresholds,
                only_sig_motifs = only_sig_motifs,
                only_sig_add_motifs = only_sig_add_motifs,
                additional_features_builder = builder,
                interaction_scale_factor = interaction_scale_factor,
                min_signal_correlation = interaction_min_signal_correlation,
                max_motif_n = max_motif_interaction_n,
                max_add_n = max_add_interaction_n,
                max_n = max_n_interactions,
                seed = seed
            )
        }
        final_model <- infer_and_save(traj_model, "iq_regression_final_model.rds")
    }


    cli::cli_alert_success("Finished IQ regression")
    cli::cli_text("\n")
    cli::cli_text("Run {.code plot_traj_model_report(traj_model)} to visualize the model features")
    cli::cli_text("Run {.code plot_prediction_scatter(traj_model)} to visualize the model predictions")
    return(final_model)
}
