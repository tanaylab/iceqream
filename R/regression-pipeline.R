#' Perform IQ regression on peak intervals
#'
#' Perform IQ regression on peak intervals using the provided ATAC-seq scores, ATAC-seq score differences, normalized intervals, motif energies, and additional features, after dividing the intervals into training and testing sets.
#'
#' @param frac_train  A numeric value indicating the fraction of intervals to use for training (default is 0.8).
#' @param filter_model A logical value indicating whether to filter the model (default is TRUE).
#' @param filter_sample_frac The fraction of samples to use for computing the r2 without each model at the filtering step. When NULL, all samples are used.
#' @param n_cores The number of cores to use for parallel processing. When NULL, the number of threads is automatically determined as 80% of the available cores. See \code{\link{prego::set_parallel}} for more details.
#' @param add_sequences_features Add CG content and dinuceotide content to the additional features.
#' @param train_idxs A vector of indices to use for training. If NULL, the training set is randomly selected.
#' @param test_idxs A vector of indices to use for testing. If NULL, the testing set is the complement of the training set.
#' @param output_dir A directory to save intermediate results. If not NULL, the train and test indices are saved to a CSV file, together with the models at each step (before filtering, after filtering, and after adding interactions). The models are saved as RDS files. If the directory exists, the files are overwritten.
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
    peak_intervals,
    atac_scores = NULL,
    atac_diff = NULL,
    normalize_bins = TRUE,
    norm_intervals = NULL,
    motif_energies = NULL,
    additional_features = NULL,
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
    r2_threshold = 0.0005,
    bits_threshold = 1.75,
    filter_sample_frac = 0.1,
    include_interactions = FALSE,
    interaction_threshold = 0.001,
    max_motif_interaction_n = NULL,
    max_add_interaction_n = NULL,
    n_cores = NULL,
    output_dir = NULL,
    ...) {
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

    infer_and_save <- function(traj_model, file = NULL) {
        traj_model_all <- infer_trajectory_motifs(traj_model, peak_intervals[test_idxs, ],
            additional_features = additional_features[test_idxs, ],
            atac_scores = atac_scores[test_idxs, ]
        )
        cli::cli_text("\n")
        cli::cli_text("Number of motifs: {.val {length(traj_model_all@motif_models)}}")
        cli::cli_text("R^2 train: {.val {round(cor(traj_model_all@diff_score[traj_model_all@type == 'train'], traj_model_all@predicted_diff_score[traj_model_all@type == 'train'], use = 'pairwise.complete.obs')^2, digits = 3)}}")
        cli::cli_text("R^2 test: {.val {round(cor(traj_model_all@diff_score[traj_model_all@type == 'test'], traj_model_all@predicted_diff_score[traj_model_all@type == 'test'], use = 'pairwise.complete.obs')^2, digits = 3)}}")
        cli::cli_text("\n")
        if (!is.null(file) && !is.null(output_dir)) {
            out_file <- file.path(output_dir, file)
            cli::cli_alert("Saving the model to {.val {out_file}}")
            readr::write_rds(traj_model_all, out_file)
        }
        traj_model_all
    }

    cli::cli_alert("Regressing on train set")

    traj_model <- regress_trajectory_motifs(
        peak_intervals = peak_intervals[train_idxs, ],
        atac_scores = atac_scores[train_idxs, ],
        atac_diff = atac_diff[train_idxs, ],
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
        ...
    )

    final_model <- infer_and_save(traj_model, "iq_regression_model.rds")

    if (filter_model) {
        cli::cli_alert("Filtering the model")
        traj_model <- filter_traj_model(traj_model, r2_threshold = r2_threshold, bits_threshold = bits_threshold, sample_frac = filter_sample_frac, seed = seed)
        final_model <- infer_and_save(traj_model, "iq_regression_filtered_model.rds")
    }

    if (include_interactions) {
        traj_model <- add_interactions(traj_model, interaction_threshold = interaction_threshold, max_motif_n = max_motif_interaction_n, max_add_n = max_add_interaction_n, seed = seed)
        final_model <- infer_and_save(traj_model, "iq_regression_final_model.rds")
    }


    cli::cli_alert_success("Finished IQ regression")
    cli::cli_text("\n")
    cli::cli_text("Run {.code plot_traj_model_report(traj_model)} to visualize the model features")
    cli::cli_text("Run {.code plot_prediction_scatter(traj_model)} to visualize the model predictions")
    return(final_model)
}
