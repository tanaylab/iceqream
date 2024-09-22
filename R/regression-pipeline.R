#' Perform IQ regression on peak intervals
#'
#' Perform IQ regression on peak intervals using the provided ATAC-seq scores, ATAC-seq score differences, normalized intervals, motif energies, and additional features, after dividing the intervals into training and testing sets.
#'
#' @param frac_train  A numeric value indicating the fraction of intervals to use for training (default is 0.8).
#' @param filter_model A logical value indicating whether to filter the model (default is TRUE).
#'
#'
#' @inheritParams regress_trajectory_motifs
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
    max_motif_num = 30,
    traj_prego = NULL,
    peaks_size = 300,
    bin_start = 1,
    bin_end = NULL,
    seed = 60427,
    frac_train = 0.8,
    filter_model = TRUE,
    ...) {
    n_intervals <- nrow(peak_intervals)
    train_idxs <- sample(1:n_intervals, frac_train * n_intervals)
    test_idxs <- setdiff(1:n_intervals, train_idxs)

    cli::cli_alert_info("Training on {.val {length(train_idxs)}} intervals ({scales::percent(frac_train)}) and testing on {.val {length(test_idxs)}} intervals ({scales::percent(1 - frac_train)})")

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
        ...
    )

    if (filter_model) {
        cli::cli_alert("Filtering the model")
        traj_model <- filter_traj_model(traj_model, ...)
    }

    cli::cli_alert("Infering trajectory motifs on the test set")
    traj_model_all <- infer_trajectory_motifs(traj_model, peak_intervals[test_idxs, ],
        additional_features = additional_features[test_idxs, ],
        atac_scores = atac_scores[test_idxs, ]
    )

    cli_alert_success("Finished IQ regression.")
    cli::cli_text("\n")
    cli::cli_text("Number of motifs: {.val {length(traj_model_all@motif_models)}}")
    cli::cli_text("R^2 train: {.val {round(cor(traj_model_all@diff_score[traj_model_all@type == 'train'], traj_model_all@predicted_diff_score[traj_model_all@type == 'train'], use = 'pairwise.complete.obs')^2, digits = 3)}}")
    cli::cli_text("R^2 test: {.val {round(cor(traj_model_all@diff_score[traj_model_all@type == 'test'], traj_model_all@predicted_diff_score[traj_model_all@type == 'test'], use = 'pairwise.complete.obs')^2, digits = 3)}}")
    cli::cli_text("\n")
    cli::cli_text("Run {.code plot_traj_model_report(traj_model)} to visualize the model features")
    cli::cli_text("Run {.code plot_prediction_scatter(traj_model)} to visualize the model predictions")
    return(traj_model_all)
}
