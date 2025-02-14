#' Merge Trajectory Motifs
#'
#' This function merges multiple motifs (chosen manually) into a new motif in a trajectory model.
#'
#' @param traj_model The trajectory model object.
#' @param motifs A character vector specifying the names of the motifs to be merged.
#' @param new_motif_name A character string specifying the name of the new merged motif.
#' @param min_diff A numeric value specifying the minimum absolute difference in ATAC scores for a motif to be considered in the merging process. Default is 0.1.
#' @param seed An integer specifying the seed for random number generation. Default is 60427.
#'
#' @return A modified trajectory model object with the merged motif and updated attributes.
#'
#' @details This function validates the trajectory model, checks if the specified motifs exist in the model, and performs the merging process. It computes the merged motif using distillation on the motifs partial response, updates the model features, and returns the modified trajectory model object.
#'
#' @export
merge_trajectory_motifs <- function(traj_model, motifs, new_motif_name, min_diff = 0.1, seed = 60427) {
    validate_traj_model(traj_model)
    if (length(motifs) <= 1) {
        cli_abort("At least two motifs are required to merge.")
    }

    if (!all(motifs %in% names(traj_model@motif_models))) {
        cli_abort("Motif{?s} {.val {motifs[!motifs %in% names(traj_model@motif_models)]}} not found in the trajectory model.")
    }
    params <- traj_model@params
    seqs <- prego::intervals_to_seq(traj_model@peak_intervals, params$peaks_size)
    norm_seqs <- prego::intervals_to_seq(traj_model@normalization_intervals, params$peaks_size)
    atac_diff <- traj_model@diff_score
    atac_diff_n <- norm01(atac_diff)

    diff_filter <- abs(atac_diff) >= min_diff
    features <- grep(paste0("(", paste(motifs, collapse = "|"), ")(_low-energy|_high-energy|_higher-energy|_sigmoid)"), colnames(traj_model@model_features), value = TRUE)


    cli_alert_info("Merging {.val {motifs}} into a new motif named {.val {new_motif_name}}...")
    distilled <- run_prego_on_clust_residuals(
        model = traj_model@model,
        feats = traj_model@model_features[diff_filter, ],
        clust_motifs = features,
        sequences = seqs[diff_filter],
        lambda = params$lambda,
        seed = seed,
        spat_num_bins = params$spat_num_bins,
        spat_bin_size = params$spat_bin_size,
        kmer_sequence_length = params$kmer_sequence_length,
        symmetrize_spat = params$symmetrize_spat %||% 1
    )

    motif_e <- prego::compute_pwm(seqs, distilled$pssm, spat = distilled$spat, spat_min = distilled$spat_min %||% 1, spat_max = distilled$spat_max, func = "logSumExp")
    motif_norm_e <- prego::compute_pwm(norm_seqs, distilled$pssm, spat = distilled$spat, spat_min = distilled$spat_min %||% 1, spat_max = distilled$spat_max, func = "logSumExp")
    motif_e <- norm_energy_dataset(motif_e, motif_norm_e, min_energy = params$min_energy, q = params$energy_norm_quantile, norm_energy_max = params$norm_energy_max)

    motif_m <- as.matrix(motif_e)
    colnames(motif_m) <- new_motif_name
    motif_e_logist <- create_logist_features(motif_m)

    feat_mat <- traj_model@model_features
    feat_mat <- feat_mat[, !colnames(feat_mat) %in% features]
    feat_mat <- cbind(feat_mat, motif_e_logist)

    model <- glmnet::glmnet(feat_mat, atac_diff_n, binomial(link = "logit"), alpha = params$alpha, lambda = params$lambda, parallel = TRUE, seed = seed)
    model <- strip_glmnet(model)
    pred <- logist(glmnet::predict.glmnet(model, newx = feat_mat, type = "link", s = params$lambda))[, 1]
    pred <- norm01(pred)
    pred <- rescale(pred, atac_diff)
    r2 <- cor(pred, atac_diff)^2

    traj_model_new <- traj_model
    traj_model_new@model <- model
    new_motif_list <- list(distilled)
    names(new_motif_list) <- new_motif_name
    traj_model_new@motif_models <- c(traj_model@motif_models, new_motif_list)
    traj_model_new@predicted_diff_score <- pred
    traj_model_new@model_features <- feat_mat
    traj_model_new@coefs <- get_model_coefs(model)
    traj_model_new@normalized_energies <- cbind(traj_model@normalized_energies, motif_m)
    traj_model_new@features_r2 <- traj_model_new@features_r2[!names(traj_model_new@features_r2) %in% features]

    cli_alert_success("Finished merging motifs. R^2: {.val {r2}} (Before: {.val {cor(traj_model@predicted_diff_score, atac_diff)^2}})")

    return(traj_model_new)
}
