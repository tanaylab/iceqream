#' Infer trajectory motifs using a pre-trained trajectory model
#'
#' @description This function infers the motif energies for a set of peaks using a pre-trained trajectory model.
#'
#' @param traj_model A trajectory model object, as returned by \code{regress_trajectory_motifs}
#' @inheritParams infer_trajectory_motifs
#'
#' @return a `TrajectoryModel` object which contains both the original ('train') peaks and the newly inferred ('test') peaks. The field `@type` indicates whether a peak is a 'train' or 'test' peak.
#'
#' @export
infer_trajectory_motifs <- function(traj_model, peak_intervals, atac_scores = NULL, bin_start = 1, bin_end = ncol(atac_scores), additional_features = NULL) {
    validate_traj_model(traj_model)
    validate_additional_features(additional_features, peak_intervals)
    if (ncol(traj_model@additional_features) > 0) {
        if (is.null(additional_features)) {
            additional_features <- matrix(0, nrow = nrow(peak_intervals), ncol = length(traj_model@additional_features))
            colnames(additional_features) <- traj_model@additional_features
            cli_warn("No additional features were provided. Using 0 for all features. The following features are needed: {.val {traj_model@additional_features}}")
        } else {
            for (feat in colnames(traj_model@additional_features)) {
                if (!(feat %in% colnames(additional_features))) {
                    cli_warn("Additional feature {.val {feat}} is missing. Using 0 for this feature.")
                    additional_features[, feat] <- 0
                }
            }
        }
    }

    withr::local_options(list(gmax.data.size = 1e9))

    all_intervals <- bind_rows(
        traj_model@peak_intervals %>% mutate(type = "train"),
        peak_intervals %>% mutate(type = "test")
    ) %>%
        select(chrom, start, end, type)

    intervals_unique <- all_intervals %>%
        distinct(chrom, start, end) %>%
        mutate(id = 1:n())

    cli_alert_info("Extracting sequences...")
    sequences <- toupper(misha::gseq.extract(intervals_unique))

    cli_alert_info("Computing motif energies for {.val {nrow(intervals_unique)}} intervals (train: {.val {sum(all_intervals$type == 'train')}}, test: {.val {sum(all_intervals$type == 'test')}})")
    clust_energies <- plyr::llply(purrr::discard(traj_model@motif_models, is.null), function(x) {
        prego::compute_pwm(sequences, x$pssm, spat = x$spat, spat_min = x$spat_min %||% 1, spat_max = x$spat_max)
    }, .parallel = TRUE)

    clust_energies <- do.call(cbind, clust_energies)

    clust_energies <- apply(clust_energies, 2, norm_energy, min_energy = -10, q = traj_model@params$energy_norm_quantile)

    idxs <- peak_intervals %>%
        select(chrom, start, end) %>%
        left_join(intervals_unique, by = c("chrom", "start", "end")) %>%
        pull(id)

    e_test <- clust_energies[idxs, , drop = FALSE]

    if (!is.null(additional_features)) {
        additional_features[is.na(additional_features)] <- 0
        e_test <- cbind(e_test, additional_features)
    }

    e_test_logist <- create_logist_features(e_test)
    e_test_logist <- e_test_logist[, colnames(traj_model@model_features), drop = FALSE]

    predicted_diff_score <- logist(glmnet::predict.glmnet(traj_model@model, newx = e_test_logist, type = "link", s = traj_model@params$lambda))[, 1]
    # predicted_diff_score <- (predicted_diff_score * max(traj_model@diff_score)) + min(traj_model@diff_score)

    predicted_diff_score <- rescale(predicted_diff_score, min(traj_model@diff_score), max(traj_model@diff_score))

    traj_model@model_features <- rbind(traj_model@model_features, e_test_logist[, intersect(colnames(e_test_logist), colnames(traj_model@model_features))])
    traj_model@normalized_energies <- rbind(traj_model@normalized_energies, e_test[, intersect(colnames(e_test), colnames(traj_model@normalized_energies))])
    if (!is.null(atac_scores)) {
        traj_model@diff_score <- c(traj_model@diff_score, atac_scores[, bin_end] - atac_scores[, bin_start])
    }
    traj_model@predicted_diff_score <- c(traj_model@predicted_diff_score, predicted_diff_score)
    traj_model@type <- c(traj_model@type, rep("test", nrow(e_test_logist)))
    traj_model@peak_intervals <- bind_rows(traj_model@peak_intervals, peak_intervals)
    if (!is.null(additional_features)) {
        traj_model@additional_features <- bind_rows(traj_model@additional_features, as.data.frame(additional_features))
    }


    return(traj_model)
}
