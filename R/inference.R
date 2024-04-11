#' Infer trajectory motifs using a pre-trained trajectory model
#'
#' @description This function infers the motif energies for a set of peaks using a pre-trained trajectory model.
#'
#' @param traj_model A trajectory model object, as returned by \code{regress_trajectory_motifs}
#' @param test_energies An already computed matrix of motif energies for the test peaks. An advanced option to provide the energies directly.
#' @param diff_score The difference in ATAC-seq scores between the end and start of the peak. If provided, the function will ignore the atac_scores parameter.
#' @inheritParams regress_trajectory_motifs
#'
#' @return a `TrajectoryModel` object which contains both the original ('train') peaks and the newly inferred ('test') peaks. The field `@type` indicates whether a peak is a 'train' or 'test' peak. R^2 statistics are computed at `object@params$stats`.
#'
#' @export
infer_trajectory_motifs <- function(traj_model, peak_intervals, atac_scores = NULL, bin_start = 1, bin_end = ncol(atac_scores), additional_features = NULL, test_energies = NULL, diff_score = NULL) {
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

    if (is.null(test_energies)) {
        e_test <- calc_traj_model_energies(traj_model, peak_intervals)
    } else {
        if (nrow(test_energies) != nrow(peak_intervals)) {
            cli_abort("The number of rows in test_energies should be equal to the number of peaks in peak_intervals.")
        }
        if (!all(colnames(traj_model@normalized_energies) %in% colnames(test_energies))) {
            cli_abort("The columns in test_energies should match the columns in the normalized energies of the trajectory model.")
        }
        e_test <- test_energies[, intersect(colnames(test_energies), colnames(traj_model@normalized_energies))]
    }

    if (!is.null(additional_features)) {
        additional_features[is.na(additional_features)] <- 0
        e_test <- cbind(e_test, additional_features)
    }

    e_test_logist <- create_logist_features(e_test)
    e_test_logist <- e_test_logist[, colnames(traj_model@model_features), drop = FALSE]

    pred <- predict_traj_model(traj_model, e_test_logist)

    traj_model@model_features <- as.matrix(rbind(traj_model@model_features, e_test_logist[, intersect(colnames(e_test_logist), colnames(traj_model@model_features))]))
    traj_model@normalized_energies <- as.matrix(rbind(traj_model@normalized_energies, e_test[, intersect(colnames(e_test), colnames(traj_model@normalized_energies))]))
    if (!is.null(diff_score)) {
        traj_model@diff_score <- c(traj_model@diff_score, diff_score)
    } else if (!is.null(atac_scores)) {
        atac_scores <- as.data.frame(atac_scores)
        traj_model@diff_score <- c(traj_model@diff_score, atac_scores[, bin_end] - atac_scores[, bin_start])
    }
    traj_model@predicted_diff_score <- c(traj_model@predicted_diff_score, pred)
    traj_model@type <- c(traj_model@type, rep("test", nrow(e_test_logist)))
    traj_model@peak_intervals <- bind_rows(traj_model@peak_intervals, peak_intervals)
    if (!is.null(additional_features)) {
        traj_model@additional_features <- bind_rows(traj_model@additional_features, as.data.frame(additional_features))
    }

    traj_model <- add_traj_model_stats(traj_model)

    return(traj_model)
}

calc_traj_model_energies <- function(traj_model, peak_intervals = traj_model@peak_intervals, func = "logSumExp") {
    withr::local_options(list(gmax.data.size = 1e9))

    cli_alert_info("Extracting sequences...")
    if (!is.null(traj_model@params$peaks_size)) {
        sequences <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(peak_intervals, traj_model@params$peaks_size)))
    } else {
        sequences <- toupper(misha::gseq.extract(peak_intervals))
    }
    norm_sequences <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(traj_model@normalization_intervals, traj_model@params$peaks_size)))

    cli_alert_info("Computing motif energies for {.val {nrow(peak_intervals)}} intervals")
    e_test <- infer_energies(sequences, norm_sequences, traj_model@motif_models, traj_model@params$min_energy, traj_model@params$energy_norm_quantile, traj_model@params$norm_energy_max, func)

    return(e_test)
}

infer_energies <- function(sequences, norm_sequences, motif_list, min_energy, energy_norm_quantile, norm_energy_max, func = "logSumExp") {
    ml <- purrr::discard(motif_list, is.null)
    energies <- plyr::llply(ml, function(x) {
        prego::compute_pwm(sequences, x$pssm, spat = x$spat, spat_min = x$spat_min %||% 1, spat_max = x$spat_max, func = func)
    }, .parallel = TRUE)
    names(energies) <- names(ml)
    energies <- do.call(cbind, energies)

    norm_energies <- plyr::llply(ml, function(x) {
        prego::compute_pwm(norm_sequences, x$pssm, spat = x$spat, spat_min = x$spat_min %||% 1, spat_max = x$spat_max, func = func)
    }, .parallel = TRUE)
    names(norm_energies) <- names(ml)
    norm_energies <- do.call(cbind, norm_energies)

    energies <- norm_energy_matrix(energies, norm_energies, min_energy = min_energy, q = energy_norm_quantile, norm_energy_max = norm_energy_max)
    return(energies)
}

predict_traj_model <- function(traj_model, feats) {
    pred <- logist(glmnet::predict.glmnet(traj_model@model, newx = feats, type = "link", s = traj_model@params$lambda))[, 1]
    pred_train <- logist(glmnet::predict.glmnet(traj_model@model, newx = traj_model@model_features[, colnames(feats), drop = FALSE], type = "link", s = traj_model@params$lambda))[, 1]

    min_val <- min(pred_train)
    max_val <- max(pred_train)
    normalized_pred_train <- (pred_train - min_val) / (max_val - min_val)

    pred <- (pred - min_val) / (max_val - min_val)
    # pred[pred > max_val] <- max_val
    pred <- rescale(pred, traj_model@diff_score)

    return(pred)
}
