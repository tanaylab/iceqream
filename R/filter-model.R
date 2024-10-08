#' Filter a trajectory model using punctuated regression.
#'
#' @description Run a regression without each feature and filter features that do not improve the model more
#' than \code{r2_threshold}.
#'
#' @param traj_model An instance of \code{TrajectoryModel}.
#' @param r2_threshold minimal R^2 for a feature to be included in the model.
#' @param bits_threshold minimal sum of bits for a feature to be included in the model.
#' @param seed The seed to use for reproducibility when sampling the data.
#'
#' @return An instance of \code{TrajectoryModel} with the filtered model.
#'
#' @inheritParams add_features_r2
#' @export
filter_traj_model <- function(traj_model, r2_threshold = 0.0005, bits_threshold = 1.75, sample_frac = 0.1, seed = 60427) {
    traj_model_all <- traj_model
    if (traj_model_has_test(traj_model)) {
        cli_alert_info("Using only the training set for filtering")
        traj_model <- split_traj_model_to_train_test(traj_model_all)$train
    }

    if (!is.null(sample_frac)) {
        traj_model_s <- sample_model(traj_model, sample_frac = sample_frac, seed = seed, verbose = TRUE)
    } else {
        traj_model_s <- traj_model
    }

    full_model_r2 <- cor(traj_model_s@predicted_diff_score, traj_model_s@diff_score)^2
    cli_alert_info("R^2 of the full model: {.val {full_model_r2}}")

    motif_models <- names(traj_model@motif_models)

    if (!is.null(r2_threshold)) {
        cli_alert_info("Filtering features with R^2 < {.val {r2_threshold}} and bits < {.val {bits_threshold}}")
        traj_model_s <- add_features_r2(traj_model_s, sample_frac = NULL)
        vars_r2 <- traj_model_s@features_r2
    } else {
        r2_threshold <- 0
        vars_r2 <- rep(0, length(motif_models))
        names(vars_r2) <- motif_models
    }

    vars_bits <- calc_features_bits(traj_model_s)
    if (any(vars_bits < bits_threshold)) {
        cli_alert_info("Removing the following features with bits < {.val {bits_threshold}}: {.val {motif_models[vars_bits < bits_threshold]}}")
        bit_vars_to_remove <- names(vars_bits)[vars_bits < bits_threshold]
    } else {
        bit_vars_to_remove <- c()
    }

    if (any(vars_r2 < r2_threshold)) {
        cli_alert_info("Trying to remove the following features with R^2 < {.val {r2_threshold}}: {.val {motif_models[vars_r2 < r2_threshold]}}")
        r2_vars_to_remove <- vars_r2[vars_r2 < r2_threshold]
        r2_vars_to_remove <- r2_vars_to_remove[!(names(r2_vars_to_remove) %in% bit_vars_to_remove)]
        r2_vars_to_remove <- names(sort(r2_vars_to_remove))
        vars_to_remove <- c(bit_vars_to_remove, r2_vars_to_remove[1])
        traj_model_new <- remove_motif_models_from_traj(traj_model_s, vars_to_remove, verbose = FALSE)

        # go over the rest of the features and remove them if they do not change the R^2 by more than r2_threshold
        if (length(r2_vars_to_remove) > 1) {
            r2_before <- cor(traj_model_new@predicted_diff_score, traj_model_new@diff_score)^2
            for (var in r2_vars_to_remove[2:length(r2_vars_to_remove)]) {
                traj_model_f <- remove_motif_models_from_traj(traj_model_new, var, verbose = FALSE)
                r2_after <- cor(traj_model_f@predicted_diff_score, traj_model_f@diff_score)^2
                if (r2_before - r2_after < r2_threshold) {
                    cli_alert("Removing {.val {var}} changed the R^2 by {.val {r2_before - r2_after}}")
                    traj_model_new <- traj_model_f
                    r2_before <- r2_after
                    vars_to_remove <- c(vars_to_remove, var)
                } else {
                    cli_alert("{.field Not removing} {.val {var}} (changed the R^2 by only {.val {r2_before - r2_after}}).")
                    next
                }
            }
        }
    } else {
        vars_to_remove <- c()
    }

    if (length(vars_to_remove) > 0) {
        cli_alert_info("Removed {.val {length(bit_vars_to_remove)}} features with bits < {.val {bits_threshold}}")
        cli_alert_info("Removed {.val {length(vars_to_remove) - length(bit_vars_to_remove)}} features with R^2 < {.val {r2_threshold}}")
        traj_model_new <- remove_motif_models_from_traj(traj_model_all, vars_to_remove, verbose = FALSE)
    } else {
        cli_alert_info("No features removed")
        traj_model_new <- traj_model_all
    }

    traj_model_new@features_r2 <- vars_r2[names(traj_model_new@motif_models)]
    traj_model_new@params$features_bits <- vars_bits[names(traj_model_new@motif_models)]
    traj_model_new@params$r2_threshold <- r2_threshold
    traj_model_new@params$bits_threshold <- bits_threshold

    cli_alert_success("After filtering: Number of non-zero coefficients: {.val {sum(traj_model_new@model$beta != 0)}} (out of {.val {ncol(traj_model_new@model_features)}}). R^2: {.val {cor(traj_model_new@predicted_diff_score, norm01(traj_model_new@diff_score))^2}}. Number of models: {.val {length(traj_model_new@motif_models)}}")

    return(traj_model_new)
}
