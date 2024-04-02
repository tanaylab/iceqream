#' Remove motif models from trajectory model
#'
#' This function removes specified motif models from a given trajectory model.
#' It updates the model, motif models, predicted difference score, model features,
#' coefficients, normalized energies, and features R^2 of the trajectory model.
#'
#' @param traj_model The trajectory model object to remove motif models from.
#' @param motif_models A character vector specifying the names of the motif models to remove.
#' @param verbose A logical value indicating whether to display information about the R^2 after removing the motif models. Default is TRUE.
#'
#' @return The updated trajectory model object after removing the motif models.
#'
#' @export
remove_motif_models_from_traj <- function(traj_model, motif_models, verbose = TRUE) {
    X <- traj_model@model_features
    vars <- names(traj_model@motif_models)
    if (!all(motif_models %in% vars)) {
        cli_abort("Motif{?s} {.val {motif_models[!motif_models %in% vars]}} not found in the trajectory model.")
    }
    vars_f <- vars[!(vars %in% motif_models)]

    X_f <- X[, grep(paste0("(", paste(motif_models, collapse = "|"), ")(_low-energy|_high-energy|_higher-energy|_sigmoid)"), colnames(X), invert = TRUE)]

    diff_score <- traj_model@diff_score
    y <- norm01(diff_score)
    model_f <- glmnet::glmnet(X_f, y, binomial(link = "logit"), alpha = traj_model@params$alpha, lambda = traj_model@params$lambda, parallel = FALSE, seed = traj_model@params$seed)

    pred_f <- logist(glmnet::predict.glmnet(model_f, newx = X_f, type = "link", s = traj_model@params$lambda))[, 1]
    pred_f <- norm01(pred_f)
    pred_f <- rescale(pred_f, diff_score)
    r2_f <- cor(pred_f, y)^2
    if (verbose) {
        cli_alert_info("R^2 after removing {.field {motif_models}}: {.val {r2_f}}")
    }

    traj_model@model <- model_f
    traj_model@motif_models <- traj_model@motif_models[vars_f]
    traj_model@predicted_diff_score <- pred_f
    traj_model@model_features <- X_f
    traj_model@coefs <- get_model_coefs(model_f)
    traj_model@normalized_energies <- traj_model@normalized_energies[, vars_f, drop = FALSE]
    traj_model@features_r2 <- traj_model@features_r2[vars_f[vars_f %in% names(traj_model@features_r2)]]

    return(traj_model)
}

#' Filter a trajectory model using punctuated regression.
#'
#' @description Run a regression without each feature and filter features that do not improve the model more
#' than \code{r2_threshold}.
#'
#' @param traj_model An instance of \code{TrajectoryModel}.
#' @param r2_threshold minimal R^2 for a feature to be included in the model.
#' @param bits_threshold minimal sum of bits for a feature to be included in the model.
#'
#' @return An instance of \code{TrajectoryModel} with the filtered model.
#'
#' @export
filter_traj_model <- function(traj_model, r2_threshold = 0.0005, bits_threshold = 1.75) {
    full_model <- traj_model@model
    full_model_r2 <- cor(traj_model@predicted_diff_score, traj_model@diff_score)^2

    motif_models <- names(traj_model@motif_models)

    var_stats <- plyr::llply(motif_models, function(var) {
        pssm <- traj_model@motif_models[[var]]$pssm
        traj_model_f_var <- remove_motif_models_from_traj(traj_model, var, verbose = FALSE)
        bits <- sum(prego::bits_per_pos(pssm))
        r2 <- cor(traj_model_f_var@predicted_diff_score, traj_model_f_var@diff_score)^2
        cli::cli_alert("R^2 added by {.field {var}} ({.strong {prego::consensus_from_pssm(pssm)}}): {.val {full_model_r2 - r2}}. Bits: {.val {bits}}")
        if (full_model_r2 - r2 < r2_threshold) {
            cli::cli_alert_info("Variable {.field {var}} removed due to low R^2")
        }
        if (bits < bits_threshold) {
            cli::cli_alert_info("Variable {.field {var}} removed due to low information content")
        }
        list(r2 = r2, bits = bits)
    })
    vars_r2 <- purrr::map_dbl(var_stats, ~ .x$r2)
    names(vars_r2) <- motif_models

    vars_bits <- purrr::map_dbl(var_stats, ~ .x$bits)
    names(vars_bits) <- motif_models

    f_bits <- vars_bits > bits_threshold
    f_r2 <- (full_model_r2 - vars_r2) > r2_threshold

    vars_to_remove <- motif_models[!(f_bits & f_r2)]
    if (length(vars_to_remove) > 0) {
        cli_alert_info("Removed {.val {sum(!f_bits)}} features with bits < {.val {bits_threshold}}")
        cli_alert_info("Removed {.val {sum(!f_r2)}} features with R^2 < {.val {r2_threshold}}")
        traj_model_new <- remove_motif_models_from_traj(traj_model, vars_to_remove)
    } else {
        cli_alert_info("No features removed")
        traj_model_new <- traj_model
    }

    traj_model@features_r2 <- vars_r2
    traj_model@params$features_bits <- vars_bits
    traj_model@params$r2_threshold <- r2_threshold
    traj_model@params$bits_threshold <- bits_threshold

    cli_alert_success("After filtering: Number of non-zero coefficients: {.val {sum(traj_model@model$beta != 0)}} (out of {.val {ncol(traj_model@model_features)}}). R^2: {.val {cor(traj_model@predicted_diff_score, norm01(traj_model@diff_score))^2}}")

    return(traj_model_new)
}

filter_model_using_coefs <- function(X, coefs, diff_score, alpha, lambda, seed, full_model, n_motifs, ignore_variables = NULL) {
    y <- norm01(diff_score)
    variables <- coefs$variable
    if (!is.null(ignore_variables)) {
        variables <- variables[!(variables %in% ignore_variables)]
    }
    coefs_max <- coefs %>%
        tibble::column_to_rownames("variable") %>%
        .[variables, ] %>%
        apply(1, max) %>%
        sort(decreasing = TRUE)

    vars_f <- names(coefs_max)[1:n_motifs]

    X_f <- X[, grep(paste0("(", paste(c(vars_f, ignore_variables), collapse = "|"), ").+"), colnames(X))]
    cli_alert_info("Number of features left: {.val {length(vars_f)}}")

    model_f <- glmnet::glmnet(X_f, y, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = FALSE, seed = seed)
    pred_f <- logist(glmnet::predict.glmnet(model_f, newx = X_f, type = "link", s = lambda))[, 1]
    pred_f <- norm01(pred_f)
    pred_f <- rescale(pred_f, diff_score)
    r2_f <- cor(pred_f, y)^2
    cli_alert_info("R^2 after filtering: {.val {r2_f}}")

    return(list(model = model_f, pred = pred_f, X = X_f, r2 = r2_f, vars = vars_f))
}


filter_traj_model_using_coefs <- function(traj_model, n_motifs) {
    res <- filter_model_using_coefs(traj_model@model_features, traj_model@coefs, traj_model@diff_score, traj_model@params$alpha, traj_model@params$lambda, traj_model@params$seed, traj_model@model, ignore_variables = colnames(traj_model@additional_features), n_motifs = n_motifs)

    traj_model@model <- res$model
    traj_model@predicted_diff_score <- res$pred
    traj_model@model_features <- res$X
    traj_model@coefs <- get_model_coefs(res$model)
    traj_model@normalized_energies <- traj_model@normalized_energies[, res$vars, drop = FALSE]

    cli_alert_success("After filtering: Number of non-zero coefficients: {.val {sum(traj_model@model$beta != 0)}} (out of {.val {ncol(traj_model@model_features)}}). R^2: {.val {cor(traj_model@predicted_diff_score, norm01(traj_model@diff_score))^2}}")

    return(traj_model)
}
