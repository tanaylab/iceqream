#' Relearn Trajectory Model
#'
#' This function relearns a trajectory model using the glmnet package.
#'
#' @param traj_model The trajectory model object.
#' @param verbose Logical indicating whether to display additional information.
#' @return The updated trajectory model object.
#'
#' @export
relearn_traj_model <- function(traj_model, verbose = FALSE) {
    X <- traj_model@model_features
    y <- norm01(traj_model@diff_score)

    X_train <- X[traj_model@type == "train", ]
    y_train <- y[traj_model@type == "train"]

    model <- glmnet::glmnet(X_train, y_train, binomial(link = "logit"), alpha = traj_model@params$alpha, lambda = traj_model@params$lambda, seed = traj_model@params$seed)

    pred <- logist(glmnet::predict.glmnet(model, newx = X, type = "link", s = traj_model@params$lambda))[, 1]
    pred <- norm01(pred)
    pred <- rescale(pred, traj_model@diff_score)
    r2_f <- cor(pred, y)^2
    if (verbose) {
        cli_alert_info("R^2 after relearning: {.val {r2_f}}")
    }

    traj_model@model <- model
    traj_model@predicted_diff_score <- pred
    traj_model@coefs <- get_model_coefs(model)

    return(traj_model)
}


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

    traj_model@model_features <- X_f
    traj_model@motif_models <- traj_model@motif_models[vars_f]
    traj_model@features_r2 <- traj_model@features_r2[vars_f[vars_f %in% names(traj_model@features_r2)]]
    traj_model@normalized_energies <- traj_model@normalized_energies[, vars_f, drop = FALSE]

    traj_model <- relearn_traj_model(traj_model, verbose = FALSE)
    r2_f <- cor(traj_model@predicted_diff_score, norm01(traj_model@diff_score))^2
    if (verbose) {
        cli_alert_info("R^2 after removing {.field {motif_models}}: {.val {r2_f}}")
    }

    return(traj_model)
}
