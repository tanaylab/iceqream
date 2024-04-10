#' Relearn Trajectory Model
#'
#' This function relearns a trajectory model using the glmnet package.
#'
#' @param traj_model The trajectory model object.
#' @param new_energies If TRUE - recreate the energies. Default is FALSE.
#' @param new_logist If TRUE - recreate the logistic features. Default is FALSE. If \code{new_energies} is TRUE, this is automatically set to TRUE.
#' @param lambda The lambda value to use for relearning. If NULL, the lambda value from the trajectory model is used.
#' @param verbose Logical indicating whether to display additional information.
#' @return The updated trajectory model object.
#'
#' @export
relearn_traj_model <- function(traj_model, new_energies = FALSE, new_logist = FALSE, lambda = NULL, verbose = FALSE) {
    if (new_energies) {
        traj_model@normalized_energies <- calc_traj_model_energies(traj_model)
        new_logist <- TRUE
    }

    if (new_logist) {
        X <- cbind(
            create_logist_features(traj_model@normalized_energies),
            traj_model@validate_additional_features
        )
    } else {
        X <- traj_model@model_features
    }

    if (is.null(lambda)) {
        lambda <- traj_model@params$lambda
    }

    X <- traj_model@model_features
    y <- norm01(traj_model@diff_score)

    X_train <- X[traj_model@type == "train", ]
    y_train <- y[traj_model@type == "train"]

    model <- glmnet::glmnet(X_train, y_train, binomial(link = "logit"), alpha = traj_model@params$alpha, lambda = lambda, seed = traj_model@params$seed)

    pred <- logist(glmnet::predict.glmnet(model, newx = X, type = "link", s = traj_model@params$lambda))[, 1]
    pred <- norm01(pred)
    pred <- rescale(pred, traj_model@diff_score)
    r2_f <- cor(pred, y)^2
    r2_t <- cor(pred[traj_model@type == "test"], y[traj_model@type == "test"])^2
    if (verbose) {
        cli_alert_info("R^2 all after relearning: {.val {r2_f}}")
        cli_alert_info("R^2 test after relearning: {.val {r2_t}}")
    }

    traj_model@model <- model
    traj_model@predicted_diff_score <- pred
    traj_model@coefs <- get_model_coefs(model)
    traj_model@model_features <- X

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

sample_model <- function(traj_model, sample_frac = 0.1, seed = 60427, verbose = FALSE) {
    traj_model_s <- traj_model
    idxs <- prego::sample_quantile_matched_rows(as.data.frame(traj_model_s@model_features) %>% mutate(id = 1:n()), traj_model_s@diff_score, sample_frac = sample_frac, num_quantiles = 10, seed = seed, verbose = FALSE) %>% pull(id)
    traj_model_s@type <- rep("test", nrow(traj_model_s@model_features))
    traj_model_s@type[idxs] <- "train"

    traj_model_s <- relearn_traj_model(traj_model_s, verbose = FALSE)
    if (verbose) {
        cli_alert_info("Using {.val {length(idxs)}} samples for filtering")
    }
    return(traj_model_s)
}

#' Add Features R^2
#'
#' This function adds the added R-squared values of each feature to the trajectory model.
#'
#' @param traj_model The trajectory model object.
#' @param sample_frac The fraction of samples to use for computing the r2 without each model. When NULL, all samples are used.
#'
#' @return The trajectory model object with the added R-squared values.
#'
#' @export
add_features_r2 <- function(traj_model, sample_frac = 0.1, seed = 60427) {
    traj_model_full <- traj_model
    if (!is.null(sample_frac)) {
        traj_model <- sample_model(traj_model, sample_frac = sample_frac, seed = seed, verbose = TRUE)
    }

    motif_models <- names(traj_model@motif_models)
    full_model_r2 <- cor(traj_model@predicted_diff_score, traj_model@diff_score)^2
    var_stats <- plyr::llply(motif_models, function(var) {
        pssm <- traj_model@motif_models[[var]]$pssm
        traj_model_f_var <- remove_motif_models_from_traj(traj_model, var, verbose = FALSE)
        bits <- sum(prego::bits_per_pos(pssm), na.rm = TRUE)
        r2 <- cor(traj_model_f_var@predicted_diff_score, traj_model_f_var@diff_score)^2
        cli::cli_alert("R^2 added by {.field {var}} ({.strong {prego::consensus_from_pssm(pssm)}}): {.val {full_model_r2 - r2}}. Bits: {.val {bits}}")
        list(r2 = r2, bits = bits)
    })
    vars_r2 <- purrr::map_dbl(var_stats, ~ .x$r2)
    names(vars_r2) <- motif_models

    traj_model_full@features_r2 <- full_model_r2 - vars_r2

    return(traj_model_full)
}

calc_features_bits <- function(traj_model) {
    bits <- purrr::map_dbl(traj_model@motif_models, ~ sum(prego::bits_per_pos(.x$pssm), na.rm = TRUE))
    names(bits) <- names(traj_model@motif_models)
    return(bits)
}

add_traj_model_stats <- function(traj_model) {
    obs <- traj_model@diff_score
    pred <- traj_model@predicted_diff_score
    if (sum(names(obs) == "") == 0) {
        pred <- pred[names(obs)]
    }

    names(obs) <- NULL
    names(pred) <- NULL
    train_idxs <- which(traj_model@type == "train")
    test_idxs <- which(traj_model@type == "test")
    traj_model@params$stats <- list(
        r2_train = cor(pred[train_idxs], obs[train_idxs], use = "pairwise.complete.obs")^2
    )
    if (length(test_idxs) > 0) {
        traj_model@params$stats$r2_test <- cor(pred[test_idxs], obs[test_idxs], use = "pairwise.complete.obs")^2
    }

    traj_model@params$stats$r2_all <- cor(pred, obs, use = "pairwise.complete.obs")^2
    return(traj_model)
}

traj_model_has_test <- function(traj_model) {
    return(any(traj_model@type == "test"))
}

split_traj_model_to_train_test <- function(traj_model) {
    if (!traj_model_has_test(traj_model)) {
        return(list(train = traj_model, test = NULL))
    }

    train_idxs <- traj_model@type == "train"
    test_idxs <- traj_model@type == "test"

    traj_model_train <- traj_model
    traj_model_train@model_features <- traj_model@model_features[train_idxs, ]
    traj_model_train@normalized_energies <- traj_model@normalized_energies[train_idxs, ]
    traj_model_train@diff_score <- traj_model@diff_score[train_idxs]
    traj_model_train@predicted_diff_score <- traj_model@predicted_diff_score[train_idxs]
    traj_model_train@type <- traj_model@type[train_idxs]
    stopifnot(all(traj_model_train@type == "train"))
    traj_model_train@peak_intervals <- traj_model@peak_intervals[train_idxs, ]
    if (ncol(traj_model@additional_features) > 0) {
        traj_model_train@additional_features <- traj_model@additional_features[train_idxs, ]
    }

    traj_model_test <- traj_model
    traj_model_test@model_features <- traj_model@model_features[test_idxs, ]
    traj_model_test@normalized_energies <- traj_model@normalized_energies[test_idxs, ]
    traj_model_test@diff_score <- traj_model@diff_score[test_idxs]
    traj_model_test@predicted_diff_score <- traj_model@predicted_diff_score[test_idxs]
    traj_model_test@type <- traj_model@type[test_idxs]
    stopifnot(all(traj_model_test@type == "test"))
    traj_model_test@peak_intervals <- traj_model@peak_intervals[test_idxs, ]
    if (ncol(traj_model@additional_features) > 0) {
        traj_model_test@additional_features <- traj_model@additional_features[test_idxs, ]
    }

    return(list(train = traj_model_train, test = traj_model_test))
}
