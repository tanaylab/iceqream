#' Relearn Trajectory Model
#'
#' This function relearns a trajectory model using the glmnet package, without re-inferring the motif models.
#'
#' @param traj_model The trajectory model object.
#' @param new_energies If TRUE - recreate the energies. Default is FALSE.
#' @param new_logist If TRUE - recreate the logistic features. Default is FALSE. If \code{new_energies} is TRUE, this is automatically set to TRUE.
#' @param lambda The lambda value to use for relearning. If NULL, the lambda value from the trajectory model is used.
#' @param use_additional_features A logical value indicating whether to use additional features. Default is TRUE.
#' @param use_motifs A logical value indicating whether to use motif models. Default is TRUE.
#' @param verbose Logical indicating whether to display additional information.
#' @param rescale_pred Logical indicating whether to rescale the predicted values. Default is TRUE.
#' @param relearn_model Logical indicating whether to relearn the model. Default is TRUE.
#' @param family The family to use for the glmnet model. Either "binomial" (default) or "gaussian" for linear regression.
#' @param new_interactions Logical indicating whether to add new interactions. Default is FALSE.
#' @param max_n_interactions The maximum number of interactions to add. Default is NULL.
#' @param use_cv Logical indicating whether to use cross-validation for lambda selection. Default is FALSE.
#' @param nfolds Number of folds for cross-validation. Default is 10.
#' @param interaction_scale_factor The factor to scale the interactions by. Default is 1 (no scaling).
#' @param logist_interactions Logical indicating whether to transform interactions to logistic functions. Default is FALSE.
#' @param logist_dinucs Logical indicating whether to transform dinucleotides to logistic functions. Default is FALSE.
#' @return The updated trajectory model object.
#'
#' @export
relearn_traj_model <- function(traj_model, new_energies = FALSE, new_logist = FALSE, lambda = NULL, use_additional_features = TRUE, use_motifs = TRUE, verbose = FALSE, rescale_pred = TRUE, relearn_model = TRUE, family = "binomial", new_interactions = FALSE, max_n_interactions = NULL, use_cv = FALSE, nfolds = 10, interaction_scale_factor = 1, logist_interactions = FALSE, logist_dinucs = FALSE) {
    if (verbose) {
        r2_train_before <- cor(traj_model@predicted_diff_score[traj_model@type == "train"], traj_model@diff_score[traj_model@type == "train"])^2
        r2_test_before <- cor(traj_model@predicted_diff_score[traj_model@type == "test"], traj_model@diff_score[traj_model@type == "test"])^2
    }
    if (new_energies) {
        traj_model@normalized_energies <- calc_traj_model_energies(traj_model, traj_model@peak_intervals)
        new_logist <- TRUE
    }

    if (new_logist) {
        if (use_additional_features) {
            add_feats <- traj_model@additional_features[rownames(traj_model@normalized_energies), ]
            if (logist_dinucs) {
                add_feats <- create_logist_features(add_feats)
            } else {
                dinucs <- c("TT", "CT", "GT", "AT", "TC", "CC", "GC", "AC", "TG", "CG", "GG", "AG", "TA", "CA", "GA", "AA")
                dinuc_feats <- traj_model@additional_features[, colnames(traj_model@additional_features) %in% dinucs]
                other_feats <- create_logist_features(traj_model@additional_features[, !colnames(traj_model@additional_features) %in% dinucs])
                add_feats <- cbind(other_feats, dinuc_feats)
            }
            if (use_motifs) {
                if (has_interactions(traj_model)) {
                    ftv_inter <- feat_to_variable(traj_model, add_type = TRUE) %>%
                        filter(type == "interaction") %>%
                        distinct(variable, term1, term2)
                    interactions <- create_specifc_terms(cbind(traj_model@normalized_energies, traj_model@additional_features), ftv_inter)
                    interactions <- interactions[, colnames(traj_model@interactions), drop = FALSE]
                    interactions <- interactions * interaction_scale_factor
                    traj_model@interactions <- interactions
                    traj_model@params$logist_interactions <- logist_interactions
                    if (logist_interactions) {
                        interactions_logist <- create_logist_features(interactions)
                    } else {
                        interactions_logist <- interactions
                    }
                    traj_model@model_features <- as.matrix(cbind(
                        create_logist_features(traj_model@normalized_energies),
                        add_feats,
                        interactions_logist
                    ))
                } else {
                    traj_model@model_features <- as.matrix(cbind(
                        create_logist_features(traj_model@normalized_energies),
                        add_feats
                    ))
                }
            } else {
                traj_model@model_features <- create_logist_features(traj_model@additional_features)
                traj_model@motif_models <- list()
                traj_model@normalized_energies <- matrix()
            }
        } else {
            if (!use_motifs) {
                cli_abort("No features to use for relearning.")
            }
            traj_model@additional_features <- data.frame()
            traj_model@model_features <- create_logist_features(traj_model@normalized_energies)
        }
    }

    if (new_interactions) {
        cli::cli_alert("Adding interactions")
        return(add_interactions(traj_model, max_n = max_n_interactions, force = TRUE, logist_interactions = logist_interactions, use_cv = use_cv, nfolds = nfolds, family = family, rescale_pred = rescale_pred))
    }

    X <- traj_model@model_features

    y <- norm01(traj_model@diff_score)

    X_train <- X[traj_model@type == "train", ]
    y_train <- y[traj_model@type == "train"]

    # Ensure data types are correct for glmnet
    X_train <- as.matrix(X_train)
    y_train <- as.numeric(y_train)

    lambda <- lambda %||% traj_model@params$lambda

    if (relearn_model) {
        if (family == "binomial") {
            family_param <- binomial(link = "logit")
        } else if (family == "gaussian") {
            family_param <- gaussian()
        } else {
            cli_abort("Family must be either 'binomial' or 'gaussian'")
        }

        if (use_cv) {
            cv_model <- glmnet::cv.glmnet(X_train, y_train, family = family_param, alpha = traj_model@params$alpha, nfolds = nfolds, seed = traj_model@params$seed)
            traj_model@params$cv_stats <- list(
                lambda.min = cv_model$lambda.min,
                lambda.1se = cv_model$lambda.1se,
                cvm = cv_model$cvm,
                cvsd = cv_model$cvsd,
                lambda = cv_model$lambda
            )

            lambda <- cv_model$lambda.min
            cli::cli_alert("Using lambda: {.val {lambda}}")
            model <- glmnet::glmnet(X_train, y_train,
                family = family_param,
                alpha = traj_model@params$alpha,
                lambda = lambda,
                seed = traj_model@params$seed
            )
            traj_model@params$lambda <- lambda
        } else {
            model <- glmnet::glmnet(X_train, y_train, family = family_param, alpha = traj_model@params$alpha, lambda = lambda, seed = traj_model@params$seed)
        }
    } else {
        model <- traj_model@model
    }

    model <- strip_glmnet(model)

    # Ensure X is a matrix for prediction
    X <- as.matrix(X)

    if (family == "binomial") {
        pred <- logist(glmnet::predict.glmnet(model, newx = X, type = "link", s = lambda))[, 1]
    } else {
        pred <- glmnet::predict.glmnet(model, newx = X, type = "response", s = lambda)[, 1]
    }

    if (rescale_pred) {
        pred <- norm01(pred)
        pred <- rescale(pred, traj_model@diff_score)
    }

    if (verbose) {
        r2_f <- cor(pred, y)^2
        r2_train <- cor(pred[traj_model@type == "train"], y[traj_model@type == "train"])^2
        r2_t <- cor(pred[traj_model@type == "test"], y[traj_model@type == "test"])^2
        cli_alert_info("R^2 train before relearning: {.val {r2_train_before}}, after: {.val {r2_train}}, difference: {.val {r2_train - r2_train_before}}")
        cli_alert_info("R^2 test before relearning: {.val {r2_test_before}}, after: {.val {r2_t}}, difference: {.val {r2_t - r2_test_before}}")
    }

    traj_model@model <- model
    traj_model@predicted_diff_score <- pred
    traj_model@coefs <- get_model_coefs(model)
    traj_model@model_features <- X

    traj_model <- add_traj_model_stats(traj_model)

    return(traj_model)
}

#' Add motif models to trajectory model
#'
#' This function adds specified motif models to a given trajectory model.
#' It updates the model, motif models, predicted difference score, model features,
#' coefficients, normalized energies, and features R^2 of the trajectory model.
#'
#' @param traj_model The trajectory model object to add motif models to.
#' @param new_motif_models A named list of motif models to add. Each element should contain a 'pssm' component and optionally a 'spat' component.
#' @param verbose A logical value indicating whether to display information about the R^2 after adding the motif models. Default is TRUE.
#'
#' @return The updated trajectory model object after adding the motif models.
#'
#' @export
add_motif_models_to_traj <- function(traj_model, new_motif_models, verbose = TRUE) {
    if (!all(purrr::map_lgl(new_motif_models, ~ !is.null(.x$pssm)))) {
        cli_abort("All motif models must contain a 'pssm' component.")
    }

    # Check for name conflicts
    existing_names <- names(traj_model@motif_models)
    new_names <- names(new_motif_models)
    if (any(new_names %in% existing_names)) {
        cli_abort("Motif{?s} {.val {new_names[new_names %in% existing_names]}} already exist{?s} in the trajectory model.")
    }

    traj_model_new <- traj_model

    traj_model_new_motifs <- traj_model
    traj_model_new_motifs@motif_models <- new_motif_models

    # infer energies
    traj_model_new@normalized_energies <- cbind(
        traj_model_new@normalized_energies,
        calc_traj_model_energies(traj_model_new_motifs, traj_model@peak_intervals)
    )

    traj_model_new@motif_models <- c(traj_model@motif_models, new_motif_models)

    # Create new logistic features and update model features
    X <- traj_model_new@model_features
    X_new <- create_logist_features(traj_model_new@normalized_energies)
    X <- cbind(X, X_new)
    traj_model_new@model_features <- X

    traj_model_new <- relearn_traj_model(traj_model_new, verbose = verbose)

    return(traj_model_new)
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

    ftv <- feat_to_variable(traj_model, add_types = TRUE)
    interactions <- traj_model@interactions
    if (has_interactions(traj_model)) {
        ftv <- ftv %>% filter(
            !(variable %in% motif_models), !(term1 %in% motif_models), !(term2 %in% motif_models)
        )

        feats <- ftv %>% pull(feature)
        interactions <- interactions[, intersect(colnames(interactions), ftv$variable)]
    } else {
        feats <- ftv %>%
            filter(!(variable %in% motif_models)) %>%
            pull(feature)
    }

    X_f <- X[, feats]

    traj_model@model_features <- X_f
    traj_model@motif_models <- traj_model@motif_models[vars_f]
    traj_model@features_r2 <- traj_model@features_r2[vars_f[vars_f %in% names(traj_model@features_r2)]]
    traj_model@normalized_energies <- traj_model@normalized_energies[, vars_f, drop = FALSE]
    traj_model@interactions <- interactions

    traj_model <- relearn_traj_model(traj_model, verbose = FALSE)
    r2_f <- cor(traj_model@predicted_diff_score, norm01(traj_model@diff_score))^2
    if (verbose) {
        cli_alert_info("R^2 after removing {.field {motif_models}}: {.val {r2_f}}")
    }

    return(traj_model)
}

remove_additional_feature_from_traj <- function(traj_model, feature_name, verbose = TRUE) {
    X <- traj_model@model_features
    vars <- colnames(traj_model@additional_features)

    if (!feature_name %in% vars) {
        cli_abort("Feature {.val {feature_name}} not found in the trajectory model.")
    }
    ftv <- feat_to_variable(traj_model, add_types = TRUE)

    interactions <- traj_model@interactions
    if (has_interactions(traj_model)) {
        ftv <- ftv %>% filter(
            variable != feature_name, term1 != feature_name, term2 != feature_name
        )

        feats <- ftv %>% pull(feature)
        interactions <- interactions[, intersect(colnames(interactions), ftv$variable)]
    } else {
        feats <- ftv %>%
            filter(variable != feature_name) %>%
            pull(feature)
    }

    X_f <- X[, feats]

    traj_model@model_features <- X_f
    traj_model@additional_features <- traj_model@additional_features[, vars[vars != feature_name], drop = FALSE]
    traj_model@interactions <- interactions

    traj_model <- relearn_traj_model(traj_model, verbose = FALSE)
    r2_f <- cor(traj_model@predicted_diff_score, norm01(traj_model@diff_score))^2
    if (verbose) {
        cli_alert_info("R^2 after removing {.field {feature_name}}: {.val {r2_f}}")
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
#' This function adds the added R-squared values of each motif model to the trajectory model.
#'
#' @param traj_model The trajectory model object.
#' @param sample_frac The fraction of samples to use for computing the r2 without each model. When NULL, all samples are used.
#' @param additional calculate also for additional features. Default is FALSE.
#' @param seed The seed to use for sampling the data. Default is 60427.
#'
#' @return The trajectory model object with the added R-squared values.
#'
#' @export
add_features_r2 <- function(traj_model, sample_frac = 0.1, additional = FALSE, seed = 60427) {
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
        cli::cli_alert("R^2 added by {.field {var}} ({.strong {gsub('N', '-', prego::pssm_to_kmer(pssm, pos_bits_thresh = 0.2))}}): {.val {full_model_r2 - r2}}. Bits: {.val {bits}}")
        list(r2 = r2, bits = bits)
    })
    vars_r2 <- purrr::map_dbl(var_stats, ~ .x$r2)
    names(vars_r2) <- motif_models

    if (additional) {
        add_r2 <- purrr::map_dbl(colnames(traj_model@additional_features), function(var) {
            traj_model_f_var <- remove_additional_feature_from_traj(traj_model, var, verbose = FALSE)
            r2 <- cor(traj_model_f_var@predicted_diff_score, traj_model_f_var@diff_score)^2
            cli::cli_alert("R^2 added by {.field {var}}: {.val {full_model_r2 - r2}}")
            r2
        })
        names(add_r2) <- colnames(traj_model@additional_features)
        vars_r2 <- c(vars_r2, add_r2)
    }

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
    if (!is.null(names(obs)) && sum(names(obs) == "") == 0) {
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

    if (length(obs) == length(pred)) {
        traj_model@params$stats$r2_all <- cor(pred, obs, use = "pairwise.complete.obs")^2
    }

    return(traj_model)
}

traj_model_has_test <- function(traj_model) {
    return(any(traj_model@type == "test"))
}

#' Split a trajectory model into train and test sets
#'
#' This function takes a trajectory model and splits it into train and test sets.
#'
#' @param traj_model The trajectory model to split
#' @return A list with the train and test sets
#' @export
split_traj_model_to_train_test <- function(traj_model) {
    if (!traj_model_has_test(traj_model)) {
        return(list(train = traj_model, test = NULL))
    }

    train_idxs <- traj_model@type == "train"
    test_idxs <- traj_model@type == "test"

    traj_model_train <- filter_traj_model_intervals(traj_model, train_idxs)
    stopifnot(all(traj_model_train@type == "train"))
    traj_model_test <- filter_traj_model_intervals(traj_model, test_idxs)
    stopifnot(all(traj_model_test@type == "test"))

    return(list(train = traj_model_train, test = traj_model_test))
}


#' Filter trajectory model intervals
#'
#' This function filters a trajectory model by selecting specific intervals based on their indices.
#'
#' @param traj_model The trajectory model object to filter.
#' @param idxs The indices of the intervals to select.
#'
#' @return The filtered trajectory model object.
#'
#'
#' @export
filter_traj_model_intervals <- function(traj_model, idxs) {
    traj_model@model_features <- traj_model@model_features[idxs, ]
    if (ncol(traj_model@additional_features) > 0) {
        traj_model@additional_features <- traj_model@additional_features[idxs, ]
    }
    traj_model@normalized_energies <- traj_model@normalized_energies[idxs, ]
    traj_model@diff_score <- traj_model@diff_score[idxs]
    traj_model@predicted_diff_score <- traj_model@predicted_diff_score[idxs]
    traj_model@type <- traj_model@type[idxs]
    traj_model@peak_intervals <- traj_model@peak_intervals[idxs, ]
    return(traj_model)
}

#' Match trajectory model motif names
#'
#' This function matches the motif names in a trajectory model with a given dataset.
#' This is used in order to give more 'friendly' names to the motif models.
#' The default dataset is "HOMER".
#' Note that the run might be slow.
#'
#' @param traj_model The trajectory model object.
#' @param dataset The dataset to match the motif names with. Default is the "HOMER" dataset.
#'
#' @return A named character vector mapping the motif names in the trajectory model to the matched motif names in the dataset.
#'
#' @export
match_traj_model_motif_names <- function(traj_model, dataset = prego::all_motif_datasets() %>% filter(dataset == "HOMER")) {
    motmatch <- purrr::imap_dfr(traj_model@motif_models, ~ {
        cli::cli_alert("Matching {.field {.y}}")
        m <- prego::pssm_match(.x$pssm, dataset) %>% mutate(motif1 = .y)
        cli::cli_alert("Matched with {.val {m$motif[1]}}, PSSM correlation = {.val {m$cor[1]}}")
        m
    })
    names_map <- motmatch %>%
        arrange(motif1, desc(cor)) %>%
        distinct(motif1, .keep_all = TRUE) %>%
        select(motif1, motif, cor) %>%
        select(motif1, motif) %>%
        mutate(motif = gsub("^HOMER.", "", motif)) %>%
        deframe()

    names_map <- names_map %>%
        enframe("motif", "name") %>%
        group_by(name) %>%
        mutate(i = 1:n()) %>%
        mutate(name = ifelse(i > 1, paste0(name, ".", i), name)) %>%
        select(motif, name) %>%
        deframe()

    return(names_map)
}


#' Rename motif models in a trajectory model
#'
#' This function renames motif models in a trajectory model based on a provided names_map.
#' The names_map should be a named character vector where the names are the current names of the motif models
#' and the values are the new names.
#'
#' @param traj_model The trajectory model object to modify.
#' @param names_map A named character vector specifying the mapping of current motif model names to new names.
#' Can be generated using the \code{match_traj_model_motif_names} function.
#' @return The modified trajectory model object with renamed motif models.
#'
#' @export
rename_motif_models <- function(traj_model, names_map) {
    if (!all(names(traj_model@motif_models) %in% names(names_map))) {
        cli_abort("Some of the motif models to rename are not found in the trajectory model. The format of the names_map should be a named character vector where the names are the current names of the motif models and the values are the new names.")
    }

    if (length(unique(names_map)) != length(names_map)) {
        cli_abort("The new names in the names_map should be unique.")
    }

    traj_model@params$names_map <- names_map

    add_feat_map <- colnames(traj_model@additional_features)
    names(add_feat_map) <- add_feat_map
    names_map <- c(add_feat_map, names_map)

    names(traj_model@motif_models) <- names_map[names(traj_model@motif_models)]
    colnames(traj_model@normalized_energies) <- names_map[colnames(traj_model@normalized_energies)]
    names(colnames(traj_model@normalized_energies)) <- NULL
    names(traj_model@features_r2) <- names_map[names(traj_model@features_r2)]

    ext_names_map <- purrr::map(c("low-energy", "high-energy", "higher-energy", "sigmoid"), ~ {
        m <- paste0(names_map, "_", .x)
        names(m) <- paste0(names(names_map), "_", .x)
        m
    }) %>% do.call(c, .)

    colnames(traj_model@model_features) <- ext_names_map[colnames(traj_model@model_features)]
    names(colnames(traj_model@model_features)) <- NULL
    traj_model@coefs <- traj_model@coefs %>%
        mutate(variable = names_map[variable])
    if (!is.null(traj_model@params$distilled_features)) {
        traj_model@params$distilled_features <- traj_model@params$distilled_features %>%
            mutate(distilled = ifelse(distilled %in% names(names_map), names_map[distilled], distilled))
    }
    if (!is.null(traj_model@params$features_bits)) {
        names(traj_model@params$features_bits) <- names_map[names(traj_model@params$features_bits)]
    }

    traj_model <- relearn_traj_model(traj_model, new_energies = FALSE, new_logist = FALSE, lambda = NULL, use_additional_features = TRUE, use_motifs = TRUE, verbose = FALSE)

    return(traj_model)
}


has_additional_features <- function(traj_model) {
    return(ncol(traj_model@additional_features) > 0 && nrow(traj_model@additional_features) > 0)
}

#' Adjust Energy Values in a Trajectory Model
#'
#' This function optimizes the energy normalization parameters for each motif in a trajectory model
#' by testing different combinations of quantile thresholds and minimum energy values. For each motif,
#' it selects the parameters that maximize the absolute correlation with the differential score.
#'
#' @param traj_model The trajectory model object containing the motif models and peak intervals.
#' @param q A numeric vector of quantile values to test for normalization. Default values range from 1 to 0.95.
#' @param min_energy A numeric vector of minimum energy values to test. Default values range from -20 to -5.
#' @param diff_thresh The threshold for the difference in correlation between the current and best parameters. Default is 0.01.
#' @param sequences Optional character vector of DNA sequences corresponding to the peak intervals.
#' If NULL, sequences will be extracted from the trajectory model.
#' @param norm_sequences Optional character vector of DNA sequences for normalization.
#' If NULL, sequences will be extracted from the trajectory model.
#'
#' @return The updated trajectory model object with optimized energy values for each motif.
#'
#' @details
#' The function works by:
#' 1. Extracting DNA sequences from the peak intervals if not provided
#' 2. Computing motif energies for all sequences
#' 3. Testing different combinations of normalization parameters
#' 4. For each motif, selecting the parameters that maximize correlation with the differential score
#' 5. Applying the selected parameters to normalize the energies
#'
#' @examples
#' # Create a trajectory model and adjust its energies
#' traj_model <- create_traj_model(...)
#' traj_model <- adjust_energies(traj_model)
#'
#' # Adjust energies with custom parameters
#' traj_model <- adjust_energies(
#'     traj_model,
#'     q = c(1, 0.99, 0.95),
#'     min_energy = c(-15, -10, -5)
#' )
#'
#' @inheritParams relearn_traj_model
#' @export
adjust_energies <- function(traj_model, q = c(1, 0.999, 0.995, 0.99, 0.98, 0.97, 0.96, 0.95), min_energy = c(-20, -15, -13, -10, -7, -5), diff_thresh = 0.01, sequences = NULL, norm_sequences = NULL, relearn = TRUE, lambda = NULL, use_additional_features = TRUE, use_motifs = TRUE, verbose = TRUE) {
    if (is.null(sequences) || is.null(norm_sequences)) {
        seqs <- extract_traj_model_sequences(traj_model, traj_model@peak_intervals)
        sequences <- seqs$sequences
        norm_sequences <- seqs$norm_sequences
    }

    mdb <- motifs_to_mdb(traj_model@motif_models)
    all_energies <- prego::extract_pwm(c(sequences, norm_sequences), dataset = mdb, prior = 0.01)
    energies <- all_energies[1:length(sequences), ]
    norm_energies <- all_energies[(length(sequences) + 1):nrow(all_energies), ]
    rownames(energies) <- traj_model@peak_intervals$peak_name
    rownames(norm_energies) <- traj_model@normalization_intervals$peak_name
    energies_train <- energies[traj_model@type == "train", ]
    y_train <- traj_model@diff_score[traj_model@type == "train"]

    grid_params <- expand.grid(
        q = q,
        min_energy = min_energy
    )

    norm_matrices <- purrr::map(1:nrow(grid_params), function(i) {
        q <- grid_params$q[i]
        min_energy <- grid_params$min_energy[i]
        norm_energy_matrix(energies, norm_energies, min_energy = min_energy, q = q, norm_energy_max = traj_model@params$norm_energy_max)
    })
    names(norm_matrices) <- paste0("q.", grid_params$q, "_min_energy.", grid_params$min_energy)

    n_stats <- purrr::imap_dfr(norm_matrices, ~ {
        tgs_cor(.x[traj_model@type == "train", ], as.matrix(y_train))[, 1] %>%
            enframe("motif", "cor") %>%
            mutate(param = .y)
    })

    current_param <- paste("q.", traj_model@params$energy_norm_quantile %||% 0.999, "_min_energy.", traj_model@params$min_energy %||% -7, sep = "")

    best_per_motif <- n_stats %>%
        group_by(motif) %>%
        summarise(
            diff = max(abs(cor - cor[param == current_param])),
            best_param = param[which.max(abs(cor))],
            .groups = "drop"
        )

    new_energies <- purrr::pmap(best_per_motif, function(motif, diff, best_param) {
        if (diff > diff_thresh) {
            norm_matrices[[best_param]][, motif]
        } else {
            traj_model@normalized_energies[, motif]
        }
    }) %>% do.call(cbind, .)
    colnames(new_energies) <- best_per_motif$motif
    new_energies <- new_energies[, colnames(traj_model@normalized_energies)]

    traj_model@normalized_energies <- new_energies
    if (relearn) {
        traj_model <- relearn_traj_model(traj_model, new_energies = TRUE, new_logist = TRUE, lambda = lambda, use_additional_features = use_additional_features, use_motifs = TRUE, verbose = verbose)
    }
    return(traj_model)
}

