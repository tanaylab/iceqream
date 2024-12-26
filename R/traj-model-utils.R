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
#' @return The updated trajectory model object.
#'
#' @export
relearn_traj_model <- function(traj_model, new_energies = FALSE, new_logist = FALSE, lambda = NULL, use_additional_features = TRUE, use_motifs = TRUE, verbose = FALSE) {
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
            if (use_motifs) {
                traj_model@model_features <- as.matrix(cbind(
                    create_logist_features(traj_model@normalized_energies),
                    create_logist_features(traj_model@additional_features)
                ))
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

    X <- traj_model@model_features

    if (is.null(lambda)) {
        lambda <- traj_model@params$lambda
    }

    y <- norm01(traj_model@diff_score)

    X_train <- X[traj_model@type == "train", ]
    y_train <- y[traj_model@type == "train"]

    model <- glmnet::glmnet(X_train, y_train, binomial(link = "logit"), alpha = traj_model@params$alpha, lambda = lambda, seed = traj_model@params$seed)
    model <- strip_glmnet(model)

    pred <- logist(glmnet::predict.glmnet(model, newx = X, type = "link", s = traj_model@params$lambda))[, 1]
    pred <- norm01(pred)
    pred <- rescale(pred, traj_model@diff_score)
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
    cli::cli_alert_info("Matching motif names, note that this might take a while.")
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
    names(traj_model@features_r2) <- names_map[names(traj_model@features_r2)]

    ext_names_map <- purrr::map(c("low-energy", "high-energy", "higher-energy", "sigmoid"), ~ {
        m <- paste0(names_map, "_", .x)
        names(m) <- paste0(names(names_map), "_", .x)
        m
    }) %>% do.call(c, .)

    colnames(traj_model@model_features) <- ext_names_map[colnames(traj_model@model_features)]
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
