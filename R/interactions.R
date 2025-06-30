has_interactions <- function(traj_model) {
    if (!.hasSlot(traj_model, "interactions")) {
        return(FALSE)
    }
    n_interactions(traj_model) > 0
}

n_interactions <- function(traj_model) {
    ncol(traj_model@interactions)
}

create_specifc_terms <- function(energies, terms, scale = 1) {
    term1_matrix <- energies[, terms$term1]
    term2_matrix <- energies[, terms$term2]
    inter <- term1_matrix * term2_matrix
    inter <- t(t(inter) / apply(inter, 2, max, na.rm = TRUE))
    inter <- apply(inter, 2, norm01) * scale
    colnames(inter) <- terms$variable
    return(inter)
}

create_features_terms <- function(energies, features, data, scale = 1) {
    interactions <- purrr::map_dfc(features, ~ {
        inter <- energies[, setdiff(colnames(energies), .x), drop = FALSE] * data[, .x]
        inter <- t(t(inter) / apply(inter, 2, max, na.rm = TRUE))
        colnames(inter) <- paste0(.x, ":", colnames(inter))
        inter
    })

    interactions <- apply(interactions, 2, norm01) * scale
    interactions
}


create_interaction_terms <- function(energies, motif_feats = NULL, add_feats = NULL, additional_features = NULL, max_motif_n = NULL, max_add_n = NULL, only_sig_motifs = FALSE) {
    create_interactions <- function(features, data, max_n) {
        if (is.null(features) || is.null(data)) {
            return(NULL)
        }

        features <- head(features, n = max_n %||% length(features))
        interactions <- create_features_terms(data, features, data)

        interactions
    }


    if (length(add_feats) > 0) {
        add_inter <- create_interactions(add_feats, cbind(energies, additional_features), max_add_n)

        if (!is.null(add_inter)) {
            cli::cli_alert_info("Created {.val {ncol(add_inter)}} interactions between additional features and motif features.")
        }
    } else {
        add_inter <- NULL
    }

    if (only_sig_motifs) {
        motif_inter <- create_interactions(motif_feats, energies[, motif_feats], max_motif_n)
    } else {
        motif_inter <- create_interactions(motif_feats, energies, max_motif_n)
    }

    if (!is.null(motif_inter)) {
        cli::cli_alert_info("Created {.val {ncol(motif_inter)}} interactions between motif features.")
    }

    interactions <- cbind(motif_inter, add_inter)
    if (!is.null(interactions)) {
        if (!is.null(rownames(energies))) {
            rownames(interactions) <- rownames(energies)
        }
        cli::cli_alert_info("Created {.val {ncol(interactions)}} interactions in total.")
    }

    return(interactions)
}

get_significant_interactions <- function(
    energies, y, interaction_threshold, max_motif_n = NULL, max_add_n = NULL,
    max_n = NULL,
    additional_features = NULL, lambda = 1e-5, alpha = 1, seed = 60427,
    ignore_feats = c("TT", "CT", "GT", "AT", "TC", "CC", "GC", "AC", "TG", "CG", "GG", "AG", "TA", "CA", "GA", "AA"), idxs = NULL, only_sig_motifs = FALSE) {
    if (is.null(idxs)) {
        idxs <- seq_len(nrow(energies))
    }
    glm_model_lin <- glmnet::glmnet(as.matrix(energies[idxs, ]), y[idxs], binomial(link = "logit"), alpha = alpha, lambda = lambda, seed = seed)
    glm_model_lin <- strip_glmnet(glm_model_lin)

    feats_all <- abs(stats::coef(glm_model_lin)[-1])
    names(feats_all) <- rownames(stats::coef(glm_model_lin))[-1]
    sig_feats <- names(feats_all)[feats_all > interaction_threshold]
    sig_feats <- setdiff(sig_feats, ignore_feats)

    if (length(sig_feats) == 0) {
        cli::cli_alert_warning("No significant features to consider for interactions.")
        return(NULL)
    }
    additional_features <- additional_features[, setdiff(colnames(additional_features), ignore_feats)]
    add_feats <- intersect(sig_feats, colnames(additional_features))
    add_feats <- setdiff(add_feats, ignore_feats)
    motif_feats <- setdiff(sig_feats, add_feats)

    cli::cli_alert_info("# of significant features to consider for interactions: {.val {length(sig_feats)}} (out of {.val {ncol(energies)}}) above the threshold of {.val {interaction_threshold}}. Of these, {.val {length(motif_feats)}} are motif features and {.val {length(add_feats)}} are additional features.")

    if (!is.null(additional_features)) {
        # remove the features from energies
        energies <- energies[, setdiff(colnames(energies), colnames(additional_features))]

        # remove the ignored features
        energies <- energies[, setdiff(colnames(energies), ignore_feats)]
    }

    inter <- create_interaction_terms(energies,
        motif_feats = motif_feats, add_feats = add_feats,
        additional_features = additional_features, max_motif_n = max_motif_n, max_add_n = max_add_n, only_sig_motifs = only_sig_motifs
    )

    if (!is.null(max_n) && ncol(inter) > max_n) {
        cli::cli_alert_info("Selecting top {.val {max_n}} interactions based on correlation with the signal.")
        max_n <- min(max_n, ncol(inter))
        cm <- tgs_cor(inter[idxs, ], as.matrix(y[idxs]))[, 1]
        chosen <- names(sort(abs(cm), decreasing = TRUE)[1:max_n])
        inter <- inter[, chosen]
    }

    return(inter)
}

#' Add interactions to a trajectory model
#'
#' This function adds significant interactions to a given trajectory model if they do not already exist.
#' It identifies significant interactions based on the provided threshold and updates the model features
#' with logistic features derived from these interactions. The trajectory model is then re-learned with
#' the new features.
#' @param interaction_threshold threshold for the selecting features to create interactions. IQ learns a linear model on the features and selects the features with coefficients above this threshold. Default: 0.001
#' @param max_motif_n maximum number of motifs to consider for interactions. If NULL, all motifs above the interaction_threshold will be considered. Default: NULL
#' @param max_add_n maximum number of additional features to consider for interactions. If NULL, all additional features above the interaction_threshold will be considered. Default: NULL
#' @param max_n maximum number of interactions to consider. If NULL, all interactions will be considered. If set, the interactions will be selected based on correlation with the signal in the training data. Default: (motif models + additional features) * 10
#' @param interactions A precomputed interaction matrix. If provided, the function will not compute the interactions. Default: NULL
#' @param ignore_feats A character vector of features to ignore when creating interactions. Default: dinucleotides
#' @param force If TRUE, the function will add interactions even if they already exist. Default: FALSE
#' @param logist_interactions Logical indicating whether to transform interactions to logistic functions. Default: FALSE
#' @param only_sig_motifs Logical indicating whether to only consider significant motifs for interactions. Default: FALSE
#'
#' @inheritParams regress_trajectory_motifs
#' @inheritParams relearn_traj_model
#'
#' @return The updated trajectory model with added interactions.
#' @export
add_interactions <- function(traj_model, interaction_threshold = 0.001, max_motif_n = NULL, max_add_n = NULL, max_n = NULL, lambda = 1e-5, alpha = 1, seed = 60427, interactions = NULL, ignore_feats = c("TT", "CT", "GT", "AT", "TC", "CC", "GC", "AC", "TG", "CG", "GG", "AG", "TA", "CA", "GA", "AA"), force = FALSE, logist_interactions = FALSE, use_cv = FALSE, nfolds = 10, family = "binomial", rescale_pred = FALSE, only_sig_motifs = FALSE) {
    r2_all_before <- cor(traj_model@diff_score, traj_model@predicted_diff_score)^2
    if (traj_model_has_test(traj_model)) {
        r2_train_before <- cor(traj_model@diff_score[traj_model@type == "train"], traj_model@predicted_diff_score[traj_model@type == "train"])^2
        r2_test_before <- cor(traj_model@diff_score[traj_model@type == "test"], traj_model@predicted_diff_score[traj_model@type == "test"])^2
    }

    if (is.null(max_n)) {
        n_feats <- ncol(traj_model@additional_features[, setdiff(colnames(traj_model@additional_features), ignore_feats)]) %||% 0
        max_n <- (ncol(traj_model@normalized_energies) + n_feats) * 10
        cli::cli_alert("Setting {.field max_n} (maximal number of interactions) to {.val {max_n}}.")
    }

    if (has_interactions(traj_model)) {
        if (force) {
            traj_model <- remove_interactions(traj_model)
        } else {
            cli::cli_alert_warning("Interactions already exist. Set {.field force} to TRUE to overwrite.")
            return(traj_model)
        }
    }


    if (is.null(interactions)) {
        cli::cli_alert("Adding interactions")
        interactions <- get_significant_interactions(
            cbind(traj_model@normalized_energies, traj_model@additional_features), norm01(traj_model@diff_score), interaction_threshold,
            max_motif_n = max_motif_n, max_add_n = max_add_n,
            max_n = max_n,
            additional_features = traj_model@additional_features, lambda = lambda, alpha = alpha, seed = seed, ignore_feats = ignore_feats, idxs = which(traj_model@type == "train"), only_sig_motifs = only_sig_motifs
        )
    }

    if (!is.null(interactions)) {
        traj_model@interactions <- interactions
    }

    if (logist_interactions) {
        logist_inter <- create_logist_features(interactions)
        traj_model@model_features <- cbind(traj_model@model_features, logist_inter)
        cli::cli_alert_info("Re-learning the model with the new interactions. Number of features: {.val {ncol(traj_model@model_features)}}, out of which {.val {ncol(traj_model@normalized_energies)}}*4 are motif features, {.val {ncol(traj_model@additional_features)}}*4 are additional features and {.val {ncol(traj_model@interactions)}}*4 are interactions.")
    } else {
        traj_model@model_features <- cbind(traj_model@model_features, interactions)
        cli::cli_alert_info("Re-learning the model with the new interactions. Number of features: {.val {ncol(traj_model@model_features)}}, out of which {.val {ncol(traj_model@normalized_energies)}}*4 are motif features, {.val {ncol(traj_model@additional_features)}}*4 are additional features and {.val {ncol(traj_model@interactions)}} are interactions.")
    }
    traj_model@params$logist_interactions <- logist_interactions


    traj_model <- relearn_traj_model(traj_model, new_energies = FALSE, new_logist = FALSE, use_additional_features = TRUE, use_motifs = TRUE, verbose = FALSE, logist_interactions = logist_interactions, use_cv = use_cv, nfolds = nfolds, lambda = lambda, family = family, rescale_pred = rescale_pred)
    if (traj_model_has_test(traj_model)) {
        r2_train_after <- cor(traj_model@diff_score[traj_model@type == "train"], traj_model@predicted_diff_score[traj_model@type == "train"])^2
        cli::cli_alert_info("R^2 train before: {.val {r2_train_before}}, after: {.val {r2_train_after}}, change: {.val {r2_train_after - r2_train_before}}")
        r2_test_after <- cor(traj_model@diff_score[traj_model@type == "test"], traj_model@predicted_diff_score[traj_model@type == "test"])^2
        cli::cli_alert_info("R^2 test before: {.val {r2_test_before}}, after: {.val {r2_test_after}}, change: {.val {r2_test_after - r2_test_before}}")
    } else {
        r2_after <- cor(traj_model@diff_score, traj_model@predicted_diff_score)^2
        cli::cli_alert_info("R^2 all before: {.val {r2_all_before}}, after: {.val {r2_after}}, change: {.val {r2_after - r2_all_before}}")
    }


    return(traj_model)
}

remove_interactions <- function(traj_model) {
    if (has_interactions(traj_model)) {
        if (is.null(traj_model@params$logist_interactions) || !traj_model@params$logist_interactions) {
            inter_terms <- colnames(traj_model@interactions)
        } else {
            inter_terms <- purrr::map(c("low-energy", "high-energy", "higher-energy", "sigmoid"), ~ {
                paste0(colnames(traj_model@model_features), "_", .x)
            }) %>% do.call(c, .)
        }
        traj_model@model_features <- traj_model@model_features[, !(colnames(traj_model@model_features) %in% inter_terms)]
        traj_model@interactions <- matrix(nrow = 0, ncol = 0)
        traj_model@params$logist_interactions <- FALSE
    }
    return(traj_model)
}

add_interactions_multi_traj <- function(mm, add_to_full = TRUE, max_motif_n = 30, max_inter_n = max_motif_n * 10, max_add_n = NULL, ignore_feats = c("TT", "CT", "GT", "AT", "TC", "CC", "GC", "AC", "TG", "CG", "GG", "AG", "TA", "CA", "GA", "AA")) {
    if (add_to_full) {
        models <- mm@models_full
    } else {
        models <- mm@models
    }
    energies <- mm@models_full[[1]]@normalized_energies
    y_mat <- purrr::map(models, ~ .x@diff_score) %>%
        do.call(cbind, .) %>%
        as.matrix()
    train_idxs <- which(mm@models_full[[1]]@type == "train")

    select_top <- function(mat, n) {
        n <- min(n, ncol(mat))
        cm <- tgs_cor(mat[train_idxs, , drop = FALSE], y_mat[train_idxs, , drop = FALSE])
        chosen <- apply(abs(cm), 2, function(x) names(sort(x, decreasing = TRUE)[1:n]))
        return(chosen)
    }

    max_motif_n <- max_motif_n %||% ncol(energies)
    chosen <- select_top(energies, max_motif_n)
    cli::cli_alert_info("{.val {nrow(chosen)}} motifs selected for interactions in each model.")

    if (has_additional_features(mm@models_full[[1]])) {
        additional_features <- as.matrix(mm@models_full[[1]]@additional_features)
        additional_features <- additional_features[, setdiff(colnames(additional_features), ignore_feats)]
        max_add_n <- max_add_n %||% ncol(additional_features)
        chosen_add <- select_top(additional_features, max_add_n)
        cli::cli_alert_info("{.val {nrow(chosen_add)}} additional features selected for interactions in each model.")
    } else {
        additional_features <- NULL
        chosen_add <- NULL
    }

    cli::cli_alert("Creating interactions between motifs.")

    traj_inter <- purrr::map(names(models), ~ {
        cli::cli_alert_info("Creating interactions for model {.val {.x}}.")
        motifs <- chosen[, .x]
        create_features_terms(energies[, motifs], motifs, energies[, motifs])
    })

    cli::cli_alert_info("{.val {ncol(traj_inter[[1]])}} interactions created between motifs.")

    if (!is.null(additional_features)) {
        cli::cli_alert("Creating interactions between additional features and motifs.")
        traj_inter_add <- purrr::map(names(models), ~ {
            cli::cli_alert_info("Creating interactions for model {.val {.x}}.")
            motifs <- chosen[, .x]
            add_feats <- chosen_add[, .x]
            create_features_terms(energies[, motifs], add_feats, additional_features)
        })

        cli::cli_alert_info("{.val {ncol(traj_inter_add[[1]])}} interactions created between additional features and motifs.")
        traj_inter <- purrr::map2(traj_inter, traj_inter_add, ~ cbind(.x, .y))
    }

    if (!is.null(max_inter_n)) {
        cli::cli_alert("Selecting top interactions based on correlation with the signal.")
        traj_inter <- purrr::imap(traj_inter, ~ {
            n <- min(max_inter_n, ncol(.x))
            cm <- tgs_cor(.x[train_idxs, , drop = FALSE], y_mat[train_idxs, .y, drop = FALSE])[, 1]
            chosen <- names(sort(abs(cm), decreasing = TRUE)[1:n])
            .x[, chosen]
        })
        names(traj_inter) <- names(models)

        cli::cli_alert_info("{.val {ncol(traj_inter[[1]])}} interactions selected based on correlation with the signal.")
    }

    models_inter <- purrr::imap(models, ~ {
        cli::cli_h3("Adding interactions to model {.val {.y}}.")
        add_interactions(.x, interactions = traj_inter[[.y]])
    })

    if (add_to_full) {
        mm@models_full <- models_inter
    } else {
        mm@models <- models_inter
    }

    stats_df <- purrr::imap_dfr(models_inter, ~ {
        tibble(
            model = .y,
            r2_all = .x@params$stats$r2_all,
            r2_train = .x@params$stats$r2_train,
            r2_test = .x@params$stats$r2_test,
            n_motifs = ncol(.x@normalized_energies),
            n_features = ncol(.x@model_features),
            type = "interactions"
        )
    })

    stats_no_inter <- purrr::imap_dfr(models, ~ {
        tibble(
            model = .y,
            r2_all = .x@params$stats$r2_all,
            r2_train = .x@params$stats$r2_train,
            r2_test = .x@params$stats$r2_test,
            n_features = ncol(.x@model_features),
            n_motifs = ncol(.x@normalized_energies),
            type = "no interactions"
        )
    })

    mm@stats <- dplyr::bind_rows(stats_df, stats_no_inter, mm@stats)

    return(mm)
}
