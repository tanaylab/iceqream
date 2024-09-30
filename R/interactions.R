has_interactions <- function(traj_model) {
    n_interactions(traj_model) > 0
}

n_interactions <- function(traj_model) {
    ncol(traj_model@interactions)
}

create_specifc_terms <- function(energies, terms) {
    term1_matrix <- energies[, terms$term1]
    term2_matrix <- energies[, terms$term2]
    inter <- term1_matrix * term2_matrix
    inter <- t(t(inter) / apply(inter, 2, max, na.rm = TRUE))
    inter <- apply(inter, 2, norm01) * 1
    colnames(inter) <- terms$variable
    return(inter)
}


create_interaction_terms <- function(energies, motif_feats = NULL, add_feats = NULL, additional_features = NULL, max_motif_n = NULL, max_add_n = NULL) {
    create_interactions <- function(features, data, max_n) {
        if (is.null(features) || is.null(data)) {
            return(NULL)
        }

        features <- head(features, n = max_n %||% length(features))

        interactions <- purrr::map_dfc(features, ~ {
            inter <- energies[, setdiff(colnames(energies), .x)] * data[, .x]
            inter <- t(t(inter) / apply(inter, 2, max, na.rm = TRUE))
            colnames(inter) <- paste0(.x, ":", colnames(inter))
            inter
        })

        interactions <- apply(interactions, 2, norm01) * 1
        interactions
    }


    add_inter <- create_interactions(add_feats, additional_features, max_add_n)

    if (!is.null(add_inter)) {
        cli::cli_alert_info("Created {.val {ncol(add_inter)}} interactions between additional features and motif features.")
    }


    motif_inter <- create_interactions(motif_feats, energies, max_motif_n)
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
    additional_features = NULL, lambda = 1e-5, alpha = 1, seed = 60427,
    ignore_feats = c("TT", "CT", "GT", "AT", "TC", "CC", "GC", "AC", "TG", "CG", "GG", "AG", "TA", "CA", "GA", "AA")) {
    glm_model_lin <- glmnet::glmnet(as.matrix(energies), y, binomial(link = "logit"), alpha = alpha, lambda = lambda, seed = seed)

    feats_all <- abs(stats::coef(glm_model_lin)[-1])
    names(feats_all) <- rownames(stats::coef(glm_model_lin))[-1]
    sig_feats <- names(feats_all)[feats_all > interaction_threshold]
    sig_feats <- setdiff(sig_feats, ignore_feats)

    if (length(sig_feats) == 0) {
        cli::cli_alert_warning("No significant features to consider for interactions.")
        return(NULL)
    }

    add_feats <- intersect(sig_feats, colnames(additional_features))
    motif_feats <- setdiff(sig_feats, add_feats)

    cli::cli_alert_info("# of significant features to consider for interactions: {.val {length(sig_feats)}} (out of {.val {ncol(energies)}}) above the threshold of {.val {interaction_threshold}}. Of these, {.val {length(motif_feats)}} are motif features and {.val {length(add_feats)}} are additional features.")

    if (!is.null(additional_features)) {
        # remove the features from energies
        energies <- energies[, setdiff(colnames(energies), colnames(additional_features))]

        # remove the ignored features
        energies <- energies[, setdiff(colnames(energies), ignore_feats)]
    }

    create_interaction_terms(energies,
        motif_feats = motif_feats, add_feats = add_feats,
        additional_features = additional_features, max_motif_n = max_motif_n, max_add_n = max_add_n
    )
}

#' Add interactions to a trajectory model
#'
#' This function adds significant interactions to a given trajectory model if they do not already exist.
#' It identifies significant interactions based on the provided threshold and updates the model features
#' with logistic features derived from these interactions. The trajectory model is then re-learned with
#' the new features.
#'
#' @inheritParams regress_trajectory_motifs
#'
#' @return The updated trajectory model with added interactions.
#' @export
add_interactions <- function(traj_model, interaction_threshold = 0.001, max_motif_n = NULL, max_add_n = NULL, lambda = 1e-5, alpha = 1, seed = 60427) {
    if (!has_interactions(traj_model)) {
        cli::cli_alert("Adding interactions")
        interactions <- get_significant_interactions(
            cbind(traj_model@normalized_energies, traj_model@additional_features), norm01(traj_model@diff_score), interaction_threshold,
            max_motif_n = max_motif_n, max_add_n = max_add_n,
            additional_features = traj_model@additional_features, lambda = lambda, alpha = alpha, seed = seed
        )

        if (!is.null(interactions)) {
            traj_model@interactions <- interactions
        }

        logist_inter <- create_logist_features(interactions)
        traj_model@model_features <- cbind(traj_model@model_features, logist_inter)

        cli::cli_alert_info("Re-learning the model with the new interactions. Number of features: {.val {ncol(traj_model@model_features)}}")
        cli::cli_alert_info("R^2 all before learning: {.val {cor(traj_model@diff_score, traj_model@predicted_diff_score)^2}}")
        if (traj_model_has_test(traj_model)) {
            cli::cli_alert_info("R^2 train before learning: {.val {cor(traj_model@diff_score[traj_model@type == 'train'], traj_model@predicted_diff_score[traj_model@type == 'train'])^2}}")
            cli::cli_alert_info("R^2 test before learning: {.val {cor(traj_model@diff_score[traj_model@type == 'test'], traj_model@predicted_diff_score[traj_model@type == 'test'])^2}}")
        }

        traj_model <- relearn_traj_model(traj_model, new_energies = FALSE, new_logist = FALSE, use_additional_features = TRUE, use_motifs = TRUE, verbose = FALSE)
        cli::cli_alert_info("R^2 all after learning: {.val {cor(traj_model@diff_score, traj_model@predicted_diff_score)^2}}")
        if (traj_model_has_test(traj_model)) {
            cli::cli_alert_info("R^2 train after learning: {.val {cor(traj_model@diff_score[traj_model@type == 'train'], traj_model@predicted_diff_score[traj_model@type == 'train'])^2}}")
            cli::cli_alert_info("R^2 test after learning: {.val {cor(traj_model@diff_score[traj_model@type == 'test'], traj_model@predicted_diff_score[traj_model@type == 'test'])^2}}")
        }
    } else {
        cli::cli_alert_warning("Interactions already exist.")
    }

    return(traj_model)
}
