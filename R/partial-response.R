classify_var <- function(var, traj_model) {
    case_when(
        var %in% names(traj_model@motif_models) ~ "motif",
        var %in% colnames(traj_model@additional_features) ~ "additional",
        var %in% colnames(traj_model@interactions) ~ "interaction"
    )
}

feat_to_variable <- function(traj_model, add_types = FALSE) {
    ftv <- tibble::tibble(
        feature = colnames(traj_model@model_features),
        variable = sub("_(low-energy|high-energy|higher-energy|sigmoid)$", "", feature)
    )
    if (add_types) {
        ftv <- ftv %>%
            mutate(type = classify_var(variable, traj_model))

        if (has_interactions(traj_model)) {
            ftv_nointer <- ftv %>%
                filter(type != "interaction")
            ftv_inter <- ftv %>%
                filter(type == "interaction") %>%
                separate(variable, c("term1", "term2"), sep = ":", remove = FALSE) %>%
                mutate(
                    term1_type = classify_var(term1, traj_model),
                    term2_type = classify_var(term2, traj_model)
                )
            ftv <- bind_rows(ftv_nointer, ftv_inter)
        }
    }
    return(ftv)
}

variable_to_feat <- function(glm_model, var) {
    vars <- purrr::map(var, ~ paste0(.x, c("_low-energy", "_high-energy", "_higher-energy", "_sigmoid"))) %>% do.call(c, .)
    rownames(glm_model$beta)[rownames(glm_model$beta) %in% vars]
}

compute_partial_response <- function(traj_model, vars = NULL, lambda = 1e-5, reverse = FALSE) {
    f2v <- feat_to_variable(traj_model)
    f2v_all <- f2v

    if (!is.null(vars)) {
        f2v <- f2v %>% filter(variable %in% vars)
    }

    pr <- plyr::dlply(f2v, "variable", function(x) {
        if (reverse) {
            variables <- f2v_all %>%
                filter(feature != x$feature) %>%
                pull(feature)
        } else {
            variables <- x$feature
        }
        feats <- traj_model@model_features[, variables, drop = FALSE]

        (feats %*% coef(traj_model@model, s = lambda)[variables, , drop = FALSE])[, 1]
    }) %>% do.call(cbind, .)

    return(as.data.frame(pr))
}

#' Compute Partial Response for Trajectory Model
#'
#' This function computes the partial response for a given trajectory model.
#'
#' @param traj_model A trajectory model object.
#' @param vars A vector of variables to compute the partial response for. Default is NULL (all variables).
#' @param lambda A regularization parameter. Default is 1e-5.
#' @param reverse If TRUE, the partial response is computed for all variables except the ones in `vars`. Default is FALSE.
#'
#' @return A data frame with the computed partial responses, where the column names are the variables.
#'
#'
#' @export
traj_model_variable_response <- function(traj_model, vars = NULL, lambda = 1e-5, reverse = FALSE) {
    pr <- compute_partial_response(traj_model, vars = vars, lambda = lambda, reverse = reverse)
    return(pr)
}
