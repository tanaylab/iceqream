feat_to_variable <- function(traj_model) {
    tibble::tibble(
        feature = colnames(traj_model@model_features),
        variable = sub("_(low-energy|high-energy|higher-energy|sigmoid)$", "", feature)
    )
}

variable_to_feat <- function(glm_model, var) {
    vars <- purrr::map(var, ~ paste0(.x, c("_low-energy", "_high-energy", "_higher-energy", "_sigmoid"))) %>% do.call(c, .)
    rownames(glm_model$beta)[rownames(glm_model$beta) %in% vars]
}

compute_partial_response <- function(traj_model, vars = NULL, lambda = 1e-5) {
    f2v <- feat_to_variable(traj_model)

    if (!is.null(vars)) {
        f2v <- f2v %>% filter(variable %in% vars)
    }

    pr <- plyr::dlply(f2v, "variable", function(x) {
        feats <- traj_model@model_features[, x$feature, drop = FALSE]
        (feats %*% coef(traj_model@model, s = lambda)[x$feature, , drop = FALSE])[, 1]
    }) %>% do.call(cbind, .)

    return(as.data.frame(pr))
}
