feat_to_variable <- function(traj_model) {
    tibble::tibble(
        feature = colnames(traj_model@model_features),
        variable = sub("_(low-energy|high-energy|higher-energy|sigmoid)$", "", feature)
    )
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
