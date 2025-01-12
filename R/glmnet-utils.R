strip_glmnet <- function(model) {
    essentials <- c(
        "beta", "lambda", "a0", "dev.ratio", "nulldev",
        "df", "dim", "nobs", "offset"
    )

    stripped <- structure(
        lapply(essentials, function(x) model[[x]]),
        names = essentials,
        class = class(model)
    )

    stripped$call <- NULL

    return(stripped)
}

strip_traj_model <- function(traj_model) {
    traj_model@model <- strip_glmnet(traj_model@model)
    traj_model@params$distilled_features <- NULL
    return(traj_model)
}

strip_traj_model_multi <- function(mm) {
    mm@models <- purrr::map(mm@models, strip_traj_model)
    mm@models_full <- purrr::map(mm@models_full, strip_traj_model)
    return(mm)
}
