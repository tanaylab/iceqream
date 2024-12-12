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
    stripped$family <- binomial(link = "logit")

    return(stripped)
}

strip_traj_model <- function(traj_model) {
    traj_model@params$distilled_features <- NULL
    return(traj_model)
}
