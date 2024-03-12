filter_model <- function(X, variables, y, alpha, lambda, seed, full_model, ignore_variables = NULL, r2_threshold = 0.0005) {
    if (!is.null(ignore_variables)) {
        variables <- variables[!(variables %in% ignore_variables)]
    }
    # for each variable of X calculate the r^2 of a GLM model without it
    vars_r2 <- plyr::llply(variables, function(var) {
        cli_alert("Testing variable {.field {var}}...")
        X_f <- X[, grep(var, colnames(X), value = TRUE, invert = TRUE)]
        m <- glmnet::glmnet(X_f, y, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = FALSE, seed = seed)
        pred <- logist(glmnet::predict.glmnet(m, newx = X_f, type = "link", s = lambda))[, 1]
        r2 <- cor(pred, y)^2
        cli::cli_alert_info("R^2 without {.field {var}}: {.val {r2}}")
        r2
    }, .parallel = FALSE)

    vars_r2 <- purrr::map_dbl(vars_r2, ~.x)
    names(vars_r2) <- variables

    full_model_r2 <- cor(logist(glmnet::predict.glmnet(full_model, newx = X, type = "link", s = lambda))[, 1], y)^2

    vars_f <- variables[(full_model_r2 - vars_r2) > r2_threshold]

    X_f <- X[, grep(paste0("(", paste(c(vars_f, ignore_variables), collapse = "|"), ")(_low-energy|_high-energy|_higher-energy|_sigmoid)"), colnames(X))]
    cli_alert_info("Number of features left: {.val {length(vars_f)}}")

    model_f <- glmnet::glmnet(X_f, y, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = FALSE, seed = seed)
    pred_f <- logist(glmnet::predict.glmnet(model_f, newx = X_f, type = "link", s = lambda))[, 1]
    r2_f <- cor(pred_f, y)^2
    cli_alert_info("R^2 after filtering: {.val {r2_f}}")

    # tibble(var = variables, r2 = full_model_r2 - vars_r2) %>%
    #     ggplot(aes(x = reorder(var, r2), y = r2)) +
    #     geom_col() +
    #     theme_classic() +
    #     vertical_labs() +
    #     geom_hline(yintercept = r2_threshold, color = "red", linetype = "dashed")

    return(list(model = model_f, pred = pred_f, X = X_f, vars_r2 = full_model_r2 - vars_r2, vars = vars_f))
}

filter_model_using_coefs <- function(X, coefs, y, alpha, lambda, seed, full_model, n_motifs, ignore_variables = NULL) {
    variables <- coefs$variable
    if (!is.null(ignore_variables)) {
        variables <- variables[!(variables %in% ignore_variables)]
    }
    coefs_max <- coefs %>%
        column_to_rownames("variable") %>%
        .[variables, ] %>%
        apply(1, max) %>%
        sort(decreasing = TRUE)

    vars_f <- names(coefs_max)[1:n_motifs]

    X_f <- X[, grep(paste0("(", paste(c(vars_f, ignore_variables), collapse = "|"), ").+"), colnames(X))]
    cli_alert_info("Number of features left: {.val {length(vars_f)}}")

    model_f <- glmnet::glmnet(X_f, y, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = FALSE, seed = seed)
    pred_f <- logist(glmnet::predict.glmnet(model_f, newx = X_f, type = "link", s = lambda))[, 1]
    r2_f <- cor(pred_f, y)^2
    cli_alert_info("R^2 after filtering: {.val {r2_f}}")

    # tibble(var = variables, r2 = full_model_r2 - vars_r2) %>%
    #     ggplot(aes(x = reorder(var, r2), y = r2)) +
    #     geom_col() +
    #     theme_classic() +
    #     vertical_labs() +
    #     geom_hline(yintercept = r2_threshold, color = "red", linetype = "dashed")

    return(list(model = model_f, pred = pred_f, X = X_f, r2 = r2_f, vars = vars_f))
}

#' Filter a trajectory model using punctuated regression.
#'
#' @description Run a regression without each feature and filter features that do not improve the model more
#' than \code{r2_threshold}.
#'
#' @param traj_model An instance of \code{TrajectoryModel}.
#' @param r2_threshold minimal R^2 for a feature to be included in the model.
#'
#' @return An instance of \code{TrajectoryModel} with the filtered model.
#'
#' @export
filter_traj_model <- function(traj_model, r2_threshold = 0.0005) {
    res <- filter_model(traj_model@model_features, traj_model@coefs$variable, norm01(traj_model@diff_score), traj_model@params$alpha, traj_model@params$lambda, traj_model@params$seed, traj_model@model, ignore_variables = colnames(traj_model@additional_features), r2_threshold = r2_threshold)

    traj_model@model <- res$model
    traj_model@motif_models <- traj_model@motif_models[res$vars]
    traj_model@predicted_diff_score <- res$pred
    traj_model@model_features <- res$X
    traj_model@coefs <- get_model_coefs(res$model)
    traj_model@normalized_energies <- traj_model@normalized_energies[, res$vars, drop = FALSE]
    traj_model@features_r2 <- res$vars_r2

    cli_alert_success("After filtering: Number of non-zero coefficients: {.val {sum(traj_model@model$beta != 0)}} (out of {.val {ncol(traj_model@model_features)}}). R^2: {.val {cor(traj_model@predicted_diff_score, norm01(traj_model@diff_score))^2}}")

    return(traj_model)
}

filter_traj_model_using_coefs <- function(traj_model, n_motifs) {
    res <- filter_model_using_coefs(traj_model@model_features, traj_model@coefs, norm01(traj_model@diff_score), traj_model@params$alpha, traj_model@params$lambda, traj_model@params$seed, traj_model@model, ignore_variables = colnames(traj_model@additional_features), n_motifs = n_motifs)

    traj_model@model <- res$model
    traj_model@predicted_diff_score <- res$pred
    traj_model@model_features <- res$X
    traj_model@coefs <- get_model_coefs(res$model)
    traj_model@normalized_energies <- traj_model@normalized_energies[, res$vars, drop = FALSE]

    cli_alert_success("After filtering: Number of non-zero coefficients: {.val {sum(traj_model@model$beta != 0)}} (out of {.val {ncol(traj_model@model_features)}}). R^2: {.val {cor(traj_model@predicted_diff_score, norm01(traj_model@diff_score))^2}}")

    return(traj_model)
}
