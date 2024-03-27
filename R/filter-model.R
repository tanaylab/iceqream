filter_model <- function(X, variables, y, alpha, lambda, seed, full_model, motif_models, ignore_variables = NULL, r2_threshold = 0.0005, bits_threshold = 1.75) {
    if (!is.null(ignore_variables)) {
        variables <- variables[!(variables %in% ignore_variables)]
    }

    full_model_r2 <- cor(logist(glmnet::predict.glmnet(full_model, newx = X, type = "link", s = lambda))[, 1], y)^2

    # for each variable of X calculate the r^2 of a GLM model without it
    var_stats <- plyr::llply(variables, function(var) {
        pssm <- motif_models[[var]]$pssm
        bits <- sum(prego::bits_per_pos(pssm))

        # cli_alert("Testing variable {.field {var}}...")
        X_f <- X[, grep(var, colnames(X), value = TRUE, invert = TRUE)]
        m <- glmnet::glmnet(X_f, y, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = FALSE, seed = seed)
        pred <- logist(glmnet::predict.glmnet(m, newx = X_f, type = "link", s = lambda))[, 1]
        r2 <- cor(pred, y)^2
        cli::cli_alert("R^2 added by {.field {var}} ({.strong {prego::consensus_from_pssm(pssm)}}): {.val {full_model_r2 - r2}}. Bits: {.val {bits}}")
        if (full_model_r2 - r2 < r2_threshold) {
            cli::cli_alert_info("Variable {.field {var}} removed due to low R^2")
        }
        if (bits < bits_threshold) {
            cli::cli_alert_info("Variable {.field {var}} removed due to low information content")
        }
        list(r2 = r2, bits = bits)
    }, .parallel = FALSE)

    vars_r2 <- purrr::map_dbl(var_stats, ~ .x$r2)
    names(vars_r2) <- variables

    vars_bits <- purrr::map_dbl(var_stats, ~ .x$bits)
    names(vars_bits) <- variables

    f_bits <- vars_bits > bits_threshold
    f_r2 <- (full_model_r2 - vars_r2) > r2_threshold

    cli_alert_info("Removed {.val {sum(!f_bits)}} features with bits < {.val {bits_threshold}}")
    cli_alert_info("Removed {.val {sum(!f_r2)}} features with R^2 < {.val {r2_threshold}}")

    vars_f <- variables[f_bits & f_r2]

    X_f <- X[, grep(paste0("(", paste(c(vars_f, ignore_variables), collapse = "|"), ")(_low-energy|_high-energy|_higher-energy|_sigmoid)"), colnames(X))]
    cli_alert_info("Number of features left: {.val {length(vars_f)}}")

    model_f <- glmnet::glmnet(X_f, y, binomial(link = "logit"), alpha = alpha, lambda = lambda, parallel = FALSE, seed = seed)
    pred_f <- logist(glmnet::predict.glmnet(model_f, newx = X_f, type = "link", s = lambda))[, 1]
    r2_f <- cor(pred_f, y)^2
    cli_alert_info("R^2 after filtering: {.val {r2_f}}. (Before: {.val {full_model_r2}})")

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

    return(list(model = model_f, pred = pred_f, X = X_f, r2 = r2_f, vars = vars_f))
}

#' Filter a trajectory model using punctuated regression.
#'
#' @description Run a regression without each feature and filter features that do not improve the model more
#' than \code{r2_threshold}.
#'
#' @param traj_model An instance of \code{TrajectoryModel}.
#' @param r2_threshold minimal R^2 for a feature to be included in the model.
#' @param bits_threshold minimal sum of bits for a feature to be included in the model.
#'
#' @return An instance of \code{TrajectoryModel} with the filtered model.
#'
#' @export
filter_traj_model <- function(traj_model, r2_threshold = 0.0005, bits_threshold = 1.75) {
    res <- filter_model(traj_model@model_features, traj_model@coefs$variable, norm01(traj_model@diff_score), traj_model@params$alpha, traj_model@params$lambda, traj_model@params$seed, traj_model@model, traj_model@motif_models, ignore_variables = colnames(traj_model@additional_features), r2_threshold = r2_threshold, bits_threshold = bits_threshold)

    traj_model@model <- res$model
    traj_model@motif_models <- traj_model@motif_models[res$vars]
    traj_model@predicted_diff_score <- res$pred
    traj_model@model_features <- res$X
    traj_model@coefs <- get_model_coefs(res$model)
    traj_model@normalized_energies <- traj_model@normalized_energies[, res$vars, drop = FALSE]
    traj_model@features_r2 <- res$vars_r2
    traj_model@params$r2_threshold <- r2_threshold
    traj_model@params$bits_threshold <- bits_threshold

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
