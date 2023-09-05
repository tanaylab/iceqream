#' Plot scatter plot of observed vs predicted ATAC difference
#'
#' @param traj_model Trajectory model object. Please run \code{\link{regress_trajectory_motifs}} first.
#' @return ggplot2 object
#'
#'
#' @export
plot_prediction_scatter <- function(traj_model) {
    type_str <- ""
    if (sum(traj_model@type == "test") == 0) {
        cli::cli_alert_warning("No test intervals found. Please run {.func infer_trajectory_motifs} first. Using train intervals instead.")
        plot_df <- data.frame(
            observed = traj_model@diff_score[traj_model@type == "train"],
            predicted = traj_model@predicted_diff_score[traj_model@type == "train"]
        )
        type_str <- " (train)"
    } else {
        plot_df <- data.frame(
            observed = traj_model@diff_score[traj_model@type == "test"],
            predicted = traj_model@predicted_diff_score[traj_model@type == "test"]
        )
    }

    limits <- c(min(plot_df$observed, plot_df$predicted), max(plot_df$observed, plot_df$predicted))

    p <- ggplot(plot_df, aes(x = observed, y = predicted)) +
        geom_point(alpha = 0.5) +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        labs(x = "Observed ATAC difference", y = "Predicted ATAC difference") +
        xlim(limits) +
        ylim(limits) +
        ggtitle(paste0("R^2 = ", round(summary(lm(predicted ~ observed, data = plot_df))$r.squared, 2), type_str)) +
        theme_classic()

    return(p)
}

#' Plot scatter of a variable vs ATAC difference
#'
#' @param traj_model Trajectory model object. Please run \code{\link{regress_trajectory_motifs}} first.
#' @param variable Variable to plot.
#' @param point_size Size of the points.
#'
#' @return ggplot2 object
#'
#' @export
plot_variable_vs_response <- function(traj_model, variable, point_size = 0.5) {
    validate_traj_model(traj_model)
    if (!(variable %in% traj_model_all@coefs$variable)) {
        cli_abort("Variable {.val {variable}} not found in the model.")
    }

    if (variable %in% colnames(traj_model@additional_features)) {
        plot_df <- data.frame(
            response = traj_model@diff_score,
            variable = traj_model@additional_features[[variable]]
        )
        p <- ggplot(plot_df, aes(x = variable, y = response)) +
            geom_point(alpha = 0.5) +
            labs(x = variable, y = "ATAC difference") +
            theme_classic()
    } else {
        plot_df <- purrr::map_dfr(c("early", "linear", "late"), ~ data.frame(
            response = traj_model@diff_score,
            variable = traj_model@normalized_energies[, paste0(variable, "_", .x)],
            type = .x
        ))

        p <- ggplot(aes(x = variable, y = response, color = type), data = plot_df) +
            geom_point(size = point_size, alpha = 0.5) +
            scale_color_manual(name = "", values = c(early = "red", linear = "green", late = "blue")) +
            labs(x = variable, y = "ATAC difference") +
            theme_classic()
    }

    return(p)
}

#' Plot a report of trajectory motifs
#'
#' @param traj_model Trajectory model object. Please run \code{regress_trajectory_motifs} first.
#' @param motif_num Number of motifs to plot. If NULL, all motifs will be plotted.
#' @param free_coef_axis Whether to use a free axis for the coefficient plots.
#' @param filename Filename to save the plot to. If NULL, the plot will be returned.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @param dev Device to use for saving the plot. Default: \code{grDevices::pdf}.
#' @param ... Additional arguments to pass to the device.
#'
#' @return ggplot2 object. If filename is not NULL, the plot will be saved to the file and the function will return \code{invisible(NULL)}.
#'
#'
#' @export
plot_motifs_report <- function(traj_model, motif_num = NULL, free_coef_axis = TRUE, filename = NULL, width = NULL, height = NULL, dev = grDevices::pdf, ...) {
    validate_traj_model(traj_model)

    models <- traj_model@motif_models

    sorted_vars <- traj_model@coefs %>%
        tibble::column_to_rownames("variable") %>%
        as.matrix() %>%
        abs() %>%
        apply(1, max) %>%
        sort(decreasing = TRUE) %>%
        names()

    sorted_vars <- sorted_vars[!(sorted_vars %in% traj_model@additional_features)]

    if (!is.null(motif_num)) {
        if (motif_num > length(models)) {
            cli_abort("Motif number {.val {motif_num}} is larger than the number of motifs in the model ({.val {length(models)}})")
        }
    } else {
        motif_num <- length(sorted_vars)
    }

    cli_alert_info("Plotting {.val {motif_num}} motifs")
    models <- models[sorted_vars[1:motif_num]]

    spatial_p <- purrr::imap(models, ~ prego::plot_spat_model(.x$spat))
    motifs_p <- purrr::imap(models, ~ prego::plot_pssm_logo(.x$pssm, title = .y))

    if (free_coef_axis) {
        coef_limits <- NULL
    } else {
        coef_limits <- c(min(as.matrix(traj_model@coefs[, c("early", "linear", "late")])), max(as.matrix(traj_model@coefs[, c("early", "linear", "late")])))
    }
    coefs_p <- purrr::map(names(models), ~ plot_coefs(traj_model, .x, limits = coef_limits, title = ""))

    # scatter_p <- purrr::map(names(models), ~ plot_variable_vs_response(traj_model, .x, point_size = 0.001))

    p <- patchwork::wrap_plots(
        A = patchwork::wrap_plots(motifs_p, ncol = 1),
        B = patchwork::wrap_plots(coefs_p, ncol = 1),
        C = patchwork::wrap_plots(spatial_p, ncol = 1),
        # D = patchwork::wrap_plots(scatter_p, ncol = 1),
        design = "ABC",
        widths = c(0.4, 0.2, 0.2)
    )

    if (!is.null(filename)) {
        if (is.null(width)) {
            width <- 9
        }
        if (is.null(height)) {
            height <- motif_num * 1.5
        }
        cli_alert_info("Saving plot to {.file {filename}}")
        dev(filename, width = width, height = height, ...)
        print(p)
        dev.off()
        cli_alert_success("Plot saved to {.file {filename}}")
        invisible(p)
    } else {
        return(p)
    }
}

plot_coefs <- function(traj_model, variable, limits = NULL, title = variable) {
    validate_traj_model(traj_model)
    coef_df <- traj_model@coefs %>% filter(variable == !!variable)

    if (nrow(coef_df) == 0) {
        cli_abort("Coefficient for variable {.val {variable}} not found.")
    }

    coef_df <- coef_df %>%
        gather("type", "value", -variable) %>%
        mutate(type = factor(type, levels = c("late", "linear", "early")))

    p <- ggplot(coef_df, aes(x = type, y = value)) +
        geom_col() +
        labs(x = "", y = "Coefficient") +
        coord_flip() +
        theme_classic() +
        ggtitle(title)

    if (!is.null(limits)) {
        p <- p + ylim(limits)
    }

    return(p)
}
