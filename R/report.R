#' Plot scatter plot of observed vs predicted ATAC difference
#'
#' @param traj_model Trajectory model object. Please run \code{\link{regress_trajectory_motifs}} first.
#' @param type "train" or "test". If NULL - "test" would be used if available, otherwise "train".
#' @param point_size Size of the points.
#' @param alpha Transparency of the points.
#' @return ggplot2 object
#'
#'
#' @export
plot_prediction_scatter <- function(traj_model, type = NULL, point_size = 1, alpha = 0.5) {
    r <- get_obs_pred_df(traj_model, type = type)
    plot_df <- r$df
    type_str <- r$type_str

    limits <- c(min(plot_df$observed, plot_df$predicted), max(plot_df$observed, plot_df$predicted))

    p <- ggplot(plot_df, aes(x = observed, y = predicted)) +
        geom_point(alpha = alpha, size = point_size) +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        labs(x = "Observed ATAC difference", y = "Predicted ATAC difference") +
        xlim(limits) +
        ylim(limits) +
        ggtitle(paste0("R^2 = ", round(summary(lm(predicted ~ observed, data = plot_df))$r.squared, 2), type_str)) +
        theme_classic()

    return(p)
}

get_obs_pred_df <- function(traj_model, type = NULL) {
    if (!is.null(type) && !(type %in% c("train", "test"))) {
        cli_abort("{.field type} must be one of {.val train} or {.val test}")
    }

    type_str <- ""
    if (sum(traj_model@type == "test") == 0) {
        if (type == "test") {
            cli::cli_alert_warning("No test intervals found. Please run {.func infer_trajectory_motifs} first. Using train intervals instead.")
        }
        type <- "train"
    } else if (is.null(type)) {
        type <- "test"
    }

    if (type == "train") {
        type_str <- " (train)"
    }

    plot_df <- data.frame(
        observed = traj_model@diff_score[traj_model@type == type],
        predicted = traj_model@predicted_diff_score[traj_model@type == type]
    )

    return(list(df = plot_df, type_str = type_str))
}

#' Plot a boxplot of observed vs predicted ATAC differences
#'
#' This function takes a trajectory model and plots a boxplot of the observed vs predicted ATAC differences.
#' The plot is grouped by the observed ATAC difference, which is divided into n_groups quantiles.
#'
#' @param traj_model A trajectory model object
#' @param n_groups The number of groups to divide the observed ATAC difference into
#'
#' @return A ggplot object representing the boxplot
#'
#' @export
plot_prediction_boxplot <- function(traj_model, n_groups = 5) {
    r <- get_obs_pred_df(traj_model)
    plot_df <- r$df
    type_str <- r$type_str

    plot_df <- plot_df %>%
        mutate(obs_group = cut(observed, quantile(observed, seq(0, 1, length.out = n_groups + 1)), include.lowest = TRUE)) %>%
        gather(key = "type", value = "value", observed, predicted)

    p <- plot_df %>%
        ggplot(aes(x = obs_group, y = value, fill = type)) +
        geom_boxplot() +
        labs(x = "Observed ATAC difference", y = "ATAC difference") +
        scale_fill_manual(name = "", values = c(observed = "blue", predicted = "darkred")) +
        theme_classic() +
        theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))


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
    if (!(variable %in% traj_model@coefs$variable)) {
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

#' Plot Partial Response
#'
#' This function plots the partial response for a given trajectory model and motif.
#'
#' @param traj_model The trajectory model object.
#' @param motif The motif to plot the partial response for.
#' @param ylim The y-axis limits for the plot (optional).
#' @param xlab The label for the x-axis (optional, default is "Energy").
#' @param ylab The label for the y-axis (optional, default is "Partial response").
#' @param pointsize The size of the points in the scatter plot (optional, default is 3).
#'
#' @return A ggplot object representing the partial response plot.
#'
#' @export
plot_partial_response <- function(traj_model, motif, ylim = NULL, xlab = "Energy", ylab = "Partial response", pointsize = 3) {
    validate_traj_model(traj_model)
    pr <- compute_partial_response(traj_model)
    p <- tibble(
        e = traj_model@normalized_energies[, motif],
        pr = pr[, motif]
    ) %>%
        ggplot(aes(x = e, y = pr)) +
        scattermore::geom_scattermore(pointsize = pointsize) +
        labs(x = xlab, y = ylab) +
        theme_classic() +
        theme(aspect.ratio = 1)
    if (!is.null(ylim)) {
        p <- p + ylim(ylim)
    }
    return(p)
}

#' Plot Motif Energy vs Response Boxplot
#'
#' This function plots a boxplot of the ATAC difference (response) against the energy levels of a given motif.
#' If \code{subtitle=NULL}, the subtitle will be set to a kologorov-smirnov test between the lowest (0-3) and highest (9-10) energy levels.
#'
#' @param traj_model The trajectory model object.
#' @param motif The motif to plot the energy levels for.
#' @param xlab The label for the x-axis (default is the motif name followed by "energy").
#' @param ylab The label for the y-axis (default is "ATAC difference").
#' @param ylim The limits for the y-axis (default is -0.5 to 0.5).
#' @param fill The fill color for the boxplot (default is "lightblue1").
#' @param outliers Whether to plot outliers (default is TRUE).
#' @param title The title for the plot.
#' @param subtitle The subtitle for the plot (default is a kologorov-smirnov test between the lowest and highest energy levels). The color of the subtitle will be set to "darkred" if the p-value is less than 0.01.
#'
#' @return A ggplot object representing the boxplot.
#'
#'
#' @export
plot_motif_energy_vs_response_boxplot <- function(traj_model, motif, xlab = paste(motif, "energy"), ylab = "ATAC difference", ylim = c(-0.5, 0.5), fill = "lightblue1", title = "", subtitle = NULL, outliers = TRUE) {
    validate_traj_model(traj_model)
    if (!(motif %in% colnames(traj_model@normalized_energies))) {
        cli_abort("Motif {.val {motif}} not found in the model.")
    }

    plot_df <- tibble(
        observed = traj_model@diff_score,
        e = traj_model@normalized_energies[, motif]
    ) %>%
        mutate(energy = cut(e, c(0, 3, 6, 7, 8, 9, 10), include.lowest = TRUE, labels = c("0-3", "3-6", "6-7", "7-8", "8-9", "9-10")))

    subtitle_color <- "black"
    if (is.null(subtitle)) {
        if (sum(plot_df$energy == "0-3", na.rm = TRUE) >= 3 && sum(plot_df$energy == "9-10", na.rm = TRUE) >= 3) {
            ks <- stats::ks.test(plot_df$observed[plot_df$energy == "0-3"], plot_df$observed[plot_df$energy == "9-10"])
            subtitle <- glue::glue("KS.D={round(ks$statistic, digits = 2)} (pv={round(ks$p.value, digits = 3)})")
            if (ks$p.value <= 0.01) {
                subtitle_color <- "darkred"
            }
        } else {
            subtitle <- "Not enough samples in the lowest and highest energy levels"
        }
    }


    p <- plot_df %>%
        ggplot(aes(x = energy, y = observed)) +
        geom_hline(yintercept = 0) +
        geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.5, fill = fill, fatten = 1, linewidth = 0.5, outliers = outliers) +
        xlab(xlab) +
        ylab(ylab) +
        coord_cartesian(ylim = ylim) +
        ggtitle(title, subtitle = subtitle) +
        theme_classic() +
        theme(plot.subtitle = ggtext::element_markdown(color = subtitle_color)) +
        theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    return(p)
}



#' Plot a report of trajectory motifs
#'
#' @param traj_model Trajectory model object. Please run \code{regress_trajectory_motifs} first.
#' @param filename Filename to save the plot to. If NULL, the plot will be returned.
#' @param motif_num Number of motifs to plot. If NULL, all motifs will be plotted.
#' @param free_coef_axis Whether to use a free axis for the coefficient plots.
#' @param spatial_freqs Pre-computed spatial frequencies to plot. Use \code{compute_traj_model_spatial_freq} to compute.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @param dev Device to use for saving the plot. Default: \code{grDevices::pdf}.
#' @param ... Additional arguments to pass to the device.
#' @param title Title of the plot.
#' @param motif_titles Titles for the motifs. If NULL, the motif names will be used.
#' @param sort_motifs Whether to sort the motifs by the absolute value of the coefficients / r2 values.
#' @param names_map a named vector to map the names of the motifs to new names.
#' @param boxp_ylim ylimits for the boxplot of energy vs response.
#'
#' @return ggplot2 object. If filename is not NULL, the plot will be saved to the file and the function will return \code{invisible(NULL)}.
#'
#'
#' @export
plot_traj_model_report <- function(traj_model, filename = NULL, motif_num = NULL, free_coef_axis = TRUE, spatial_freqs = NULL, width = NULL, height = NULL, dev = grDevices::pdf, title = NULL, motif_titles = NULL, sort_motifs = TRUE, names_map = NULL, boxp_ylim = c(-0.5, 0.5), ...) {
    validate_traj_model(traj_model)
    models <- traj_model@motif_models
    has_features_r2 <- length(traj_model@features_r2) > 0 && all(names(models) %in% names(traj_model@features_r2))

    if (sort_motifs) {
        if (has_features_r2) {
            sorted_vars <- names(sort(traj_model@features_r2, decreasing = TRUE))
            sorted_vars <- sorted_vars[sorted_vars %in% names(models)]
        } else {
            sorted_vars <- traj_model@coefs %>%
                tibble::column_to_rownames("variable") %>%
                as.matrix() %>%
                abs() %>%
                apply(1, max) %>%
                sort(decreasing = TRUE) %>%
                names()
        }
        sorted_vars <- sorted_vars[!(sorted_vars %in% colnames(traj_model@additional_features))]
    } else {
        sorted_vars <- names(models)
    }

    if (!is.null(motif_num)) {
        if (motif_num > length(models)) {
            cli_abort("Motif number {.val {motif_num}} is larger than the number of motifs in the model ({.val {length(models)}})")
        }
    } else {
        motif_num <- length(sorted_vars)
    }

    cli_alert_info("Plotting {.val {motif_num}} motifs")
    models <- models[sorted_vars[1:motif_num]]

    plot_spat_model <- function(spat, title = NULL) {
        if (is.null(spat)) {
            plot_fake_spat()
        } else {
            prego::plot_spat_model(spat, title = title)
        }
    }

    if (has_features_r2) {
        spatial_p <- purrr::imap(models, ~ plot_spat_model(.x$spat, title = paste0("R^2=", round(traj_model@features_r2[.y], 6))))
    } else {
        spatial_p <- purrr::imap(models, ~ plot_spat_model(.x$spat))
    }

    if (is.null(names_map)) {
        names_map <- names(models)
        names(names_map) <- names(models)
    }


    if (is.null(motif_titles)) {
        motif_titles <- names_map[sorted_vars]
    } else {
        if (length(motif_titles) != length(models)) {
            cli_abort("Length of text must be equal to the number of motifs")
        }
    }

    motifs_p <- purrr::map2(models, motif_titles, ~ prego::plot_pssm_logo(.x$pssm, title = .y) +
        theme(plot.title = ggtext::element_markdown()))

    if (free_coef_axis) {
        coef_limits <- NULL
    } else {
        coef_limits <- c(min(as.matrix(traj_model@coefs[, c("low-energy", "sigmoid", "high-energy", "higher-energy")])), max(as.matrix(traj_model@coefs[, c("low-energy", "sigmoid", "high-energy", "higher-energy")])))
    }
    coefs_p <- purrr::map(names(models), ~ plot_coefs(traj_model, .x, limits = coef_limits, title = ""))

    if (is.null(spatial_freqs)) {
        spatial_freqs <- compute_traj_model_spatial_freq(traj_model, size = 1000, pwm_threshold = 7, top_q = 0.1, bottom_q = 0.1)
    }

    spat_freq_p <- purrr::map(names(models), ~ plot_motif_spatial_freq(spatial_freqs, .x, smooth = 10) + ggtitle(names_map[.x]))

    if ("atac_freq" %in% colnames(spatial_freqs)) {
        atac_spat_freq_p <- purrr::map(names(models), ~ plot_motif_spatial_freq(spatial_freqs, .x, smooth = 10, plot_atac = TRUE))
    } else {
        atac_spat_freq_p <- purrr::map(names(models), ~ ggplot() +
            theme_void())
    }

    pr <- compute_partial_response(traj_model)
    e_vs_pr_p <- purrr::map(names(models), ~ plot_e_vs_pr(.x, pr, traj_model))

    e_vs_r_boxp_p <- purrr::map(names(models), ~ plot_motif_energy_vs_response_boxplot(traj_model, .x, ylim = boxp_ylim, xlab = paste(names_map[.x], "energy"), outliers = FALSE))

    # scatter_p <- purrr::map(names(models), ~ plot_variable_vs_response(traj_model, .x, point_size = 0.001))

    p <- patchwork::wrap_plots(
        A = patchwork::wrap_plots(motifs_p, ncol = 1),
        B = patchwork::wrap_plots(e_vs_pr_p, ncol = 1),
        C = patchwork::wrap_plots(coefs_p, ncol = 1),
        D = patchwork::wrap_plots(spatial_p, ncol = 1),
        E = patchwork::wrap_plots(spat_freq_p, ncol = 1),
        F = patchwork::wrap_plots(e_vs_r_boxp_p, ncol = 1),
        G = patchwork::wrap_plots(atac_spat_freq_p, ncol = 1),
        design = "ABCDEFG",
        widths = c(0.5, 0.1, 0.1, 0.1, 0.3, 0.12, 0.3)
    )

    if (!is.null(title)) {
        p <- p + patchwork::plot_annotation(title = title)
    }

    if (!is.null(filename)) {
        if (is.null(width)) {
            width <- 27
        }
        if (is.null(height)) {
            height <- motif_num * 1.8
        }
        cli_alert_info("Saving plot...")
        dev(filename, width = width, height = height, ...)
        print(p)
        dev.off()
        cli_alert_success("Plot saved to {.file {filename}}")
        invisible(p)
    } else {
        return(p)
    }
}

plot_e_vs_pr <- function(motif, pr, traj_model, ylim = NULL) {
    p <- tibble(
        e = traj_model@normalized_energies[, motif],
        pr = pr[, motif]
    ) %>%
        ggplot(aes(x = e, y = pr)) +
        scattermore::geom_scattermore(pointsize = 3) +
        labs(x = "Energy", y = "Partial response") +
        theme_classic() +
        theme(aspect.ratio = 1)
    if (!is.null(ylim)) {
        p <- p + ylim(ylim)
    }
    return(p)
}

plot_coefs <- function(traj_model, variable, limits = NULL, title = variable, color = TRUE) {
    validate_traj_model(traj_model)
    coef_df <- traj_model@coefs %>% filter(variable == !!variable)

    if (nrow(coef_df) == 0) {
        cli_abort("Coefficient for variable {.val {variable}} not found.")
    }

    coef_df <- coef_df %>%
        gather("type", "value", -variable) %>%
        mutate(type = factor(type, levels = c("low-energy", "sigmoid", "high-energy", "higher-energy")))

    cli_alert_info("Plotting motif {.val {variable}}")
    if (color) {
        p <- ggplot(coef_df, aes(x = type, y = value, fill = type)) +
            geom_col() +
            scale_fill_manual(name = "", values = c("low-energy" = "blue", "high-energy" = "orange", "sigmoid" = "purple", "higher-energy" = "brown", "early" = "green", "linear" = "black")) +
            guides(fill = "none")
    } else {
        p <- ggplot(coef_df, aes(x = type, y = value)) +
            geom_col()
    }

    p <- p +
        labs(x = "", y = "Coefficient") +
        coord_flip() +
        theme_classic() +
        ggtitle(title)

    if (!is.null(limits)) {
        p <- p + ylim(limits)
    }

    return(p)
}

#' Plot trajectory model clusters report
#'
#' @param traj_model A trajectory model object
#' @param dir A directory to save the report files
#' @param k Number of clusters to split the heatmap into
#' @param spatial_freqs A vector of spatial frequencies to plot. Use \code{compute_traj_model_spatial_freq} to compute.
#'
#' @return None. A file called \code{heatmap.png} with the heatmap and a file called \code{motifs_report.pdf} with the motifs report will be saved to the directory.
#'
#' @export
plot_traj_model_clusters_report <- function(traj_model, dir, k = 10, spatial_freqs = NULL) {
    validate_traj_model(traj_model)
    e_mat <- traj_model@normalized_energies
    e_mat <- e_mat[, setdiff(colnames(e_mat), colnames(traj_model@additional_features))]
    cm <- tgs_cor(e_mat, pairwise.complete.obs = TRUE)
    hc <- hclust(as.dist(1 - cm), method = "complete")

    hm <- ComplexHeatmap::Heatmap(cm, name = "features", cluster_rows = hc, cluster_columns = hc, col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), split = k, column_split = k)

    if (dir.exists(dir)) {
        unlink(dir, recursive = TRUE)
    }

    dir.create(dir, showWarnings = FALSE, recursive = TRUE)

    png(file.path(dir, "heatmap.png"), width = 2000, height = 1000)
    if (length(traj_model@features_r2) > 0) {
        ha <- ComplexHeatmap::rowAnnotation(r2 = ComplexHeatmap::anno_points(traj_model@features_r2[rownames(cm)]), width = grid::unit(3, "cm"))
        hm <- ComplexHeatmap::draw(hm + ha, heatmap_legend_side = "left")
    } else {
        hm <- ComplexHeatmap::draw(hm, heatmap_legend_side = "left")
    }

    dev.off()

    clust_df <- ComplexHeatmap::row_order(hm) %>%
        purrr::imap_dfr(~ tibble(ord = .x, clust = .y)) %>%
        mutate(motif = rownames(cm)[ord])

    traj_model_clust <- traj_model
    traj_model_clust@motif_models <- traj_model@motif_models[clust_df$motif]

    cluster_colors <- chameleon::distinct_colors(n = length(unique(clust_df$clust)))$name
    titles <- glue::glue("<span style='color:{cluster_colors[clust_df$clust]}; font-weight: bold;'>C{clust_df$clust}: {clust_df$motif}</span>")

    plot_traj_model_report(traj_model_clust, spatial_freqs = spatial_freqs, filename = file.path(dir, "motifs_report.pdf"), motif_titles = titles, sort_motifs = FALSE)


    invisible(clust_df)
}
