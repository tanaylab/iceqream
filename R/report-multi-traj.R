#' Plot multi-trajectory model report
#'
#' This function plots a report for a multi-trajectory model.
#'
#' @param multi_traj The multi-trajectory model object.
#' @param filename Optional. The name of the output file. If not provided, the plot will be displayed on the screen.
#' @param width Optional. The width of the output plot in inches.
#' @param height Optional. The height of the output plot in inches.
#' @param dev Optional. The graphics device to use for output. Default is 'pdf'.
#' @param title Optional. The title of the plot.
#' @param use_full Optional. Whether to use the models with all the motifs or not. Default is FALSE.
#' @param names_map a named vector to map the names of the motifs to new names.
#' @param ... Additional arguments to be passed to the plotting function.
#'
#' @return A ggplot2 object (invsibly if a filename is provided).
#'
#'
#' @export
plot_multi_traj_model_report <- function(multi_traj, filename = NULL, width = NULL, height = NULL, dev = grDevices::pdf, title = NULL, use_full = FALSE, names_map = NULL, ...) {
    models <- multi_traj@motif_models

    motif_num <- length(models)

    if (use_full) {
        partial_resp_list <- purrr::map(multi_traj@models_full, compute_partial_response)
        traj_models <- multi_traj@models_full
    } else {
        partial_resp_list <- purrr::map(multi_traj@models, compute_partial_response)
        traj_models <- multi_traj@models
    }

    if (is.null(names_map)) {
        names_map <- names(models)
        names(names_map) <- names(models)
    }

    pr_mat <- partial_resp_list %>%
        purrr::imap_dfr(~ enframe(matrixStats::colMaxs(as.matrix(.x)), "var", "max_pr") %>% mutate(model = .y)) %>%
        select(var, max_pr, model) %>%
        filter(var %in% names(models)) %>%
        pivot_wider(names_from = c(model), values_from = max_pr) %>%
        tibble::column_to_rownames("var") %>%
        as.matrix()
    pr_mat[is.na(pr_mat)] <- 0
    cm <- tgs_cor(t(pr_mat), pairwise.complete.obs = TRUE)
    hc <- hclust(as.dist(1 - cm), method = "ward.D2")
    models <- models[rownames(pr_mat)]
    models <- models[hc$order]


    spatial_p <- purrr::imap(models, ~ {
        if (is.null(.x$spat)) {
            return(plot_fake_spat() + ggtitle("Spatial model", subtitle = ""))
        }
        prego::plot_spat_model(.x$spat) + ggtitle("Spatial model", subtitle = "")
    })

    motifs_p <- purrr::imap(models, ~ prego::plot_pssm_logo(.x$pssm) + ggtitle(names_map[.y], subtitle = ""))

    pr_p <- purrr::map(names(models), function(motif) {
        purrr::map(names(traj_models), function(model) {
            pr <- partial_resp_list[[model]]
            if (!(motif %in% names(pr))) {
                return(plot_fake_pr() + ggtitle(model))
            }
            if (motif %in% names(traj_models[[model]]@features_r2)) {
                subtitle <- round(traj_models[[model]]@features_r2[motif], 6)
            } else {
                subtitle <- ggplot2::waiver()
            }
            plot_e_vs_pr(motif, pr, traj_models[[model]]) + ggtitle(model, subtitle = subtitle)
        }) %>%
            purrr::discard(is.null) %>%
            patchwork::wrap_plots(nrow = 1)
    })
    names(pr_p) <- names(models)

    n_traj_models <- length(traj_models)

    p <- patchwork::wrap_plots(
        A = patchwork::wrap_plots(motifs_p, ncol = 1),
        B = patchwork::wrap_plots(spatial_p, ncol = 1),
        C = patchwork::wrap_plots(pr_p, ncol = 1),
        design = "ABC",
        widths = c(3, 0.8, n_traj_models)
    )

    if (!is.null(title)) {
        p <- p + patchwork::plot_annotation(title = title)
    }

    if (!is.null(filename)) {
        if (is.null(width)) {
            width <- 3 + 3 * n_traj_models
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

plot_fake_spat <- function() {
    prego::plot_spat_model(tibble(bin = seq(0, 400, 20), spat_factor = 1 / 21))
}

plot_fake_pr <- function() {
    ggplot() +
        labs(x = "Energy", y = "Partial response") +
        theme_classic() +
        theme(aspect.ratio = 1) +
        xlim(0, 10) +
        ylim(0, 1)
}
