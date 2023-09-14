compute_pssm_spatial_freq <- function(pssm, intervals = NULL, size = NULL, pwm_threshold = 7, sequences = NULL, atac_track = NULL) {
    if (is.null(sequences)) {
        if (is.null(intervals)) {
            cli_abort("Either {.field {intervals}} or {.field {sequences}} must be provided.")
        }

        if (!is.null(size)) {
            intervals <- misha.ext::gintervals.normalize(intervals, size)
        }

        withr::local_options(list(gmax.data.size = 1e9))

        sequences <- toupper(misha::gseq.extract(intervals))
    }

    local_pwm <- prego::compute_local_pwm(sequences, pssm)

    local_pwm_n <- norm_energy(local_pwm, min_energy = -10, q = 1)

    freqs <- local_pwm_n >= pwm_threshold

    spat_freq <- colMeans(freqs, na.rm = TRUE)

    n <- nrow(freqs)
    n_hits <- sum(freqs, na.rm = TRUE)

    if (!is.null(atac_track)) {
        atac <- gextract(atac_track, iterator = 1, intervals = intervals, colnames = "v")
        atac_mat <- atac %>%
            arrange(intervalID) %>%
            mutate(pos = start - intervals$start[intervalID] + 1) %>%
            select(intervalID, pos, v) %>%
            tidyr::spread(pos, v) %>%
            select(-intervalID) %>%
            as.matrix()
        atac_freq <- colMeans(atac_mat, na.rm = TRUE)
        return(
            tibble::tibble(pos = 1:length(spat_freq), freq = spat_freq, n = n, n_hits = n_hits, atac_freq = atac_freq)
        )
    }

    return(
        tibble::tibble(pos = 1:length(spat_freq), freq = spat_freq, n = n, n_hits = n_hits)
    )
}

#' Compute spatial frequency of motifs in trajectory model
#'
#' This function computes the spatial frequency of motifs in a trajectory model, using the top and bottom 10% of peaks based on diff_score.
#'
#' @param traj_model A trajectory model object
#' @param size The size of the region to compute spatial frequency for
#' @param pwm_threshold The threshold for the PWM score
#' @param top_q The proportion of top peaks to select
#' @param bottom_q The proportion of bottom peaks to select
#'
#' @return A data frame with the spatial frequency of each motif
#'
#' @export
compute_traj_model_spatial_freq <- function(traj_model, size, pwm_threshold = 7, top_q = 0.1, bottom_q = 0.1, atac_track = NULL, parallel = TRUE) {
    intervals <- traj_model@peak_intervals

    # select top and bottom 10% of peaks using diff_score
    intervals <- intervals %>%
        mutate(diff_score = traj_model@diff_score) %>%
        arrange(desc(diff_score)) %>%
        mutate(type = ifelse(row_number() <= nrow(intervals) * top_q, "top", ifelse(row_number() > nrow(intervals) * (1 - bottom_q), "bottom", "middle"))) %>%
        filter(type != "middle") %>%
        select(-diff_score)

    spatial_freqs <- plyr::ldply(names(traj_model@motif_models), function(motif) {
        cli::cli_alert("Computing spatial frequency for {.val {motif}}")
        bind_rows(
            compute_pssm_spatial_freq(
                traj_model@motif_models[[motif]]$pssm,
                intervals = intervals %>% filter(type == "top"),
                size = size,
                pwm_threshold = pwm_threshold,
                atac_track = atac_track,
            ) %>%
                mutate(type = "top"),
            compute_pssm_spatial_freq(
                traj_model@motif_models[[motif]]$pssm,
                intervals = intervals %>% filter(type == "bottom"),
                size = size,
                pwm_threshold = pwm_threshold,
                atac_track = atac_track
            ) %>%
                mutate(type = "bottom")
        ) %>% mutate(motif = motif)
    }, .parallel = parallel)

    return(spatial_freqs)
}

#' Plot motif spatial frequency
#'
#' This function plots the spatial frequency of a given motif across genomic positions.
#'
#' @param spatial_freqs A data frame containing the spatial frequency of motifs and/or ATAC-seq peaks. Output of \code{compute_traj_model_spatial_freq}.
#' @param motif A string specifying the motif to plot.
#' @param smooth An integer specifying the window size for smoothing the frequency values.
#' @param plot_atac A logical indicating whether to plot the ATAC-seq frequency instead of the motif frequency.
#'
#' @return A ggplot object showing the spatial frequency of the given motif.
#'
#'
#' @export
plot_motif_spatial_freq <- function(spatial_freqs, motif, smooth = 10, plot_atac = FALSE) {
    if (plot_atac) {
        if (!(("atac_freq" %in% colnames(spatial_freqs)))) {
            cli_abort("atac_freq column is missing from spatial_freqs")
        }
        spatial_freqs <- spatial_freqs %>%
            mutate(freq = atac_freq)
    }

    spatial_freqs <- spatial_freqs %>%
        filter(motif == !!motif) %>%
        mutate(freq_roll = zoo::rollmean(freq, smooth, fill = NA, align = "center"))

    p <- ggplot(spatial_freqs, aes(x = pos, y = freq_roll, color = type)) +
        geom_line() +
        scale_color_manual(name = "", values = c(bottom = "red", top = "blue")) +
        theme_classic() +
        theme(legend.position = c(0.9, 0.9)) +
        xlab("Position") +
        ggtitle(motif)

    if (plot_atac) {
        p <- p + ylab("ATAC (mean)")
    } else {
        p <- p + ylab("Motif frequency")
        n_top <- spatial_freqs %>%
            filter(type == "top") %>%
            slice(1) %>%
            pull(n)
        n_hits_top <- spatial_freqs %>%
            filter(type == "top") %>%
            slice(1) %>%
            pull(n_hits)
        n_bottom <- spatial_freqs %>%
            filter(type == "bottom") %>%
            slice(1) %>%
            pull(n)
        n_bottom_hits <- spatial_freqs %>%
            filter(type == "bottom") %>%
            slice(1) %>%
            pull(n_hits)

        # annotate with the N_hits + %
        p <- p + annotate("text", x = 1, y = 0.9 * max(spatial_freqs$freq_roll, na.rm = TRUE), label = glue::glue("N = {n_hits_top} ({round(n_hits_top / n_top, 2)})"), hjust = 0, size = 3, color = "blue")
        p <- p + annotate("text", x = 1, y = 0.75 * max(spatial_freqs$freq_roll, na.rm = TRUE), label = glue::glue("N = {n_bottom_hits} ({round(n_bottom_hits / n_bottom, 2)})"), hjust = 0, size = 3, color = "red")
    }

    return(p)
}

plot_all_motif_spatial_freq <- function(traj_model, size, pwm_threshold = 7, top_q = 0.1, bottom_q = 0.1, smooth = 10, filename = NULL, width = NULL, height = NULL, dev = grDevices::pdf, ...) {
    spatial_freqs <- compute_traj_model_spatial_freq(traj_model, size, pwm_threshold, top_q, bottom_q)

    motifs <- names(traj_model@motif_models)

    p <- purrr::map(motifs, ~ plot_motif_spatial_freq(spatial_freqs, .x, smooth = smooth)) %>%
        patchwork::wrap_plots(ncol = 1)

    if (!is.null(filename)) {
        if (is.null(width)) {
            width <- 10
        }
        if (is.null(height)) {
            height <- length(motifs) * 3
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
