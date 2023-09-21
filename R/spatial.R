compute_spat_pwm <- function(pssm, intervals = NULL, size = NULL, sequences = NULL, bidirect = TRUE, min_energy = -10, q = 1) {
    if (is.null(sequences)) {
        if (is.null(intervals)) {
            cli_abort("Either {.field {intervals}} or {.field {sequences}} must be provided.")
        }

        if (!is.null(size)) {
            intervals <- misha.ext::gintervals.normalize(intervals %>% select(any_of(c("chrom", "start", "end", "strand"))), size)
        }

        withr::local_options(list(gmax.data.size = 1e9))

        sequences <- toupper(misha::gseq.extract(intervals))
    }

    local_pwm <- prego::compute_local_pwm(sequences, pssm, bidirect = bidirect)

    local_pwm_n <- norm_energy(local_pwm, min_energy = min_energy, q = q)

    return(local_pwm_n)
}

direct_sequences <- function(sequences, pssm, bidi_seqs = sequences) {
    sequences_rc <- chartr("acgtACGT", "tgcaTGCA", sequences) %>%
        stringi::stri_reverse()

    size <- nchar(bidi_seqs[1])
    middle_point <- round(size / 2)

    s_l <- substr(bidi_seqs, 1, middle_point)
    s_r <- substr(bidi_seqs, (middle_point + 1), size)

    l_pwm <- prego::compute_pwm(s_l, pssm, bidirect = TRUE)
    r_pwm <- prego::compute_pwm(s_r, pssm, bidirect = TRUE)

    # for each sequence - if r > l reverse complement
    new_sequences <- ifelse(r_pwm > l_pwm, sequences_rc, sequences)

    return(new_sequences)
}

#' Compute motif directional hits
#'
#' This function computes the number of hits of a given motif in a set of intervals, separated by strand and direction.
#'
#' @param pssm A position-specific scoring matrix (PSSM) object.
#' @param intervals A data frame containing the genomic intervals to search for motif hits.
#' @param size The size of the intervals.
#' @param pwm_threshold The threshold for the PSSM score to consider a hit.
#'
#' @return A data frame containing the number of hits of the motif in each interval, separated by strand and direction. Column names are 'ml' (minus, left), 'mr' (minus, right), 'pl' (plus, left), 'pr' (plus, right).
#'
#' @export
compute_motif_directional_hits <- function(pssm, intervals, size, pwm_threshold = 7) {
    pwm_minus <- compute_spat_pwm(pssm, intervals %>% mutate(strand = -1), size, bidirect = TRUE)
    pwm_plus <- compute_spat_pwm(pssm, intervals %>% mutate(strand = 1), size, bidirect = TRUE)
    freqs_plus <- pwm_plus >= pwm_threshold
    freqs_minus <- pwm_minus >= pwm_threshold

    middle_point <- round(size / 2)

    # count # of hits at the left and at the right, for plus and minus.
    hits <- matrix(0, nrow = 4, ncol = nrow(intervals))
    hits[1, ] <- rowSums(freqs_minus[, (middle_point + 1):size], na.rm = TRUE)
    hits[2, ] <- rowSums(freqs_minus[, 1:middle_point], na.rm = TRUE)
    hits[3, ] <- rowSums(freqs_plus[, 1:middle_point], na.rm = TRUE)
    hits[4, ] <- rowSums(freqs_plus[, (middle_point + 1):size], na.rm = TRUE)

    hits <- t(hits)
    colnames(hits) <- c("ml", "mr", "pl", "pr")

    return(as.data.frame(hits))
}

annotate_intervals_diff_score <- function(intervals, diff_score, top_q = 0.1, bottom_q = 0.1) {
    intervals <- intervals %>%
        mutate(diff_score = !!diff_score) %>%
        arrange(desc(diff_score)) %>%
        mutate(type = ifelse(row_number() <= nrow(intervals) * top_q, "top", ifelse(row_number() > nrow(intervals) * (1 - bottom_q), "bottom", "middle"))) %>%
        filter(type != "middle") %>%
        select(-diff_score)
    return(intervals)
}

compute_motif_pwm_threshold <- function(intervals, pssm, q_thresh = 0.95, diff_score = NULL, top_q = 0.1, bottom_q = 0.1, size = 1000) {
    if (!("type" %in% colnames(intervals))) {
        intervals <- annotate_intervals_diff_score(intervals, diff_score = diff_score, top_q = top_q, bottom_q = bottom_q)
    }
    pwms <- compute_spat_pwm(pssm, intervals, size, bidirect = TRUE)

    top_pwms <- pwms[intervals$type == "top", ]
    bottom_pwms <- pwms[intervals$type == "bottom", ]

    if (mean(top_pwms, na.rm = TRUE) > mean(bottom_pwms, na.rm = TRUE)) {
        thresh <- quantile(bottom_pwms, q_thresh, na.rm = TRUE)
    } else {
        thresh <- quantile(top_pwms, q_thresh, na.rm = TRUE)
    }

    return(thresh)
}

compute_traj_model_directional_hits <- function(traj_model, size, pwm_quantile = 0.999, top_q = 0.1, bottom_q = 0.1, parallel = TRUE) {
    intervals <- traj_model@peak_intervals

    intervals <- annotate_intervals_diff_score(intervals, traj_model@diff_score, top_q, bottom_q)

    hits <- plyr::llply(names(traj_model@motif_models), function(motif) {
        pssm <- traj_model@motif_models[[motif]]$pssm
        cli::cli_alert("Computing directional hits for {.val {motif}}")
        pwm_threshold <- compute_motif_pwm_threshold(intervals, pssm, diff_score = traj_model@diff_score, q_thresh = pwm_quantile, size = size)
        cli::cli_alert_info("PWM threshold for {.val {motif}}: {.val {pwm_threshold}}")
        top <- compute_motif_directional_hits(
            pssm = pssm,
            intervals = intervals %>% filter(type == "top"),
            size = size,
            pwm_threshold = pwm_threshold
        )
        bottom <- compute_motif_directional_hits(
            pssm = pssm,
            intervals = intervals %>% filter(type == "bottom"),
            size = size,
            pwm_threshold = pwm_threshold
        )

        colnames(top) <- paste0(colnames(top), ".", motif)
        colnames(bottom) <- paste0(colnames(bottom), ".", motif)
        list(
            top = top,
            bottom = bottom,
            threshold = pwm_threshold
        )
    }, .parallel = parallel)

    top_mat <- do.call(cbind, lapply(hits, function(x) x$top))
    bottom_mat <- do.call(cbind, lapply(hits, function(x) x$bottom))
    thresholds <- sapply(hits, function(x) x$threshold)

    return(list(top = as.matrix(top_mat), bottom = as.matrix(bottom_mat), threshold = thresholds))
}


compute_pssm_spatial_freq <- function(pssm, intervals = NULL, size = NULL, pwm_threshold = 7, sequences = NULL, atac_track = NULL, bidirect_size = NULL, k4me3_track = NULL, k27me3_track = NULL, k27ac_track = NULL, ...) {
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

    if (!is.null(bidirect_size)) {
        if (is.null(intervals)) {
            cli_abort("intervals must be provided providing {.field bidirect_size}.")
        }

        sequences_d <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(intervals, bidirect_size)))
    } else {
        sequences_d <- sequences
    }

    sequences <- direct_sequences(sequences, pssm, bidi_seqs = sequences_d)

    local_pwm <- prego::compute_local_pwm(sequences, pssm, bidirect = TRUE)

    local_pwm_n <- norm_energy(local_pwm, min_energy = -10, q = 1)

    freqs <- local_pwm_n >= pwm_threshold

    spat_freq <- colMeans(freqs, na.rm = TRUE)

    n <- nrow(freqs)
    n_hits <- sum(freqs, na.rm = TRUE)

    res <- tibble::tibble(pos = 1:length(spat_freq), freq = spat_freq, n = n, n_hits = n_hits)

    if (!is.null(atac_track)) {
        if (is.null(intervals)) {
            cli_abort("Intervals must be provided when computing ATAC-seq frequency.")
        }
        if (is.null(size)) {
            size <- intervals$end[1] - intervals$start[1]
        }
        # align the intervals to the maximum in every sequence
        max_pwms <- apply(local_pwm_n, 1, which.max)
        atac_intervals <- intervals %>%
            mutate(start = start + max_pwms) %>%
            misha.ext::gintervals.normalize(size) %>%
            select(chrom, start, end)
        atac <- gextract(atac_track, iterator = 1, intervals = atac_intervals, colnames = "v")
        atac_mat <- atac %>%
            arrange(intervalID) %>%
            mutate(pos = start - atac_intervals$start[intervalID] + 1) %>%
            select(intervalID, pos, v) %>%
            tidyr::spread(pos, v) %>%
            select(-intervalID) %>%
            as.matrix()
        atac_freq <- colMeans(atac_mat, na.rm = TRUE)
        res$atac_freq <- atac_freq
    }

    dinucs <- calc_sequences_dinuc_dist(sequences)

    res <- res %>%
        left_join(dinucs, by = "pos")

    # add additional features
    if (!is.null(k4me3_track)) {
        res <- res %>%
            mutate(k4me3 = calc_track_pos_data(k4me3_track, atac_intervals))
    }

    if (!is.null(k27me3_track)) {
        res <- res %>%
            mutate(k27me3 = calc_track_pos_data(k27me3_track, atac_intervals))
    }

    if (!is.null(k27ac_track)) {
        res <- res %>%
            mutate(k27ac = calc_track_pos_data(k27ac_track, atac_intervals))
    }

    return(res)
}

calc_track_pos_data <- function(track, intervals, threshold = 7) {
    withr::local_options(list(gmultitasking = FALSE, gmax.data.size = 1e7))
    gvtrack.create("track", track, "global.percentile.max")
    chip_data <- gextract(c("-log2(1-track)"), iterator = 1, intervals = intervals, colnames = "track") %>%
        arrange(intervalID) %>%
        mutate(pos = start - intervals$start[intervalID] + 1)

    track_mat <- chip_data %>%
        select(intervalID, pos, v = track) %>%
        tidyr::spread(pos, v) %>%
        select(-intervalID) %>%
        as.matrix()

    track_mat <- track_mat >= threshold
    track <- colMeans(track_mat, na.rm = TRUE)
    return(track)
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
#' @param bidirect_size Size of the intervals to use for deciding the directionality of the sequence
#' @param k4me3_track name of k4me3 track
#' @param k27me3_track name of k27me3 track
#' @param k27ac_track name of k27ac track
#'
#' @return A data frame with the spatial frequency of each motif
#'
#' @export
compute_traj_model_spatial_freq <- function(traj_model, size, pwm_threshold = 7, top_q = 0.1, bottom_q = 0.1, atac_track = NULL, parallel = TRUE, bidirect_size = NULL, k4me3_track = NULL, k27me3_track = NULL, k27ac_track = NULL) {
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
                bidirect_size = bidirect_size,
                k4me3_track = k4me3_track,
                k27me3_track = k27me3_track,
                k27ac_track = k27ac_track
            ) %>%
                mutate(type = "top"),
            compute_pssm_spatial_freq(
                traj_model@motif_models[[motif]]$pssm,
                intervals = intervals %>% filter(type == "bottom"),
                size = size,
                pwm_threshold = pwm_threshold,
                atac_track = atac_track,
                bidirect_size = bidirect_size,
                k4me3_track = k4me3_track,
                k27me3_track = k27me3_track,
                k27ac_track = k27ac_track
            ) %>%
                mutate(type = "bottom")
        ) %>% mutate(motif = motif)
    }, .parallel = parallel)

    return(spatial_freqs)
}

plot_epi_spatial_freq <- function(spatial_freqs, motif, mark, smooth = 10) {
    spatial_freqs <- spatial_freqs %>%
        filter(motif == !!motif)

    spatial_freqs$mark <- spatial_freqs[, mark]

    spatial_freqs <- spatial_freqs %>%
        mutate(p = zoo::rollmean(mark, smooth, fill = NA, align = "center"))


    p <- ggplot(spatial_freqs, aes(x = pos, y = p, color = type)) +
        geom_line() +
        scale_color_manual(name = "", values = c(top = "red", bottom = "blue")) +
        theme_classic() +
        theme(legend.position = "none") +
        xlab("Position") +
        ylab("Frequency") +
        ggtitle(mark)

    return(p)
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
        scale_color_manual(name = "", values = c(top = "red", bottom = "blue")) +
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
        p <- p + annotate("text", x = 1, y = 0.9 * max(spatial_freqs$freq_roll, na.rm = TRUE), label = glue::glue("N = {n_hits_top} ({round(n_hits_top / n_top, 2)})"), hjust = 0, size = 3, color = "red")
        p <- p + annotate("text", x = 1, y = 0.75 * max(spatial_freqs$freq_roll, na.rm = TRUE), label = glue::glue("N = {n_bottom_hits} ({round(n_bottom_hits / n_bottom, 2)})"), hjust = 0, size = 3, color = "blue")
    }

    # dinucs plot
    # spatial_freqs %>%
    #     select(pos, AA:TT, type) %>%
    #     gather("dinuc", "p", -pos, -type) %>%
    #     mutate(p = zoo::rollmean(p, smooth, fill = NA, align = "center")) %>%
    #     as_tibble() %>%
    #     ggplot(aes(x = pos, y = p, color = type)) +
    #     geom_line() +
    #     facet_wrap(~dinuc) +
    #     theme_classic() +
    #     scale_color_manual(name = "", values = c(top = "red", bottom = "blue"))

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

compute_fold_enrichment <- function(mat) {
    mat[is.na(mat)] <- 0

    mat <- mat > 0

    # Get the probabilities of having a 1 in each column
    p_col <- colSums(mat) / nrow(mat)

    # Get the probability of having a 1 in both columns using matrix multiplication
    p_both <- (t(mat) %*% mat) / nrow(mat)

    # Get the outer product of the column probabilities to get the denominator for the fold enrichment calculation
    p_first_second <- tcrossprod(p_col)

    # Calculate fold enrichment
    fe_matrix <- log2(p_both / p_first_second)

    return(fe_matrix)
}
