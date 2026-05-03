compute_spat_pwm <- function(pssm, intervals = NULL, size = NULL, sequences = NULL, bidirect = TRUE, min_energy = -7, q = 1) {
    if (is.null(sequences)) {
        if (is.null(intervals)) {
            cli_abort("Either {.field {intervals}} or {.field {sequences}} must be provided.")
        }

        if (!is.null(size)) {
            intervals <- gintervals.normalize(intervals %>% select(any_of(c("chrom", "start", "end", "strand"))), size)
        }

        sequences <- prego::intervals_to_seq(intervals)
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


compute_pssm_spatial_freq <- function(pssm, intervals = NULL, size = NULL, pwm_threshold = 7, sequences = NULL, atac_track = NULL, k4me3_track = NULL, k27me3_track = NULL, k27ac_track = NULL, orient_to_intervals = NULL, align_to_max = TRUE, ...) {
    if (!is.null(orient_to_intervals)) {
        if (is.null(intervals)) {
            cli_abort("Intervals must be provided when orienting to intervals.")
        }

        intervals <- gintervals.neighbors(intervals, orient_to_intervals, mindist = 1, maxneighbors = 1) %>%
            mutate(strand = ifelse(dist < 0, -1, 1)) %>%
            select(chrom, start, end, strand)
    }

    if (is.null(sequences)) {
        if (is.null(intervals)) {
            cli_abort("Either {.field {intervals}} or {.field {sequences}} must be provided.")
        }

        if (!is.null(size)) {
            intervals <- gintervals.normalize(intervals, size)
        }

        sequences <- prego::intervals_to_seq(intervals)

        seq_l <- stringr::str_length(sequences)
        if (any(seq_l != size)) {
            cli_abort("Some sequences are at the edge of the chromosome and therefore cannot be extended to the size {.val {size}}. Please make sure the 'end' column is at least {.val {ceiling(size/2)}}bp from the chromosome end.")
        }
    }

    local_pwm <- prego::compute_local_pwm(sequences, pssm, bidirect = TRUE)

    local_pwm_n <- norm_energy(local_pwm, min_energy = -10, q = 1)

    freqs <- local_pwm_n >= pwm_threshold

    spat_freq <- colMeans(freqs, na.rm = TRUE)

    n <- nrow(freqs)
    n_hits <- sum(freqs, na.rm = TRUE)

    res <- tibble::tibble(pos = seq_along(spat_freq), freq = spat_freq, n = n, n_hits = n_hits)

    if (!is.null(atac_track)) {
        if (is.null(intervals)) {
            cli_abort("Intervals must be provided when computing ATAC-seq frequency.")
        }
        if (is.null(size)) {
            size <- intervals$end[1] - intervals$start[1]
        }

        if (align_to_max) {
            # take only intervals with an occurence of the motif
            pwm_maxs <- apply(local_pwm_n, 1, max, na.rm = TRUE)
            atac_intervals <- intervals[pwm_maxs >= pwm_threshold, ]

            # align the intervals to the maximum in every sequence
            max_pwms <- apply(local_pwm_n[pwm_maxs >= pwm_threshold, ], 1, which.max)
            atac_intervals <- atac_intervals %>%
                mutate(start = start + max_pwms) %>%
                gintervals.normalize(size) %>%
                select(chrom, start, end)
        } else {
            atac_intervals <- gintervals.normalize(intervals, size) %>%
                select(chrom, start, end)
        }

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

    dinucs <- prego::calc_sequences_dinuc_dist(sequences)

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

calc_track_pos_data <- function(track, intervals, threshold = 7, direction = NULL) {
    withr::local_options(list(gmultitasking = FALSE, gmax.data.size = 1e7))
    gvtrack.create("track", track, "global.percentile.max")
    on.exit(gvtrack.rm("track"), add = TRUE)
    chip_data <- gextract(c("-log2(1-track)"), iterator = 1, intervals = intervals, colnames = "track") %>%
        arrange(intervalID) %>%
        mutate(pos = start - intervals$start[intervalID] + 1)

    if (!is.null(direction)) {
        # reverse the positions if the direction is negative
        chip_data <- chip_data %>%
            mutate(pos = ifelse(direction[intervalID] == -1, intervals$end[intervalID] - start + 1, pos))
    }

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
#' @param pwm_q_threshold The genomic quantile to use for the PWM threshold. Would be used if \code{pwm_threshold} is NULL.
#' @param top_q The proportion of top peaks to select
#' @param bottom_q The proportion of bottom peaks to select
#' @param bidirect_size Size of the intervals to use for deciding the directionality of the sequence
#' @param k4me3_track name of k4me3 track
#' @param k27me3_track name of k27me3 track
#' @param k27ac_track name of k27ac track
#' @param orient_to_intervals A data frame containing the intervals to orient the sequences to
#' @param align_to_max A logical indicating whether to align the sequences to the maximum in each sequence at the epigenetic tracks and atac signal
#' @param atac_track name of ATAC-seq marginal track
#' @param parallel Whether to use parallel processing
#' @param motifs A vector of motif names to compute spatial frequency for. Default is all motifs in the trajectory model.
#'
#' @return A data frame with the spatial frequency of each motif
#'
#' @export
compute_traj_model_spatial_freq <- function(traj_model, size, pwm_threshold = 7, pwm_q_threshold = 0.99, top_q = 0.1, bottom_q = 0.1, atac_track = NULL, parallel = TRUE, bidirect_size = NULL, k4me3_track = NULL, k27me3_track = NULL, k27ac_track = NULL, orient_to_intervals = NULL, align_to_max = TRUE, motifs = names(traj_model@motif_models)) {
    intervals <- traj_model@peak_intervals

    # select top and bottom 10% of peaks using diff_score
    intervals <- intervals %>%
        mutate(diff_score = traj_model@diff_score) %>%
        arrange(desc(diff_score)) %>%
        mutate(type = ifelse(row_number() <= nrow(intervals) * top_q, "top", ifelse(row_number() > nrow(intervals) * (1 - bottom_q), "bottom", "middle"))) %>%
        filter(type != "middle") %>%
        select(-diff_score)

    spatial_freqs <- plyr::ldply(motifs, function(motif) {
        if (is.null(pwm_threshold)) {
            local_pwm_r <- gextract.local_pwm(traj_model@normalization_intervals, traj_model@motif_models[[motif]]$pssm, bidirect = TRUE)
            local_pwm_rn <- norm_energy(local_pwm_r, min_energy = -10, q = 1)

            pwm_threshold <- quantile(local_pwm_rn, pwm_q_threshold, na.rm = TRUE)
            cli::cli_alert("PWM threshold for {.field {motif}}: {.val {pwm_threshold}}")
        }
        cli::cli_alert("Computing spatial frequency for {.val {motif}}")
        bind_rows(
            compute_pssm_spatial_freq(
                traj_model@motif_models[[motif]]$pssm,
                intervals = intervals %>% filter(type == "top"),
                size = size,
                pwm_threshold = pwm_threshold,
                pwm_q_threshold = pwm_q_threshold,
                atac_track = atac_track,
                k4me3_track = k4me3_track,
                k27me3_track = k27me3_track,
                k27ac_track = k27ac_track,
                orient_to_intervals = orient_to_intervals,
                align_to_max = align_to_max
            ) %>%
                mutate(type = "top"),
            compute_pssm_spatial_freq(
                traj_model@motif_models[[motif]]$pssm,
                intervals = intervals %>% filter(type == "bottom"),
                size = size,
                pwm_threshold = pwm_threshold,
                pwm_q_threshold = pwm_q_threshold,
                atac_track = atac_track,
                k4me3_track = k4me3_track,
                k27me3_track = k27me3_track,
                k27ac_track = k27ac_track,
                orient_to_intervals = orient_to_intervals,
                align_to_max = align_to_max
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
#' @param linewidth The width of the line in the plot.
#'
#' @return A ggplot object showing the spatial frequency of the given motif.
#'
#'
#' @export
plot_motif_spatial_freq <- function(spatial_freqs, motif, smooth = 10, plot_atac = FALSE, linewidth = 1) {
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
        geom_line(linewidth = linewidth) +
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
