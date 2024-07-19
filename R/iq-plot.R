#' Plot IQ Locus
#'
#' This function generates a plot of an IQ locus, including DNA sequence, motif responses, and ATAC-seq tracks.
#'
#' @param interval A genomic interval to plot.
#' @param pbm_list A list of PBM objects.
#' @param atac_tracks A vector of ATAC-seq track names to plot.
#' @param width Width of the plot in bp.
#' @param ext_width Width of the top plot in bp.
#' @param scale_cex Numeric, scaling factor for letter sizes.
#' @param T_emax Numeric, threshold for maximum energy.
#' @param bits_threshold Numeric, threshold for trimming PSSMs.
#' @param order_motifs Logical, whether to order motifs by maximum response.
#' @param atac_names Character vector, names for ATAC-seq tracks.
#' @param atac_colors Named vector, colors for ATAC-seq tracks.
#' @param atac_sizes Named vector, line sizes for ATAC-seq tracks.
#' @param line_thresh Numeric, threshold for drawing vertical lines.
#' @param score Numeric, optional score to display in the plot title.
#' @param norm_q Numeric, quantile for normalization.
#' @param tracks_q Data frame, pre-computed quantiles for tracks. Should have a 'type' column and a 'q' column.
#' @param iterator Numeric, iterator for gextract.
#' @param norm_intervals Genomic intervals for normalization.
#' @param atac_smooth Numeric, smoothing factor for ATAC-seq tracks.
#' @param ext_atac_smooth Numeric, smoothing factor for the top plot.
#' @param tn5bias_track Character, name of the Tn5 bias track.
#' @param normalize_tn5bias Logical, whether to normalize the TN5 bias track.
#' @param filename Character, output filename (if NULL, plot is not saved).
#' @param dev Function, device to use for plotting.
#' @param plot_width Numeric, width of the output plot.
#' @param plot_height Numeric, height of the output plot.
#'
#' @return A ggplot object containing the IQ locus plot.
#'
#' @export
plot_iq_locus <- function(interval, pbm_list, atac_tracks,
                          width = 500, ext_width = 5e4, T_emax = 8,
                          bits_threshold = NULL, order_motifs = TRUE, atac_names = atac_tracks, atac_colors = NULL,
                          atac_sizes = NULL,
                          line_thresh = 0.9,
                          score = NULL,
                          norm_q = 0.995,
                          tracks_q = NULL,
                          iterator = 20,
                          norm_intervals = gintervals.all(),
                          atac_smooth = 2,
                          ext_atac_smooth = atac_smooth * 5,
                          tn5bias_track = "tn5bias",
                          normalize_tn5bias = TRUE,
                          tss_intervals = "intervs.global.tss",
                          exon_intervals = "intervs.global.exon",
                          scale_cex = 500,
                          filename = NULL, dev = grDevices::pdf, plot_width = 15, plot_height = 8) {
    pbm_list <- preprocess_pbm_list(pbm_list, bits_threshold)
    interval <- preprocess_interval(interval, width)
    dna <- prego::intervals_to_seq(interval)

    energy_response_data <- compute_energy_response(pbm_list, dna, T_emax)
    e_mat <- energy_response_data$e_mat
    r_mat <- energy_response_data$r_mat

    max_e <- compute_max_energy(e_mat, max(purrr::map_dbl(pbm_list, ~ nrow(.x@pssm))))

    dna_df <- prepare_dna_data(dna, max_e, r_mat, e_mat, scale_cex, pbm_list, eps = 0.05)

    if (order_motifs) {
        dna_df <- order_motifs_data(dna_df, r_mat)
    }

    atac_data <- prepare_atac_data(atac_names, atac_tracks, interval, iterator, atac_smooth, norm_q, norm_intervals, tracks_q, tn5bias_track, normalize_tn5bias, ext_width = ext_width)

    diffs_df <- compute_diffs(dna_df, line_thresh, width)

    atac_colors <- prepare_atac_colors(atac_colors, atac_data)
    atac_sizes <- prepare_atac_sizes(atac_sizes, atac_names)

    p_atac <- plot_atac(atac_data, diffs_df, atac_colors, atac_sizes, width, atac_smooth)
    p_atac_ext <- plot_atac_ext(atac_data, atac_colors, atac_sizes, ext_width, width, ext_atac_smooth, interval = interval, tss_intervals = tss_intervals, exon_intervals = exon_intervals)

    title <- sprintf("%s:%d-%d", interval$chrom, interval$start, interval$end)
    if (!is.null(score)) {
        title <- paste0(title, ", ", sprintf("Score: %.2f", score))
    }
    # p_atac <- p_atac + labs(title = title)
    p_atac_ext <- p_atac_ext + labs(title = title)

    p_dna <- plot_dna(dna_df, diffs_df)
    p_logos <- plot_logos(pbm_list[rev(levels(dna_df$motif))])
    p_logos_rc <- plot_logos(pbm_list[rev(levels(dna_df$motif))], rc = TRUE)

    design <- "
        ##5
        ##1
        342
    "
    p <- p_atac + p_dna + p_logos + p_logos_rc + p_atac_ext + patchwork::plot_layout(design = design, heights = c(0.1, 0.55, 0.35), width = c(0.05, 0.05, 0.9), guides = "collect")


    if (!is.null(filename)) {
        save_plot(p, filename, dev, plot_width, plot_height)
    }

    return(p)
}

plot_logos <- function(pbm_list, rc = FALSE) {
    if (rc) {
        pssms <- purrr::map(pbm_list, ~ t(pssm_rc(.x@pssm)))
    } else {
        pssms <- purrr::map(pbm_list, ~ t(.x@pssm))
    }

    ggseqlogo::ggseqlogo(pssms, ncol = 1, method = "probability") +
        theme_void() +
        theme(
            strip.text.x = element_blank(),
            axis.text.y = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        ) +
        scale_y_discrete(expand = c(0, 0))
}

plot_atac <- function(atac_data, diffs_df, atac_colors, atac_sizes, l, atac_smooth) {
    p <- atac_data %>%
        group_by(type) %>%
        mutate(atac = zoo::rollmean(atac_n, atac_smooth, fill = "extend", na.rm = TRUE)) %>%
        ungroup() %>%
        filter(pos >= 1, pos <= l + 1) %>%
        ggplot(aes(x = pos, y = atac, color = type, linewidth = type)) +
        geom_hline(yintercept = 0, color = "black") +
        geom_segment(data = diffs_df, inherit.aes = FALSE, aes(x = pos, xend = pos, y = 1, yend = 0), size = 0.5, linetype = "dashed", color = "gray") +
        geom_segment(data = diffs_df, inherit.aes = FALSE, aes(x = pos, xend = lpos, y = 0, yend = -0.3), size = 0.5) +
        geom_segment(data = diffs_df, inherit.aes = FALSE, aes(x = lpos, xend = lpos, y = -0.3, yend = -0.5), size = 0.5) +
        geom_line() +
        scale_linewidth_manual(values = atac_sizes) +
        scale_color_manual(name = "", values = atac_colors) +
        guides(linewidth = "none") +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            axis.line.x = element_blank(),
            legend.position = "right"
        ) +
        scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.5, 1)) +
        labs(x = NULL, y = NULL, title = "") +
        scale_x_continuous(expand = c(0, 0)) +
        theme(plot.margin = margin(0, 0, 0, 0))

    return(p)
}

plot_atac_ext <- function(atac_data, atac_colors, atac_sizes, l_ext, l, atac_smooth, interval, tss_intervals, exon_intervals) {
    ext_interval <- gintervals.normalize(interval, l_ext)

    tss_data <- gintervals.neighbors1(tss_intervals, ext_interval) %>%
        filter(dist == 0) %>%
        mutate(ext_pos = start - start1 + 1) %>%
        filter(ext_pos > 0, ext_pos <= l_ext) %>%
        # add < or > to geneSymbol according to strand
        mutate(geneSymbol = ifelse(strand == 1, paste0(geneSymbol, "->"), paste0("<-", geneSymbol)))

    exon_data <- gintervals.neighbors1(exon_intervals, ext_interval) %>%
        filter(dist == 0) %>%
        mutate(start = start - start1 + 1) %>%
        mutate(end = end - start1 + 1) %>%
        filter(start > 0, start <= l_ext, end > 0, end <= l_ext)

    genes <- unique(exon_data$geneSymbol)
    n_genes <- length(genes)

    p <- atac_data %>%
        group_by(type) %>%
        mutate(atac = zoo::rollmean(atac_n, atac_smooth, fill = "extend", na.rm = TRUE)) %>%
        ungroup() %>%
        filter(ext_pos >= 1, ext_pos <= l_ext + 1) %>%
        ggplot(aes(x = ext_pos, y = atac, color = type)) +
        geom_rect(aes(xmin = 1 + (atac_data$ext_pos[1] - atac_data$pos[1]), xmax = l + (atac_data$ext_pos[1] - atac_data$pos[1]), ymin = 0, ymax = 1), fill = "transparent", color = "darkgray") +
        geom_line() +
        scale_linewidth_manual(values = atac_sizes) +
        scale_color_manual(name = "", values = atac_colors) +
        # TSS
        geom_segment(
            data = tss_data,
            aes(x = ext_pos, xend = ext_pos, y = 0, yend = 1.2),
            color = "black", linetype = "dashed"
        ) +
        geom_text(
            data = tss_data,
            inherit.aes = FALSE,
            aes(x = ext_pos, y = 1.4, label = geneSymbol),
            vjust = 0.5,
            size = 3
            # vjust = 1.5, hjust = 0.5, size = 3
        ) +
        geom_rect(
            data = exon_data,
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = 1.1, ymax = 1.2, fill = geneSymbol)
        ) +
        ggsci::scale_fill_aaas()

    p <- p +
        guides(linewidth = "none") +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            # axis.line.x = element_blank(),
            legend.position = "none"
        ) +
        scale_y_continuous(
            expand = c(0, 0),
            limits = c(0, 1.6),
            breaks = c(0, 0.5, 1)
        ) +
        # ylim(0, 1) +
        labs(x = NULL, y = NULL, title = "") +
        scale_x_continuous(expand = c(0, 0)) +
        theme(plot.margin = margin(0, 0, 0, 0))

    return(p)
}


plot_dna <- function(dna_df, diffs_df) {
    dna_df %>%
        ggplot(aes(x = letter_pos, y = motif, label = nuc, size = size, color = response)) +
        geom_segment(
            data = diffs_df,
            inherit.aes = FALSE,
            aes(
                x = letter_pos - 0.003,
                xend = letter_pos - 0.003,
                y = 0.7,
                yend = length(unique(dna_df$motif)) + 0.3
            ),
            size = 0.2,
            linetype = "dashed",
            color = "gray"
        ) +
        geom_text() +
        scale_size_identity() +
        scale_color_gradient2(name = "Response", low = "blue", mid = "darkgray", high = "red", midpoint = 0) +
        theme_minimal() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            legend.position = "right"
        ) +
        labs(x = NULL, y = NULL, title = "") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme(plot.margin = margin(0, 0, 0, 0))
}


save_plot <- function(p, filename, dev, plot_width, plot_height) {
    dev(filename, width = plot_width, height = plot_height)
    print(p)
    dev.off()
}

preprocess_pbm_list <- function(pbm_list, bits_threshold) {
    if (!is.null(bits_threshold)) {
        pbm_list <- pbm_list.trim_pssm(pbm_list, bits_threshold)
    }
    return(pbm_list)
}

preprocess_interval <- function(interval, width) {
    gintervals.normalize(interval, width)
}

compute_energy_response <- function(pbm_list, dna, T_emax) {
    cli_alert("Computing energies and responses for {.val {length(pbm_list)}} PBM models")

    energies <- pbm_list.compute_local(pbm_list, dna)
    resp <- pbm_list.compute_local(pbm_list, dna, response = TRUE)

    e_mat <- do.call(rbind, energies) %>% as.matrix()
    r_mat <- do.call(rbind, resp) %>% as.matrix()

    rownames(e_mat) <- names(pbm_list)
    rownames(r_mat) <- names(pbm_list)
    colnames(e_mat) <- colnames(r_mat) <- 1:nchar(dna)

    tf_maxs <- matrixStats::rowMaxs(e_mat, na.rm = TRUE)
    f <- tf_maxs >= T_emax
    cli_alert("{.val {sum(f)}} PBM models have at least one position with energy >= {.val {T_emax}}")

    list(e_mat = e_mat[f, , drop = FALSE], r_mat = r_mat[f, , drop = FALSE])
}

compute_max_energy <- function(e_mat, max_motif_l) {
    e_mat_exp <- 2**e_mat
    max_e <- apply(e_mat_exp, 2, max) / 10
    max_e[is.na(max_e)] <- 0
    max_e <- pmax(zoo::rollmean(max_e, max_motif_l, fill = "extend", na.rm = TRUE, align = "right"), 0)
    sqrt(max_e)
}

prepare_dna_data <- function(dna, max_e, r_mat, e_mat, scale_cex, pbm_list, eps) {
    l <- nchar(dna)
    pssm_l <- purrr::map_dbl(pbm_list, ~ nrow(.x@pssm))

    tot_cex <- sum(max_e + eps, na.rm = TRUE)
    max_e <- (max_e + eps) / tot_cex
    e_min <- quantile(max_e, 0.1, na.rm = TRUE)

    dna_df <- tibble(
        pos = 1:l,
        nuc = strsplit(dna, "")[[1]],
        e = max_e,
        size = e * scale_cex * 1.1
    ) %>%
        mutate(nuc = ifelse(max_e == e_min, "*", nuc)) %>%
        left_join(
            t(r_mat) %>%
                as.data.frame() %>%
                rownames_to_column("pos") %>%
                mutate(pos = as.numeric(pos)),
            by = join_by(pos)
        ) %>%
        gather("motif", "response", -pos, -nuc, -size, -e) %>%
        left_join(
            t(e_mat) %>%
                as.data.frame() %>%
                rownames_to_column("pos") %>%
                mutate(pos = as.numeric(pos)) %>%
                gather("motif", "energy", -pos),
            by = join_by(pos, motif)
        ) %>%
        mutate(response = ifelse(abs(energy) >= 5, response, 0)) %>%
        arrange(motif, pos) %>%
        group_by(motif) %>%
        mutate(response = zoo::rollmean(response, pssm_l[motif[1]], fill = "extend", na.rm = TRUE, align = "right")) %>%
        mutate(
            cumw = cumsum(c(1, size[-length(size)])),
            letter_pos = cumw,
            letter_pos = letter_pos / max(letter_pos)
        ) %>%
        ungroup()

    return(dna_df)
}

order_motifs_data <- function(dna_df, r_mat) {
    motif_ord <- names(sort(matrixStats::rowMaxs(r_mat, na.rm = TRUE)))
    dna_df %>% mutate(motif = factor(motif, levels = motif_ord))
}

#' Compute tracks quantiles
#'
#' This function computes the quantiles for a set of ATAC-seq tracks.
#'
#' @param atac_names A character vector specifying the names of the ATAC-seq tracks.
#' @param atac_tracks A list of numeric vectors representing the ATAC-seq tracks.
#' @param iterator An iterator object used for computing quantiles.
#' @param norm_q A numeric vector specifying the quantiles to compute.
#' @param norm_intervals A numeric vector specifying the intervals for computing quantiles.
#' @param tn5bias_track A numeric vector representing the TN5 bias track.
#' @param normalize_tn5bias Logical, whether to normalize the TN5 bias track.
#'
#' @return A data frame with columns "type" and "q", where "type" represents the name of the ATAC-seq track and "q" represents the computed quantiles.
#'
#'
#' @export
compute_tracks_q <- function(atac_names, atac_tracks, iterator = 20, norm_q = 0.995, norm_intervals = gintervals.all(), tn5bias_track = "tn5bias", normalize_tn5bias = TRUE) {
    gvtrack.create("bias", tn5bias_track, func = "sum")
    purrr::walk2(atac_names, atac_tracks, ~ gvtrack.create(.x, .y, func = "sum"))

    if (is.null(tn5bias_track) || !normalize_tn5bias) {
        exprs <- atac_names
    } else {
        exprs <- paste0(atac_names, "/bias")
    }

    purrr::map2_dfr(atac_names, exprs, ~
        tibble(
            type = .x,
            q = gquantiles(.y, iterator = iterator, percentiles = norm_q, intervals = norm_intervals)
        ))
}

prepare_atac_data <- function(atac_names, atac_tracks, interval, iterator, atac_smooth, norm_q, norm_intervals, tracks_q, tn5bias_track, normalize_tn5bias, ext_width) {
    orig_interval <- interval
    interval <- gintervals.normalize(interval, ext_width)
    cli::cli_alert("Preparing ATAC-seq data ({.val {length(atac_names)}} tracks)")
    gvtrack.create("bias", tn5bias_track, func = "sum")
    purrr::walk2(atac_names, atac_tracks, ~ gvtrack.create(.x, .y, func = "sum"))

    if (is.null(tracks_q)) {
        tracks_q <- compute_tracks_q(atac_names, atac_tracks, iterator, norm_q, norm_intervals, tn5bias_track, normalize_tn5bias)
    }

    if (is.null(tn5bias_track) || !normalize_tn5bias) {
        exprs <- atac_names
    } else {
        exprs <- paste0(atac_names, "/bias")
    }

    atac_data <- misha::gextract(exprs, gintervals.expand(interval, iterator * atac_smooth), iterator = iterator, colnames = atac_names) %>%
        select(-intervalID) %>%
        mutate(pos = start - orig_interval$start + 1) %>%
        mutate(ext_pos = start - interval$start + 1) %>%
        gather("type", "atac", -pos, -ext_pos, -chrom, -start, -end) %>%
        mutate(atac = ifelse(is.na(atac), 0, atac)) %>%
        left_join(tracks_q, by = "type") %>%
        mutate(atac_n = pmin(1, atac / q))


    return(atac_data)
}

compute_diffs <- function(dna_df, line_thresh, width) {
    dna_df %>%
        filter(motif == motif[1]) %>%
        arrange(pos) %>%
        mutate(
            d = abs(size - lag(size)),
            d = ifelse(is.na(d), 0, d)
        ) %>%
        filter(d > 0) %>%
        filter(d >= quantile(d, line_thresh)) %>%
        mutate(lpos = letter_pos * width - 0.003) %>%
        mutate(cumw = cumw / max(cumw)) %>%
        mutate(grp = cut(cumw, seq(0, 1, 0.02))) %>%
        arrange(grp, desc(d)) %>%
        distinct(grp, .keep_all = TRUE) %>%
        arrange(pos)
}

prepare_atac_colors <- function(atac_colors, atac_data) {
    if (is.null(atac_colors)) {
        atac_colors <- chameleon::distinct_colors(length(unique(atac_data$type)))$name
        names(atac_colors) <- unique(atac_data$type)
    } else if (is.null(names(atac_colors))) {
        names(atac_colors) <- unique(atac_data$type)
    }
    return(atac_colors)
}


prepare_atac_sizes <- function(atac_sizes, atac_names) {
    if (is.null(atac_sizes)) {
        atac_sizes <- rep(1, length(atac_names))
        names(atac_sizes) <- atac_names
    } else if (is.null(names(atac_sizes))) {
        names(atac_sizes) <- atac_names
    }
    return(atac_sizes)
}
