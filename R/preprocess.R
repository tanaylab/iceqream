import_marginal_tracks <- function(files, track_prefix, binsize = 20, overwrite = FALSE, overwrite_marginal = FALSE, parallel = getOption("prego.parallel", TRUE)) {
    gdir.create(track_prefix, showWarnings = FALSE)
    ct_names <- tools::file_path_sans_ext(basename(files))
    tracks <- paste0(track_prefix, ".", ct_names)
    plyr::l_ply(seq_along(files), function(i) {
        track <- tracks[i]
        if (gtrack.exists(track)) {
            if (!overwrite) {
                cli::cli_alert("Track {.val {track}} already exists. Skipping. Set {.field overwrite} to {.val TRUE} to overwrite.")
                return()
            }
            gtrack.rm(track, force = TRUE)
        }
        cli::cli_alert("Importing {.val {ct_names[i]}} to {.val {track}}")
        gtrack.import(track = tracks[i], file = files[i], binsize = binsize, description = paste("Marginal signal of", ct_names[i]))
    }, .parallel = parallel)

    marginal_track <- paste0(track_prefix, ".marginal")
    if (gtrack.exists(marginal_track)) {
        if (!overwrite_marginal) {
            cli::cli_alert("Marginal track {.val {marginal_track}} already exists. Skipping. Set {.field overwrite_marginal} to {.val TRUE} to overwrite.")
            return(c(tracks, marginal_track))
        }
        gtrack.rm(marginal_track, force = TRUE)
    }
    cli::cli_alert("Creating marginal track {.val {marginal_track}}")
    psum <- misha.ext::psum
    expr <- paste0("psum(", paste0(tracks, collapse = ", "), ", na.rm = TRUE)")
    gtrack.create(expr = expr, track = marginal_track, description = paste("Marginal signal of", paste(ct_names, collapse = ", ")), iterator = binsize)

    gdb.reload()

    invisible(c(tracks, marginal_track))
}

load_peaks <- function(peaks, peaks_size = NULL) {
    if (is.data.frame(peaks)) {
        cli::cli_alert("Using peaks data frame")
    } else if (is.character(peaks)) {
        cli::cli_alert("Reading peaks from {.val {peaks}}")
        if (tools::file_ext(peaks) == "bed") {
            cli::cli_alert("Reading peaks from BED file")
            peaks <- readr::read_tsv(peaks, show_col_types = FALSE, col_names = FALSE)
            colnames(peaks)[1:3] <- c("chrom", "start", "end")
        } else {
            peaks <- readr::read_tsv(peaks, show_col_types = FALSE)
        }
    } else {
        cli::cli_abort("Peaks must be a data frame or a path to a file")
    }

    if (!all(c("chrom", "start", "end") %in% colnames(peaks))) {
        cli::cli_abort("Peaks must have 'chrom', 'start', and 'end' columns")
    }
    cli::cli_alert_info("# Peaks: {.val {nrow(peaks)}}")

    if (!is.null(peaks_size)) {
        cli::cli_alert("Normalizing peaks to {.val {peaks_size}}bp")
        peaks <- misha.ext::gintervals.normalize(peaks, peaks_size)
    }

    if (nrow(peaks) != nrow(gintervals.canonic(peaks))) {
        cli::cli_alert_warning("There are overlapping peaks. Removing duplicates.")
        peaks <- gintervals.canonic(peaks)
        cli::cli_alert_info("# Peaks after removing duplicates: {.val {nrow(peaks)}}")
    }
    return(peaks)
}

#' Preprocess and normalize ATAC-seq data
#'
#' @description
#' This function processes ATAC-seq data by importing signal tracks, normalizing peaks,
#' and performing multiple normalization steps including regional normalization,
#' constitutive peak normalization, and probability-based normalization.
#'
#' @param project_name Character string. The prefix used for track names and project identification.
#' @param files Optional character vector. Paths to input ATAC-seq signal files, can be in bigWig or tsv format, see \code{misha::gtrack.import} for more details. Required if tracks don't exist.
#' @param cell_types Optional character vector. Names of cell types to process. If NULL, derived from track names.
#' @param peaks Data frame or file path. Peak intervals with required columns 'chrom', 'start', and 'end'.
#' @param anchor_cell_type Optional character. Cell type to use as reference for normalization. If NULL, the mean of all cell types is used.
#' @param figures_dir Optional character. Directory path to save normalization plots.
#' @param peaks_size Numeric. Size to normalize peaks to in base pairs. Default: 500
#' @param binsize Numeric. Bin size for signal track import in base pairs. Default: 20
#' @param overwrite_tracks Logical. Whether to overwrite existing individual cell type tracks. Default: FALSE
#' @param overwrite_marginal Logical. Whether to overwrite existing marginal track. Default: FALSE
#' @param window_size Numeric. Window size for regional normalization in base pairs. Default: 2e4
#' @param minimal_quantile Numeric. Minimum quantile for regional normalization. Default: 0.1
#' @param const_threshold Numeric. Log2 threshold for identifying constitutive peaks. Default: -16
#' @param const_norm_quant Numeric. Quantile for constitutive peak normalization. Default: 1
#' @param const_scaling_quant Numeric. Scaling quantile for constitutive normalization. Default: 1
#' @param const_quantile Numeric. Quantile for probability normalization threshold. Default: 0.9
#' @param prob1_thresh Optional numeric. Threshold for probability=1 in normalization.
#'   If NULL, calculated from const_quantile.
#' @param add_tss_dist Logical. Whether to add TSS distance to peaks. Default: TRUE
#' @param tss_intervals Character. Name of TSS intervals track. Default: "intervs.global.tss"
#'
#' @return A list containing:
#'   \itemize{
#'     \item atac: Raw ATAC signal matrix
#'     \item atac_norm: Region-normalized signal matrix
#'     \item atac_norm_const: Constitutive peak-normalized signal matrix
#'     \item atac_norm_prob: Probability-normalized signal matrix
#'     \item peaks: Data frame of peak information
#'     \item params: List of parameters used for normalization
#'   }
#'
#' @details
#' The function performs several normalization steps:
#' 1. Regional normalization using punctured windows around peaks
#' 2. Identification and normalization of constitutive peaks
#' 3. Conversion to probability scores
#'
#' If visualization is enabled (figures_dir is provided), the function generates
#' scatter plots showing the effects of each normalization step.
#'
#' @examples
#' \dontrun{
#' # Basic usage with existing tracks
#' result <- preprocess_data(
#'     project_name = "my_project",
#'     peaks = "peaks.bed",
#'     figures_dir = "figures"
#' )
#'
#' # Full preprocessing with new data
#' result <- preprocess_data(
#'     project_name = "my_project",
#'     files = c("celltype1.bw", "celltype2.bw"),
#'     peaks = "peaks.bed",
#'     anchor_cell_type = "celltype1",
#'     figures_dir = "figures",
#'     overwrite_tracks = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{normalize_regional}}, \code{\link{normalize_const}}, \code{\link{normalize_to_prob}}
#'
#' @export
preprocess_data <- function(project_name, files = NULL, cell_types = NULL, peaks = NULL, anchor_cell_type = NULL, figures_dir = NULL, peaks_size = 500, binsize = 20, overwrite_tracks = FALSE, overwrite_marginal = FALSE, window_size = 2e4, minimal_quantile = 0.1, const_threshold = -16, const_norm_quant = 1, const_scaling_quant = 1, const_quantile = 0.9, prob1_thresh = NULL, add_tss_dist = TRUE, tss_intervals = "intervs.global.tss") {
    track_prefix <- project_name

    if (!is.null(files)) {
        ct_names <- tools::file_path_sans_ext(basename(files))
    } else {
        ct_names <- cell_types %||% gsub(paste0("^", track_prefix, "\\."), "", gtrack.ls(paste0(track_prefix, "\\.")))
        cli::cli_alert("Using cell types: {.val {ct_names}}")
    }

    tracks <- paste0(track_prefix, ".", ct_names)
    marginal_track <- paste0(track_prefix, ".marginal")

    if (all(gtrack.exists(c(tracks, marginal_track))) && !overwrite_tracks && !overwrite_marginal) {
        if (!is.null(files)) {
            cli::cli_alert("All tracks exist. Set {.field overwrite_tracks} or {.field overwrite_marginal} to {.val TRUE} to re-import.")
        }
    } else {
        if (is.null(files)) {
            cli::cli_abort("{.field files} must be provided if the tracks do not exist.")
        }
        tracks <- import_marginal_tracks(files, track_prefix, binsize = binsize, overwrite = overwrite_tracks, overwrite_marginal = overwrite_marginal)
    }

    # test that the tracks were created
    missing_tracks <- tracks[!gtrack.exists(tracks)]
    if (length(missing_tracks) > 0) {
        cli::cli_abort("Failed to create the following tracks: {.val {missing_tracks}}")
    }

    if (is.null(peaks)) {
        # TODO: create peaks
        cli::cli_abort("Peaks must be provided")
    } else {
        peaks <- load_peaks(peaks, peaks_size = peaks_size)
    }

    cli::cli_alert("Loading ATAC data")
    purrr::walk2(tracks, ct_names, ~ gvtrack.create(.y, .x, func = "sum"))
    atac_data <- gextract(ct_names, iterator = peaks, intervals = peaks) %>%
        arrange(intervalID) %>%
        select(-intervalID)

    atac_mat <- misha.ext::intervs_to_mat(atac_data)
    atac_mat[is.na(atac_mat)] <- 0

    cli::cli_alert("Normalizing regional effect using a punctured window of {.val {window_size}}bp. Minimal quantile: {.val {minimal_quantile}}")
    norm_mat <- normalize_regional(peaks, atac_mat, marginal_track = marginal_track, window_size = window_size, minimal_quantile = minimal_quantile)

    legc <- as.matrix(log2(norm_mat + 1e-5))
    const_peaks <- matrixStats::rowMins(legc, na.rm = TRUE) >= const_threshold
    peaks <- peaks %>% mutate(const = const_peaks)

    cli::cli_alert_info("There are {.val {sum(const_peaks)}} constant peaks above the threshold of {.val {const_threshold}} in all cell types")

    if (!is.null(anchor_cell_type)) {
        if (!anchor_cell_type %in% ct_names) {
            cli::cli_abort("Anchor cell type {.val {anchor_cell_type}} not found in the list of cell types: {.val {ct_names}}")
        }
        cli::cli_alert("Using {.val {anchor_cell_type}} as the anchor cell type")
    }

    cli::cli_alert_info("Normalizing constant peaks using the anchor cell type: {.val {anchor_cell_type}}. Normalization quantile: {.val {const_norm_quant}}. Scaling quantile: {.val {const_scaling_quant}}")
    norm_mat_const <- normalize_const(peaks, norm_mat, anchor_cell_type, norm_quant = const_norm_quant, scaling_quant = const_scaling_quant)
    legc_const <- as.matrix(log2(norm_mat_const + 1e-5))

    prob1_thresh <- prob1_thresh %||% quantile(norm_mat_const[peaks$const, ], const_quantile)

    norm_mat_p <- normalize_to_prob(peaks, norm_mat_const, prob1_thresh = NULL, const_quantile = const_quantile)

    if (add_tss_dist) {
        peaks <- peaks %>%
            misha.ext::gintervals.neighbors1(tss_intervals) %>%
            select(chrom:const, tss_dist = dist)
    }

    peaks <- peaks %>%
        mutate(peak_name = paste0(chrom, "_", start, "_", end))
    stopifnot(all(peaks$peak_name == rownames(norm_mat_p)))

    obj <- list(
        atac = atac_mat,
        atac_norm = norm_mat,
        atac_norm_const = norm_mat_const,
        atac_norm_prob = norm_mat_p,
        peaks = peaks,
        params = list(
            window_size = window_size,
            minimal_quantile = minimal_quantile,
            const_threshold = const_threshold,
            const_norm_quant = const_norm_quant,
            const_scaling_quant = const_scaling_quant,
            const_quantile = const_quantile,
            prob1_thresh = prob1_thresh,
            anchor_cell_type = anchor_cell_type
        )
    )

    if (!is.null(figures_dir)) {
        plot_normalization_scatters(obj, anchor_cell_type, figures_dir)
    }

    return(obj)
}

#' Generate and save normalization visualization plots
#'
#' @description
#' Creates a series of scatter plots showing the progression of ATAC-seq signal normalization
#' steps, comparing each cell type against an anchor cell type. The plots are automatically
#' saved to the specified directory if provided.
#'
#' @param obj List. A preprocessed ATAC-seq data object, output from \code{\link{preprocess_data}}
#' @param anchor_cell_type Character. The reference cell type to use for comparisons
#' @param figures_dir Character or NULL. Directory path where plots should be saved.
#'   If NULL, plots are not saved to disk. Directory will be created if it doesn't exist.
#'
#' @details
#' The function generates four visualization plots showing the progression of normalization:
#'
#' 1. Raw signal comparison ("1_cell_type_scatter_before_norm.png")
#'    * Log2-transformed raw signals after basic library size normalization
#'
#' 2. Regional normalization ("2_cell_type_scatter_region_norm.png")
#'    * Shows effects of local background correction
#'
#' 3. Constitutive peak normalization ("3_cell_type_scatter_const_norm.png")
#'    * Displays signals after normalizing using constitutive peaks
#'    * Includes probability threshold line if specified
#'
#' 4. Probability normalization ("4_cell_type_scatter_prob.png")
#'    * Shows final probability-scaled signals
#'
#' Each plot compares all cell types against the anchor cell type, with points
#' colored to distinguish constitutive from variable peaks.
#'
#' @note
#' File names are automatically generated using a numerical prefix to maintain
#' proper ordering when viewing in a file explorer.
#'
#' @return None (called for side effects)
#'
#' @examples
#' \dontrun{
#' # Generate and save plots
#' plot_normalization_scatters(
#'     atac_obj,
#'     anchor_cell_type = "CD4_T",
#'     figures_dir = "analysis/normalization_plots"
#' )
#'
#' # Generate plots without saving
#' plot_normalization_scatters(
#'     atac_obj,
#'     anchor_cell_type = "CD4_T",
#'     figures_dir = NULL
#' )
#' }
#'
#' @seealso
#' \code{\link{plot_cell_type_scatter}} for the underlying plotting function,
#' \code{\link{normalize_regional}}, \code{\link{normalize_const}},
#' \code{\link{normalize_to_prob}} for the normalization methods being visualized
#'
#' @export
plot_normalization_scatters <- function(obj, anchor_cell_type, figures_dir) {
    if (!is.null(figures_dir)) {
        dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
        cli::cli_alert("Saving figures to {.val {figures_dir}}")
    }

    const_peaks <- obj$peaks$const

    legc_raw <- t(t(obj$atac) / colSums(obj$atac, na.rm = TRUE))
    legc_raw <- as.matrix(log2(legc_raw + 1e-5))
    plot_cell_type_scatter(legc_raw, anchor_cell_type, const_peaks, filename = file.path(figures_dir, "1_cell_type_scatter_before_norm.png"))

    legc <- as.matrix(log2(obj$atac_norm + 1e-5))
    plot_cell_type_scatter(legc, anchor_cell_type, const_peaks, filename = file.path(figures_dir, "2_cell_type_scatter_region_norm.png"), ylab = "ATAC signal (region normalized)")

    legc_const <- as.matrix(log2(obj$atac_norm_const + 1e-5))
    plot_cell_type_scatter(legc_const, anchor_cell_type, const_peaks, filename = file.path(figures_dir, "3_cell_type_scatter_const_norm.png"), ylab = "ATAC signal (const normalized)", prob1_thresh = obj$params$prob1_thresh)

    plot_cell_type_scatter(obj$atac_norm_prob, anchor_cell_type, const_peaks, filename = file.path(figures_dir, "4_cell_type_scatter_prob.png"), ylab = "ATAC signal (probability)")
}

#' Create comparison scatter plots showing normalization effects between two cell types
#'
#' @description
#' Creates a series of scatter plots comparing ATAC-seq signals between two cell types
#' at each stage of the normalization process. The plots show raw signals, regional
#' normalization, constitutive peak normalization, and probability normalization in
#' a single combined figure.
#'
#' @param obj List. A preprocessed ATAC-seq data object, output from \code{\link{preprocess_data}}
#' @param cell_type1 Character. Name of the first cell type to plot (x-axis)
#' @param cell_type2 Character. Name of the second cell type to plot (y-axis)
#' @param peaks_f Logical vector or NULL. Optional filter for selecting specific peaks.
#'   If NULL, all peaks are used. Default: NULL
#'
#' @details
#' The function creates four scatter plots showing:
#' 1. Raw ATAC signal (log2 transformed)
#' 2. Region-normalized signal
#' 3. Constitutive peak-normalized signal
#' 4. Probability-normalized signal
#'
#' Each plot includes:
#' * Points colored by peak type (constitutive vs variable)
#' * A diagonal reference line
#'
#' @return A ggplot2 object containing four scatter plots arranged horizontally using
#'   patchwork, showing the progression of normalization effects
#'
#' @examples
#' \dontrun{
#' # Basic usage comparing two cell types
#' p <- plot_cell_type_normalization_scatter(
#'     atac_obj,
#'     cell_type1 = "CD4_T",
#'     cell_type2 = "CD8_T"
#' )
#'
#' # Plot specific peaks using a filter
#' p <- plot_cell_type_normalization_scatter(
#'     atac_obj,
#'     cell_type1 = "CD4_T",
#'     cell_type2 = "CD8_T",
#'     peaks_f = atac_obj$peaks$tss_dist < 1000
#' )
#' }
#'
#' @seealso
#' \code{\link{normalize_regional}}, \code{\link{normalize_const}},
#' \code{\link{normalize_to_prob}} for the normalization methods being visualized
#'
#' @export
plot_cell_type_normalization_scatter <- function(obj, cell_type1, cell_type2, peaks_f = NULL) {
    if (is.null(peaks_f)) {
        peaks_f <- rep(TRUE, nrow(obj$peaks))
    }
    obj$peaks <- obj$peaks[peaks_f, ]
    df_raw <- obj$atac[peaks_f, c(cell_type1, cell_type2)]
    df_raw <- t(t(df_raw) / colSums(df_raw, na.rm = TRUE))
    df_raw[is.na(df_raw)] <- 0
    df_norm <- obj$atac_norm[peaks_f, c(cell_type1, cell_type2)]
    df_norm_const <- obj$atac_norm_const[peaks_f, c(cell_type1, cell_type2)]
    df_norm_prob <- obj$atac_norm_prob[peaks_f, c(cell_type1, cell_type2)]

    plot_df <- function(df, title) {
        as.data.frame(df) %>%
            mutate(const = ifelse(obj$peaks$const, "const", "variable")) %>%
            ggplot(aes(x = !!sym(cell_type1), y = !!sym(cell_type2), color = const)) +
            scattermore::geom_scattermore(pointsize = 2) +
            geom_abline(linetype = "dashed", color = "darkblue") +
            scale_color_manual(name = "Peak Type", values = c("variable" = "black", "const" = "red")) +
            ggtitle(title) +
            theme_classic() +
            theme(aspect.ratio = 1)
    }

    p_raw <- plot_df(log2(df_raw + 1e-5), "ATAC signal (log2)")
    p_norm <- plot_df(log2(df_norm + 1e-5), "Region normalized")
    p_norm_const <- plot_df(log2(df_norm_const + 1e-5), "Const normalized")
    p_norm_prob <- plot_df(df_norm_prob, "Probability")

    p_all <- patchwork::wrap_plots(p_raw, p_norm, p_norm_const, p_norm_prob, nrow = 1, guides = "collect")

    return(p_all)
}

plot_cell_type_scatter <- function(mat, anchor_cell_type, const_peaks = NULL, filename = NULL, width = NULL, height = NULL, ylab = "ATAC signal (log2)", prob1_thresh = NULL) {
    p <- mat %>%
        misha.ext::mat_to_intervs() %>%
        mutate(const = ifelse(const_peaks, "const", "variable")) %>%
        gather("type", "val", -(chrom:end), -!!sym(anchor_cell_type), -const) %>%
        ggplot(aes(x = !!sym(anchor_cell_type), y = val, color = const)) +
        scattermore::geom_scattermore(pointsize = 2) +
        geom_abline(linetype = "dashed", color = "darkblue") +
        scale_color_manual(name = "Peak Type", values = c("variable" = "black", "const" = "red")) +
        ylab(ylab) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        facet_wrap(~type)

    if (!is.null(prob1_thresh)) {
        p <- p + geom_vline(xintercept = log2(prob1_thresh), linetype = "dotted", color = "red")
        p <- p + geom_hline(yintercept = log2(prob1_thresh), linetype = "dotted", color = "red")
    }

    if (!is.null(filename)) {
        cli::cli_alert_info("Saving plot to {.val {filename}}")
        width <- width %||% length(unique(p$data$type))
        height <- height %||% length(unique(p$data$type))
        png(filename, width = width, height = height, units = "in", res = 300)
        print(p)
        dev.off()
    }

    return(p)
}

#' Normalize regional ATAC data using punctured windows
#'
#' @description
#' Performs regional normalization of ATAC-seq data by calculating marginal coverage in windows around peaks.
#' The function uses a "punctured window" approach, where the signal in a window around each peak
#' (excluding the peak itself) is used to normalize the peak signal. This helps account for local
#' background variation and accessibility biases that are generally due to GC content or mappability.
#'
#' @param peaks A data frame containing peak intervals with columns 'chrom', 'start', and 'end'
#' @param mat Numeric matrix. The raw ATAC signal matrix to normalize, with peaks as rows and cell types as columns
#' @param marginal_track Character. Name of the track containing marginal (summed) signal across all cell types
#' @param window_size Numeric. Size of the window around each peak to use for normalization (in base pairs). Default: 2e4
#' @param minimal_quantile Numeric. Minimum quantile of the punctured window coverage to use,
#'   prevents division by very small values. Default: 0.1
#'
#' @details
#' The normalization process follows these steps:
#' 1. Creates virtual tracks for peak signal and window signal
#' 2. Calculates punctured window coverage (window minus peak)
#' 3. Normalizes by the punctured coverage while enforcing a minimum based on minimal_quantile
#' 4. Performs final column-wise normalization
#'
#'
#' @return A normalized matrix with the same dimensions as the input, where each value
#'   has been adjusted for local background signal
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' norm_mat <- normalize_regional(peaks_df, raw_mat, "marginal_track")
#'
#' # With custom window size and quantile
#' norm_mat <- normalize_regional(
#'     peaks_df,
#'     raw_mat,
#'     "marginal_track",
#'     window_size = 2e4,
#'     minimal_quantile = 0.05
#' )
#' }
#'
#' @seealso
#' \code{\link{normalize_const}}, \code{\link{normalize_to_prob}} for subsequent normalization steps
#'
#' @export
normalize_regional <- function(peaks, mat, marginal_track, window_size = 2e4, minimal_quantile = 0.1) {
    gvtrack.create("marginal", marginal_track, func = "sum")
    gvtrack.create("marginal_20k", marginal_track, func = "sum")
    gvtrack.iterator("marginal_20k", sshift = -window_size / 2, eshift = window_size / 2)

    peaks_metadata <- misha.ext::gextract.left_join(
        c("marginal", "marginal_20k"),
        intervals = peaks,
        iterator = peaks
    ) %>%
        mutate(marginal_20k_punc = marginal_20k - marginal) %>%
        mutate(norm_f = marginal_20k_punc / quantile(marginal_20k_punc, minimal_quantile)) %>%
        mutate(norm_f = ifelse(norm_f < 1, 1, norm_f)) %>%
        mutate(norm_f = 1 / norm_f)

    mat_s <- mat * peaks_metadata$norm_f
    egc_s <- t(t(mat_s) / colSums(mat_s))
    return(egc_s)
}

#' Normalize ATAC data using constitutive peaks
#'
#' @description
#' Normalizes ATAC-seq data using constitutive peaks (peaks that are consistently active
#' across cell types) as reference points.
#'
#' @param peaks A data frame with peak intervals and a logical 'const' column indicating constitutive peaks
#' @param mat Numeric matrix. The matrix to normalize (typically output from \code{normalize_regional})
#' @param anchor_cell_type Character or NULL. The cell type to use as reference for normalization.
#'   If NULL, uses the mean of all cell types as reference. Default: NULL
#' @param norm_quant Numeric. Quantile used for primary normalization step. Default: 1
#' @param scaling_quant Numeric. Quantile used for final scaling step. Default: 0.9
#'
#' @details
#' The normalization process:
#' 1. Calculates total signal in constitutive peaks for each cell type
#' 2. Normalizes using either a specific cell type or mean of all cell types as reference
#' 3. Applies quantile-based normalization to standardize signal ranges
#' 4. Caps values at 1 and applies final scaling
#'
#' This approach assumes that constitutive peaks are highly accesible in all cell types and can be used as a reference.
#'
#' @return A normalized matrix with the same dimensions as the input, where cell type-specific
#'   biases have been corrected using constitutive peaks as reference
#'
#'
#' @examples
#' \dontrun{
#' # Using mean of all cell types as reference
#' norm_const_mat <- normalize_const(peaks_df, norm_mat)
#'
#' # Using specific anchor cell type
#' norm_const_mat <- normalize_const(
#'     peaks_df,
#'     norm_mat,
#'     anchor_cell_type = "CD4_T",
#'     norm_quant = 0.95,
#'     scaling_quant = 0.9
#' )
#' }
#'
#' @seealso
#' \code{\link{normalize_regional}} for prior normalization step,
#' \code{\link{normalize_to_prob}} for subsequent probability normalization
#'
#' @export
normalize_const <- function(peaks, mat, anchor_cell_type = NULL, norm_quant = 1, scaling_quant = 0.9) {
    if (!has_name(peaks, "const")) {
        cli_abort("No field named {.field const} in peaks.")
    }

    const_peaks <- peaks$const

    const_peak_score_all <- colSums(mat[const_peaks, ])
    if (is.null(anchor_cell_type)) {
        cli::cli_alert("Using the mean of all cell types as the anchor cell type")
        const_peak_score_norm_type <- mean(const_peak_score_all, na.rm = TRUE)
    } else {
        const_peak_score_norm_type <- const_peak_score_all[anchor_cell_type]
    }

    const_peak_score_norm <- const_peak_score_all / const_peak_score_norm_type

    egc_norm <- t(t(mat) / const_peak_score_norm)
    egc_norm <- egc_norm / quantile(as.vector(egc_norm), norm_quant)
    egc_norm[egc_norm > 1] <- 1
    if (is.null(anchor_cell_type)) {
        egc_norm <- egc_norm * max(mat) / scaling_quant
    } else {
        egc_norm <- egc_norm * max(mat[, anchor_cell_type]) / scaling_quant
    }

    return(egc_norm)
}

#' Convert normalized ATAC data to probability-like values
#'
#' @description
#' Transforms normalized ATAC-seq signals into probability-like values between 0 and 1,
#' using either a specified threshold or a quantile of constitutive peak values to
#' determine the maximum signal level.
#'
#' @param peaks A data frame with peak intervals and a logical 'const' column for constitutive peaks
#' @param mat Numeric matrix. The matrix to normalize (typically output from \code{normalize_const})
#' @param prob1_thresh Numeric or NULL. Signal threshold for probability = 1.
#'   If NULL, calculated from constitutive peaks using const_quantile. Default: NULL
#' @param const_quantile Numeric. Quantile of constitutive peak values to use as threshold
#'   when prob1_thresh is NULL. Default: 0.9
#'
#' @details
#' The normalization process:
#' 1. Determines threshold for probability = 1 (either directly or from constitutive peaks)
#' 2. Divides all values by this threshold
#' 3. Caps values at 1
#' 4. Applies min-max normalization to ensure full [0,1] range utilization
#'
#'
#' @return A matrix with the same dimensions as the input, containing probability-like
#'   values between 0 and 1
#'
#' @section Warning:
#' If prob1_thresh is NULL, the peaks data frame must include a logical 'const' column.
#' The function will stop with an error if this column is missing.
#'
#' @examples
#' \dontrun{
#' # Using constitutive peak quantile
#' prob_mat <- normalize_to_prob(peaks_df, norm_const_mat)
#'
#' # Using specific threshold
#' prob_mat <- normalize_to_prob(
#'     peaks_df,
#'     norm_const_mat,
#'     prob1_thresh = 0.75,
#'     const_quantile = 0.95
#' )
#' }
#'
#' @seealso
#' \code{\link{normalize_regional}}, \code{\link{normalize_const}} for prior normalization steps
#'
#' @export
normalize_to_prob <- function(peaks, mat, prob1_thresh = NULL, const_quantile = 0.9) {
    if (is.null(prob1_thresh)) {
        if (!has_name(peaks, "const")) {
            cli_abort("No field named {.field const} in peaks.")
        }
        cli::cli_alert("Using a quantile of {.val {const_quantile}} of the constitutive peaks as the p=1 threshold")
        prob1_thresh <- quantile(mat[peaks$const, ], const_quantile)
    }
    cli::cli_alert("Using {.val {log2(prob1_thresh)}} (in log scale) as the p=1 threshold")

    egc <- mat / prob1_thresh
    egc[egc > 1] <- 1
    egc <- apply(egc, 2, norm01)

    return(egc)
}
