#' Compute Motif Energies
#'
#' This function computes motif energies for given peak intervals using a specified motif database.
#' Optionally, it can normalize the computed motif energies using additional normalization intervals.
#'
#' @param peak_intervals A data frame containing the peak intervals with columns: chrom, start, and end. If an additional column 'peak_name' is present, it will be used as row names in the output matrix.
#' @param db A motif database to use for extracting motif energies. Default is `iceqream::motif_db`.
#' @param normalization_intervals A data frame containing intervals for normalization. Default is `peak_intervals`.
#' @param normalize A logical value indicating whether to normalize the motif energies. Default is TRUE.
#' @param energy_norm_quantile A numeric value for the quantile used in normalization. Default is 1.
#' @param db_quantiles A matrix of motif energies quantiles to use for normalization. Rows are motifs and columns are quantiles. If not NULL, this would be used instead of computing the quantiles from the normalization intervals. A precomputed matrix for mouse gastrulation can be found in `iceqream::mouse_db_quantiles`.
#'
#' @return A matrix of motif energies.
#'
#' @inheritParams norm_energy_matrix
#' @inheritParams prego::gextract_pwm
#'
#'
#' @export
compute_motif_energies <- function(peak_intervals, db = iceqream::motif_db, normalization_intervals = peak_intervals, prior = 0.01, normalize = TRUE, energy_norm_quantile = 1, norm_energy_max = 10, min_energy = -7, db_quantiles = NULL) {
    n_motifs <- length(unique(db$motif))
    n_peaks <- nrow(peak_intervals)
    cli_alert("Computing motif energies for {.val {n_peaks}} intervals using {.val {n_motifs}} motifs. This might take a while. You can set the number of cores to use with {.code prego::set_parallel()}")
    motif_energies <- prego::gextract_pwm(peak_intervals %>% select(chrom, start, end), dataset = db, prior = prior) %>%
        select(-chrom, -start, -end) %>%
        as.matrix()

    if (normalize) {
        if (!is.null(db_quantiles)) {
            cli_alert_info("Using pre-computed quantiles for normalization")
            motif_energies <- normalize_with_db_quantiles(
                motif_energies,
                db_quantiles,
                energy_norm_quantile,
                min_energy,
                norm_energy_max
            )
        } else {
            missing_intervals <- normalization_intervals %>% anti_join(peak_intervals, by = c("chrom", "start", "end"))
            if (nrow(missing_intervals) > 0) {
                cli::cli_alert_info("Extracting {.val {nrow(missing_intervals)}} intervals to normalize the motif energies")
                missing_energies <- prego::gextract_pwm(missing_intervals %>% select(chrom, start, end), dataset = db, prior = prior) %>%
                    select(-chrom, -start, -end) %>%
                    as.matrix()
                norm_energies <- rbind(motif_energies, missing_energies)
            } else {
                norm_energies <- motif_energies
            }

            cli_alert_info("Normalizing {.val {nrow(motif_energies)}} motif energies using {.val {nrow(norm_energies)}} intervals")
            motif_energies <- norm_energy_matrix(motif_energies, norm_energies, min_energy = min_energy, q = energy_norm_quantile, norm_energy_max = norm_energy_max)
        }
    }

    if ("peak_name" %in% colnames(peak_intervals)) {
        rownames(motif_energies) <- peak_intervals$peak_name
    }

    return(motif_energies)
}

validate_motif_energies <- function(motif_energies, peak_intervals, pssm_db = iceqream::motif_db) {
    if (nrow(motif_energies) != nrow(peak_intervals)) {
        cli_abort("Motif energies must have the same number of rows as peak_intervals. peak_intervals number of rows {.val {nrow(peak_intervals)}} != motif_energies number of rows {.val {nrow(motif_energies)}}", call = parent.frame(1))
    }

    missing_motifs <- colnames(motif_energies)[!(colnames(motif_energies) %in% pssm_db$motif)]
    if (length(missing_motifs) > 0) {
        cli_abort("Motifs {.val {missing_motifs}} are missing from the motif database ({.field pssm_db} parameter).", call = parent.frame(1))
    }
}
