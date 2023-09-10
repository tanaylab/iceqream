#' Learn 'prego' models for ATAC difference of a trajectory
#'
#' @param peak_intervals A data frame, indicating the genomic positions ('chrom', 'start', 'end') of each peak, with an additional column named "const" indicating whether the peak is constitutive. Optionally, a column named "cluster" can be added with indication of the cluster of each peak.
#' @param atac_diff A numeric vector, indicating the ATAC difference of each peak
#' @param n_motifs Number of motifs to learn
#' @param min_diff Minimum ATAC difference to include a peak in the training
#' @param energy_norm_quantile quantile of the energy used for normalization. Default: 1
#' @param min_energy Minimum energy value after normalization (default: -10)
#' @param sample_fraction Fraction of peaks to sample for training. Default: 0.1
#' @param sequences A character vector of sequences to learn the motifs on. If NULL, the sequences of the peaks are used.
#' @param seed Random seed
#'
#' @export
learn_traj_prego <- function(peak_intervals, atac_diff, n_motifs, min_diff = 0.2, energy_norm_quantile = 1, min_energy = -10, sample_fraction = 0.1, sequences = NULL, seed = NULL) {
    withr::local_options(list(gmax.data.size = 1e9))
    if (length(atac_diff) != nrow(peak_intervals)) {
        cli_abort("Length of {.field {atac_diff}} must be equal to the number of rows of {.field {peak_intervals}}. Current lengths: {.val {length(atac_diff)}} and {.val {nrow(peak_intervals)}}")
    }

    if (is.null(sequences)) {
        sequences <- toupper(misha::gseq.extract(peak_intervals))
    }

    peaks_df <- peak_intervals %>%
        select(chrom, start, end, const) %>%
        mutate(id = seq_len(n())) %>%
        mutate(score = atac_diff) %>%
        filter(abs(score) >= min_diff)

    if (has_name(peak_intervals, "const")) {
        peaks_df <- peaks_df %>%
            filter(!const) %>%
            select(-const)
    }

    if (nrow(peaks_df) < 200) {
        cli_warn("Not enough peaks to run prego. Please consider increasing {.field {min_diff}} (current value: {.val {min_prego_diff}})")
        return(NULL)
    }

    seqs <- toupper(misha::gseq.extract(peaks_df))

    cli_alert_info("Inferring {.val {n_motifs}} prego motifs...")
    reg <- prego::regress_pwm(seqs, peaks_df$score, motif_num = n_motifs, multi_kmers = TRUE, internal_num_folds = 1, screen_db = FALSE, match_with_db = FALSE, seed = seed, sample_for_kmers = TRUE, sample_frac = sample_fraction)

    prego_e <- reg$predict_multi(sequences)
    prego_e <- apply(prego_e, 2, norm_energy, min_energy = min_energy, q = energy_norm_quantile)

    prego_models <- prego::export_multi_regression(reg)$models
    names(prego_models) <- colnames(prego_e)

    prego_pssm <- purrr::imap_dfr(prego_models, ~ .x$pssm %>% mutate(motif = .y)) %>%
        select(motif, pos, everything())

    return(
        list(
            energies = prego_e,
            pssm = prego_pssm,
            models = prego_models
        )
    )
}
