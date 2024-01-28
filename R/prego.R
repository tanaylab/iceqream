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
#' @param peaks_size size of the peaks to extract sequences from. Default: 300bp
#' @param additional_features A matrix of additional features to filter out before learning the motifs (e.g. CpG content, dinucleotide content, etc.)
#' @param ... Additional arguments to be passed to \code{prego::regress_pwm}
#'
#' @export
learn_traj_prego <- function(peak_intervals, atac_diff, n_motifs, min_diff = 0.2, energy_norm_quantile = 1, min_energy = -10, sample_fraction = 0.1, sequences = NULL, seed = NULL, peaks_size = 300, additional_features = NULL, ...) {
    withr::local_options(list(gmax.data.size = 1e9))
    if (length(atac_diff) != nrow(peak_intervals)) {
        cli_abort("Length of {.field {atac_diff}} must be equal to the number of rows of {.field {peak_intervals}}. Current lengths: {.val {length(atac_diff)}} and {.val {nrow(peak_intervals)}}")
    }

    if (is.null(sequences)) {
        cli::cli_alert_info("Extracting sequences...")
        sequences <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(peak_intervals, peaks_size)))
    }

    peaks_df <- peak_intervals %>%
        select(chrom, start, end, const) %>%
        mutate(id = seq_len(n()))

    if (!is.null(additional_features)) {
        if (nrow(additional_features) != nrow(peaks_df)) {
            cli_abort("Number of rows of {.field additional_features} must be equal to the number of rows of {.field peak_intervals}. Current lengths: {.val {nrow(additional_features)}} and {.val {nrow(peak_intervals)}}")
        }

        rownames(additional_features) <- peaks_df$id
    }

    peaks_df <- peaks_df %>%
        mutate(score = atac_diff) %>%
        filter(abs(score) >= min_diff)

    if (has_name(peak_intervals, "const")) {
        peaks_df <- peaks_df %>%
            filter(!const) %>%
            select(-const)
    }

    if (nrow(peaks_df) < 200) {
        cli_warn("Not enough peaks to run prego. Please consider increasing {.field {min_diff}} (current value: {.val {min_diff}})")
        return(NULL)
    }

    if (!is.null(additional_features)) {
        additional_features <- additional_features[peaks_df$id, ]
        additional_features[is.na(additional_features)] <- 0
        cli_alert_info("Learning a model using only additional features...")
        glm_feats <- glmnet::glmnet(as.matrix(additional_features), norm01(peaks_df$score), binomial(link = "logit"), alpha = 1, lambda = 1e-5, parallel = TRUE, seed = seed)
        pred <- logist(glmnet::predict.glmnet(glm_feats, newx = as.matrix(additional_features), type = "link", s = 1e-5))[, 1]
        score <- norm01(peaks_df$score) - pred
    } else {
        score <- peaks_df$score
    }

    seqs <- toupper(misha::gseq.extract(misha.ext::gintervals.normalize(peaks_df, peaks_size)))

    cli_alert_info("Inferring {.val {n_motifs}} prego motifs...")
    reg <- prego::regress_pwm(seqs, score, motif_num = n_motifs, multi_kmers = TRUE, internal_num_folds = 1, screen_db = FALSE, match_with_db = FALSE, seed = seed, sample_for_kmers = TRUE, sample_frac = sample_fraction, ...)

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
