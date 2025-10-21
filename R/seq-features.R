#' Create Sequence Features
#'
#' This function generates sequence features from given intervals. It calculates
#' dinucleotide frequencies and GC content, normalizes these features, and returns
#' them as a matrix.
#'
#' @param intervals A data frame containing interval information with columns `start` and `end`.
#' @param size The size of the sequences to extract. If NULL, the size is calculated from the first interval.
#' @param normalize A logical value indicating whether to normalize the features to the range 0-10.
#' @param norm_quant The quantile to use for normalization. Values below this quantile and above 1-quantile are truncated, and the rest are linearly scaled to 0-10.
#' @param norm_intervals Optional. A data frame of intervals to use for computing normalization quantiles. If NULL, quantiles are computed from the input intervals.
#' @return A matrix of normalized sequence features including GC content and dinucleotide frequencies.
#'
#' @export
create_sequence_features <- function(intervals, size = NULL, normalize = TRUE, norm_quant = 0.05, norm_intervals = NULL) {
    if (is.null(size)) {
        size <- intervals$end[1] - intervals$start[1]
    }
    seqs <- prego::intervals_to_seq(intervals, size)
    dinucs <- prego::calc_sequences_dinucs(seqs)
    dinucs <- dinucs / (size - 1)
    gc_content <- stringr::str_count(seqs, "G|C") / nchar(seqs)
    seq_feats <- cbind(gc_content, dinucs)
    if (normalize) {
        if (!is.null(norm_intervals)) {
            # Compute normalization parameters from reference intervals
            norm_seqs <- prego::intervals_to_seq(norm_intervals, size)
            norm_dinucs <- prego::calc_sequences_dinucs(norm_seqs) / (size - 1)
            norm_gc_content <- stringr::str_count(norm_seqs, "G|C") / nchar(norm_seqs)
            norm_feats <- cbind(norm_gc_content, norm_dinucs)

            # Apply normalization using reference quantiles
            for (i in 1:ncol(seq_feats)) {
                q_low <- quantile(norm_feats[, i], norm_quant, na.rm = TRUE)
                q_high <- quantile(norm_feats[, i], 1 - norm_quant, na.rm = TRUE)
                vals <- pmin(pmax(seq_feats[, i], q_low), q_high)
                rng <- q_high - q_low
                if (!is.finite(rng) || rng <= 0) {
                    seq_feats[, i] <- 0
                } else {
                    seq_feats[, i] <- ((vals - q_low) / rng) * 10
                }
            }
        } else {
            # Compute quantiles from the data itself (original behavior)
            seq_feats <- apply(seq_feats, 2, function(x) normalize_feature(x, norm_quant))
        }
    }
    return(seq_feats)
}
