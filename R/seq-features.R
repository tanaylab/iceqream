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
#' @return A matrix of normalized sequence features including GC content and dinucleotide frequencies.
#'
#' @export
create_sequence_features <- function(intervals, size = NULL, normalize = TRUE, norm_quant = 0.05) {
    if (is.null(size)) {
        size <- intervals$end[1] - intervals$start[1]
    }
    seqs <- prego::intervals_to_seq(misha.ext::gintervals.normalize(intervals, size))
    dinucs <- prego::calc_sequences_dinucs(seqs)
    dinucs <- dinucs / (size - 1)
    gc_content <- stringr::str_count(seqs, "G|C") / nchar(seqs)
    seq_feats <- cbind(gc_content, dinucs)
    if (normalize) {
        seq_feats <- apply(seq_feats, 2, function(x) normalize_feature(x, norm_quant))
    }
    return(seq_feats)
}
