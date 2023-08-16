#' Generate a set of non-overlapping peaks from a set of intervals
#'
#' @description Given a set of intervals, usually generated from a threshold on the marginal ATAC signal, the function generates a set of non-overlapping peaks of a given size. The function can also take a marginal track and use it to generate peaks with the highest marginal signal. \cr
#' The function is iterative, where in each step the intervals are centered on the highest marginal signal, expanded to the desired size, and then merged to remove overlaps. The process is repeated until no more overlaps are present. \cr
#'
#' @param intervals A data frame with columns \code{chrom}, \code{start}, and \code{end}.
#' @param size The desired size of the peaks.
#' @param marginal_track A marginal track to use for generating peaks (or any other misha track expression).
#' @param iterator An iterator to use for extracting the marginal track.
#' @param keep_marginal Whether to keep the marginal track in the output. If TRUE, a column named \code{marginal} will be added to the output.
#'
#' @return A data frame with columns \code{chrom}, \code{start}, \code{end}, and optionally \code{marginal}.
#'
#' @examples
#' \dontrun{
#' intervals <- gscreen("atac_marginal >= 300")
#' peaks <- canonize_peaks(intervals, 300, marginal_track = "atac_marginal")
#' head(peaks)
#' }
#'
#' @export
canonize_peaks <- function(intervals, size, marginal_track = NULL, iterator = NULL, keep_marginal = TRUE) {
    intervals <- intervals %>%
        distinct(chrom, start, end)

    normalize <- function(intervals) {
        if (!is.null(marginal_track)) {
            intervals <- misha::gextract(marginal_track, intervals = intervals, iterator = iterator, colnames = "marginal") %>%
                arrange(desc(marginal)) %>%
                distinct(intervalID, .keep_all = TRUE) %>%
                select(-intervalID) %>%
                misha.ext::gintervals.normalize(size)
        } else {
            intervals <- intervals %>%
                misha.ext::gintervals.normalize(size)
        }
        return(intervals)
    }

    canonize <- function(intervals) {
        intervals <- intervals %>%
            misha::gintervals.force_range() %>%
            misha::gintervals.canonic()
        return(intervals)
    }

    intervals <- normalize(intervals)
    cli_alert("Number of peaks: {.val {nrow(intervals)}}")

    while (nrow(intervals) != nrow(canonize(intervals))) {
        intervals <- canonize(intervals)
        intervals <- normalize(intervals)
        cli_alert("Number of peaks: {.val {nrow(intervals)}}")
    }

    intervals <- intervals %>%
        arrange(chrom, start, end)

    if (!keep_marginal) {
        intervals <- intervals %>%
            select(-marginal)
    }

    cli_alert_success("Generated {.val {nrow(intervals)}} non-overlapping peaks of size {.val {size}}")


    return(intervals)
}
