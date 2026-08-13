# Internal replacements for the deprecated misha.ext::intervs_to_mat() /
# misha.ext::mat_to_intervs(). Those now forward to misha::gintervals.to_mat() /
# misha::gintervals.from_mat(), which label rows "<chrom>:<start>-<end>" and
# return a classed `intervs_mat` - neither is what iceqream expects.

intervs_to_mat <- function(intervs) {
    if (is.null(intervs) || nrow(intervs) == 0) {
        return(matrix(nrow = 0, ncol = 0))
    }

    mat <- misha::gintervals.to_mat(intervs)

    # ponytail: return a plain matrix keyed by "<chrom>_<start>_<end>".
    # iceqream matches rownames against `peaks$peak_name` (see preprocess_data)
    # and treats these as plain matrices - `[.intervs_mat` does not drop single
    # columns and carries an `intervals` attribute that t() silently invalidates.
    class(mat) <- NULL
    attr(mat, "intervals") <- NULL
    rownames(mat) <- paste0(intervs$chrom, "_", intervs$start, "_", intervs$end)
    return(mat)
}

mat_to_intervs <- function(mat) {
    rn <- rownames(mat)
    if (is.null(rn)) {
        stop("mat_to_intervs: mat has no rownames")
    }
    parts <- strsplit(rn, "_")
    intervals <- data.frame(
        chrom = vapply(parts, function(x) paste(x[1:(length(x) - 2)], collapse = "_"), character(1)),
        start = as.numeric(vapply(parts, function(x) x[length(x) - 1], character(1))),
        end = as.numeric(vapply(parts, function(x) x[length(x)], character(1))),
        stringsAsFactors = FALSE
    )
    res <- cbind(intervals, as.data.frame(mat))
    rownames(res) <- NULL
    return(res)
}
