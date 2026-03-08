# Internal compatibility functions to remove misha.ext dependency

intervs_to_mat <- function(intervs) {
    if (is.null(intervs) || nrow(intervs) == 0) {
        return(matrix(nrow = 0, ncol = 0))
    }

    # Usually gextract returns chrom, start, end, intervalID as first 4 columns
    # We drop them and convert the rest to a matrix
    mat <- as.matrix(intervs[, -(1:4), drop = FALSE])
    return(mat)
}

mat_to_intervs <- function(mat, intervals = NULL) {
    # Takes a matrix and an interval data frame (usually first 3 or 4 columns)
    # and combines them.
    # If intervals is NULL, reconstruct from rownames (expected format: chrom_start_end)
    if (is.null(intervals)) {
        rn <- rownames(mat)
        if (is.null(rn)) {
            stop("mat_to_intervs: intervals not provided and mat has no rownames")
        }
        parts <- strsplit(rn, "_")
        intervals <- data.frame(
            chrom = vapply(parts, function(x) paste(x[1:(length(x) - 2)], collapse = "_"), character(1)),
            start = as.numeric(vapply(parts, function(x) x[length(x) - 1], character(1))),
            end = as.numeric(vapply(parts, function(x) x[length(x)], character(1))),
            stringsAsFactors = FALSE
        )
    }
    # If intervals has more than 4 columns, we usually only want the first 3 or 4.
    cols <- min(ncol(intervals), 4)
    res <- cbind(intervals[, 1:cols, drop = FALSE], as.data.frame(mat))
    return(res)
}

gintervals.expand <- function(inv, expansion = 100) {
    if (is.character(inv)) {
        inv <- misha::gintervals.load(inv)
    }
    inv$start <- inv$start - expansion
    inv$end <- inv$end + expansion
    inv <- as.data.frame(inv)
    inv <- misha::gintervals.force_range(inv)
    return(inv)
}

gintervals.centers <- function(intervals) {
    if (is.null(intervals) || nrow(intervals) == 0) {
        return(intervals)
    }

    mid <- floor((intervals$start + intervals$end) / 2)
    intervals$start <- mid
    intervals$end <- mid + 1

    return(intervals)
}

gextract.left_join <- function(expr, colnames = NULL, intervals = NULL, iterator = NULL, ...) {
    # Implements gextract + left_join logic using base misha functions.
    # Extracts track expressions with gextract, then merges back onto
    # the original intervals so that all input rows are preserved (left join).
    extracted <- misha::gextract(expr, intervals = intervals, iterator = iterator, colnames = colnames)
    # Merge back onto original intervals to preserve all rows (left join)
    res <- merge(intervals, extracted, by = c("chrom", "start", "end"), all.x = TRUE)
    # Restore original row order
    res <- res[order(match(paste(res$chrom, res$start, res$end),
                           paste(intervals$chrom, intervals$start, intervals$end))), ]
    rownames(res) <- NULL
    return(res)
}
