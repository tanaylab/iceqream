#' IQFeatureGroup class
#'
#' This class represents a group of related IQ features, such as all dinucleotide features.
#'
#' @slot features A list of IQSeqFeature objects.
#' @slot compute_func A function to compute all features in the group at once.
#' @slot size The size of the sequences to use for feature computation (numeric).
#'
#' @export
IQFeatureGroup <- setClass(
    "IQFeatureGroup",
    slots = list(
        features = "list",
        compute_func = "function",
        size = "numeric"
    )
)

#' Initialize an IQFeatureGroup object
#'
#' @param .Object The IQFeatureGroup object to initialize.
#' @param features A list of IQSeqFeature objects.
#' @param compute_func A function to compute all features in the group at once.
#' @param size The size of the sequences to use for feature computation.
#'
#' @return An initialized IQFeatureGroup object.
#'
#' @export
setMethod(
    "initialize", "IQFeatureGroup",
    function(.Object, features, compute_func, size) {
        .Object@features <- features
        .Object@compute_func <- compute_func
        .Object@size <- size
        .Object
    }
)

#' Show method for IQFeatureGroup
#'
#' This method defines how an IQFeatureGroup object should be displayed.
#'
#' @param object An IQFeatureGroup object
#'
#' @export
setMethod("show", signature = "IQFeatureGroup", definition = function(object) {
    cli::cli({
        cli::cli_text("An {.cls IQFeatureGroup} object with {.val {length(object@features)}} features")
        cli::cli_text("Feature names: {.val {names(object@features)}}")
        cli::cli_text("Sequence size: {.val {object@size}}")
    })
})

#' Compute features for an IQFeatureGroup
#'
#' This function computes all features in the group for a given set of sequences or intervals.
#'
#' @param group An IQFeatureGroup object.
#' @param sequences A vector of DNA sequences. Optional if intervals are provided.
#' @param intervals A data frame of genomic intervals with columns 'chrom', 'start', and 'end'. Optional if sequences are provided.
#'
#' @return A matrix of computed and normalized feature values.
#'
#' @export
iq_feature_group.compute <- function(group, sequences = NULL, intervals = NULL) {
    if (is.null(sequences) && is.null(intervals)) {
        cli::cli_abort("Either sequences or intervals must be provided.")
    }

    if (is.null(sequences)) {
        if (!all(c("chrom", "start", "end") %in% colnames(intervals))) {
            cli::cli_abort("Intervals must have 'chrom', 'start', and 'end' columns.")
        }
        cli::cli_alert_info("Extracting sequences from {.val {nrow(intervals)}} intervals...")
        sequences <- prego::intervals_to_seq(gintervals.normalize(intervals, group@size))
    }

    # Ensure all sequences are of the correct size
    sequences <- ensure_sequence_size(sequences, group@size)

    # Compute all features at once
    raw_values <- group@compute_func(sequences)

    # Normalize each feature
    for (feature_name in names(group@features)) {
        feature <- group@features[[feature_name]]
        raw_values[, feature_name] <- iq_seq_feature.normalize(feature, raw_values[, feature_name])
    }

    return(raw_values)
}
