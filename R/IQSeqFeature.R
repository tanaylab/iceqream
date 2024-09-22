#' IQSeqFeature class
#'
#' This class represents a sequence feature for IQ models, inheriting from IQFeature.
#'
#' @slot name The name of the IQ sequence feature (character).
#' @slot coefs The coefficients of the IQ sequence feature (numeric).
#' @slot compute_func The function to compute the feature (function).
#' @slot min_value The minimum value for normalization (numeric).
#' @slot max_value The maximum value for normalization (numeric).
#' @slot quantile The quantile to use for normalization (numeric).
#' @slot size The size of the sequences to use for feature computation (numeric).
#'
#' @exportClass IQSeqFeature
IQSeqFeature <- setClass(
    "IQSeqFeature",
    slots = list(
        compute_func = "function",
        min_value = "numeric",
        max_value = "numeric",
        quantile = "numeric",
        size = "numeric"
    ),
    contains = "IQFeature"
)

#' Initialize an IQSeqFeature object
#'
#' @param .Object The IQSeqFeature object to initialize.
#' @param name The name of the feature.
#' @param coefs The coefficients of the feature.
#' @param compute_func The function to compute the feature.
#' @param min_value The minimum value for normalization.
#' @param max_value The maximum value for normalization.
#' @param quantile The quantile to use for normalization.
#' @param size The size of the sequences to use for feature computation.
#'
#' @return An initialized IQSeqFeature object.
#'
#' @exportMethod initialize
setMethod(
    "initialize", "IQSeqFeature",
    function(.Object, name, coefs, compute_func, min_value, max_value, quantile = 0.99, size) {
        .Object@name <- name
        .Object@coefs <- coefs
        .Object@compute_func <- compute_func
        .Object@min_value <- min_value
        .Object@max_value <- max_value
        .Object@quantile <- quantile
        .Object@size <- size
        .Object
    }
)

#' Show method for IQSeqFeature
#'
#' This method defines how an IQSeqFeature object should be displayed.
#'
#' @param object An IQSeqFeature object
#'
#' @rdname IQSeqFeature-class
#' @exportMethod show
setMethod("show", signature = "IQSeqFeature", definition = function(object) {
    cli::cli({
        cli::cli_text("An {.cls IQSeqFeature} object named {.val {object@name}}")
        if (length(object@coefs) > 0) {
            cli::cli_text("Contains {.val {length(object@coefs)}} coefficients: {.val {names(object@coefs)}}")
        } else {
            cli::cli_text("Contains no coefficients")
        }
        cli::cli_text("Minimum value: {.val {object@min_value}}")
        cli::cli_text("Maximum value: {.val {object@max_value}}")
        cli::cli_text("Normalization quantile: {.val {object@quantile}}")
        cli::cli_text("Sequence size: {.val {object@size}}")
    })
})

#' Compute the IQSeqFeature
#'
#' This function computes the IQ sequence feature for a given set of sequences.
#'
#' @param iq An IQSeqFeature object.
#' @param sequences A vector of DNA sequences.
#'
#' @return A vector of computed and normalized IQ sequence feature values.
#'
#' @export
iq_seq_feature.compute <- function(iq, sequences) {
    # Ensure sequences are of the correct size
    sequences <- ensure_sequence_size(sequences, iq@size)

    raw_values <- iq@compute_func(sequences)
    iq_seq_feature.normalize(iq, raw_values)
}

#' @noRd
ensure_sequence_size <- function(sequences, size) {
    seq_lengths <- nchar(sequences)
    if (any(seq_lengths != size)) {
        cli::cli_warn("Some sequences are not of the correct size {.val {size}}. ")
    }
    return(sequences)
}

iq_seq_feature.normalize <- function(iq, values) {
    vals <- values - iq@min_value
    vals <- vals / iq@max_value
    vals[vals > 1] <- 1
    return(vals * 10) # Scaling to 0-10 range
}

#' Create an IQSeqFeature from a trajectory model
#'
#' @param traj_model A trajectory model object.
#' @param feature_name The name of the feature to create.
#' @param compute_func The function to compute the feature.
#' @param quantile The quantile to use for normalization (default is 0.99).
#'
#' @return An IQSeqFeature object.
#'
#' @export
create_iq_seq_feature <- function(traj_model, feature_name, compute_func, quantile = 0.99) {
    # Get the size from the trajectory model
    size <- traj_model@params$peaks_size

    # Compute feature values on normalization intervals
    norm_sequences <- prego::intervals_to_seq(misha.ext::gintervals.normalize(traj_model@normalization_intervals, size))
    feature_values <- compute_func(norm_sequences)

    # Calculate min and max values
    min_value <- min(feature_values, na.rm = TRUE)
    max_value <- quantile(feature_values, quantile, na.rm = TRUE)

    # Get coefficients
    f2v <- feat_to_variable(traj_model)
    variables <- f2v %>%
        filter(variable == feature_name) %>%
        pull(feature)
    coefs <- coef(traj_model@model, s = traj_model@params$lambda)[variables, , drop = TRUE]
    names(coefs) <- gsub(paste0("^", feature_name, "_"), "", names(coefs))

    # Create and return the IQSeqFeature object
    IQSeqFeature(
        name = feature_name,
        coefs = coefs,
        compute_func = compute_func,
        min_value = min_value,
        max_value = max_value,
        quantile = quantile,
        size = size
    )
}

#' Create Dinucleotide Features
#'
#' This function creates IQSeqFeature objects for all possible dinucleotides.
#'
#' @param traj_model A trajectory model object.
#' @param quantile The quantile to use for normalization (default is 0.99).
#' @param size The size of the sequences to use for feature computation.
#' @param dinucleotides A vector of dinucleotides to create features for (default is all possible dinucleotides).
#' @param feat_names A vector of feature names to use for each dinucleotide (default is the dinucleotide itself).
#'
#' @return A list of IQSeqFeature objects, one for each dinucleotide.
#'
#' @export
create_dinuc_features <- function(traj_model, quantile = 0.99, size = traj_model@params$peaks_size, dinucleotides = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"), feat_names = dinucleotides) {
    if (length(feat_names) != length(dinucleotides)) {
        cli::cli_abort("The number of feature names should be equal to the number of dinucleotides.")
    }
    dinuc_to_name <- setNames(feat_names, dinucleotides)
    # Create a function to compute dinucleotide distributions
    compute_dinuc_dist <- function(sequences) {
        prego::calc_sequences_dinucs(sequences)
    }

    # Compute dinucleotide distribution for normalization intervals
    norm_sequences <- prego::intervals_to_seq(gintervals.normalize(traj_model@normalization_intervals, size))
    dinuc_dist <- compute_dinuc_dist(norm_sequences)

    # Create a list to store IQSeqFeature objects
    dinuc_features <- list()

    for (dinuc in dinucleotides) {
        # Create a function to compute the specific dinucleotide frequency
        compute_func <- function(sequences) {
            dist <- compute_dinuc_dist(sequences)
            return(dist[, dinuc, drop = TRUE])
        }

        # Calculate min and max values for this dinucleotide
        feature_values <- dinuc_dist[, dinuc, drop = TRUE]
        min_value <- min(feature_values, na.rm = TRUE)
        max_value <- quantile(feature_values, quantile, na.rm = TRUE)

        # Get coefficients for this dinucleotide
        feature_name <- dinuc_to_name[[dinuc]]
        f2v <- feat_to_variable(traj_model)
        variables <- f2v %>%
            filter(variable == feature_name) %>%
            pull(feature)
        coefs <- coef(traj_model@model, s = traj_model@params$lambda)[variables, , drop = TRUE]
        names(coefs) <- gsub(paste0("^", feature_name, "_"), "", names(coefs))

        # Create and store the IQSeqFeature object
        dinuc_features[[feature_name]] <- IQSeqFeature(
            name = feature_name,
            coefs = coefs,
            compute_func = compute_func,
            min_value = min_value,
            max_value = max_value,
            quantile = quantile,
            size = size
        )
    }

    return(dinuc_features)
}

#' Create Dinucleotide Feature Group
#'
#' This function creates an IQFeatureGroup object for all dinucleotide features.
#'
#' @param traj_model A trajectory model object.
#' @param quantile The quantile to use for normalization (default is 0.99).
#' @param size The size of the sequences to use for feature computation.
#' @param dinucleotides A vector of dinucleotides to create features for (default is all possible dinucleotides).
#'
#' @return An IQFeatureGroup object containing all dinucleotide features.
#'
#' @export
create_dinuc_feature_group <- function(traj_model, quantile = 0.99, size = traj_model@params$peaks_size, dinucleotides = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")) {
    # Create individual dinucleotide features
    dinuc_features <- create_dinuc_features(traj_model, quantile, size, dinucleotides)

    # Create the compute function for all dinucleotides at once
    compute_func <- function(sequences) {
        prego::calc_sequences_dinucs(sequences)
    }

    # Create and return the IQFeatureGroup object
    IQFeatureGroup(
        features = dinuc_features,
        compute_func = compute_func,
        size = size
    )
}
