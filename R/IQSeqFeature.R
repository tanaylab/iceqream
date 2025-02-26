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
    if (!is.null(traj_model@params$feats_peaks_size)) {
        size <- traj_model@params$feats_peaks_size
    } else {
        size <- traj_model@params$peaks_size
    }

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
#' This function creates IQSeqFeature objects for all possible dinucleotides and GC content.
#'
#' @param traj_model A trajectory model object.
#' @param quantile The quantile to use for normalization (default is 0.99).
#' @param size The size of the sequences to use for feature computation.
#' @param dinucleotides A vector of dinucleotides to create features for (default is all possible dinucleotides).
#' @param feat_names A vector of feature names to use for each dinucleotide (default is the dinucleotide itself).
#' @param include_gc Logical indicating whether to include GC content feature (default is TRUE).
#'
#' @return A list of IQSeqFeature objects, one for each dinucleotide and GC content if included.
#'
#' @export
create_dinuc_features <- function(traj_model, quantile = 0.99, size = NULL,
                                  dinucleotides = c(
                                      "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                                      "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"
                                  ),
                                  feat_names = dinucleotides, include_gc = TRUE) {
    if (is.null(size)) {
        if (!is.null(traj_model@params$feats_peaks_size)) {
            size <- traj_model@params$feats_peaks_size
        } else {
            size <- traj_model@params$peaks_size
        }
    }
    if (length(feat_names) != length(dinucleotides)) {
        cli::cli_abort("The number of feature names should be equal to the number of dinucleotides.")
    }
    dinuc_to_name <- setNames(feat_names, dinucleotides)
    # Create a function to compute dinucleotide distributions
    compute_dinuc_dist <- function(sequences) {
        prego::calc_sequences_dinucs(sequences)
    }

    # Compute dinucleotide distribution for normalization intervals
    norm_sequences <- prego::intervals_to_seq(traj_model@normalization_intervals, size)
    dinuc_dist <- compute_dinuc_dist(norm_sequences)

    # Compute GC content for normalization intervals if needed
    if (include_gc) {
        gc_content <- stringr::str_count(norm_sequences, "G|C") / nchar(norm_sequences)
    }

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

    # Add GC content feature if requested
    if (include_gc) {
        # Create a function to compute GC content
        compute_gc_func <- function(sequences) {
            stringr::str_count(sequences, "G|C") / nchar(sequences)
        }

        # Calculate min and max values for GC content
        min_gc <- min(gc_content, na.rm = TRUE)
        max_gc <- quantile(gc_content, quantile, na.rm = TRUE)

        # Get coefficients for GC content
        feature_name <- "gc_content"
        f2v <- feat_to_variable(traj_model)
        variables <- f2v %>%
            filter(variable == feature_name) %>%
            pull(feature)

        # If GC content coefficients exist in the model, use them
        if (length(variables) > 0) {
            coefs <- coef(traj_model@model, s = traj_model@params$lambda)[variables, , drop = TRUE]
            names(coefs) <- gsub(paste0("^", feature_name, "_"), "", names(coefs))
        } else {
            # Otherwise create an empty coefficient vector
            coefs <- numeric(0)
        }

        # Create and store the GC content IQSeqFeature
        dinuc_features[[feature_name]] <- IQSeqFeature(
            name = feature_name,
            coefs = coefs,
            compute_func = compute_gc_func,
            min_value = min_gc,
            max_value = max_gc,
            quantile = quantile,
            size = size
        )
    }

    return(dinuc_features)
}

#' Create Dinucleotide Feature Group
#'
#' This function creates an IQFeatureGroup object for all dinucleotide features and GC content.
#'
#' @param traj_model A trajectory model object.
#' @param quantile The quantile to use for normalization (default is 0.99).
#' @param size The size of the sequences to use for feature computation.
#' @param dinucleotides A vector of dinucleotides to create features for (default is all possible dinucleotides).
#' @param include_gc Logical indicating whether to include GC content feature (default is TRUE).
#'
#' @return An IQFeatureGroup object containing all dinucleotide features and GC content if included.
#'
#' @export
create_dinuc_feature_group <- function(traj_model, quantile = 0.99, size = NULL,
                                       dinucleotides = c(
                                           "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                                           "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"
                                       ),
                                       include_gc = TRUE) {
    if (is.null(size)) {
        if (!is.null(traj_model@params$feats_peaks_size)) {
            size <- traj_model@params$feats_peaks_size
        } else {
            size <- traj_model@params$peaks_size
        }
    }
    # Create individual dinucleotide features
    dinuc_features <- create_dinuc_features(traj_model, quantile, size, dinucleotides, include_gc = include_gc)

    # Create the compute function for all features at once
    compute_func <- function(sequences) {
        # Get dinucleotide frequencies
        dinuc_dist <- prego::calc_sequences_dinucs(sequences)

        # Add GC content if included
        if (include_gc && "gc_content" %in% names(dinuc_features)) {
            gc_content <- stringr::str_count(sequences, "G|C") / nchar(sequences)
            result <- cbind(dinuc_dist, gc_content)
            colnames(result)[ncol(result)] <- "gc_content"
            return(result)
        }

        return(dinuc_dist)
    }

    # Create and return the IQFeatureGroup object
    IQFeatureGroup(
        features = dinuc_features,
        compute_func = compute_func,
        size = size
    )
}

#' Create CG Content Feature
#'
#' This function creates an IQSeqFeature that calculates -log2(1-CG_fraction) as a measure
#' of CG content in DNA sequences. The feature is normalized to be between 0-10.
#'
#' @param traj_model A trajectory model object.
#' @param quantile The quantile to use for normalization (default is 0.99).
#'
#' @return An IQSeqFeature object for CG content.
#'
#' @export
create_cg_content_feature <- function(traj_model, quantile = 0.99) {
    # Get the size from the trajectory model
    if (!is.null(traj_model@params$feats_peaks_size)) {
        size <- traj_model@params$feats_peaks_size
    } else {
        size <- traj_model@params$peaks_size
    }

    # Create function to compute -log2(1-CG_fraction)
    compute_cg_content <- function(sequences) {
        # Extract all CG dinucleotides
        cg_counts <- stringr::str_count(sequences, "CG")

        # Calculate total possible positions for CG dinucleotides
        total_positions <- nchar(sequences) - 1

        # Calculate CG fraction
        cg_fraction <- cg_counts / total_positions

        # Apply -log2(1-CG_fraction) transformation
        # Handle edge cases where CG_fraction = 1 (would result in -Inf)
        result <- -log2(pmax(1 - cg_fraction, .Machine$double.eps))

        return(result)
    }

    # Compute feature values on normalization intervals
    norm_sequences <- prego::intervals_to_seq(misha.ext::gintervals.normalize(traj_model@normalization_intervals, size))
    feature_values <- compute_cg_content(norm_sequences)

    # Calculate min and max values
    min_value <- min(feature_values, na.rm = TRUE)
    max_value <- quantile(feature_values, quantile, na.rm = TRUE)

    # Get coefficients if they exist in the model
    feature_name <- "cg_cont"
    f2v <- feat_to_variable(traj_model)
    variables <- f2v %>%
        filter(variable == feature_name) %>%
        pull(feature)

    if (length(variables) > 0) {
        coefs <- coef(traj_model@model, s = traj_model@params$lambda)[variables, , drop = TRUE]
        names(coefs) <- gsub(paste0("^", feature_name, "_"), "", names(coefs))
    } else {
        # Create an empty coefficient vector if not present in model
        coefs <- numeric(0)
        warning("No coefficients found for 'cg_cont' feature in the model")
    }

    # Create and return the IQSeqFeature object
    IQSeqFeature(
        name = feature_name,
        coefs = coefs,
        compute_func = compute_cg_content,
        min_value = min_value,
        max_value = max_value,
        quantile = quantile,
        size = size
    )
}
