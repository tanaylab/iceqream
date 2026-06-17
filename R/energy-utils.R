#' Normalize Values to the 0-1 Range
#'
#' This function scales and translates the values in the input vector `x`
#' such that the minimum value becomes 0 and the maximum value becomes 1.
#' The normalization is done by first subtracting the minimum value
#' from each element (translation) and then dividing each element by
#' the new maximum value (scaling).
#'
#' @param x A numeric vector that needs to be normalized to the 0-1 range. IF `x` is a matrix, the normalization is applied column-wise.
#'
#' @return A numeric vector with values normalized to the 0-1 range. If `x` is a matrix, the function returns a matrix with columns normalized to the 0-1 range.
#'
#' @examples
#' x <- rnorm(100)
#' normed_x <- norm01(x)
#'
#' @export
norm01 <- function(x) {
    if (is.matrix(x)) {
        x <- apply(x, 2, norm01)
        return(x)
    }
    x <- x - min(x, na.rm = TRUE)
    max_x <- max(x, na.rm = TRUE)
    if (max_x == 0) {
        return(x)
    }
    x / max_x
}


#' Normalize a vector to a specified quantile
#'
#' This function takes a vector \code{x} and normalizes it to a specified quantile.
#' The normalization is done by subtracting the minimum value of \code{x} and dividing
#' the result by the specified quantile of \code{x}. Values greater than 1 are capped
#' at 1.
#'
#' @param x A numeric vector to be normalized
#' @param quant The quantile to normalize to (default is 0.99)
#' @return A normalized vector
#' @examples
#' x <- rnorm(100)
#' normed_x <- norm0q(x)
#'
#' @export
norm0q <- function(x, quant = 0.99) {
    x <- x - min(x, na.rm = TRUE)
    x <- x / quantile(x, quant, na.rm = TRUE)
    x[x > 1] <- 1
    return(x)
}

#' Rescale Values Based on Original Scale Factors
#'
#' This function applies a reverse operation of normalization for a numeric vector `x`
#' based on the scale factors derived from an original vector `orig_x`. It effectively
#' attempts to map the values in `x` back to their original scale by applying
#' the inverse operations of scaling and translation based on the minimum and
#' maximum values found in `orig_x`. The process involves scaling `x` by the
#' maximum value (after subtracting the minimum value) of `orig_x` and then
#' adding the minimum value of `orig_x` to each element.
#'
#' @param x A numeric vector that needs to be rescaled to its original range.
#' @param orig_x The original numeric vector from which the scale factors are derived.
#'
#' @return A numeric vector with values rescaled to their original range based
#' on the scale factors from `orig_x`.
#'
#' @examples
#' orig_x <- rnorm(100)
#' rescaled_x <- rescale(norm01(orig_x), orig_x)
#'
#' @export
rescale <- function(x, orig_x) {
    norm_factors <- min(orig_x, na.rm = TRUE)
    norm_factors <- c(norm_factors, max(orig_x - norm_factors, na.rm = TRUE))
    x <- x * norm_factors[2]
    x <- x + norm_factors[1]
    return(x)
}


#' Normalize Energy Values of a Vector
#'
#' This function normalizes the energy values of a given vector. The normalization process
#' involves transforming the energy values using the exponential function, followed by a
#' logarithm with base 2. The function then scales the values to lie between a user-specified
#' minimum energy and the q-th quantile of the initial energy values.
#'
#' @param x A numeric vector containing the energy values to be normalized.
#' @param min_energy A numeric value representing the minimum energy value after normalization.
#' Default is set to -7.
#' @param q A numeric value between 0 and 1, representing the quantile of the energy values
#' to be used as the maximum after normalization. Default is set to 1 (max).
#'
#' @return A numeric vector of normalized energy values.
#'
#' @examples
#' # Generate random energy values
#' x <- runif(n = 100, min = -11, max = 0)
#' # Normalize the energy values
#' norm_energy(x)
#'
#' @export
norm_energy <- function(x, min_energy = -7, q = 1) {
    x <- exp(1)^x
    y <- log2(x / quantile(x, q, na.rm = TRUE))
    y[y > 0] <- 0
    y[y < min_energy] <- min_energy
    y <- y - min_energy
    return(y)
}

#' Normalize energy values using a reference dataset
#' @param x Numeric vector of energy values to normalize
#' @param dataset_x Numeric vector of reference energy values
#' @param min_energy Minimum energy value (default: -7)
#' @param q Quantile for normalization (default: 1)
#' @param norm_energy_max Maximum normalized energy value (default: 10)
#' @return Normalized energy values
#' @export
norm_energy_dataset <- function(x, dataset_x, min_energy = -7, q = 1, norm_energy_max = 10) {
    # Convert input from natural log to log2
    x_log2 <- x / log(2)
    dataset_x_log2 <- dataset_x / log(2)

    # Calculate the reference value in log2 space
    log2_max_x <- log2(quantile(2^dataset_x_log2, q, na.rm = TRUE))

    # Process dataset_x and x in log2 space
    process <- function(z) pmax(pmin(z - log2_max_x, 0), min_energy)
    dataset_y <- process(dataset_x_log2)
    y <- process(x_log2)

    # Scale y, handling potential division by zero
    range_dataset_y <- diff(range(dataset_y, na.rm = TRUE))
    if (range_dataset_y > 0) {
        y <- (y - min(dataset_y, na.rm = TRUE)) / range_dataset_y * norm_energy_max
    } else {
        y <- rep(0, length(y))
    }

    return(y)
}


#' Normalize Energy Matrix
#'
#' This function normalizes an energy matrix by applying logarithmic transformation and scaling.
#'
#' @param x The input matrix to be normalized.
#' @param dataset_x The reference dataset matrix used for normalization.
#' @param min_energy The minimum energy value to be assigned after normalization. Default is -7.
#' @param q The quantile value used for calculating the maximum value in the reference dataset. Default is 1.
#' @param norm_energy_max The maximum value to which the normalized energy values are scaled. Default is 10.
#'
#' @return A normalized energy matrix with the same dimensions as the input matrix.
#'
#' @examples
#' data <- matrix(rnorm(100), nrow = 10)
#' norm_energy_matrix(data, data)
#'
#' @export
norm_energy_matrix <- function(x, dataset_x = x, min_energy = -7, q = 1, norm_energy_max = 10) {
    # Coerce vectors/data.frames to matrices for safe subsetting
    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }
    if (!is.matrix(dataset_x)) {
        dataset_x <- as.matrix(dataset_x)
    }

    # Fill missing column names to allow alignment
    if (is.null(colnames(x)) && !is.null(colnames(dataset_x))) {
        colnames(x) <- colnames(dataset_x)
    }
    if (is.null(colnames(dataset_x)) && !is.null(colnames(x))) {
        colnames(dataset_x) <- colnames(x)
    }
    if (is.null(colnames(x)) && is.null(colnames(dataset_x))) {
        if (ncol(x) != ncol(dataset_x)) {
            stop("Input matrices have different column counts and no column names.")
        }
        default_cols <- paste0("V", seq_len(ncol(x)))
        colnames(x) <- default_cols
        colnames(dataset_x) <- default_cols
    }

    # Check for missing columns
    not_in_x <- setdiff(colnames(dataset_x), colnames(x))
    if (length(not_in_x) > 0) {
        stop(paste(
            "The following columns are missing in the input matrix:",
            paste(not_in_x, collapse = ", ")
        ))
    }

    # Align dataset_x columns with x
    dataset_x <- dataset_x[, colnames(x), drop = FALSE]

    # Convert input from natural log to log2
    x_log2 <- x / log(2)
    dataset_x_log2 <- dataset_x / log(2)

    # Calculate the reference values in log2 space
    log2_max_x <- log2(matrixStats::colQuantiles(2^dataset_x_log2, probs = q, na.rm = TRUE))

    # Process dataset_x and x in log2 space
    process <- function(z) pmax(pmin(sweep(z, 2, log2_max_x, `-`), 0), min_energy)
    dataset_y <- process(dataset_x_log2)
    y <- process(x_log2)

    # Calculate column-wise min and range for dataset_y
    range_dataset_y <- matrixStats::colRanges(dataset_y, na.rm = TRUE)
    range <- range_dataset_y[, 2] - range_dataset_y[, 1]
    min_dataset_y <- range_dataset_y[, 1]

    y <- t((t(y) - min_dataset_y) / (range)) * norm_energy_max
    y[is.na(y)] <- 0

    colnames(y) <- colnames(x)
    rownames(y) <- rownames(x)
    return(y)
}

#' Normalize motif energies using pre-computed quantiles
#'
#' @description
#' Fast-path normalization used by [compute_motif_energies()] when a
#' pre-computed `db_quantiles` matrix is supplied (e.g.
#' [mouse_db_quantiles]).
#'
#' @section Relationship to [norm_energy_matrix] (no downstream impact):
#' This function rescales the clamped energies to `[0, norm_energy_max]`
#' using the **fixed theoretical range** `-min_energy`
#' (`(y - min_energy) / (-min_energy) * norm_energy_max`), whereas
#' [norm_energy_matrix] (the `db_quantiles = NULL` path) rescales by the
#' **observed empirical range** of the background. The two formulas
#' therefore map the same input to different normalized values (pinned in
#' `test-energy-utils.R`). However, this denominator difference is a pure
#' affine rescale (correlation = 1) and has **no effect on a trained
#' model**: the only consumer of `compute_motif_energies()` output is the
#' correlation-based motif *selection* in [regress_trajectory_motifs()] /
#' [regress_trajectory_motifs_manifold()], which is scale-invariant. Both
#' regression paths then re-extract and re-normalize the selected motifs'
#' energies with [norm_energy_matrix] (observed range) during distillation
#' (`distill_motifs` / [distill_traj_model_multi]), so the model's
#' `@normalized_energies` are observed-range and match what inference
#' (`calc_traj_model_energies` / [pbm.normalize_energies]) recomputes -
#' `predict()` / [create_iq_model()] round-trip exactly whether or not
#' `db_quantiles` was used (verified on genome models, max|diff| = 0).
#'
#' The practically meaningful difference between the two paths is **not**
#' this denominator but the **reference quantile**: `db_quantiles` carries
#' a stable genome-wide reference, while the default path takes the
#' per-call observed quantile of `normalization_intervals`. That changes
#' *which* motifs pass the correlation-based selection threshold (and hence
#' the resulting model), which is the intended tradeoff of supplying
#' `db_quantiles` (no large normalization background required). It is a
#' selection choice, not a numerical bug, and does not need a fix.
#'
#' @param motif_energies Matrix of motif energies to normalize
#' @param db_quantiles Matrix of pre-computed quantiles for normalization
#' @param energy_norm_quantile Quantile to use for normalization
#' @param min_energy Minimum energy value
#' @param norm_energy_max Maximum normalized energy value
#'
#' @return Normalized motif energies matrix
#' @keywords internal
normalize_with_db_quantiles <- function(motif_energies, db_quantiles, energy_norm_quantile, min_energy, norm_energy_max) {
    # Check if energy_norm_quantile exists in db_quantiles columns
    if (!as.character(energy_norm_quantile) %in% colnames(db_quantiles)) {
        cli_abort("The specified energy_norm_quantile {.val {energy_norm_quantile}} is not present in db_quantiles.")
    }

    # Check if all motifs in motif_energies exist in db_quantiles
    missing_motifs <- setdiff(colnames(motif_energies), rownames(db_quantiles))
    if (length(missing_motifs) > 0) {
        cli_abort("The following motifs are missing from db_quantiles: {.val {missing_motifs}}")
    }

    # Ensure db_quantiles has the same motifs as motif_energies
    db_quantiles <- db_quantiles[colnames(motif_energies), , drop = FALSE]

    # Convert motif_energies to log2 space (as db_quantiles are in natural log)
    x_log2 <- motif_energies / log(2)

    # Get the reference values (quantiles) in log2 space
    log2_max_x <- db_quantiles[, as.character(energy_norm_quantile)] / log(2)

    # Process motif_energies in log2 space
    y <- pmax(pmin(sweep(x_log2, 2, log2_max_x, `-`), 0), min_energy)

    normalized_energies <- (y - min_energy) / (-min_energy) * norm_energy_max

    # Preserve dimension names
    colnames(normalized_energies) <- colnames(motif_energies)
    rownames(normalized_energies) <- rownames(motif_energies)

    return(normalized_energies)
}


#' Logistic Function
#'
#' Calculates the logistic function value given parameters and input.
#'
#' @param x Numeric vector, the values at which the logistic function will be evaluated.
#' @param x_0 Numeric, the x-value of the sigmoid's midpoint. Default is 0.
#' @param L Numeric, the maximum value of the sigmoid. Default is 1.
#' @param k Numeric, the steepness or slope of the sigmoid. Default is 1.
#'
#' @return A numeric vector of logistic function values.
#'
#' @examples
#' x_vals <- seq(0, 10, by = 0.1)
#' logist(x_vals, x_0 = 0, L = 2, k = 0.5)
#'
#' @export
logist <- function(x, x_0 = 0, L = 1, k = 1) {
    L / (1 + exp(-k * (x - x_0)))
}


#' Create Logistic Features
#'
#' This function takes a matrix or dataframe of features, removes columns that are entirely NA,
#' and then applies three logistic transformations to each column. Each transformed set of features
#' is appended with suffixes "_early", "_linear", or "_late" to differentiate between them.
#' The resulting matrix combines all transformed features.
#'
#' @param features A matrix or dataframe where each column is a feature to be transformed.
#'
#' @return A matrix containing the transformed features with columns named according to
#' the transformation applied (i.e., "_early", "_linear", or "_late").
#'
#' @seealso \code{\link{logist}} for the logistic transformation function.
#'
#' @examples
#' sample_features <- matrix(rnorm(100), ncol = 5)
#' create_logist_features(sample_features)
#'
#' @export
create_logist_features <- function(features) {
    if (is.null(nrow(features))) {
        features <- t(t(features))
    }
    # remove features that are all NA
    features <- features[, colSums(is.na(features)) != nrow(features), drop = FALSE]

    # No columns left (e.g. a model with no additional features, or an
    # all-dinucleotide additional-features matrix): return an empty matrix that
    # still carries the right number of rows so downstream cbind() is a no-op.
    if (ncol(features) == 0) {
        return(matrix(numeric(0), nrow = nrow(features), ncol = 0))
    }

    if (is.null(colnames(features))) {
        colnames(features) <- paste0("V", seq_len(ncol(features)))
    }

    features1 <- logist(features, x_0 = 0, L = 2, k = 0.5) - 1
    colnames(features1) <- paste0(colnames(features), "_low-energy")
    # features2 <- features
    # colnames(features2) <- paste0(colnames(features), "_linear")
    features3 <- logist(features, x_0 = 10, L = 2, k = 0.50)
    colnames(features3) <- paste0(colnames(features), "_high-energy")
    features4 <- logist(features - 5, x_0 = 0, L = 1, k = 1)
    colnames(features4) <- paste0(colnames(features), "_sigmoid")
    features5 <- logist(features, x_0 = 10, L = 2, k = 1)
    colnames(features5) <- paste0(colnames(features), "_higher-energy")
    # features6 <- logist(features, x_0 = 0, L = 2, k = 1) - 1
    # colnames(features6) <- paste0(colnames(features), "_early-2")

    features <- as.matrix(cbind(features1, features3, features4, features5))
    features <- features[, order(colnames(features))]

    return(features)
}

#' Homogenize a list of PSSM models
#'
#' This function adjusts each PSSM model in the input list so that the GC content is not higher
#' than the AT content. If a model's GC content is higher than its AT content, the function applies
#' a reverse complement to the model using the prego pssm_rc function.
#'
#' @param models A list of PSSM models. Each model should be a list with a pssm element, which
#'               is a data frame containing columns 'A', 'C', 'G', 'T', and 'pos'.
#' @return A list of homogenized PSSM models.
#'
#' @examples
#' # Create simulated data
#' pssm1 <- data.frame(
#'     pos = 1:4,
#'     A = c(0.1, 0.2, 0.3, 0.4),
#'     C = c(0.3, 0.3, 0.2, 0.1),
#'     G = c(0.3, 0.3, 0.3, 0.3),
#'     T = c(0.3, 0.2, 0.2, 0.2)
#' )
#' pssm2 <- data.frame(
#'     pos = 1:4,
#'     A = c(0.1, 0.2, 0.3, 0.4),
#'     C = c(0.1, 0.1, 0.1, 0.1),
#'     G = c(0.2, 0.2, 0.2, 0.2),
#'     T = c(0.6, 0.5, 0.4, 0.3)
#' )
#'
#' models <- list(list(pssm = pssm1), list(pssm = pssm2))
#'
#' # Homogenize the models
#' homogenized_models <- homogenize_pssm_models(models)
#'
#' @export
homogenize_pssm_models <- function(models) {
    models <- models %>%
        purrr::map(homogenize_model)

    return(models)
}

homogenize_model <- function(model) {
    # if T is higher than A reverse complement the model
    if (sum(model$pssm$T) > sum(model$pssm$A)) {
        model$pssm <- prego::pssm_rc(model$pssm)
    }
    return(model)
}

inverse_logist_by_type <- function(y, type) {
    switch(type,
        "low-energy" = {
            # Reverse of: logist(x, x_0 = 0, L = 2, k = 0.5) - 1
            -2 * log(2 / (y + 1) - 1)
        },
        "high-energy" = {
            # Reverse of: logist(x, x_0 = 10, L = 2, k = 0.50)
            10 - 2 * log(2 / y - 1)
        },
        "sigmoid" = {
            # Reverse of: logist(x - 5, x_0 = 0, L = 1, k = 1)
            -log(1 / y - 1) + 5
        },
        "higher-energy" = {
            # Reverse of: logist(x, x_0 = 10, L = 2, k = 1)
            10 - log(2 / y - 1)
        },
        stop("Unknown transformation type")
    )
}
