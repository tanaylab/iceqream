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
#' # Generate random values
#' x <- rnorm(100)
#' range(x)
#' normed_x <- norm01(x)
#' range(normed_x) # This should show values between 0 and 1
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
#' range(normed_x) # This should show values between 0 and 1
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
#' # Generate random values and normalize
#' orig_x <- rnorm(100)
#' normed_x <- norm01(orig_x)
#' # Rescale normalized values back to original range
#' rescaled_x <- rescale(normed_x, orig_x)
#' range(rescaled_x) # This should closely match the range of orig_x
#' range(orig_x)
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

norm_energy_dataset_old <- function(x, dataset_x, min_energy = -7, q = 1, norm_energy_max = 10) {
    dataset_x <- exp(1)^dataset_x
    max_x <- quantile(dataset_x, q, na.rm = TRUE)
    x <- exp(1)^x
    y <- log2(x / max_x)
    y[y > 0] <- 0
    y[y < min_energy] <- min_energy
    y <- y - min_energy

    # y <- norm01(y) * norm_energy_max
    # Process dataset_x in the same way
    dataset_y <- log2(dataset_x / max_x)
    dataset_y[dataset_y > 0] <- 0
    dataset_y[dataset_y < min_energy] <- min_energy
    dataset_y <- dataset_y - min_energy

    # Calculate min and max of dataset_y for scaling
    min_dataset_y <- min(dataset_y, na.rm = TRUE)
    max_dataset_y <- max(dataset_y, na.rm = TRUE)

    # Ensure y is not less than min_dataset_y
    y <- pmax(y, min_dataset_y)

    # Scale y using both min and max of dataset_y, ensuring non-negative output
    y <- (y - min_dataset_y) / (max_dataset_y - min_dataset_y) * norm_energy_max

    # Ensure output is non-negative (in case of numerical precision issues)
    y <- pmax(y, 0)
    return(y)
}

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



norm_energy_matrix_old <- function(x, dataset_x, min_energy = -7, q = 1, norm_energy_max = 10) {
    not_in_x <- colnames(dataset_x)[!(colnames(dataset_x) %in% colnames(x))]
    if (length(not_in_x) > 0) {
        cli_abort("The following columns are missing in the input matrix: {.val {not_in_x}}")
    }
    dataset_x <- dataset_x[, colnames(x)]
    dataset_x <- exp(1)^dataset_x
    max_x <- matrixStats::colQuantiles(dataset_x, probs = q, na.rm = TRUE)
    x <- exp(1)^x
    y <- log2(t(t(x) / max_x))
    y[y > 0] <- 0
    y[y < min_energy] <- min_energy
    y <- y - min_energy

    y <- norm01(y) * norm_energy_max
    colnames(y) <- colnames(x)
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
#' # Example usage:
#' data <- matrix(rnorm(100), nrow = 10)
#' normalized_data <- norm_energy_matrix(data, data, min_energy = -7, q = 1, norm_energy_max = 10)
#'
#' @export
norm_energy_matrix <- function(x, dataset_x, min_energy = -7, q = 1, norm_energy_max = 10) {
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

    colnames(y) <- colnames(x)
    rownames(y) <- rownames(x)
    return(y)
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
#'
#' # Calculate the features for each scenario
#' features_low_energy <- logist(x_vals, x_0 = 0, L = 2, k = 0.5) - 1
#' features_high_energy <- logist(x_vals, x_0 = 10, L = 2, k = 0.5)
#' features_sigmoid <- logist(x_vals - 5, x_0 = 0, L = 1, k = 1)
#' features_higher_energy <- logist(x_vals, x_0 = 10, L = 2, k = 1)
#' features_early2 <- logist(x_vals, x_0 = 0, L = 2, k = 1) - 1
#'
#' # Base plot setup
#' plot(x_vals, features_low_energy * 10,
#'     type = "l", col = "blue",
#'     main = "Variations of the Logistic Function",
#'     xlab = "x", ylab = "y", ylim = c(0, 10), lwd = 2
#' )
#'
#' # Adding other variations
#' lines(x_vals, features_high_energy * 10, col = "orange", lwd = 2)
#' lines(x_vals, features_sigmoid * 10, col = "purple", lwd = 2)
#' lines(x_vals, features_higher_energy * 10, col = "brown", lwd = 2)
#' lines(x_vals, features_early2 * 10, col = "green", lwd = 2)
#' lines(x_vals, x_vals, col = "black", lwd = 2, lty = 2)
#'
#' legend("bottomright",
#'     legend = c("Low Energy", "High Energy", "Sigmoid", "Higher Energy", "Early 2", "Linear"),
#'     col = c("blue", "orange", "purple", "brown", "green", "black"),
#'     lty = 1,
#'     lwd = 2
#' )
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
#' # Create a sample matrix
#' sample_features <- matrix(rnorm(100), ncol = 5)
#' transformed_features <- create_logist_features(sample_features)
#' head(transformed_features)
#'
#' @export
create_logist_features <- function(features) {
    # remove features that are all NA
    features <- features[, colSums(is.na(features)) != nrow(features), drop = FALSE]

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
