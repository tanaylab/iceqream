#' Normalize Values to the \[0, 1\] Range
#'
#' This function scales and translates the values in the input vector `x`
#' such that the minimum value becomes 0 and the maximum value becomes 1.
#' The normalization is done by first subtracting the minimum value
#' from each element (translation) and then dividing each element by
#' the new maximum value (scaling).
#'
#' @param x A numeric vector that needs to be normalized to the [0, 1] range.
#'
#' @return A numeric vector with values normalized to the [0, 1] range.
#'
#' @examples
#' # Generate random values
#' x <- rnorm(100)
#' range(x)
#' # Normalize the values to [0, 1] range
#' normed_x <- norm01(x)
#' range(normed_x) # This should show values between 0 and 1
#'
#' @export
norm01 <- function(x) {
    x <- x - min(x, na.rm = TRUE)
    x / max(x, na.rm = TRUE)
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
#' Default is set to -10.
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
norm_energy <- function(x, min_energy = -10, q = 1) {
    x <- exp(1)^x
    y <- log2(x / quantile(x, q, na.rm = TRUE))
    y[y > 0] <- 0
    y[y < min_energy] <- min_energy
    y <- y - min_energy
    return(y)
}


norm_energy_dataset <- function(x, ds_x, min_energy = -10, q = 1) {
    x_both <- c(x, ds_x)
    x_idx <- seq_along(x)
    y_both <- norm_energy(x_both, min_energy, q)
    y <- y_both[x_idx]
    return(y)
}

#' Rescale numeric values to a specified range
#'
#' This function rescales numeric values from their original range to a specified
#' new range, which defaults to between -1 and 1.
#'
#' @param x A numeric vector of values to be rescaled.
#' @param new_min The minimum value of the desired output range. Defaults to -1.
#' @param new_max The maximum value of the desired output range. Defaults to 1.
#'
#' @return A numeric vector of rescaled values.
#' @examples
#' values <- c(2, 3, 5, 6, 8, 10)
#' rescale(values)
#' rescale(values, 0, 100)
#' rescale(values, 5, 50)
#'
#' @export
rescale <- function(x, new_min = -1, new_max = 1) {
    old_min <- min(x)
    old_max <- max(x)

    return(((x - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min)
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
