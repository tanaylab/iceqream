norm01 <- function(x) {
    x <- x - min(x, na.rm = TRUE)
    x / max(x, na.rm = TRUE)
}

#' Normalize energy of a vector
#'
#' @param x A vector of energy values
#' @param min_energy Minimum energy value
#' @param q Quantile of the energies to consider as maximum
#'
#' @examples
#' x <- runif(n = 100, min = -11, max = 0)
#' norm_energy(x)
#'
#' @export
norm_energy <- function(x, min_energy = -10, q = 0.99) {
    x <- exp(1)^x
    y <- log2(x / quantile(x, q, na.rm = TRUE))
    y[y > 1] <- 1
    y[y < min_energy] <- min_energy
    y <- y - min(y, na.rm = TRUE)
    y
}

norm_energy_dataset <- function(x, ds_x, min_energy = -10, q = 0.99) {
    x_both <- c(x, ds_x)
    x_idx <- 1:length(x)
    y_both <- norm_energy(x_both, min_energy, q)
    y <- y_both[x_idx]
    return(y)
}

logist <- function(x, x_0 = 0, L = 1, k = 1) {
    L / (1 + exp(-k * (x - x_0)))
}

create_logist_features <- function(features) {
    # remove features that are all NA
    features <- features[, colSums(is.na(features)) != nrow(features)]

    features1 <- logist(features, x_0 = 0, L = 2, k = 1) - 1
    colnames(features1) <- paste0(colnames(features), "_early")
    features2 <- logist(features, x_0 = 0, L = 2, k = 0.5) - 1
    colnames(features2) <- paste0(colnames(features), "_linear")
    features3 <- logist(features, x_0 = 0, L = 2, k = 0.25) - 1
    colnames(features3) <- paste0(colnames(features), "_late")

    features <- as.matrix(cbind(features1, features2, features3))

    return(features)
}
