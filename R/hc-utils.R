get_correlation_based_clusters <- function(hc, cor_matrix, threshold = 0.7) {
    # Verify we have colnames
    if (is.null(colnames(cor_matrix))) {
        stop("cor_matrix must have column names")
    }

    # Get the merge matrix from hclust object
    merge <- hc$merge
    n <- nrow(merge) + 1 # number of original observations
    names <- colnames(cor_matrix)

    # Function to get all leaves under a node
    get_leaves <- function(node) {
        if (node < 0) {
            return(-node)
        } # Negative numbers in merge matrix are leaves

        # Recursive call for internal nodes
        left_leaves <- get_leaves(merge[node, 1])
        right_leaves <- get_leaves(merge[node, 2])
        return(c(left_leaves, right_leaves))
    }

    # Function to calculate average intra-correlation for a group of indices
    calc_avg_correlation <- function(indices) {
        if (length(indices) == 1) {
            return(1)
        }

        # Get the correlation submatrix
        cor_submatrix <- cor_matrix[indices, indices]

        # Calculate mean correlation (excluding diagonal)
        mean_cor <- (sum(cor_submatrix) - length(indices)) /
            (length(indices)^2 - length(indices))
        return(mean_cor)
    }

    # Vector to store cluster assignments and correlations
    clusters <- rep(NA, n)
    intra_correlations <- rep(NA, n)
    next_cluster <- 1

    # Process each merge point from bottom to top
    for (i in nrow(merge):1) {
        # Get leaves under current merge
        current_leaves <- get_leaves(i)

        # Calculate average intra-correlation
        avg_cor <- calc_avg_correlation(current_leaves)

        # If correlation is above threshold and no cluster assignment yet
        if (avg_cor >= threshold && any(is.na(clusters[current_leaves]))) {
            # Assign cluster number and correlation to all unassigned leaves in this group
            unassigned <- current_leaves[is.na(clusters[current_leaves])]
            clusters[unassigned] <- next_cluster
            intra_correlations[unassigned] <- avg_cor
            next_cluster <- next_cluster + 1
        }
    }

    # Assign remaining unclustered points to individual clusters
    unclustered <- which(is.na(clusters))
    if (length(unclustered) > 0) {
        for (i in unclustered) {
            clusters[i] <- next_cluster
            intra_correlations[i] <- 1 # Single element clusters have intra-correlation of 1
            next_cluster <- next_cluster + 1
        }
    }

    # Create and return the data frame
    result_df <- data.frame(
        name = names,
        cluster = clusters,
        intra_cor = intra_correlations
    )

    # Sort by cluster and then by name
    result_df <- result_df[order(result_df$cluster, result_df$name), ]

    return(result_df)
}
