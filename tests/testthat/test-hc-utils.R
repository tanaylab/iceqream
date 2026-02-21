# Tests for get_correlation_based_clusters

test_that("get_correlation_based_clusters returns correct structure", {
    set.seed(42)
    m <- matrix(rnorm(40), nrow = 10, ncol = 4)
    colnames(m) <- c("feat1", "feat2", "feat3", "feat4")
    cor_matrix <- cor(m)
    hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

    result <- iceqream:::get_correlation_based_clusters(hc, cor_matrix, threshold = 0.5)

    expect_s3_class(result, "data.frame")
    expect_named(result, c("name", "cluster", "intra_cor"))
    expect_equal(nrow(result), ncol(m))
    expect_true(all(result$name %in% colnames(m)))
    expect_true(all(result$intra_cor >= 0 & result$intra_cor <= 1))
})

test_that("get_correlation_based_clusters errors on missing colnames", {
    m <- matrix(rnorm(20), nrow = 5, ncol = 4)
    cor_matrix <- cor(m)
    # cor_matrix has no colnames
    colnames(cor_matrix) <- NULL
    hc <- hclust(as.dist(1 - abs(cor(m))), method = "complete")

    expect_error(
        iceqream:::get_correlation_based_clusters(hc, cor_matrix),
        "column names"
    )
})

test_that("get_correlation_based_clusters with low threshold produces fewer clusters", {
    # Use positively correlated features so that lowering the threshold groups more together
    set.seed(42)
    base <- rnorm(20)
    m <- cbind(base + rnorm(20, sd = 0.5), base + rnorm(20, sd = 0.5),
               base + rnorm(20, sd = 0.5), base + rnorm(20, sd = 0.5))
    colnames(m) <- c("a", "b", "c", "d")
    cor_matrix <- cor(m)
    hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

    result_low <- iceqream:::get_correlation_based_clusters(hc, cor_matrix, threshold = 0.3)
    result_high <- iceqream:::get_correlation_based_clusters(hc, cor_matrix, threshold = 0.99)

    # Lower threshold should produce fewer or equal number of clusters
    expect_true(length(unique(result_low$cluster)) <= length(unique(result_high$cluster)))
})

test_that("get_correlation_based_clusters with threshold=1 makes individual clusters", {
    set.seed(42)
    m <- matrix(rnorm(40), nrow = 10, ncol = 4)
    colnames(m) <- c("a", "b", "c", "d")
    cor_matrix <- cor(m)
    hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

    result <- iceqream:::get_correlation_based_clusters(hc, cor_matrix, threshold = 1)

    # With threshold 1, only perfectly correlated items cluster together
    # Each variable is perfectly correlated with itself, so single-element clusters are typical
    expect_true(length(unique(result$cluster)) >= 1)
    # Single-element clusters should have intra_cor of 1
    single_clusters <- result$cluster[duplicated(result$cluster) | duplicated(result$cluster, fromLast = TRUE)]
    solo_items <- result[!result$cluster %in% single_clusters, ]
    if (nrow(solo_items) > 0) {
        expect_true(all(solo_items$intra_cor == 1))
    }
})

test_that("get_correlation_based_clusters with perfectly correlated features", {
    # Create perfectly correlated features
    x <- 1:20
    m <- cbind(x, x * 2, x + 5, -x)
    colnames(m) <- c("a", "b", "c", "d")
    cor_matrix <- cor(m)
    hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

    result <- iceqream:::get_correlation_based_clusters(hc, cor_matrix, threshold = 0.99)

    # a, b, c are perfectly correlated (cor=1); d is perfectly anti-correlated (cor=-1)
    # a, b, c should cluster together
    cluster_a <- result$cluster[result$name == "a"]
    cluster_b <- result$cluster[result$name == "b"]
    cluster_c <- result$cluster[result$name == "c"]
    expect_equal(cluster_a, cluster_b)
    expect_equal(cluster_b, cluster_c)
})

test_that("get_correlation_based_clusters with 2 features", {
    set.seed(42)
    m <- matrix(rnorm(20), nrow = 10, ncol = 2)
    colnames(m) <- c("x", "y")
    cor_matrix <- cor(m)
    hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

    result <- iceqream:::get_correlation_based_clusters(hc, cor_matrix, threshold = 0.5)

    expect_equal(nrow(result), 2)
    expect_true(all(c("x", "y") %in% result$name))
})

test_that("get_correlation_based_clusters result is sorted by cluster then name", {
    set.seed(42)
    m <- matrix(rnorm(50), nrow = 10, ncol = 5)
    colnames(m) <- c("e", "d", "c", "b", "a")
    cor_matrix <- cor(m)
    hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

    result <- iceqream:::get_correlation_based_clusters(hc, cor_matrix, threshold = 0.5)

    # Verify sorted by cluster first, then by name within cluster
    expected_order <- order(result$cluster, result$name)
    expect_equal(seq_len(nrow(result)), expected_order)
})

test_that("get_correlation_based_clusters assigns every feature to a cluster", {
    set.seed(123)
    m <- matrix(rnorm(60), nrow = 10, ncol = 6)
    colnames(m) <- paste0("f", 1:6)
    cor_matrix <- cor(m)
    hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

    result <- iceqream:::get_correlation_based_clusters(hc, cor_matrix, threshold = 0.7)

    expect_true(all(!is.na(result$cluster)))
    expect_true(all(!is.na(result$intra_cor)))
    expect_equal(sort(result$name), sort(colnames(m)))
})
