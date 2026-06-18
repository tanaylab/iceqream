# Tests for the de-novo prego learning helpers in R/prego.R that don't need a
# genome / the prego engine.

test_that("regress_out_additional_features handles a single additional feature", {
    # glmnet needs >= 2 columns; learn_traj_prego() used to crash with
    # "x should be a matrix with 2 or more columns" when a user passed a single
    # additional-feature column. The helper now falls back to a logistic glm.
    set.seed(1)
    score <- rnorm(300)
    af1 <- data.frame(gc = rnorm(300))
    res <- suppressWarnings(
        iceqream:::regress_out_additional_features(af1, score)
    )
    expect_length(res, 300L)
    expect_true(all(is.finite(res)))
})

test_that("regress_out_additional_features works with multiple additional features", {
    set.seed(2)
    score <- rnorm(300)
    af2 <- data.frame(gc = rnorm(300), cpg = rnorm(300))
    res <- suppressWarnings(
        iceqream:::regress_out_additional_features(af2, score)
    )
    expect_length(res, 300L)
    expect_true(all(is.finite(res)))
})

test_that("regress_out_additional_features residual removes the feature's signal", {
    # A score that is a noisy linear function of the feature should leave a
    # residual that is far less correlated with the feature than the input was.
    set.seed(3)
    gc <- rnorm(300)
    score <- 3 * gc + rnorm(300, sd = 0.5)
    af1 <- data.frame(gc = gc)
    res <- suppressWarnings(
        iceqream:::regress_out_additional_features(af1, score)
    )
    cor_before <- abs(cor(iceqream:::norm01(score), gc))
    cor_after <- abs(cor(res, gc))
    expect_lt(cor_after, cor_before)
})
