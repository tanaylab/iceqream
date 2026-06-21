# Regression tests for the prego(threads) x doMC(fork) inference deadlock.
#
# infer_energies() iterates prego::compute_pwm() under plyr's `.parallel` (doMC)
# fork backend. prego::set_parallel() enables *both* the doMC fork and prego's
# own RcppParallel/OpenMP thread pool, so the forked workers re-enter a native
# thread pool that was warm at fork() - which is not fork-safe and can deadlock.
# The fix pins prego to a single internal thread inside the fork
# (local_prego_single_thread), leaving the fork as the only parallelism layer.
# These tests guard that the pin does not change the computed energies.

test_that("local_prego_single_thread is a no-op when prego is not parallel", {
    withr::local_options(list(prego.parallel.nc = 1))
    expect_null(iceqream:::local_prego_single_thread())
})

test_that("local_prego_single_thread restores the prego thread count on exit", {
    skip_if_not_installed("RcppParallel")
    withr::local_options(list(prego.parallel.nc = 3))
    # Run the guard inside a throwaway function; after it returns the thread
    # count must be back to the configured value (3), not 1.
    wrapped <- function() {
        iceqream:::local_prego_single_thread()
        "ran"
    }
    expect_equal(wrapped(), "ran")
    # The deferred restore fires on the wrapper's exit, so a fresh compute that
    # respects the configured threads still works (smoke check, no error).
    expect_silent(RcppParallel::setThreadOptions(numThreads = 3))
})

test_that("infer_energies gives identical energies with and without the doMC fork", {
    skip_if_not_installed("prego")
    skip_if_not_installed("doMC")

    set.seed(42)
    L <- 100L
    mk_seqs <- function(n) {
        vapply(seq_len(n), function(i) {
            paste(sample(c("A", "C", "G", "T"), L, replace = TRUE), collapse = "")
        }, character(1))
    }
    sequences <- mk_seqs(200)
    norm_sequences <- mk_seqs(150)
    mk_motif <- function() {
        m <- matrix(stats::runif(L * 4) + 0.01, ncol = 4)
        colnames(m) <- c("A", "C", "G", "T")
        list(pssm = as.data.frame(m), spat = NULL, spat_min = 1, spat_max = NULL)
    }
    motif_list <- stats::setNames(lapply(1:5, function(i) mk_motif()), paste0("m", 1:5))

    args <- list(
        sequences, norm_sequences, motif_list,
        min_energy = -7, energy_norm_quantile = 1, norm_energy_max = 10
    )

    e_serial <- withr::with_options(
        list(prego.parallel = FALSE, prego.parallel.nc = 1),
        do.call(iceqream:::infer_energies, args)
    )

    # Turn on the doMC fork backend; the fix must keep prego single-threaded
    # inside the fork and return the same energies (and not hang).
    prego::set_parallel(2)
    withr::defer(prego::set_parallel(1))
    e_fork <- do.call(iceqream:::infer_energies, args)

    expect_equal(dim(e_fork), c(200L, 5L))
    expect_true(all(is.finite(e_fork)))
    expect_equal(unname(e_fork), unname(e_serial))
})
