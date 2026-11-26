test_that("distill_motifs passes each prego worker only its cluster's columns", {
    set.seed(1)
    n <- 200
    feats <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, paste0("m", 1:6)))
    feats[, 2] <- feats[, 1] + rnorm(n, sd = 0.01) # m1 and m2 fall in one cluster
    glm_model <- list(beta = matrix(seq(0.6, 0.1, length.out = 6), 6, 1, dimnames = list(colnames(feats), "s0")))
    seen <- list()
    local_mocked_bindings(
        run_prego_on_clust_residuals = function(model, feats, clust_motifs, ...) {
            seen[[length(seen) + 1]] <<- list(cols = colnames(feats), clust = clust_motifs)
            list()
        },
        infer_energies_new = function(seqs, norm_seqs, motifs, ...) {
            matrix(0, length(seqs), length(motifs), dimnames = list(NULL, names(motifs)))
        }
    )
    withr::local_options(prego.parallel = FALSE)
    seqs <- rep(strrep("A", 50), n)
    suppressMessages(distill_motifs(feats,
        target_number = 5, glm_model = glm_model, y = rnorm(n), seqs = seqs,
        norm_seqs = seqs, diff_filter = rep(TRUE, n), pssm_db = data.frame()
    ))
    expect_length(seen, 5)
    for (s in seen) expect_identical(s$cols, s$clust)
    expect_true(any(lengths(lapply(seen, `[[`, "cols")) == 2))
})
