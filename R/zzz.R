.onLoad <- function(libname, pkgname) {
    # this is here in order to avoid R CMD check warnings
    globals <- c(".", "A", "atac_freq", "C", "chrom", "clust", "const", "cor", "diff_score", "dist", "e", "end", "feat", "feature", "freq", "freq_roll", "G", "intervalID", "max_beta", "metacell", "motif", "n_hits", "name", "new_clust", "obs_group", "observed", "ord", "pos", "predicted", "response", "s1", "spread", "start", "strand", "tmp_A", "tmp_C", "tmp_G", "tmp_T", "track", "type", "v", "value", "variable")
    utils::suppressForeignCheck(globals)
    utils::globalVariables(globals)
}
