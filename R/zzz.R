.onLoad <- function(libname, pkgname) {
    # this is here in order to avoid R CMD check warnings
    globals <- c(".", "A", "atac_freq", "C", "chrom", "clust", "const", "cor", "diff_score", "dist", "e", "end", "feat", "feature", "freq", "freq_roll", "G", "intervalID", "max_beta", "metacell", "motif", "n_hits", "name", "new_clust", "obs_group", "observed", "ord", "pos", "predicted", "response", "s1", "spread", "start", "strand", "tmp_A", "tmp_C", "tmp_G", "tmp_T", "track", "type", "v", "value", "variable", "annot", "atac", "atac_n", "clust_i", "clust_name", "cumw", "d", "distilled", "energy", "ext_pos", "Gene", "geneSymbol", "grp", "i", "id", "intra_cor", "letter_pos", "lpos", "max_pr", "model", "motif1", "nuc", "rename", "size", "start1", "xmax", "xmin", "motif_db", "cluster", "term1", "term2", "param", "marginal_20k", "marginal", "marginal_20k_punc", "norm_f", "val", "punc", ".misha")
    utils::suppressForeignCheck(globals)
    utils::globalVariables(globals)
}
