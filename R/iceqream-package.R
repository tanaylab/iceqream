#' @keywords internal
"_PACKAGE"

#' @importFrom cli cli_abort cli_warn cli_alert cli_alert_info cli_ul cli_alert_success cli_h1 cli_h2 cli_h3 cli_progress_bar cli_progress_update cli_progress_done cli_progress_along cli_alert_warning
#' @importFrom tgstat tgs_cor tgs_dist tgs_knn
#' @importFrom dplyr select mutate filter slice left_join right_join inner_join anti_join arrange desc  any_of one_of pull group_by slice_max ungroup n_distinct distinct n everything count bind_rows summarise n summarize row_number starts_with add_count case_when lead last sample_n join_by rename
#' @importFrom tidyr gather pivot_wider pivot_longer unite separate
#' @importFrom glue glue
#' @importFrom zoo rollmean
#' @importFrom tibble tibble deframe enframe column_to_rownames rownames_to_column
#' @importFrom stats ks.test lm
#' @importFrom purrr `%||%` map_dfr imap_dfr map
#' @importFrom Matrix t colSums rowSums rowMeans colMeans
#' @importFrom methods new slotNames show initialize
#' @importFrom stats binomial coef cutree hclust quantile
#' @importFrom utils head
#' @importFrom grDevices png pdf dev.off
#' @importFrom misha gextract gvtrack.create gintervals.all
#' @importFrom misha.ext gintervals.neighbors1 gintervals.normalize gintervals.expand
#' @import ggplot2
## usethis namespace: start
#' @importFrom ComplexHeatmap %v%
#' @importFrom rlang has_name
#' @importFrom stats hclust as.dist lag setNames var
#' @importFrom tgstat tgs_cor tgs_matrix_tapply tgs_dist
#' @importFrom prego gextract.local_pwm
## usethis namespace: end
NULL
