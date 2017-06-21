#' Deprecated functions
#'
#' @rdname deprecated
#' @name deprecated
#' @keywords internal
#' @usage NULL
NULL



#' @rdname deprecated
#' @export
correlation_heatmap <- function(...) {
    .Deprecated("plot_correlation_heatmap")
    plot_correlation_heatmap(...)
}



#' @rdname deprecated
#' @export
create_local_project <- function() {
    .Deprecated()
}



#' @rdname deprecated
#' @export
create_new_project <- function() {
    .Deprecated()
}



#' @rdname deprecated
#' @export
gene_level_annotations <- function() {
}



#' @rdname deprecated
#' @export
import_file <- function() {
    .Deprecated()
}



#' @rdname deprecated
#' @export
import_metadata <- function() {
    .Deprecated()
}



#' @rdname deprecated
#' @export
import_summary <- function() {
    .Deprecated()
}



#' @rdname deprecated
#' @export
load_bcbio_run <- function(...) {
    .Deprecated("load_run")
    load_run(...)
}



#' @rdname deprecated
#' @export
load_bcbio_run_from_yaml <- function() {
    .Deprecated("load_run")
}



#' @rdname deprecated
#' @export
load_run_S4 <- function(...) {
    .Deprecated("load_run")
    load_run(...)
}



#' @rdname deprecated
#' @export
pool_counts <- function(...) {
    .Deprecated("aggregate_replicates")
    aggregate_replicates(...)
}



#' @rdname deprecated
#' @export
pool_dds <- function(...) {
    .Deprecated("aggregate_replicates")
    aggregate_replicates(...)
}



#' @rdname deprecated
#' @export
plot_deseq_pca <- function(...) {
    .Deprecated("plot_pca")
    plot_pca(...)
}



#' @rdname deprecated
#' @export
pool_lane_split_counts <- function(...) {
    .Deprecated("aggregate_replicates")
    aggregate_replicates(...)
}



#' @rdname deprecated
#' @export
pool_lane_split_dds <- function(...) {
    .Deprecated("aggregate_replicates")
    aggregate_replicates(...)
}



#' @rdname deprecated
#' @export
read_metadata <- function() {
    .Deprecated("metadata")
}



select_interesting_groups_coldata <- function() {
    .Deprecated()
}



#' @rdname deprecated
#' @export
tmm_normalize <- function() {
    .Deprecated("tmm")
}
