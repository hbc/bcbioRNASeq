#' Deprecated functions
#'
#' @rdname deprecated
#' @name deprecated
#' @keywords internal
#'
#' @param ... Passthrough parameters.


#' @rdname deprecated
#' @export
correlation_heatmap <- function(...) {
    .Deprecated("plot_correlation_heatmap")
    plot_correlation_heatmap(...)
}

.create_local_project <- function() {
    .Deprecated("create_new_project")
}

#' @rdname deprecated
#' @export
import_file <- function() {
    .Deprecated("read_bcbio_file")
}

import_metadata <- function() {
    .Deprecated("read_bcbio_metadata")
}

import_summary <- function() {
    .Deprecated("read_bcbio_metrics")
}

load_bcbio_run <- function() {
    .Deprecated("load_run")
}

load_bcbio_run_from_yaml <- function() {
    .Deprecated("load_run")
}

load_run_S4 <- function(...) {
    .Deprecated("load_run")
    load_run(...)
}

plot_deseq_pca <- function() {
    .Deprecated("plot_pca")
}

pool_lane_split_counts <- function() {
    .Deprecated("pool_counts")
}

pool_lane_split_dds <- function() {
    .Deprecated("pool_dds")
}

read_metadata <- function() {
    .Deprecated(".metadata")
}
