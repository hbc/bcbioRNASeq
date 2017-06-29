#' Deprecated functions
#'
#' @rdname deprecated
#' @keywords internal
#'
#' @param ... Additional parameters.



#' @rdname deprecated
#' @export
correlation_heatmap <- function(...) {
    .Deprecated("plot_correlation_heatmap")
    plot_correlation_heatmap(...)
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
load_run_S4 <- function(...) {  # nolint
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



#' @rdname deprecated
#' @export
tmm_normalize <- function() {
    .Deprecated("tmm")
}
