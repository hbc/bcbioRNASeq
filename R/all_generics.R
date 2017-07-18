#' S4 generics
#'
#' @rdname all_generics
#'
#' @param object Primary object.
#' @param value Value to assign.
#' @param ... Additional arguments.



#' @rdname all_generics
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



#' @rdname all_generics
#' @export
setGeneric("bcbio", function(object, ...) {
    standardGeneric("bcbio")
})



#' @rdname all_generics
#' @export
setGeneric("bcbio<-", function(object, ..., value) {
    standardGeneric("bcbio<-")
})



#' @rdname all_generics
#' @export
setGeneric("melt_log10", function(object, ...) {
    standardGeneric("melt_log10")
})



#' @rdname all_generics
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})



#' @rdname all_generics
#' @export
setGeneric("metrics", function(object) {
    standardGeneric("metrics")
})



#' @rdname all_generics
#' @export
setGeneric("plot_correlation_heatmap", function(object, ...) {
    standardGeneric("plot_correlation_heatmap")
})



#' @rdname all_generics
#' @export
setGeneric("plot_deg_heatmap", function(object, counts, ...) {
    standardGeneric("plot_deg_heatmap")
})



#' @rdname all_generics
#' @export
setGeneric("plot_dispersion", function(object) {
    standardGeneric("plot_dispersion")
})



#' @rdname all_generics
#' @export
setGeneric("plot_gender_markers", function(object, ...) {
    standardGeneric("plot_gender_markers")
})



#' @rdname all_generics
#' @export
setGeneric("plot_gene", function(object, ...) {
    standardGeneric("plot_gene")
})



#' @rdname all_generics
#' @export
setGeneric("plot_gene_heatmap", function(object, ...) {
    standardGeneric("plot_gene_heatmap")
})



#' @rdname all_generics
#' @export
setGeneric("plot_ma", function(object, ...) {
    standardGeneric("plot_ma")
})



#' @rdname all_generics
#' @export
setGeneric("plot_mean_sd", function(object, ...) {
    standardGeneric("plot_mean_sd")
})



#' @rdname all_generics
#' @export
setGeneric("plot_mirna_counts", function(object, ...) {
    standardGeneric("plot_mirna_counts")
})



#' @rdname all_generics
#' @export
setGeneric("plot_pattern", function(object, ...) {
    standardGeneric("plot_pattern")
})



#' @rdname all_generics
#' @export
setGeneric("plot_pca", function(object, ...) {
    standardGeneric("plot_pca")
})



#' @rdname all_generics
#' @export
setGeneric("plot_pca_covariates", function(object, ...) {
    standardGeneric("plot_pca_covariates")
})



#' @rdname all_generics
#' @export
setGeneric("plot_size_distribution", function(object, ...) {
    standardGeneric("plot_size_distribution")
})



#' @rdname all_generics
#' @export
setGeneric("plot_total_reads", function(object, ...) {
    standardGeneric("plot_total_reads")
})



#' @rdname all_generics
#' @export
setGeneric("plot_volcano", function(object, ...) {
    standardGeneric("plot_volcano")
})



#' @rdname all_generics
#' @export
setGeneric("sample_dirs", function(object) {
    standardGeneric("sample_dirs")
})



#' @rdname all_generics
#' @export
setGeneric("tmm", function(object) {
    standardGeneric("tmm")
})



#' @rdname all_generics
#' @export
setGeneric("tpm", function(object) {
    standardGeneric("tpm")
})



#' @rdname all_generics
#' @export
setGeneric("txi", function(object) {
    standardGeneric("txi")
})
