#' S4 Generics
#'
#' @rdname all_generics
#' @name all_generics
#'
#' @param object Object.
#' @param value Value to assign.
#' @param ... Additional arguments.
NULL



#' Quality Control Plots
#'
#' @rdname qc_plots
#' @name qc_plots
#'
#' @param pass_limit Threshold to plot pass color marker.
#' @param warn_limit Threshold to plot warning color marker.
#' @param interesting_group *Optional*. Category to use to group samples (color
#'   and shape). If unset, this is automatically determined by the metadata set
#'   inside the [bcbioRNADataSet].
#' @param min_counts Numeric value for filtering the counts matrix before
#'   plotting.
#' @param normalized Count normalization method. See [counts()] documentation
#'   for more information.
#'
#' @return [ggplot].
NULL



#' @rdname aggregate_replicates
#' @inherit all_generics
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



#' @rdname alpha_summary
#' @inherit all_generics
#' @export
setGeneric("alpha_summary", function(object, ...) {
    standardGeneric("alpha_summary")
})



#' @rdname bcbio
#' @inherit all_generics
#' @export
setGeneric("bcbio", function(object, ...) {
    standardGeneric("bcbio")
})



#' @rdname bcbio
#' @inherit all_generics
#' @export
setGeneric("bcbio<-", function(object, ..., value) {
    standardGeneric("bcbio<-")
})



#' @rdname melt_log10
#' @inherit all_generics
#' @export
setGeneric("melt_log10", function(object, ...) {
    standardGeneric("melt_log10")
})



#' @rdname metadata_table
#' @inherit all_generics
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})



#' @rdname metrics
#' @inherit all_generics
#' @export
setGeneric("metrics", function(object) {
    standardGeneric("metrics")
})



#' @rdname plot_correlation_heatmap
#' @inherit all_generics
#' @export
setGeneric("plot_correlation_heatmap", function(object, ...) {
    standardGeneric("plot_correlation_heatmap")
})



#' @rdname plot_count_density
#' @inherit all_generics
#' @export
setGeneric("plot_count_density", function(object, ...) {
    standardGeneric("plot_count_density")
})



#' @rdname plot_counts_per_gene
#' @inherit all_generics
#' @export
setGeneric("plot_counts_per_gene", function(object, ...) {
    standardGeneric("plot_counts_per_gene")
})



#' @rdname plot_deg_heatmap
#' @inherit all_generics
#' @export
setGeneric("plot_deg_heatmap", function(object, counts, ...) {
    standardGeneric("plot_deg_heatmap")
})



#' @rdname plot_dispersion
#' @inherit all_generics
#' @export
setGeneric("plot_dispersion", function(object) {
    standardGeneric("plot_dispersion")
})



#' @rdname plot_exonic_mapping_rate
#' @inherit all_generics
#' @export
setGeneric("plot_exonic_mapping_rate", function(object, ...) {
    standardGeneric("plot_exonic_mapping_rate")
})



#' @rdname plot_gender_markers
#' @inherit all_generics
#' @export
setGeneric("plot_gender_markers", function(object, ...) {
    standardGeneric("plot_gender_markers")
})



#' @rdname plot_gene
#' @inherit all_generics
#' @export
setGeneric("plot_gene", function(object, ...) {
    standardGeneric("plot_gene")
})



#' @rdname plot_gene_detection_saturation
#' @inherit all_generics
#' @export
setGeneric("plot_gene_detection_saturation", function(object, ...) {
    standardGeneric("plot_gene_detection_saturation")
})



#' @rdname plot_gene_heatmap
#' @inherit all_generics
#' @export
setGeneric("plot_gene_heatmap", function(object, ...) {
    standardGeneric("plot_gene_heatmap")
})



#' @rdname plot_genes_detected
#' @inherit all_generics
#' @export
setGeneric("plot_genes_detected", function(object, ...) {
    standardGeneric("plot_genes_detected")
})



#' @rdname plot_ma
#' @inherit all_generics
#' @export
setGeneric("plot_ma", function(object, ...) {
    standardGeneric("plot_ma")
})



#' @rdname plot_mapped_reads
#' @inherit all_generics
#' @export
setGeneric("plot_mapped_reads", function(object, ...) {
    standardGeneric("plot_mapped_reads")
})



#' @rdname plot_mapping_rate
#' @inherit all_generics
#' @export
setGeneric("plot_mapping_rate", function(object, ...) {
    standardGeneric("plot_mapping_rate")
})



#' @rdname plot_mean_sd
#' @inherit all_generics
#' @export
setGeneric("plot_mean_sd", function(object, ...) {
    standardGeneric("plot_mean_sd")
})



#' @rdname plot_mirna_counts
#' @inherit all_generics
#' @export
setGeneric("plot_mirna_counts", function(object, ...) {
    standardGeneric("plot_mirna_counts")
})



#' @rdname plot_pattern
#' @inherit all_generics
#' @export
setGeneric("plot_pattern", function(object, ...) {
    standardGeneric("plot_pattern")
})



#' @rdname plot_pca
#' @inherit all_generics
#' @export
setGeneric("plot_pca", function(object, ...) {
    standardGeneric("plot_pca")
})



#' @rdname plot_pca_covariates
#' @inherit all_generics
#' @export
setGeneric("plot_pca_covariates", function(object, ...) {
    standardGeneric("plot_pca_covariates")
})



#' @rdname plot_size_distribution
#' @inherit all_generics
#' @export
setGeneric("plot_size_distribution", function(object, ...) {
    standardGeneric("plot_size_distribution")
})



#' @rdname plot_total_reads
#' @inherit all_generics
#' @export
setGeneric("plot_total_reads", function(object, ...) {
    standardGeneric("plot_total_reads")
})



#' @rdname plot_volcano
#' @inherit all_generics
#' @export
setGeneric("plot_volcano", function(object, ...) {
    standardGeneric("plot_volcano")
})



#' @rdname sample_dirs
#' @inherit all_generics
#' @export
setGeneric("sample_dirs", function(object) {
    standardGeneric("sample_dirs")
})



#' @rdname tmm
#' @inherit all_generics
#' @export
setGeneric("tmm", function(object) {
    standardGeneric("tmm")
})



#' @rdname tpm
#' @inherit all_generics
#' @export
setGeneric("tpm", function(object) {
    standardGeneric("tpm")
})



#' @rdname txi
#' @inherit all_generics
#' @export
setGeneric("txi", function(object) {
    standardGeneric("txi")
})
