#' Plot Gene Detection Saturation
#'
#' @rdname plot_gene_detection_saturation
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#' plot_gene_detection_saturation(bcb)
.plot_gene_detection_saturation <- function(
    metrics,
    counts,
    interesting_group,
    min_counts = 1L) {
    if (is.null(metrics)) return(NULL)
    ggplot(
        metrics,
        aes_(x = ~mapped_reads / 1e6L,
             y = colSums(counts >= min_counts),
             color = as.name(interesting_group))) +
        geom_point(size = 3L) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(title = "gene detection saturation",
             x = "mapped reads (million)",
             y = "gene count",
             color = "",
             shape = "")
}



#' @rdname plot_gene_detection_saturation
#' @export
setMethod("plot_gene_detection_saturation", "bcbioRNADataSet", function(
    object, normalized = "tmm", ...) {
    .plot_gene_detection_saturation(
        metrics = metrics(object),
        counts = counts(object, normalized = normalized),
        interesting_group = .interesting_group(object),
        ...)
})
