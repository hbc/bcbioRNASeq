#' Plot Counts Per Gene
#'
#' @rdname plot_counts_per_gene
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#' plot_counts_per_gene(bcb)



#' @rdname plot_counts_per_gene
.plot_counts_per_gene <- function(
    melted,
    interesting_group = "sample_name",
    flip = TRUE) {
    p <- ggplot(melted,
                aes_(x = ~sample_name,
                     y = ~counts,
                     color = as.name(interesting_group))) +
        geom_boxplot(outlier.shape = NA) +
        labs(title = "counts per gene",
             x = "sample",
             y = "log10 counts per gene")
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plot_counts_per_gene
#' @export
setMethod("plot_counts_per_gene", "bcbioRNADataSet", function(
    object, normalized = "tmm", ...) {
    .plot_counts_per_gene(
        melted = melt_log10(object, normalized = normalized),
        interesting_group = .interesting_group(object),
        ...)
})
