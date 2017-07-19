#' Plot Genes Detected
#'
#' @rdname plot_genes_detected
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#' plot_genes_detected(bcb, pass_limit = NULL)
.plot_genes_detected <- function(
    metrics,
    counts,
    interesting_group,
    pass_limit = 20000L,
    min_counts = 1L,
    flip = TRUE) {
    if (is.null(metrics)) return(NULL)
    p <- ggplot(metrics,
                aes_(x = ~sample_name,
                     y = colSums(counts >= min_counts),
                     fill = as.name(interesting_group))) +
        geom_bar(stat = "identity") +
        labs(title = "genes detected",
             x = "sample",
             y = "gene count",
             fill = "")
    if (!is.null(pass_limit)) {
        p <- p + geom_hline(alpha = 0.75,
                            color = pass_color,
                            size = 2L,
                            yintercept = pass_limit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plot_genes_detected
#' @export
setMethod("plot_genes_detected", "bcbioRNADataSet", function(object, ...) {
    .plot_genes_detected(
        metrics = metrics(object),
        counts = assay(object),
        interesting_group = .interesting_group(object),
        ...)
})
