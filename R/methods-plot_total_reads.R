#' Plot Total Reads
#'
#' @rdname plot_total_reads
#' @family Quality Control Metrics Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit plot_metrics
#'
#' @examples
#' data(bcb)
#' plot_total_reads(bcb)
.plot_total_reads <- function(
    metrics,
    interesting_group,
    pass_limit = 20L,
    warn_limit = 10L,
    flip = TRUE) {
    p <- ggplot(metrics,
                aes_(x = ~sample_name,
                     y = ~total_reads / 1e6L,
                     fill = as.name(interesting_group))) +
        labs(title = "total reads",
             x = "sample",
             y = "total reads (million)",
             fill = "") +
        geom_bar(stat = "identity")
    if (!is.null(pass_limit)) {
        p <- p +
            geom_hline(alpha = 0.75,
                       color = pass_color,
                       size = 2L,
                       yintercept = pass_limit)
    }
    if (!is.null(warn_limit)) {
        p <- p +
            geom_hline(alpha = 0.75,
                       color = warn_color,
                       size = 2L,
                       yintercept = warn_limit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plot_total_reads
#' @export
setMethod("plot_total_reads", "bcbioRNADataSet", function(object, ...) {
    .plot_total_reads(
        metrics(object),
        interesting_group = .interesting_group(object),
        ...)
})
