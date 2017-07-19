#' Plot Mapped Reads
#'
#' @rdname plot_mapped_reads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#' plot_mapped_reads(bcb)
.plot_mapped_reads <- function(
    metrics,
    interesting_group,
    pass_limit = 20L,
    warn_limit = 10L,
    flip = TRUE) {
    if (is.null(metrics)) return(NULL)
    p <- ggplot(metrics,
                aes_(x = ~sample_name,
                     y = ~mapped_reads / 1e6L,
                     fill = as.name(interesting_group))) +
        geom_bar(stat = "identity") +
        labs(title = "mapped reads",
             x = "sample",
             y = "mapped reads (million)",
             fill = "")
    if (!is.null(pass_limit)) {
        p <- p +
            geom_hline(alpha = 0.75,
                       color = pass_color,
                       size = 2L,
                       yintercept = pass_limit)
    }
    if (is.null(warn_limit)) {
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



#' @rdname plot_mapped_reads
#' @export
setMethod("plot_mapped_reads", "bcbioRNADataSet", function(object, ...) {
    .plot_mapped_reads(
        metrics(object),
        interesting_group = .interesting_group(object),
        ...)
})
