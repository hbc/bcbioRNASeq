#' Plot Total Reads
#'
#' @rdname plot_total_reads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plot_total_reads(bcb)
#'
#' # data.frame
#' metrics <- metrics(bcb)
#' plot_total_reads(metrics)



#' @rdname plot_total_reads
.plot_total_reads <- function(
    metrics,
    interesting_group,
    pass_limit = 20L,
    warn_limit = 10L,
    flip = TRUE) {
    if (is.null(metrics)) return(NULL)
    p <- ggplot(metrics,
                aes_(x = ~sample_name,
                     y = ~total_reads / 1e6L,
                     fill = as.name(interesting_group))) +
        geom_bar(stat = "identity") +
        labs(title = "total reads",
             x = "sample",
             y = "total reads (million)")
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



#' @rdname plot_total_reads
#' @export
setMethod("plot_total_reads", "data.frame", function(
    object, interesting_group = "sample_name", ...) {
    .plot_total_reads(object, interesting_group, ...)
})
