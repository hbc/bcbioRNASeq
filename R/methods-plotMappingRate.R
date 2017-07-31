#' Plot Mapping Rate
#'
#' @rdname plot_mapping_rate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#' plot_mapping_rate(bcb)
.plot_mapping_rate <- function(
    metrics,
    interesting_group = "sample_name",
    pass_limit = 90L,
    warn_limit = 70L,
    flip = TRUE) {
    if (is.null(metrics)) return(NULL)
    p <- ggplot(metrics,
                aes_(x = ~sample_name,
                     y = ~mapped_reads / total_reads * 100L,
                     fill = as.name(interesting_group))) +
        geom_bar(stat = "identity") +
        ylim(0L, 100L) +
        labs(title = "mapping rate",
             x = "sample",
             y = "mapping rate (%)")
    if (!is.null(pass_limit)) {
        p <- p +
            geom_hline(alpha = qc_line_alpha,
                       color = qc_pass_color,
                       size = qc_line_size,
                       yintercept = pass_limit)
    }
    if (!is.null(warn_limit)) {
        p <- p +
            geom_hline(alpha = qc_line_alpha,
                       color = qc_warn_color,
                       size = qc_line_size,
                       yintercept = warn_limit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plot_mapping_rate
#' @export
setMethod("plot_mapping_rate", "bcbioRNADataSet", function(object, ...) {
    .plot_mapping_rate(
        metrics(object),
        interesting_group = .interesting_group(object),
        ...)
})
