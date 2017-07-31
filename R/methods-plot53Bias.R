#' Plot 5'->3' Bias
#'
#' @rdname plot_5_3_bias
#' @author Michael Steinbaugh
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plot_5_3_bias(bcb)
#'
#' # data.frame
#' metrics <- metrics(bcb)
#' plot_5_3_bias(metrics)



#' @rdname plot_5_3_bias
.plot_5_3_bias <- function(
    metrics,
    interesting_group = "sample_name",
    warn_limit = 2L,
    flip = TRUE) {
    if (is.null(metrics)) return(NULL)
    p <- ggplot(metrics,
                aes_(x = ~sample_name,
                     y = ~x5_3_bias,
                     fill = as.name(interesting_group))) +
        geom_bar(stat = "identity") +
        labs(title = "5'->3' bias",
             x = "sample",
             y = "5'->3' bias")
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



#' @rdname plot_5_3_bias
#' @export
setMethod("plot_5_3_bias", "bcbioRNADataSet", function(object, ...) {
    .plot_5_3_bias(
        metrics(object),
        interesting_group = .interesting_group(object),
        ...)
})



#' @rdname plot_5_3_bias
#' @export
setMethod("plot_5_3_bias", "data.frame", function(object, ...) {
    .plot_5_3_bias(object, ...)
})
