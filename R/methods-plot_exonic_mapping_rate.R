#' Plot Exonic Mapping Rate
#'
#' @rdname plot_exonic_mapping_rate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#' plot_exonic_mapping_rate(bcb)



#' @rdname plot_exonic_mapping_rate
#' @export
setMethod("plot_exonic_mapping_rate", "bcbioRNADataSet", function(
    object, ...) {
    .plot_exonic_mapping_rate(
        metrics = metrics(object),
        interesting_group = .interesting_group(object),
        ...)
})



#' @rdname plot_exonic_mapping_rate
.plot_exonic_mapping_rate <- function(
    metrics,
    interesting_group,
    pass_limit = 60L,
    flip = TRUE) {
    if (is.null(metrics)) return(NULL)
    p <- ggplot(metrics,
                aes_(x = ~sample_name,
                     y = ~exonic_rate * 100L,
                     fill = as.name(interesting_group))) +
        geom_bar(stat = "identity") +
        labs(title = "exonic mapping rate",
             x = "sample",
             y = "exonic mapping rate (%)",
             fill = "") +
        ylim(0L, 100L)
    if (!is.null(pass_limit)) {
        p <- p +
            geom_hline(alpha = 0.75,
                       color = pass_color,
                       size = 2L,
                       yintercept = pass_limit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}
