#' Plot rRNA Mapping Rate
#'
#' @rdname plot_rrna_mapping_rate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plot_rrna_mapping_rate(bcb)
#'
#' # data.frame
#' metrics <- metrics(bcb)
#' plot_rrna_mapping_rate(metrics)



#' @rdname plot_rrna_mapping_rate
.plot_rrna_mapping_rate <- function(
    metrics,
    interesting_group = "sample_name",
    warn_limit = 10L,
    flip = TRUE) {
    p <- ggplot(metrics,
                aes_(x = ~sample_name,
                     y = ~r_rna_rate * 100L,
                     fill = as.name(interesting_group))) +
        geom_bar(stat = "identity") +
        geom_hline(alpha = qc_line_alpha,
                   color = qc_warn_color,
                   size = qc_line_size,
                   yintercept = warn_limit) +
        labs(title = "rrna mapping rate",
             x = "sample",
             y = "rRNA mapping rate (%)")
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plot_rrna_mapping_rate
#' @export
setMethod("plot_rrna_mapping_rate", "bcbioRNADataSet", function(object, ...) {
    .plot_rrna_mapping_rate(
        metrics(object),
        interesting_group = .interesting_group(object),
        ...)
})



#' @rdname plot_rrna_mapping_rate
#' @export
setMethod("plot_rrna_mapping_rate", "data.frame", function(object, ...) {
    .plot_rrna_mapping_rate(object, ...)
})
