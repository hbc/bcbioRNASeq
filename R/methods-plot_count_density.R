#' Plot Count Density
#'
#' @rdname plot_count_density
#' @family Quality Control Metrics Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit plot_metrics
#'
#' @examples
#' data(bcb)
#' plot_count_density(bcb)
.plot_count_density <- function(melted, interesting_group) {
    ggplot(melted,
        aes_(x = ~counts,
             group = ~sample_name)) +
        geom_density() +
        labs(title = "count density",
             x = "log10 counts per gene",
             y = "density")
}



#' @rdname plot_count_density
#' @export
setMethod("plot_count_density", "bcbioRNADataSet", function(
    object, normalized = "tmm") {
    .plot_count_density(
        melted = melt_log10(object, normalized = normalized),
        interesting_group = .interesting_group(object))
})
