#' Plot Count Density
#'
#' @rdname plot_count_density
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qc_plots
#'
#' @examples
#' data(bcb)
#' plot_count_density(bcb)



#' @rdname plot_count_density
.plot_count_density <- function(
    melted,
    interesting_group = "sample_name") {
    ggplot(melted,
        aes_(x = ~counts,
             group = ~sample_name)) +
        geom_density() +
        labs(title = "count density",
             x = "log10 counts per gene")
}



#' @rdname plot_count_density
#' @export
setMethod("plot_count_density", "bcbioRNADataSet", function(
    object, normalized = "tmm") {
    .plot_count_density(
        melted = melt_log10(object, normalized = normalized),
        interesting_group = .interesting_group(object))
})
