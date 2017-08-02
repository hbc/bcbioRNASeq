#' Plot Count Density
#'
#' @rdname plotCountDensity
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qcPlots
#'
#' @examples
#' data(bcb)
#' plotCountDensity(bcb)



#' @rdname plotCountDensity
.plotCountDensity <- function(
    melted,
    interestingGroup = "sampleName") {
    ggplot(melted,
        aes_(x = ~counts,
             group = ~sampleName,
             color = ~sampleName)) +
        geom_density() +
        labs(title = "count density",
             x = "log10 counts per gene")
}



#' @rdname plotCountDensity
#' @export
setMethod("plotCountDensity", "bcbioRNADataSet", function(
    object, normalized = "tmm") {
    .plotCountDensity(
        melted = meltLog10(object, normalized = normalized),
        interestingGroup = .interestingGroup(object))
})
