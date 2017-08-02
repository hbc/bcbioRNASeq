#' Plot Counts Per Gene
#'
#' @rdname plotCountsPerGene
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qcPlots
#'
#' @examples
#' data(bcb)
#' plotCountsPerGene(bcb)



#' @rdname plotCountsPerGene
.plotCountsPerGene <- function(
    melted,
    interestingGroup = "sampleName",
    flip = TRUE) {
    p <- ggplot(melted,
                aes_(x = ~sampleName,
                     y = ~counts,
                     color = as.name(interestingGroup))) +
        geom_boxplot(outlier.shape = NA) +
        labs(title = "counts per gene",
             x = "sample",
             y = "log10 counts per gene")
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plotCountsPerGene
#' @export
setMethod("plotCountsPerGene", "bcbioRNADataSet", function(
    object, normalized = "tmm", ...) {
    .plotCountsPerGene(
        melted = meltLog10(object, normalized = normalized),
        interestingGroup = .interestingGroup(object),
        ...)
})
