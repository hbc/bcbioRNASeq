#' Plot Counts Per Gene
#'
#' @rdname plotCountsPerGene
#' @name plotCountsPerGene
#'
#' @examples
#' data(bcb)
#' plotCountsPerGene(bcb)
NULL



# Constructors ====
.plotCountsPerGene <- function(
    object,
    interestingGroup = "sampleName",
    flip = TRUE) {
    p <- ggplot(object,
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



# Methods ====
#' @rdname plotCountsPerGene
#' @export
setMethod("plotCountsPerGene", "bcbioRNADataSet", function(
    object, normalized = "tmm", ...) {
    .plotCountsPerGene(
        meltLog10(object, normalized = normalized),
        interestingGroup = .interestingGroup(object),
        ...)
})



#' @rdname plotCountsPerGene
#' @export
setMethod("plotCountsPerGene", "data.frame", .plotCountsPerGene)
