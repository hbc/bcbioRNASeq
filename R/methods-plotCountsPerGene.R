#' Plot Counts Per Gene
#'
#' @rdname plotCountsPerGene
#' @name plotCountsPerGene
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotCountsPerGene(bcb)
#'
#' # data.frame
#' meltLog10(bcb) %>% plotCountsPerGene
NULL



# Constructors ====
.plotCountsPerGene <- function(
    object,
    interestingGroup = "sampleName",
    flip = TRUE) {
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~counts,
                     fill = as.name(interestingGroup))) +
        geom_boxplot(color = lineColor, outlier.shape = NA) +
        labs(title = "counts per gene",
             x = "sample",
             y = "log10 counts per gene") +
        scale_fill_viridis(discrete = TRUE)
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotCountsPerGene
#' @export
setMethod("plotCountsPerGene", "bcbioRNADataSet", function(
    object,
    normalized = "tmm",
    flip = TRUE) {
    .plotCountsPerGene(
        meltLog10(object, normalized = normalized),
        interestingGroup = .interestingGroup(object),
        flip = flip)
})



#' @rdname plotCountsPerGene
#' @export
setMethod("plotCountsPerGene", "data.frame", .plotCountsPerGene)
