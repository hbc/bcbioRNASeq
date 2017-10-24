#' Plot Counts Per Gene
#'
#' @rdname plotCountsPerGene
#' @name plotCountsPerGene
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#' @inheritParams plotCountDensity
#'
#' @examples
#' plotCountsPerGene(bcb)
#'
#' \dontrun{
#' plotCountsPerGene(
#'     bcb,
#'     interestingGroups = "group",
#'     fill = NULL)
#' }
#'
#' # data.frame
#' \dontrun{
#' meltLog10(bcb, normalized = "tmm") %>% plotCountsPerGene()
#' }
NULL



# Constructors ====
#' @importFrom ggplot2 aes_string geom_boxplot ggplot labs
#' @importFrom viridis scale_fill_viridis
.plotCountsPerGene <- function(
    object,
    interestingGroups = "sampleName",
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    .checkInterestingGroups(object, interestingGroups)
    p <- ggplot(
        object,
        mapping = aes_string(
            x = "sampleName",
            y = "counts",
            fill = interestingGroups)
    ) +
        geom_boxplot(color = lineColor, outlier.shape = NA) +
        labs(title = "counts per gene",
             x = "sample",
             y = "log10 counts per gene")
    if (!is.null(fill)) {
        p <- p + fill
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotCountsPerGene
#' @importFrom S4Vectors metadata
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod(
    "plotCountsPerGene",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = "tmm",
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }
        .plotCountsPerGene(
            meltLog10(object, normalized = normalized),
            interestingGroups = interestingGroups,
            fill = fill,
            flip = flip)
    })



#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("data.frame"),
    .plotCountsPerGene)
