#' Plot Counts Per Gene
#'
#' @rdname plotCountsPerGene
#' @name plotCountsPerGene
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNASeq
#' plotCountsPerGene(bcb)
#'
#' \dontrun{
#' plotCountsPerGene(bcb, interestingGroups = "group")
#'
#' # data.frame
#' meltLog10(bcb, normalized = "tmm") %>%
#'     plotCountsPerGene()
#' }
NULL



# Constructors ====
.plotCountsPerGene <- function(
    object,
    interestingGroups = "sampleName",
    flip = TRUE) {
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
             y = "log10 counts per gene") +
        scale_fill_viridis(discrete = TRUE)
    if (isTRUE(flip)) {
        p <- p +
            coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("bcbioRNASeqANY"),
    function(
        object,
        interestingGroups,
        normalized = "tmm",
        flip = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- interestingGroups(object)[[1]]
        }
        .plotCountsPerGene(
            meltLog10(object, normalized = normalized),
            interestingGroups = interestingGroups,
            flip = flip)
    })



#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("data.frame"),
    .plotCountsPerGene)
