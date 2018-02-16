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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotCountsPerGene(bcb)
#' plotCountsPerGene(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- meltLog10(bcb, normalized = "tmm")
#' plotCountsPerGene(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_string geom_boxplot ggplot guides labs
.plotCountsPerGene <- function(
    object,
    interestingGroups = "sampleName",
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    if (isTRUE(title)) {
        title <- "counts per gene"
    } else if (!is.character(title)) {
        title <- NULL
    }

    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = "counts",
            fill = "interestingGroups")
    ) +
        geom_boxplot(color = lineColor, outlier.shape = NA) +
        labs(
            title = title,
            x = "sample",
            y = "log10 counts per gene",
            fill = paste(interestingGroups, collapse = ":\n"))

    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(fill = FALSE)
    }

    p
}



# Methods ======================================================================
#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = "tmm",
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotCountsPerGene(
            meltLog10(object, normalized = normalized),
            interestingGroups = interestingGroups,
            fill = fill,
            flip = flip,
            title = title)
    })



#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("data.frame"),
    .plotCountsPerGene)
