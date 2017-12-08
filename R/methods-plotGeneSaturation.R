#' Plot Gene Detection Saturation
#'
#' @rdname plotGeneSaturation
#' @name plotGeneSaturation
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#' @inheritParams plotCountDensity
#' @inheritParams plotGene
#'
#' @param counts Object containing a count matrix.
#' @param minCounts Numeric value for filtering the counts matrix before
#'   plotting.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotGeneSaturation(bcb)
#' plotGeneSaturation(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     color = NULL)
#'
#' # data.frame, matrix
#' df <- metrics(bcb)
#' counts <- counts(bcb, normalized = "rlog")
#' plotGeneSaturation(df, counts = counts)
NULL



# Constructors ====
#' @importFrom basejump uniteInterestingGroups
#' @importFrom ggplot2 aes_ geom_point geom_smooth ggplot labs
#' @importFrom viridis scale_color_viridis
.plotGeneSaturation <- function(
    object,
    counts,
    interestingGroups = "sampleName",
    minCounts = 0,
    color = viridis::scale_color_viridis(discrete = TRUE)) {
    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_(
            x = ~mappedReads / 1e6,
            y = colSums(counts > minCounts),
            color = ~interestingGroups)
    ) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(title = "gene saturation",
             x = "mapped reads (million)",
             y = "genes",
             color = paste(interestingGroups, collapse = ":\n"))
    if (!is.null(color)) {
        p <- p + color
    }
    p
}



# Methods ====
#' @rdname plotGeneSaturation
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGeneSaturation",
    signature(object = "bcbioRNASeq",
              counts = "missing"),
    function(
        object,
        interestingGroups,
        normalized = "tmm",
        minCounts = 0,
        color = viridis::scale_color_viridis(discrete = TRUE)) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        .plotGeneSaturation(
            metrics(object),
            counts = counts(object, normalized = normalized),
            interestingGroups = interestingGroups,
            minCounts = minCounts,
            color = color)
    })



#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature(object = "data.frame",
              counts = "matrix"),
    .plotGeneSaturation)
