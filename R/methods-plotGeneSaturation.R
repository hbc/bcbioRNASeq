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
#' plotGeneSaturation(bcb)
#'
#' \dontrun{
#' plotGeneSaturation(bcb, interestingGroups = "group")
#' }
#'
#' # data.frame, matrix
#' \dontrun{
#' plotGeneSaturation(metrics(bcb), assay(rld))
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_color_viridis
.plotGeneSaturation <- function(
    object,
    counts,
    interestingGroups = "sampleName",
    minCounts = 0,
    color = scale_color_viridis(discrete = TRUE)) {
    interestingGroups <- .checkInterestingGroups(object, interestingGroups)
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~mappedReads / 1e6,
            y = colSums(counts > minCounts),
            color = as.name(interestingGroups))
    ) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(title = "gene saturation",
             x = "mapped reads (million)",
             y = "genes")
    if (!is.null(color)) {
        p <- p + color
    }
    p
}



# Methods ====
#' @rdname plotGeneSaturation
#' @importFrom S4Vectors metadata
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
        color = scale_color_viridis(discrete = TRUE)) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
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
