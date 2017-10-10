#' Plot Gene Detection Saturation
#'
#' @rdname plotGeneSaturation
#' @name plotGeneSaturation
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' data(bcb, rld)
#'
#' # bcbioRNASeq
#' plotGeneSaturation(bcb)
#'
#' \dontrun{
#' plotGeneSaturation(bcb, interestingGroups = "group")
#'
#' # data.frame, matrix
#' plotGeneSaturation(metrics(bcb), assay(rld))
#' }
NULL



# Constructors ====
.plotGeneSaturation <- function(
    object,
    counts,
    interestingGroups = "sampleName",
    minCounts = 0) {
    ggplot(
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
             y = "genes") +
        scale_color_viridis(discrete = TRUE)
}



# Methods ====
#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature(object = "bcbioRNASeqANY",
              counts = "missing"),
    function(
        object,
        interestingGroups,
        normalized = "tmm",
        minCounts = 0) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <- interestingGroups(object)[[1]]
        }
        .plotGeneSaturation(
            metrics(object),
            counts = counts(object, normalized = normalized),
            interestingGroups = interestingGroups(object)[[1]],
            minCounts = minCounts)
    })



#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature(object = "data.frame",
              counts = "matrix"),
    .plotGeneSaturation)
