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
#' plotGeneSaturation(bcb, interestingGroup = "group")
#'
#' \dontrun{
#' # data.frame, matrix
#' plotGeneSaturation(metrics(bcb), assay(rld))
#' }
NULL



# Constructors ====
.plotGeneSaturation <- function(
    object,
    counts,
    interestingGroup = "sampleName",
    minCounts = 0L) {
    ggplot(object,
           aes_(x = ~mappedReads / 1e6L,
                y = colSums(counts > minCounts),
                color = as.name(interestingGroup))) +
        geom_point(size = 3L) +
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
        interestingGroup,
        normalized = "tmm",
        minCounts = 0L) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroup)) {
            interestingGroup <- .interestingGroup(object)
        }
        .plotGeneSaturation(
            metrics(object),
            counts = counts(object, normalized = normalized),
            interestingGroup = .interestingGroup(object),
            minCounts = minCounts)
    })



#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature(object = "data.frame",
              counts = "matrix"),
    .plotGeneSaturation)
