#' Plot Gene Detection Saturation
#'
#' @rdname plotGeneDetectionSaturation
#' @name plotGeneDetectionSaturation
#'
#' @examples
#' data(bcb)
#' plotGeneDetectionSaturation(bcb)
NULL



# Constructors ====
.plotGeneDetectionSaturation <- function(
    object,
    counts,
    interestingGroup = "sampleName",
    minCounts = 0L) {
    if (is.null(object)) return(NULL)
    ggplot(object,
           aes_(x = ~mappedReads / 1e6L,
                y = colSums(counts > minCounts),
                color = as.name(interestingGroup))) +
        geom_point(size = 3L) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(title = "gene detection saturation",
             x = "mapped reads (million)",
             y = "gene count")
}



# Methods ====
#' @rdname plotGeneDetectionSaturation
#' @export
setMethod("plotGeneDetectionSaturation", "bcbioRNADataSet", function(
    object, normalized = "tmm", ...) {
    .plotGeneDetectionSaturation(
        metrics(object),
        counts = counts(object, normalized = normalized),
        interestingGroup = .interestingGroup(object),
        ...)
})



#' @rdname plotGeneDetectionSaturation
#' @export
setMethod("plotGeneDetectionSaturation",
          "data.frame",
          .plotGeneDetectionSaturation)
