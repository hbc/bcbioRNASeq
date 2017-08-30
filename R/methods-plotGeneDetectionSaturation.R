#' Plot Gene Detection Saturation
#'
#' @rdname plotGeneSaturation
#' @name plotGeneSaturation
#'
#' @examples
#' data(bcb, rld)
#'
#' # bcbioRNADataSet
#' plotGeneSaturation(bcb)
#'
#' # data.frame, matrix
#' plotGeneSaturation(metrics(bcb), assay(rld))
NULL



# Constructors ====
.plotGeneSaturation <- function(
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
    signature(object = "bcbioRNADataSet",
              counts = "missing"),
    function(
        object, normalized = "tmm", ...) {
        .plotGeneSaturation(
            metrics(object),
            counts = counts(object, normalized = normalized),
            interestingGroup = .interestingGroup(object),
            ...)
    })



#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature(object = "data.frame",
              counts = "matrix"),
    .plotGeneSaturation)
