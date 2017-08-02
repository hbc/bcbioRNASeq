#' Plot Ribosomal RNA (rRNA) Mapping Rate
#'
#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotRRNAMappingRate(bcb)
#'
#' # data.frame
#' metrics <- metrics(bcb)
#' plotRRNAMappingRate(metrics)
NULL



# Constructors ====
.plotRRNAMappingRate <- function(
    object,
    interestingGroup = "sampleName",
    warnLimit = 10L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~rrnaRate * 100L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        geom_hline(alpha = qcLineAlpha,
                   color = qcWarnColor,
                   size = qcLineSize,
                   yintercept = warnLimit) +
        labs(title = "rrna mapping rate",
             x = "sample",
             y = "rRNA mapping rate (%)")
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotRRNAMappingRate
#' @export
setMethod("plotRRNAMappingRate", "bcbioRNADataSet", function(object, ...) {
    .plotRRNAMappingRate(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        ...)
})



#' @rdname plotRRNAMappingRate
#' @export
setMethod("plotRRNAMappingRate", "data.frame", .plotRRNAMappingRate)
