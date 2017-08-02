#' Plot Exonic Mapping Rate
#'
#' @rdname plotExonicMappingRate
#' @name plotExonicMappingRate
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotExonicMappingRate(bcb)
#'
#' # data.frame
#' metrics <- metrics(bcb)
#' plotExonicMappingRate(metrics)
NULL



# Constructors ====
.plotExonicMappingRate <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 60L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~exonicRate * 100L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "exonic mapping rate",
             x = "sample",
             y = "exonic mapping rate (%)") +
        ylim(0L, 100L)
    if (!is.null(passLimit)) {
        p <- p +
            geom_hline(alpha = qcLineAlpha,
                       color = qcPassColor,
                       size = qcLineSize,
                       yintercept = passLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotExonicMappingRate
#' @export
setMethod("plotExonicMappingRate", "bcbioRNADataSet", function(object, ...) {
    .plotExonicMappingRate(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        ...)
})



#' @rdname plotExonicMappingRate
#' @export
setMethod("plotExonicMappingRate", "data.frame", .plotExonicMappingRate)
