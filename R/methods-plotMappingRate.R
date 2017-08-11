#' Plot Mapping Rate
#'
#' @rdname plotMappingRate
#' @name plotMappingRate
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotMappingRate(bcb)
#'
#' # data.frame
#' metrics(bcb) %>% plotMappingRate
NULL



# Constructors ====
.plotMappingRate <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 90L,
    warnLimit = 70L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~mappedReads / totalReads * 100L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        ylim(0L, 100L) +
        labs(title = "mapping rate",
             x = "sample",
             y = "mapping rate (%)")
    if (!is.null(passLimit)) {
        p <- p +
            geom_hline(alpha = qcLineAlpha,
                       color = qcPassColor,
                       size = qcLineSize,
                       yintercept = passLimit)
    }
    if (!is.null(warnLimit)) {
        p <- p +
            geom_hline(alpha = qcLineAlpha,
                       color = qcWarnColor,
                       size = qcLineSize,
                       yintercept = warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotMappingRate
#' @export
setMethod("plotMappingRate", "bcbioRNADataSet", function(object, ...) {
    .plotMappingRate(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        ...)
})



#' @rdname plotMappingRate
#' @export
setMethod("plotMappingRate", "data.frame", .plotMappingRate)
