#' Plot Mapped Reads
#'
#' @rdname plotMappedReads
#' @name plotMappedReads
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotMappedReads(bcb)
#'
#' # data.frame
#' metrics <- metrics(bcb)
#' plotMappedReads(metrics)
NULL



# Constructors ====
.plotMappedReads <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 20L,
    warnLimit = 10L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~mappedReads / 1e6L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "mapped reads",
             x = "sample",
             y = "mapped reads (million)")
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
#' @rdname plotMappedReads
#' @export
setMethod("plotMappedReads", "bcbioRNADataSet", function(object, ...) {
    .plotMappedReads(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        ...)
})



#' @rdname plotMappedReads
#' @export
setMethod("plotMappedReads", "data.frame", .plotMappedReads)
