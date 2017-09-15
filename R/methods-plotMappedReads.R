#' Plot Mapped Reads
#'
#' @rdname plotMappedReads
#' @name plotMappedReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotMappedReads(bcb)
#'
#' # data.frame
#' metrics(bcb) %>% plotMappedReads
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
             y = "mapped reads (million)") +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
    if (!is.null(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotMappedReads
#' @export
setMethod("plotMappedReads", "bcbioRNADataSet", function(
    object,
    passLimit = 20L,
    warnLimit = 10L,
    flip = TRUE) {
    .plotMappedReads(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        passLimit = passLimit,
        warnLimit = warnLimit,
        flip = flip)
})



#' @rdname plotMappedReads
#' @export
setMethod("plotMappedReads", "data.frame", .plotMappedReads)
