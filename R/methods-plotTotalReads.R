#' Plot Total Reads
#'
#' @rdname plotTotalReads
#' @name plotTotalReads
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotTotalReads(bcb)
#'
#' # data.frame
#' metrics(bcb) %>% plotTotalReads
NULL



# Constructors ====
.plotTotalReads <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 20L,
    warnLimit = 10L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~totalReads / 1e6L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "total reads",
             x = "sample",
             y = "total reads (million)") +
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
#' @rdname plotTotalReads
#' @export
setMethod("plotTotalReads", "bcbioRNADataSet", function(
    object,
    passLimit = 20L,
    warnLimit = 10L,
    flip = TRUE) {
    .plotTotalReads(
        metrics(object),
        passLimit = passLimit,
        warnLimit = warnLimit,
        flip = flip,
        # Automatic
        interestingGroup = .interestingGroup(object))
})



#' @rdname plotTotalReads
#' @export
setMethod("plotTotalReads", "data.frame", .plotTotalReads)
