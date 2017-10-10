#' Plot Total Reads
#'
#' @rdname plotTotalReads
#' @name plotTotalReads
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNASeq
#' plotTotalReads(bcb)
#'
#' # data.frame
#' metrics(bcb) %>%
#'     plotTotalReads()
NULL



# Constructors ====
.plotTotalReads <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 20,
    warnLimit = 10,
    flip = TRUE) {
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = ~totalReads / 1e6,
            fill = as.name(interestingGroup))
    ) +
        geom_bar(stat = "identity") +
        labs(title = "total reads",
             x = "sample",
             y = "total reads (million)") +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(passLimit)) {
        p <- p +
            qcPassLine(passLimit)
    }
    if (!is.null(warnLimit)) {
        p <- p +
            qcWarnLine(warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p +
            coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotTotalReads
#' @export
setMethod("plotTotalReads", "bcbioRNASeqANY", function(
    object,
    interestingGroup,
    passLimit = 20,
    warnLimit = 10,
    flip = TRUE) {
    if (is.null(metrics(object))) {
        return(NULL)
    }
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    .plotTotalReads(
        metrics(object),
        interestingGroup = interestingGroup,
        passLimit = passLimit,
        warnLimit = warnLimit,
        flip = flip)
})



#' @rdname plotTotalReads
#' @export
setMethod("plotTotalReads", "data.frame", .plotTotalReads)
