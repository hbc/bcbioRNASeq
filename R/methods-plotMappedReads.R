#' Plot Mapped Reads
#'
#' @rdname plotMappedReads
#' @name plotMappedReads
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' plotMappedReads(bcb)
#'
#' \dontrun{
#' plotMappedReads(
#'     bcb,
#'     interestingGroups = "group",
#'     fill = NULL)
#' }
#'
#' # data.frame
#' \dontrun{
#' metrics(bcb) %>% plotMappedReads()
#' }
NULL



# Constructors ====
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot labs
#' @importFrom viridis scale_fill_viridis
.plotMappedReads <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 20,
    warnLimit = 10,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    interestingGroups <- .checkInterestingGroups(object, interestingGroups)
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = ~mappedReads / 1e6,
            fill = as.name(interestingGroups))
    ) +
        geom_bar(stat = "identity") +
        labs(title = "mapped reads",
             x = "sample",
             y = "mapped reads (million)")
    if (!is.null(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
    if (!is.null(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }
    if (!is.null(fill)) {
        p <- p + fill
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotMappedReads
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotMappedReads",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 20,
        warnLimit = 10,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }
        .plotMappedReads(
            metrics(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip)
    })



#' @rdname plotMappedReads
#' @export
setMethod(
    "plotMappedReads",
    signature("data.frame"),
    .plotMappedReads)
