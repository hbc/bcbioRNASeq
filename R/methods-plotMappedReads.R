#' Plot Mapped Reads
#'
#' @rdname plotMappedReads
#' @name plotMappedReads
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
#' plotMappedReads(bcb)
#'
#' \dontrun{
#' plotMappedReads(bcb, interestingGroups = "group")
#'
#' # data.frame
#' metrics(bcb) %>%
#'     plotMappedReads()
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotMappedReads <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 20,
    warnLimit = 10,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
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
             y = "mapped reads (million)") +
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
#' @rdname plotMappedReads
#' @export
setMethod(
    "plotMappedReads",
    signature("bcbioRNASeqANY"),
    function(
        object,
        interestingGroups,
        passLimit = 20,
        warnLimit = 10,
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
            flip = flip)
    })



#' @rdname plotMappedReads
#' @export
setMethod(
    "plotMappedReads",
    signature("data.frame"),
    .plotMappedReads)
