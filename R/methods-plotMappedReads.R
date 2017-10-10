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
#' plotMappedReads(bcb, interestingGroup = "group")
#'
#' \dontrun{
#' # data.frame
#' metrics(bcb) %>%
#'     plotMappedReads()
#' }
NULL



# Constructors ====
.plotMappedReads <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 20,
    warnLimit = 10,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = ~mappedReads / 1e6,
            fill = as.name(interestingGroup))
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
        .plotMappedReads(
            metrics(object),
            interestingGroup = interestingGroup,
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
