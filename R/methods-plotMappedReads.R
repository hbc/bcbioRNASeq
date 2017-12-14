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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotMappedReads(bcb)
#' plotMappedReads(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotMappedReads(df)
NULL



# Constructors =================================================================
#' @importFrom basejump uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs
#' @importFrom viridis scale_fill_viridis
.plotMappedReads <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 20,
    warnLimit = 10,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    if (isTRUE(title)) {
        title <- "mapped reads"
    } else if (!is.character(title)) {
        title <- NULL
    }

    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_(
            x = ~sampleName,
            y = ~mappedReads / 1e6,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(title = title,
             x = "sample",
             y = "mapped reads (million)",
             fill = paste(interestingGroups, collapse = ":\n"))

    if (!is.null(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
    if (!is.null(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }

    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    if (interestingGroups == "sampleName") {
        p <- p + guides(fill = FALSE)
    }

    p
}



# Methods ======================================================================
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
        fill = viridis::scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE) {
        if (is.null(metrics(object))) return(NULL)
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        .plotMappedReads(
            metrics(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip,
            title = title)
    })



#' @rdname plotMappedReads
#' @export
setMethod(
    "plotMappedReads",
    signature("data.frame"),
    .plotMappedReads)
