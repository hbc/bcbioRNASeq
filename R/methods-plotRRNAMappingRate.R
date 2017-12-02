#' Plot Ribosomal RNA (rRNA) Mapping Rate
#'
#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' # bcbioRNASeq
#' bcb <- examples[["bcb"]]
#' plotRRNAMappingRate(bcb)
#' plotRRNAMappingRate(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL,
#'     warnLimit = NULL)
#'
#' # data.frame
#' \dontrun{
#' metrics <- examples[["metrics"]]
#' plotRRNAMappingRate(metrics)
#' }
NULL



# Constructors ====
#' @importFrom basejump uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot labs
#' @importFrom viridis scale_fill_viridis
.plotRRNAMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 10,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    # Fix for camel variant mismatch (e.g. rRnaRate). This is safe to remove
    # in a future update.
    if (!"rrnaRate" %in% colnames(object)) {
        # grep match the outdated camel variant
        col <- grep(
            x = colnames(object),
            pattern = "rrnarate",
            ignore.case = TRUE,
            value = TRUE)
        object[["rrnaRate"]] <- object[[col]]
        object[[col]] <- NULL
    }
    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_(
            x = ~sampleName,
            y = ~rrnaRate * 100,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(title = "rRNA mapping rate",
             x = "sample",
             y = "rRNA mapping rate (%)",
             fill = paste(interestingGroups, collapse = ":\n"))
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
#' @rdname plotRRNAMappingRate
#' @importFrom S4Vectors metadata
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 10,
        fill = viridis::scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        .plotRRNAMappingRate(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip)
    })



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("data.frame"),
    .plotRRNAMappingRate)
