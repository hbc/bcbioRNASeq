#' Plot 5'->3' Bias
#'
#' @rdname plot53Bias
#' @name plot53Bias
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @inherit plotTotalReads
#'
#' @examples
#' # bcbioRNASeq
#' bcb <- examples[["bcb"]]
#' plot53Bias(bcb)
#' plot53Bias(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' \dontrun{
#' metrics <- examples[["metrics"]]
#' plot53Bias(metrics)
#' }
NULL



# Constructors ====
#' @importFrom basejump uniteInterestingGroups
#' @importFrom ggplot2 aes_string coord_flip geom_bar ggplot labs
#' @importFrom viridis scale_fill_viridis
.plot53Bias <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 2,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = "x53Bias",
            fill = "interestingGroups")
    ) +
        geom_bar(stat = "identity") +
        labs(title = "5'->3' bias",
             x = "sample",
             y = "5'->3' bias",
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
#' @rdname plot53Bias
#' @importFrom S4Vectors metadata
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod(
    "plot53Bias",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 2,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        .plot53Bias(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip)
    })



#' @rdname plot53Bias
#' @export
setMethod(
    "plot53Bias",
    signature("data.frame"),
    .plot53Bias)
