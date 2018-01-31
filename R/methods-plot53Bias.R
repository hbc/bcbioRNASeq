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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plot53Bias(bcb)
#' plot53Bias(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plot53Bias(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom dplyr rename
#' @importFrom ggplot2 aes_string coord_flip geom_bar ggplot guides labs
#' @importFrom viridis scale_fill_viridis
.plot53Bias <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 2L,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    if (isTRUE(title)) {
        title <- "5'->3' bias"
    } else if (!is.character(title)) {
        title <- NULL
    }

    metrics <- uniteInterestingGroups(object, interestingGroups)

    # Legacy code: make sure `x53Bias` is camel sanitized to `x5x3Bias`.
    # The internal camel method has been updated in basejump 0.1.11.
    if ("x53Bias" %in% colnames(metrics)) {
        metrics <- rename(metrics, "x5x3Bias" = "x53Bias")
    }

    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "sampleName",
            y = "x5x3Bias",
            fill = "interestingGroups")
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = "5'->3' bias",
            x = "sample",
            y = "5'->3' bias",
            fill = paste(interestingGroups, collapse = ":\n"))

    if (is.numeric(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }

    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(fill = FALSE)
    }

    p
}



# Methods ======================================================================
#' @rdname plot53Bias
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod(
    "plot53Bias",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 2L,
        fill = viridis::scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE) {
        if (is.null(metrics(object))) return(NULL)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plot53Bias(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip,
            title = title)
    })



#' @rdname plot53Bias
#' @export
setMethod(
    "plot53Bias",
    signature("data.frame"),
    .plot53Bias)
