#' Plot Total Reads
#'
#' @rdname plotTotalReads
#' @name plotTotalReads
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams AllGenerics
#'
#' @param interestingGroups Category to use to group samples. In the plotting
#'   functions, this will define color and shape, where applicable. If unset,
#'   this is automatically determined by the metadata set inside the
#'   [bcbioRNASeq] object. When set to `NULL`, this will default to
#'   `sampleName`.
#' @param passLimit Threshold to plot pass color marker.
#' @param warnLimit Threshold to plot warning color marker.
#' @param fill Desired ggplot fill scale. Defaults to
#'   [viridis::scale_fill_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [ggplot2::scale_fill_manual()].
#' @param flip Flip x and y axes.
#'
#' @return [ggplot].
#'
#' @examples
#' # bcbioRNASeq
#' bcb <- examples[["bcb"]]
#' plotTotalReads(bcb)
#' plotTotalReads(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' metrics <- examples[["metrics"]]
#' plotTotalReads(metrics)
NULL



# Constructors ====
#' @importFrom basejump uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot labs
#' @importFrom rlang !!! syms
#' @importFrom tidyr unite
#' @importFrom viridis scale_fill_viridis
.plotTotalReads <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 20,
    warnLimit = 10,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_(
            x = ~sampleName,
            y = ~totalReads / 1e6,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(title = "total reads",
             x = "sample",
             y = "total reads (million)",
             fill = paste(interestingGroups, collapse = ":\n"))
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
#' @rdname plotTotalReads
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod(
    "plotTotalReads",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 20,
        warnLimit = 10,
        fill = viridis::scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        .plotTotalReads(
            metrics(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip)
    })



#' @rdname plotTotalReads
#' @export
setMethod(
    "plotTotalReads",
    signature("data.frame"),
    .plotTotalReads)
