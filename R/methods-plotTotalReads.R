#' Plot Total Reads
#'
#' @rdname plotTotalReads
#' @name plotTotalReads
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @importFrom viridis scale_fill_viridis
#'
#' @inheritParams AllGenerics
#'
#' @param counts Object containing a count matrix.
#' @param flip Flip X and Y axes.
#' @param interestingGroups Category to use to group samples (color and shape).
#'   If unset, this is automatically determined by the metadata set inside the
#'   [bcbioRNASeq] object.
#' @param minCounts Numeric value for filtering the counts matrix before
#'   plotting.
#' @param normalized Count normalization method. See [counts()] documentation
#'   for more information.
#' @param passLimit Threshold to plot pass color marker.
#' @param warnLimit Threshold to plot warning color marker.
#'
#' @param fill Desired ggplot color scale. Defaults to
#'   [viridis::scale_color_viridis(discrete = TRUE)]. Must supply discrete
#'   values. When set to `NULL`, the default ggplot2 color palette will be used.
#'   If manual color definitions are desired, we recommend using
#'   [ggplot2::scale_fill_manual()].
#'
#' @return [ggplot].
#'
#' @examples
#' plotTotalReads(bcb)
#'
#' \dontrun{
#' plotTotalReads(
#'     bcb,
#'     interestingGroups = "group",
#'     fill = NULL)
#'
#' # data.frame
#' plotTotalReads(metrics(bcb))
#' }
NULL



# Constructors ====
.plotTotalReads <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 20,
    warnLimit = 10,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = ~totalReads / 1e6,
            fill = as.name(interestingGroups))
    ) +
        geom_bar(stat = "identity") +
        labs(title = "total reads",
             x = "sample",
             y = "total reads (million)")
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
#' @export
setMethod(
    "plotTotalReads",
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
