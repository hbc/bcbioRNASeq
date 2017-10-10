#' Plot Genes Detected
#'
#' @rdname plotGenesDetected
#' @name plotGenesDetected
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' data(bcb, dds)
#'
#' # bcbioRNASeq
#' plotGenesDetected(
#'     bcb,
#'     passLimit = NULL,
#'     warnLimit = NULL)
#'
#' \dontrun{
#' plotGenesDetected(bcb, interestingGroups = "group")
#'
#' # data.frame, DESeqDataSet
#' plotGenesDetected(
#'     metrics(bcb),
#'     counts = dds)
#'
#' # data.frame, matrix
#' plotGenesDetected(
#'     metrics(bcb),
#'     counts = assay(dds))
#' }
NULL



# Constructors ====
.plotGenesDetected <- function(
    object,
    counts,
    interestingGroups = "sampleName",
    passLimit = 20000,
    warnLimit = 15000,
    minCounts = 0,
    flip = TRUE) {
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = colSums(counts > minCounts),
            fill = as.name(interestingGroups))
    ) +
        geom_bar(stat = "identity") +
        labs(title = "genes detected",
             x = "sample",
             y = "gene count") +
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
#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(object = "bcbioRNASeqANY",
              counts = "missing"),
    function(
        object,
        interestingGroups,
        passLimit = 20000,
        warnLimit = 15000,
        minCounts = 0,
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }
        .plotGenesDetected(
            metrics(object),
            counts = assay(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            minCounts = minCounts,
            flip = flip)
    })



#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(object = "data.frame",
              counts = "DESeqDataSet"),
    function(
        object,
        counts,
        interestingGroups = "sampleName",
        passLimit = 20000,
        warnLimit = 15000,
        minCounts = 0,
        flip = TRUE) {
        .plotGenesDetected(
            object,
            counts = assay(counts),
            interestingGroups = "sampleName",
            passLimit = passLimit,
            warnLimit = warnLimit,
            minCounts = minCounts,
            flip = flip)
    })



#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(object = "data.frame",
              counts = "matrix"),
    .plotGenesDetected)
