#' Plot Genes Detected
#'
#' @rdname plotGenesDetected
#' @name plotGenesDetected
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#' @inheritParams plotGeneSaturation
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotGenesDetected(
#'     bcb,
#'     passLimit = NULL,
#'     warnLimit = NULL)
#' plotGenesDetected(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL,
#'     passLimit = NULL,
#'     warnLimit = NULL)
#'
#' # data.frame, DESeqDataSet
#' df <- metrics(bcb)
#' dds <- bcbio(bcb, "DESeqDataSet")
#' plotGenesDetected(df, counts = dds)
#'
#' # data.frame, matrix
#' counts <- counts(bcb, normalized = TRUE)
#' plotGenesDetected(df, counts = counts)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot labs
#' @importFrom viridis scale_fill_viridis
.plotGenesDetected <- function(
    object,
    counts,
    interestingGroups = "sampleName",
    passLimit = 20000,
    warnLimit = 15000,
    minCounts = 0,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_(
            x = ~sampleName,
            y = colSums(counts > minCounts),
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = "genes detected",
            x = "sample",
            y = "gene count",
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



# Methods ======================================================================
#' @rdname plotGenesDetected
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGenesDetected",
    signature(
        object = "bcbioRNASeq",
        counts = "missing"),
    function(
        object,
        interestingGroups,
        passLimit = 20000,
        warnLimit = 15000,
        minCounts = 0,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotGenesDetected(
            metrics(object),
            counts = assay(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            minCounts = minCounts,
            fill = fill,
            flip = flip)
    })



#' @rdname plotGenesDetected
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGenesDetected",
    signature(
        object = "data.frame",
        counts = "DESeqDataSet"),
    function(
        object,
        counts,
        interestingGroups = "sampleName",
        passLimit = 20000,
        warnLimit = 15000,
        minCounts = 0,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        .plotGenesDetected(
            object,
            counts = assay(counts),
            interestingGroups = "sampleName",
            passLimit = passLimit,
            warnLimit = warnLimit,
            minCounts = minCounts,
            fill = fill,
            flip = flip)
    })



#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(
        object = "data.frame",
        counts = "matrix"),
    .plotGenesDetected)
