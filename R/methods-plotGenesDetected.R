#' Plot Genes Detected
#'
#' @rdname plotGenesDetected
#' @name plotGenesDetected
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' plotGenesDetected(bcb, passLimit = NULL, warnLimit = NULL)
#'
#' \dontrun{
#' plotGenesDetected(bcb, interestingGroups = "group")
#' }
#'
#' # data.frame, DESeqDataSet
#' \dontrun{
#' plotGenesDetected(metrics(bcb), counts = dds)
#' }
#'
#' # data.frame, matrix
#' \dontrun{
#' plotGenesDetected(metrics(bcb), counts = assay(dds))
#' }
NULL



# Constructors ====
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
             y = "gene count")
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
#' @rdname plotGenesDetected
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGenesDetected",
    signature(object = "bcbioRNASeq",
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
            fill = fill,
            flip = flip)
    })



#' @rdname plotGenesDetected
#' @importFrom viridis scale_color_viridis
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
    signature(object = "data.frame",
              counts = "matrix"),
    .plotGenesDetected)
