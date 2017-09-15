#' Plot Genes Detected
#'
#' @rdname plotGenesDetected
#' @name plotGenesDetected
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @examples
#' data(bcb, dds)
#'
#' # bcbioRNADataSet
#' plotGenesDetected(bcb)
#' plotGenesDetected(bcb, passLimit = NULL, warnLimit = NULL)
#'
#' # data.frame, DESeqDataSet
#' plotGenesDetected(metrics(bcb), dds,
#'                   passLimit = NULL, warnLimit = NULL)
#'
#' # data.frame, matrix
#' plotGenesDetected(metrics(bcb), assay(dds),
#'                   passLimit = NULL, warnLimit = NULL)
NULL



# Constructors ====
.plotGenesDetected <- function(
    object,
    counts,
    interestingGroup = "sampleName",
    passLimit = 20000L,
    warnLimit = 15000L,
    minCounts = 0L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = colSums(counts > minCounts),
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "genes detected",
             x = "sample",
             y = "gene count") +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
    if (!is.null(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(object = "bcbioRNADataSet",
              counts = "missing"),
    function(
        object,
        passLimit = 20000L,
        warnLimit = 15000L,
        minCounts = 0L,
        flip = TRUE) {
        .plotGenesDetected(
            metrics(object),
            counts = assay(object),
            interestingGroup = .interestingGroup(object),
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
        passLimit = 20000L,
        warnLimit = 15000L,
        minCounts = 0L,
        flip = TRUE) {
        .plotGenesDetected(
            object,
            counts = assay(counts),
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
