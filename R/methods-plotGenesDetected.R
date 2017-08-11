#' Plot Genes Detected
#'
#' @rdname plotGenesDetected
#' @name plotGenesDetected
#'
#' @examples
#' data(bcb, dds)
#'
#' # bcbioRNADataSet
#' plotGenesDetected(bcb, passLimit = NULL, warnLimit = NULL)
#'
#' # data.frame + DESeqDataSet
#' plotGenesDetected(metrics(bcb), dds,
#'                   passLimit = NULL, warnLimit = NULL)
#'
#' # data.frame + matrix
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
             y = "gene count")
    if (!is.null(passLimit)) {
        p <- p + geom_hline(alpha = qcLineAlpha,
                            color = qcPassColor,
                            size = qcLineSize,
                            yintercept = passLimit)
    }
    if (!is.null(warnLimit)) {
        p <- p + geom_hline(alpha = qcLineAlpha,
                            color = qcWarnColor,
                            size = qcLineSize,
                            yintercept = warnLimit)
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
    function(object, ...) {
        .plotGenesDetected(
            metrics(object),
            counts = assay(object),
            interestingGroup = .interestingGroup(object),
            ...)
    })



#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(object = "data.frame",
              counts = "DESeqDataSet"),
    function(object, counts, ...) {
        .plotGenesDetected(
            object,
            counts = assay(counts),
            ...)
    })



#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(object = "data.frame",
              counts = "matrix"),
    .plotGenesDetected)
