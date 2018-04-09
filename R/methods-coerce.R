#' Methods for Coercing an Object to a Class
#'
#' @name coerce
#' @aliases as
#' @family S4 Class Definition
#' @author Michael Steinbaugh
#'
#' @return Object of new class.
#'
#' @seealso
#' - [methods::as()].
#' - [methods::coerce()].
#'
#' @examples
#' # DESeqDataSet ====
#' x <- as(bcb_small, "DESeqDataSet")
#' class(x)
#' show(x)
#'
#' # RangedSummarizedExperiment ====
#' x <- as(bcb_small, "RangedSummarizedExperiment")
#' slotNames(x)
#' show(x)
#'
#' # SummarizedExperiment ====
#' # Coerce to RangedSummarizedExperiment first.
#' # Otherwise, rowData will be NULL.
#' x <- as(bcb_small, "RangedSummarizedExperiment")
#' x <- as(x, "SummarizedExperiment")
#' class(x)
#' slotNames(x)
#' show(x)
#'
#' # list ====
#' x <- as(bcb_small, "list")
#' class(x)
#' names(x)
NULL



# Methods ======================================================================
#' @rdname coerce
#' @name coerce-bcbioRNASeq-DESeqDataSet
setAs(
    from = "bcbioRNASeq",
    to = "DESeqDataSet",
    function(from) {
        validObject(from)
        txi <- .regenerateTximportList(from)
        dds <- DESeqDataSetFromTximport(
            txi = txi,
            colData = colData(from),
            # Use an empty design formula
            design = ~ 1  # nolint
        )
        # Suppress warning about empty design formula
        to <- suppressWarnings(DESeq(dds))
        validObject(to)
        to
    }
)



#' @rdname coerce
#' @name coerce-bcbioRNASeq-list
setAs(
    from = "bcbioRNASeq",
    to = "list",
    function(from) {
        flatFiles(from)
    }
)



#' @rdname coerce
#' @name coerce-bcbioRNASeq-SummarizedExperiment
setAs(
    from = "bcbioRNASeq",
    to = "SummarizedExperiment",
    function(from) {
        # Otherwise rowData will be NULL
        rse <- as(from, "RangedSummarizedExperiment")
        se <- as(rse, "SummarizedExperiment")
        se
    }
)
