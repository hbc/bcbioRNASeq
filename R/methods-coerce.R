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
#' dds <- as(bcb_small, "DESeqDataSet")
#' class(dds)
#' show(dds)
#'
#' # RangedSummarizedExperiment ====
#' rse <- as(bcb_small, "RangedSummarizedExperiment")
#' slotNames(rse)
#' show(rse)
#'
#' # SummarizedExperiment ====
#' # Coerced to RangedSummarizedExperiment first.
#' # Otherwise, rowData will be NULL.
#' se <- as(rse, "SummarizedExperiment")
#' class(se)
#' slotNames(se)
#' show(se)
#'
#' # list ====
#' list <- as(bcb_small, "list")
#' class(list)
#' names(list)
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
