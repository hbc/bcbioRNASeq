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
#' names(S4Vectors::mcols(x))
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
        if (metadata(from)[["level"]] != "genes") {
            stop("Gene-level counts are required")
        }
        # Creating `DESeqDataSet` from `RangedSummarizedExperiment` is
        # preferable to `DESeqDataSetFromTximport` method because `rowRanges`
        # are defined, with richer metadata
        rse <- as(from, "RangedSummarizedExperiment")
        # Integer counts are required
        assay(rse) <- round(assay(rse))
        # Prepare using an empty design formula
        dds <- DESeqDataSet(
            se = rse,
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
