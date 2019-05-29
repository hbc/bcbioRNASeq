#' Methods for coercing an object to a class
#'
#' @name coerce
#' @aliases as
#' @family S4 Object
#' @author Michael Steinbaugh
#' @importFrom methods coerce
#' @exportMethod coerce
#'
#' @return Object of new class.
#'
#' @seealso
#' - [methods::as()].
#' - [methods::canCoerce()].
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
#' x <- as(bcb_small, "RangedSummarizedExperiment")
#' x <- as(x, "SummarizedExperiment")
#' class(x)
#' slotNames(x)
#' show(x)
NULL



#' @rdname coerce
#' @name coerce,bcbioRNASeq,DESeqDataSet-method
setAs(
    from = "bcbioRNASeq",
    to = "DESeqDataSet",
    function(from) {
        validObject(from)
        if (metadata(from)[["level"]] != "genes") {
            stop("Gene-level counts are required")  # nocov
        }
        message(paste(
            "Coercing bcbioRNASeq to DESeqDataSet with DESeq2",
            packageVersion("DESeq2")
        ))
        # Creating `DESeqDataSet` from `RangedSummarizedExperiment` is
        # preferable to `DESeqDataSetFromTximport` method because `rowRanges`
        # are defined, with richer metadata
        rse <- as(from, "RangedSummarizedExperiment")
        # Integer counts are required
        assay(rse) <- round(assay(rse), digits = 0L)
        # Prepare using an empty design formula
        to <- DESeqDataSet(se = rse, design = ~ 1L)
        interestingGroups(to) <- interestingGroups(from)
        validObject(to)
        to
    }
)
