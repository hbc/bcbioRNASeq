#' Coerce Object
#'
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @return Object of new class.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#'
#' # DESeqDataSet ====
#' dds <- as(bcb_small, "DESeqDataSet")
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
#' slotNames(se)
#' show(se)
#'
#' # list ====
#' list <- as(bcb_small, "list")
#' names(list)
NULL



# Methods ======================================================================
#' @rdname coerce
#' @name coerce-bcbioRNASeq-DESeqDataSet
#' @importFrom DESeq2 DESeq DESeqDataSetFromTximport
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
        to <- lapply(slotNames(from), function(slot) {
            if (.hasSlot(from, slot)) {
                slot(from, slot)
            } else {
                NULL
            }
        })
        names(to) <- slotNames(from)
        to
    }
)
