#' Coerce Object
#'
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @return Object of new class.
#'
#' @examples
#' # DESeqDataSet ====
#' dds <- as(bcb_small, "DESeqDataSet")
#' show(dds)
#'
#' # RangedSummarizedExperiment ====
#' rse <- as(bcb_small, "RangedSummarizedExperiment")
#' show(rse)
#'
#' # SummarizedExperiment ====
#' se <- as(bcb_small, "SummarizedExperiment")
#' show(se)
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
        # Regenerate tximport list
        txi <- list(
            "abundance" = assays(from)[["tpm"]],
            "counts" = assays(from)[["raw"]],
            "length" = assays(from)[["length"]],
            "countsFromAbundance" = metadata(from)[["countsFromAbundance"]]
        )
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
