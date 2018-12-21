#' @name as
#' @aliases coerce
#' @author Michael Steinbaugh
#' @exportMethod coerce
#' @importFrom methods coerce
#' @inherit methods::as
#' @inheritParams coerce
#'
#' @section bcbioRNASeq to DESeqDataSet:
#'
#' 1. Coerces to `RangedSummarizedExperiment`.
#' 2. Rounds raw counts to `integer matrix`.
#' 3. Subsets `colData` to include only clean factor columns.
#'    See [`sampleData()`][basejump::sampleData] for more details.
#' 4. Simplifies `metadata` to include only relevant information and updates
#'    `sessionInfo` slot.
#'
#' Note that gene-level counts are required. Alternatively
#' [tximport::summarizeToGene()] can be called to convert transcript-level
#' counts to gene-level. By default, we're using length-scaled TPM, so a
#' corresponding average transcript length matrix isn't necessary. The average
#' transcript length matrix is only necessary when raw counts matrix isn't
#' scaled during tximport call (see `countsFromAbundance` in tximport).
#'
#' @section bcbioRNASeq to DESeqTransform:
#'
#' 1. Coerces to `DESeqDataSet`.
#' 2. Calls [DESeq2::DESeq()].
#' 3. Calls [DESeq2::varianceStabilizingTransformation()].
#'
#' @seealso
#' - `tximport::tximport()`
#' - `DESeq2::DESeqDataSetFromTximport()`.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq to DESeqDataSet ====
#' x <- as(bcb, "DESeqDataSet")
#' names(S4Vectors::mcols(x))
#' class(x)
#' show(x)
#'
#' ## bcbioRNASeq to RangedSummarizedExperiment ====
#' x <- as(bcb, "RangedSummarizedExperiment")
#' slotNames(x)
#' show(x)
#'
#' ## bcbioRNASeq to SummarizedExperiment ====
#' ## Coerce to RangedSummarizedExperiment first.
#' x <- as(bcb, "RangedSummarizedExperiment")
#' x <- as(x, "SummarizedExperiment")
#' class(x)
#' slotNames(x)
#' show(x)
NULL



#' @rdname as
#' @name coerce,bcbioRNASeq,DESeqDataSet-method
setAs(
    from = "bcbioRNASeq",
    to = "DESeqDataSet",
    function(from) {
        validObject(from)
        se <- as(from, "RangedSummarizedExperiment")
        colData(se) <- sampleData(from, clean = TRUE)
        validObject(se)
        to <- .new.DESeqDataSet(se = se)
        interestingGroups(to) <- interestingGroups(from)
        validObject(to)
        to
    }
)



#' @rdname as
#' @name coerce,bcbioRNASeq,DESeqTransform-method
setAs(
    from = "bcbioRNASeq",
    to = "DESeqTransform",
    function(from) {
        validObject(from)
        dds <- as(from, "DESeqDataSet")
        # Expect warning about empty design formula.
        dds <- suppressWarnings(DESeq(dds))
        validObject(dds)
        message("Applying variance stabilizing transformation.")
        dt <- varianceStabilizingTransformation(dds)
        validObject(dt)
        dt
    }
)
