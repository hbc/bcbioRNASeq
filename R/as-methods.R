#' @name as
#' @aliases coerce
#' @inherit methods::as
#' @importFrom methods coerce
#' @exportMethod coerce
#' @author Michael Steinbaugh
#'
#' @section bcbioRNASeq to DESeqDataSet:
#'
#' 1. Coerces to `RangedSummarizedExperiment`.
#' 2. Rounds raw counts to `integer matrix`.
#' 3. Subsets [colData()] to include only clean factor columns. See
#'    [sampleData()] for more details.
#' 4. Simplifies [metadata()] to include only relevant information and
#'    updates `sessionInfo`.
#'
#' @section bcbioRNASeq to DESeqTransform:
#'
#' 1. Coerces to `DESeqDataSet`.
#' 2. Calls [DESeq2::DESeq2()].
#' 3. Calls [DESeq2::varianceStabilizingTransformation()].
#'
#' @section DESeqAnalysis:
#'
#' Supported coercion methods will extract any of these internal objects:
#'
#' - `DESeqDataSet`.
#' - `DESeqTransform`.
#' - `DESeqResults`. Extracts the first results slotted. Note that this
#'   corresponds to results containing log2 fold change (LFC) values that
#'   *have not been shrunken* using [DESeq2::lfcShrink()].
#'
#' @examples
#' ## bcbioRNASeq to DESeqDataSet ====
#' x <- as(bcb_small, "DESeqDataSet")
#' names(S4Vectors::mcols(x))
#' class(x)
#' show(x)
#'
#' ## bcbioRNASeq to RangedSummarizedExperiment ====
#' x <- as(bcb_small, "RangedSummarizedExperiment")
#' slotNames(x)
#' show(x)
#'
#' ## bcbioRNASeq to SummarizedExperiment ====
#' ## Coerce to RangedSummarizedExperiment first.
#' x <- as(bcb_small, "RangedSummarizedExperiment")
#' x <- as(x, "SummarizedExperiment")
#' class(x)
#' slotNames(x)
#' show(x)
#'
#' ## DESeqAnalysis ====
#' dds <- as(deseq_small, "DESeqDataSet")
#' print(dds)
#' dt <- as(deseq_small, "DESeqTransform")
#' print(dt)
#' ## Pulls the first results slotted.
#' res <- as(deseq_small, "DESeqResults")
#' contrastName(res)
#' summary(res)
NULL



# FIXME Add this to documentation as a section.
# FIXME Minimize the metadata we're passing through.
# FIXME Ensure sessionInfo is updated.

# Note that gene-level counts are required. Mention `summarizeToGene()`.
# By default, we're using length-scaled TPM, so a corresponding average
# transcript length matrix isn't necessary. The average transcript length matrix
# is only necessary when raw counts matrix isn't scaled during tximport call
# (see `countsFromAbundance`).
# @seealso `tximport::tximport`, `DESeq2::DESeqDataSetFromTximport()`.

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



#' @rdname as
#' @name coerce,DESeqAnalysis,DESeqDataSet-method
setAs(
    from = "DESeqAnalysis",
    to = "DESeqDataSet",
    function(from) {
        validObject(from)
        to <- slot(from, "data")
        validObject(to)
        to
    }
)



#' @rdname as
#' @name coerce,DESeqAnalysis,DESeqTransform-method
setAs(
    from = "DESeqAnalysis",
    to = "DESeqTransform",
    function(from) {
        validObject(from)
        to <- slot(from, "transform")
        validObject(to)
        to
    }
)



#' @rdname as
#' @name coerce,DESeqAnalysis,DESeqResults-method
setAs(
    from = "DESeqAnalysis",
    to = "DESeqResults",
    function(from) {
        validObject(from)
        to <- slot(from, "results")[[1L]]
        validObject(to)
        to
    }
)
