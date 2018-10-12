#' Methods for Coercing an Object to a Class
#'
#' @name coerce
#' @aliases as
#' @family S4 Functions
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
# Note that gene-level counts are required. Mention `summarizeToGene()`.
# By default, we're using length-scaled TPM, so a corresponding average
# transcript length matrix isn't necessary. The average transcript length matrix
# is only necessary when raw counts matrix isn't scaled during tximport call
# (see `countsFromAbundance`).
# @seealso `tximport::tximport`, `DESeq2::DESeqDataSetFromTximport()`.
#' @rdname coerce
#' @name coerce,bcbioRNASeq,DESeqDataSet-method
setAs(
    from = "bcbioRNASeq",
    to = "DESeqDataSet",
    function(from) {
        validObject(from)
        se <- as(from, "RangedSummarizedExperiment")
        # Don't include the metrics columns inside bcbioRNASeq object.
        colData(se) <- sampleData(from, clean = TRUE)
        # FIXME Minimize the metadata we're passing through.
        # FIXME Ensure sessionInfo is updated.
        validObject(se)
        to <- .new.DESeqDataSet(se = se)
        interestingGroups(to) <- interestingGroups(from)
        validObject(to)
        to
    }
)



#' @rdname coerce
#' @name coerce,bcbioRNASeq,DESeqTransform-method
setAs(
    from = "bcbioRNASeq",
    to = "DESeqTransform",
    function(from) {
        validObject(from)
        # First, generate a DESeqDataSet.
        dds <- as(from, "DESeqDataSet")
        validObject(dds)
        # Coerce using a new on-the-fly calculation.
        # We're using VST here because it's faster than rlog.
        message("Applying variance stabilizing transformation.")
        dt <- varianceStabilizingTransformation(dds)
        validObject(dt)
        dt
    }
)



#' @rdname coerce
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



#' @rdname coerce
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



#' @rdname coerce
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
