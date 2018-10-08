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
#' # bcbioRNASeq to DESeqDataSet ====
#' x <- as(bcb_small, "DESeqDataSet")
#' names(S4Vectors::mcols(x))
#' class(x)
#' show(x)
#'
#' # bcbioRNASeq to RangedSummarizedExperiment ====
#' x <- as(bcb_small, "RangedSummarizedExperiment")
#' slotNames(x)
#' show(x)
#'
#' # bcbioRNASeq to SummarizedExperiment ====
#' # Coerce to RangedSummarizedExperiment first.
#' x <- as(bcb_small, "RangedSummarizedExperiment")
#' x <- as(x, "SummarizedExperiment")
#' class(x)
#' slotNames(x)
#' show(x)
#'
#' # DESeqAnalysis ====
#' dds <- as(deseq_small, "DESeqDataSet")
#' print(dds)
#' dt <- as(deseq_small, "DESeqTransform")
#' print(dt)
#' # Pulls the first results slotted.
#' res <- as(deseq_small, "DESeqResults")
#' contrastName(res)
#' summary(res)
NULL



# Refer to `DESeq2::DESeqDataSetFromTximport()`.
#' @rdname coerce
#' @name coerce,bcbioRNASeq,DESeqDataSet-method
setAs(
    from = "bcbioRNASeq",
    to = "DESeqDataSet",
    function(from) {
        validObject(from)
        if (metadata(from)[["level"]] != "genes") {
            # Consider adding summarize to gene support here.
            # Consult the tximport vignette if we decide to add this.
            stop("Gene-level counts are required.")
        }
        message(paste0(
            "Coercing bcbioRNASeq to DESeqDataSet with DESeq2 ",
            packageVersion("DESeq2"), "..."
        ))
        to <- .new.DESeqDataSet(se = as(from, "RangedSummarizedExperiment"))
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
        message("Applying variance stabilizing transformation...")
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
