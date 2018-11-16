#' @name markdown
#' @inherit basejump::markdown
#' @author Michael Steinbaugh
#'
#' @inheritParams params
#' @inheritParams basejump::params
#'
#' @examples
#' data(deseq)
#' x <- DESeqResultsTables(deseq)
#' markdown(x)
NULL



#' @importFrom basejump markdown
#' @aliases NULL
#' @export
basejump::markdown



# bcbioRNASeq ==================================================================
markdown.bcbioRNASeq <-  # nolint
    function(object) {
        rse <- as(object, "RangedSummarizedExperiment")
        sampleData(rse) <- sampleData(object, clean = TRUE)
        markdown(rse)
    }



#' @describeIn markdown Sample metadata table.
#' @export
setMethod(
    f = "markdown",
    signature = signature("bcbioRNASeq"),
    definition = markdown.bcbioRNASeq
)
