#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' @name runGSEA
#' @inheritParams general
#' @author Michael Steinbaugh
#'
#' @param gmtFile `string`. MSigDB GMT file path.
NULL



# FIXME Need to point to gmtFile.
runGSEA.DESeqAnalysis <-  # nolint
    function(object, gmtFile) {
        stop("Draft update")
    }



#' @rdname runGSEA
#' @export
setMethod(
    f = "runGSEA",
    signature = signature("DESeqAnalysis"),
    definition = runGSEA.DESeqAnalysis
)
