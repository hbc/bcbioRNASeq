setClassUnion("missingOrNULL", c("missing", "NULL"))



#' @rdname bcbioRNASeq
#' @aliases NULL
#' @exportClass bcbioRNASeq
#' @usage NULL
setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment"
)



#' @rdname DESeqAnalysis
#' @aliases NULL
#' @exportClass DESeqAnalysis
#' @usage NULL
setClass(
    Class = "DESeqAnalysis",
    slots = list(
        DESeqDataSet = "DESeqDataSet",
        DESeqTransform = "DESeqTransform",
        DESeqResults = "list",
        lfcShrink = "list"
    )
)
