#' @rdname bcbioRNASeq
#' @aliases NULL
#' @exportClass bcbioRNASeq
#' @usage NULL
bcbioRNASeq <- setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment"
)

setClassUnion("missingOrNULL", c("missing", "NULL"))
