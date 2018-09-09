setClassUnion("missingOrNULL", c("missing", "NULL"))



#' @rdname bcbioRNASeq
#' @aliases NULL
#' @exportClass bcbioRNASeq
#' @usage NULL
bcbioRNASeq <- setClass(
    Class = "bcbioRNASeq",
    contains = "RangedSummarizedExperiment"
)



#' DESeq2 Differential Expression Analysis
#'
#' Class containing all elements generated during differential expression
#' analysis with DESeq2.
#'
#' @section Slots:
#'
#'     - `data`: `DESeqDataSet`.
#'     - `transform`: `DESeqTransform`.
#'     - `results`: `list` containing one or more `DESeqResults`.
#'
#' @export
#'
#' @examples
#' x <- new(
#'     "DESeqAnalysis",
#'     data = dds_small,
#'     transform = vst_small,
#'     results = list(res_small)
#' )
#' class(x)
#' slotNames(x)
setClass(
    Class = "DESeqAnalysis",
    slots = list(
        data = "DESeqDataSet",
        transform = "DESeqTransform",
        results = "list"
    )
)
