#' Column Data
#'
#' Improved assignment method support for [bcbioRNASeq] object.
#'
#' @rdname colData
#' @name colData
#'
#' @inheritParams AllGenerics
#' @inherit SummarizedExperiment::colData
#'
#' @examples
#' bcb <- examples[["bcb"]]
#' cd <- colData(bcb)
#' cd[["age"]] <- factor(c(14, 30, 14, 30))
#' colData(bcb) <- cd
#' colData(bcb)
NULL



# Methods ====
#' @rdname colData
#' @export
#' @seealso
#' `getMethod(
#'     "colData<-",
#'     signature(x = "SummarizedExperiment",
#'               value = "DataFrame"))`
setMethod(
    "colData<-",
    signature(x = "bcbioRNASeq", value = "DataFrame"),
    function(x, ..., value) {
        if (nrow(value) != ncol(x)) {
            stop("nrow of supplied 'colData' must equal ncol of object")
        }
        if (!is.null(bcbio(x, "DESeqDataSet"))) {
            colData(bcbio(x, "DESeqDataSet")) <- value
        }
        if (!is.null(assays(x)[["rlog"]])) {
            colData(assays(x)[["rlog"]]) <- value
        }
        if (!is.null(assays(x)[["vst"]])) {
            colData(assays(x)[["vst"]]) <- value
        }
        BiocGenerics:::replaceSlots(x, colData = value, check = FALSE)
    })
