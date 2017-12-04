#' Column Data
#'
#' @description
#' Improved assignment method support for [bcbioRNASeq] object.
#'
#' This method support will also update the `colData` inside the `bcbio` and
#' `assays` slots.
#'
#' @rdname colData
#' @name colData
#'
#' @inheritParams AllGenerics
#'
#' @seealso
#' `help("colData", "SummarizedExperiment")`
#'
#' @examples
#' bcb <- examples[["bcb"]]
#' cd <- colData(bcb)
#' cd[["age"]] <- factor(c(14, 30, 14, 30))
#' colData(bcb) <- cd
#' colData(bcb)
#'
#' # These internal objects will also get updated
#' bcbio(bcb, "DESeqDataSet") %>% colData()
#' assays(bcb)[["rlog"]] %>% colData()
#' assays(bcb)[["vst"]] %>% colData()
NULL



# Methods ====
#' @rdname colData
#' @export
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
        slot(x, "colData") <- value
        x
    })
